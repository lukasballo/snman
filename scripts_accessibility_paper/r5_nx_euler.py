import warnings
warnings.filterwarnings('ignore')

import geopandas
import r5py
import shapely
import pandas as pd
import geopandas as gpd
import numpy as np
import datetime
import copy

import os
import snman
from snman import osmnx_customized as oxc
from snman.constants import *

PERIMETER = '_accessibility_debug_euler'

# Set these paths according to your own setup
data_directory = os.path.join(os.sep, 'cluster', 'home', 'lballo', 'data')
inputs_path = os.path.join(data_directory, 'inputs')
outputs_path = os.path.join(data_directory, 'outputs', PERIMETER)

CRS_internal = 2056      # for Zurich

# --------------------------------------------------------------------

REBUILD_INPUTS = True
STATE = 'before'
N_ORIGIN_CELLS = 100
SAMPLE_RANGE = [0,10]

# --------------------------------------------------------------------

# read transport network and timetable for r5
if 1:
    r5_transit_network = r5py.TransportNetwork(
        os.path.join(outputs_path, 'before_oneway_links_default.osm.pbf'),
        #os.path.join(inputs_path, 'switzerland', 'switzerland', 'gtfs', 'vbz_2024.zip'),
        os.path.join(inputs_path, 'switzerland', 'switzerland', 'gtfs', 'gtfs_fp2024_2024-06-27_mod.zip'),
    )

if REBUILD_INPUTS:
    G_modes = {}
    L_modes = {}

if REBUILD_INPUTS:
    # Load small street graph
    G = snman.io.load_street_graph(
        os.path.join(data_directory, outputs_path, 'G_edges.gpkg'),
        os.path.join(data_directory, outputs_path, 'G_nodes.gpkg'),
        crs=CRS_internal
    )

    if STATE == 'before':
        ln_desc = KEY_LANES_DESCRIPTION
    else:
        ln_desc = KEY_LANES_DESCRIPTION_AFTER

    # create mode graphs for other modes
    for mode in [MODE_CYCLING, MODE_FOOT]:
        print(f'Make street and lane graphs for {mode}')
        G_modes[mode] = copy.deepcopy(G)
        snman.street_graph.filter_lanes_by_modes(
            G_modes[mode], [mode], lane_description_key=ln_desc
        )
        L_modes[mode] = snman.lane_graph.create_lane_graph(
            G_modes[mode], lanes_attribute=ln_desc
        )

if REBUILD_INPUTS:
    # load the new car street graph created from MATSim output
    G_modes[MODE_PRIVATE_CARS] = snman.io.load_street_graph(
        os.path.join(outputs_path, f'tt_{STATE}_edges.gpkg'),
        os.path.join(outputs_path, f'tt_{STATE}_nodes.gpkg'),
        crs=CRS_internal,
        unstringify_attributes={'lanes': snman.space_allocation.space_allocation_from_string}
    )

    # we must avoid that the start/end point of trips are matched onto nodes that not accessible by car
    # therefore, we reduce the graph only to those links and nodes that are accessible by car
    snman.street_graph.filter_lanes_by_modes(
        G_modes[MODE_PRIVATE_CARS], [MODE_PRIVATE_CARS], lane_description_key='lanes'
    )

    L_modes[MODE_PRIVATE_CARS] = snman.lane_graph.create_lane_graph(
        G_modes[MODE_PRIVATE_CARS], lanes_attribute='lanes', cast_attributes={'matsim_tt_cars': 'car_18:00'}
    )

# --------------------------------------------------------------------

# Import Statent grid
if REBUILD_INPUTS:
    statent = gpd.read_parquet(
        os.path.join(inputs_path, 'switzerland', 'switzerland', 'statent', 'STATENT_2021_ebc_10_reduced_fields.gzip')
    )
    statent.rename(columns={'RELI': 'id'}, inplace=True)

if REBUILD_INPUTS:
    # read statpop that has been pre-joined with statent
    statpop_with_statent_ids = gpd.read_parquet(
        os.path.join(inputs_path, 'switzerland', 'switzerland', 'statpop', 'statpop2017_with_statent_reduced_columns.gzip')
    )

perimeter = snman.io.import_geofile_to_gdf(
    os.path.join(inputs_path, 'perimeters', 'perimeters.shp'), index='id'
).geometry['ebc_zrh_v01'].simplify(500, preserve_topology=True)
statent_in_perimeter = statent[statent.within(perimeter)]

# --------------------------------------------------------------------

cells = list(statent_in_perimeter['id'].sample(n=N_ORIGIN_CELLS, random_state=0))[SAMPLE_RANGE[0]:SAMPLE_RANGE[1]]

def miniwrapper(args):
    cell, i, N = args
    print(f'calculating cell {cell} ({i}/{N})')
    return snman.accessibility.calculate_accessibility_for_statent_cell(
        r5_transit_network,
        G_modes,
        L_modes,
        statpop_with_statent_ids, statent, cell,
        datetime.datetime(2024, 2, 22, 18, 00),
        distance_limit=70*1000,
        population_sample=0.05,
        destinations_sample=0.1
    )

accessibility = pd.concat(
    map(
        miniwrapper,
        zip(
            cells,
            range(len(cells)),
            [len(cells)] * len(cells)
        )
    ),
    ignore_index=True
)


print('exporting file')
snman.io.export_gdf(
    accessibility,
    os.path.join(outputs_path, f'accessibility_{STATE}.gpkg')
)
# --------------------------------------------------------------------
