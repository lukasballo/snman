import r5py
import shapely
import pandas as pd
import geopandas as gpd
import numpy as np
import datetime

import os
import snman
from snman import osmnx_customized as oxc
from snman.constants import *

import warnings
warnings.filterwarnings('ignore')

PERIMETER = '_accessibility_debug_euler'

# Set these paths according to your own setup
data_directory = os.path.join(os.sep, 'cluster', 'home', 'lballo', 'data')
inputs_path = os.path.join(data_directory, 'inputs')
outputs_path = os.path.join(data_directory, 'outputs', PERIMETER)

CRS_internal = 2056      # for Zurich

# --------------------------------------------------------------------

statent = gpd.read_parquet(
    os.path.join(inputs_path, 'switzerland', 'switzerland', 'statent', 'STATENT_2021_ebc_10_reduced_fields.gzip')
)
statent.rename(columns={'RELI': 'id'}, inplace=True)

statpop = gpd.read_parquet(
    os.path.join(inputs_path, 'switzerland', 'switzerland', 'statpop', 'statpop2017_with_statent_reduced_columns.gzip')
)

# --------------------------------------------------------------------

print('read default transport network')
transport_network_default = r5py.TransportNetwork(
    os.path.join(outputs_path, 'before_oneway_links_default.osm.pbf'),
    os.path.join(inputs_path, 'switzerland', 'switzerland', 'gtfs', 'gtfs_fp2024_2024-06-27_mod.zip')
    #os.path.join(inputs_path, 'switzerland', 'switzerland', 'gtfs', 'vbz_2024.zip')
)

print('read transport network cycling')
transport_network_cycling = r5py.TransportNetwork(
    os.path.join(outputs_path, 'before_oneway_links_cycling.osm.pbf'),
)

print('read transport network cycling after')
transport_network_cycling_after = r5py.TransportNetwork(
    os.path.join(outputs_path, 'after_oneway_links_cycling.osm.pbf'),
)

print('read transport network car')
transport_network_car = r5py.TransportNetwork(
    os.path.join(outputs_path, 'tt.osm.pbf'),
)

print('read perimeter')
perimeter = snman.io.import_geofile_to_gdf(
    os.path.join(inputs_path, 'perimeters', 'perimeters.shp'), index='id'
).geometry['ebc_zrh_v01'].simplify(500, preserve_topology=True)
statent_in_perimeter = statent[statent.within(perimeter)]

# --------------------------------------------------------------------

print('travel time matrix before')
cells = list(statent_in_perimeter['id'].sample(n=10, random_state=0))

def miniwrapper(cell):
    print(f'calculating cell {cell}')
    return snman.accessibility.calculate_accessibility_for_statent_cell(
        transport_network_default,
        transport_network_car,
        transport_network_cycling,
        statpop, statent, cell,
        distance_limit=50*1000,
        population_sample=0.05,
        destinations_sample=0.02
    )

accessibility_before = pd.concat(
    map(miniwrapper, cells),
    ignore_index=True
)

print('exporting file')
snman.io.export_gdf(
    gpd.GeoDataFrame(accessibility_before),
    os.path.join(outputs_path, 'accessibility_before.gpkg')
)

# --------------------------------------------------------------------

print('travel time matrix after')
cells = list(statent_in_perimeter['id'].sample(n=10, random_state=0))

def miniwrapper(cell):
    print(f'calculating cell {cell}')
    return snman.accessibility.calculate_accessibility_for_statent_cell(
        transport_network_default,
        transport_network_car,
        transport_network_cycling_after,
        statpop, statent, cell,
        distance_limit=50*1000,
        population_sample=0.05,
        destinations_sample=0.02
    )

accessibility_after = pd.concat(
    map(miniwrapper, cells),
    ignore_index=True
)

print('exporting file')
snman.io.export_gdf(
    gpd.GeoDataFrame(accessibility_after),
    os.path.join(outputs_path, 'accessibility_before.gpkg')
)

# -------------------------------------------------------------------

print('merge files and export')

accessibility_before_after = pd.merge(
    accessibility_before, accessibility_after[['accessibility']], left_index=True, right_index=True, suffixes=['_before', '_after']
)

accessibility_before_after['accessibility_after-before'] =(
        accessibility_before_after['accessibility_after'] - accessibility_before_after['accessibility_before'])

snman.io.export_gdf(
    gpd.GeoDataFrame(accessibility_before_after),
    os.path.join(outputs_path, 'accessibility_before_after.gpkg')
)
