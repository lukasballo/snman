# import osmnx_ebc as ox
# import matplotlib.pyplot as plt
# import fiona
import osmnx_ebc as ox
import networkx as nx
import math
import config
import re
import pandas as pd

print('Starting...')

# --------------------------------------------------------
# Variables

osm_tags = ['bridge', 'tunnel', 'oneway', 'ref', 'name',
                    'highway', 'maxspeed', 'service', 'access', 'area',
                    'landuse', 'width', 'est_width', 'junction', 'surface',
                    'lanes', 'lanes:forward', 'lanes:backward',
                    'cycleway', 'cycleway:both', 'cycleway:left', 'cycleway:right',
                    'sidewalk', 'sidewalk:left', 'sidewalk:right']
ox.utils.config(useful_tags_way=osm_tags)

custom_filter = [
    (
        f'["highway"]["area"!~"yes"]["access"!~"private"]'
        f'["highway"!~"abandoned|bridleway|bus_guideway|construction|corridor|cycleway|elevator|'
        f'escalator|footway|path|pedestrian|planned|platform|proposed|raceway|service|'
        f'steps|track"]'
        f'["motor_vehicle"!~"no"]["motorcar"!~"no"]'
        f'["service"!~"alley|driveway|emergency_access|parking|parking_aisle|private"]'
    ),
    (
        f'["bicycle"="designated"]'
    )
]

# --------------------------------------------------------
# Get data from OSM server
street_graph = ox.graph_from_place(
    'Seebach, Kreis 11, Zurich, Switzerland',
    custom_filter=custom_filter,
    simplify=True,
    simplify_strict=False,
    retain_all=True,
    one_edge_per_direction=False,
)

# Consolidate intersections
street_graph = ox.simplification.consolidate_intersections(
    street_graph, tolerance=0.00015, rebuild_graph=True, dead_ends=True, reconnect_edges=True
)

# Convert to GeoDataFrame
nodes, edges = ox.graph_to_gdfs(street_graph)

# Function for reverse-engineering the lane configuration
def generate_lanes_description(row):
    lanes_list = []
    row = row.fillna(-1)

    n_motorized_lanes = int(row['lanes'])
    n_motorized_lanes_forward = int(row['lanes:forward'])
    n_motorized_lanes_backward = int(row['lanes:backward'])

    if row['highway'] == 'footway' or row['highway'] == 'path':
        lanes_list.append('c-')
    else:
        if row['oneway'] == 1:
            if row['cycleway:left'] == 'lane' or row['cycleway:both'] == 'lane':
                lanes_list.append('c>')
            if n_motorized_lanes >= 1:
                lanes_list.extend(['m>'] * n_motorized_lanes)
            else:
                lanes_list.extend(['m>'])
            if row['cycleway:right'] == 'lane' or row['cycleway:both'] == 'lane' or row['cycleway'] == 'lane':
                lanes_list.append('c>')
        else:
            if row['cycleway:left'] == 'lane' or row['cycleway:both'] or row['cycleway'] == 'lane':
                lanes_list.append('c<')
            if n_motorized_lanes >= 2:
                if n_motorized_lanes_backward >= 0 and n_motorized_lanes_forward >= 0:
                    lanes_list.extend(['m<'] * n_motorized_lanes_backward)
                    lanes_list.extend(['m>'] * n_motorized_lanes_forward)
                else:
                    lanes_list.extend(['m<'] * math.floor(n_motorized_lanes / 2))
                    lanes_list.extend(['m>'] * math.ceil(n_motorized_lanes / 2))
            else:
                lanes_list.append('m-')
            if row['cycleway:right'] == 'lane' or row['cycleway:both'] or row['cycleway'] == 'lane':
                lanes_list.append('c>')

    return ' | '.join(lanes_list)

# Run the reverse-engineering of lanes on all edges
edges['ln_desc'] = edges.apply(generate_lanes_description, axis=1)

# Calculate lane counts and overall width
def lane_stats(row):
    lanes_description = row['ln_desc']
    lanes_list = lanes_description.split(' | ')
    n_lanes_cycling = lanes_list.count('c>') + lanes_list.count('c<') + lanes_list.count('c-') * 2
    n_lanes_motorized = lanes_list.count('m>') + lanes_list.count('m<') + lanes_list.count('m-') * 1.5
    width_cycling_m = n_lanes_cycling * config.lane_width_cycling_m
    width_motorized_m = n_lanes_motorized * config.lane_width_motorized_m
    width_total_m = width_cycling_m + width_motorized_m
    proportion_cycling = width_cycling_m / width_total_m

    res = {
        'n_ln_cyc': n_lanes_cycling,
        'n_ln_mot': n_lanes_motorized,
        'w_cyc_m': width_cycling_m,
        'w_mot_m': width_motorized_m,
        'w_tot_m': width_total_m,
        'prop_cyc': proportion_cycling
    }

    return res


def rowfunc(row):
    res = pd.Series(
        lane_stats(row).values()
    )
    #print(res)
    return res

edges[[
    'n_ln_cyc',
    'n_ln_mot',
    'w_cyc_m',
    'w_mot_m',
    'w_tot_m',
    'prop_cyc'
]] = edges.apply(rowfunc, axis=1)

# Export SHP
edges.to_file(config.data_path + 'edges.shp')

# Make manual corrections on the network in QGIS
# *** do in QGIS ***


print('Done!')
