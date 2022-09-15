# import osmnx as ox
# import matplotlib.pyplot as plt
# import fiona

from snman import osmnx as ox, config

import snman

# Steps in this file:
# - Acquisition
# - Simplification
# - Reconstruct the lane composition
# - Merge parallel links (pending)
# - Export to SHP for manual corrections

print('Starting...')

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

# Get data from OSM server
street_graph = ox.graph_from_place(
    'Oerlikon, Zurich, Switzerland',
    custom_filter=custom_filter,
    simplify=True,
    simplify_strict=False,
    retain_all=True,
    one_edge_per_direction=False,
)

# Consolidate intersections
street_graph = ox.simplification.consolidate_intersections(
    # TODO: Work with a projected CRS in consolidating intersections
    # TODO: Fix problems with "edge u-v is not in the graph"
    street_graph, tolerance=0.00020, rebuild_graph=True, dead_ends=True, reconnect_edges=True
)

#street_graph = ox.utils_graph.get_undirected(street_graph)

# Generate lanes
snman.generate_lanes(street_graph)

# Merge parallel edges
snman.merge_parallel_edges(street_graph)

# Merge consecutive edges
snman.merge_consecutive_edges(street_graph)

# Merge parallel edges again
snman.merge_parallel_edges(street_graph)

# Add lane stats
snman.generate_lane_stats(street_graph)

# Export SHP to make manual corrections on the network in QGIS
snman.export_to_shp(street_graph)

print('Done!')
