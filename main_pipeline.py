# import osmnx as ox
# import matplotlib.pyplot as plt
# import fiona

from snman import osmnx as ox, config

import snman

print('Starting...')

# Variables
osm_tags = ['bridge', 'tunnel', 'layer', 'oneway', 'ref', 'name',
                    'highway', 'maxspeed', 'service', 'access', 'area',
                    'landuse', 'width', 'est_width', 'junction', 'surface',
                    'lanes', 'lanes:forward', 'lanes:backward',
                    'cycleway', 'cycleway:both', 'cycleway:left', 'cycleway:right',
                    'sidewalk', 'sidewalk:left', 'sidewalk:right']
ox.utils.config(useful_tags_way=osm_tags)

# TODO: Include highway=living_street and service roads
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

print('Get data from OSM server')
street_graph = ox.graph_from_place(
    #'Zurich, Zurich, Switzerland',
    #'Seebach, Zurich, Switzerland',
    'Altstetten, Zurich, Switzerland',
    custom_filter=custom_filter,
    simplify=True,
    simplify_strict=False,
    retain_all=True,
    one_edge_per_direction=False,
)

print('Consolidate intersections')
street_graph = ox.simplification.consolidate_intersections(
    # TODO: Work with a projected CRS in consolidating intersections (right now, this only works with EPSG:4326)
    street_graph, tolerance=0.00020, rebuild_graph=True, dead_ends=True, reconnect_edges=True
)

print('Convert CRS of street graph to 2056')
snman.convert_crs_of_street_graph(street_graph, str(config.crs))

print('Generate lanes')
snman.generate_lanes(street_graph)

print('Normalize edge directions, enforce direction from lower to higher node id')
snman.normalize_edge_directions(street_graph)

print('Convert into an undirected graph')
street_graph = ox.utils_graph.get_undirected(street_graph)

print('Identify hierarchy')
snman.add_hierarchy(street_graph)

print('Merge parallel and consecutive edges, repeat a few times')
# TODO: Merge parallel lanes even if there is a one-sided intersection (Example: Bottom station of Seilbahn Rigiblick)
for i in range(5):
    snman.merge_parallel_edges(street_graph)
    snman.merge_consecutive_edges(street_graph)
    pass

print('Add public transport')
pt_network = snman.import_shp_to_gdf("C:/DATA/CLOUD STORAGE/polybox/Data/stadt_zuerich_open_data/Linien_des_offentlichen_Verkehrs_-OGD/ZVV_LINIEN_GEN_L.shp")
snman.match_pt(street_graph, pt_network)

print('Add lane stats')
snman.generate_lane_stats(street_graph)

print('Export SHP to make manual corrections on the network in QGIS')
snman.export_streetgraph_to_shp(street_graph)

print('Done!')
