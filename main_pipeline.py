# import osmnx as ox
# import matplotlib.pyplot as plt
# import fiona
import geopandas as gpd
import copy
import shapely as shp
import networkx as nx

from snman import osmnx as ox

import snman

print('Starting...')

# Constants
INTERSECTION_TOLERANCE = 10

osm_tags = ['bridge', 'tunnel', 'layer', 'oneway', 'ref', 'name',
                    'highway', 'maxspeed', 'service', 'access', 'area',
                    'landuse', 'width', 'est_width', 'junction', 'surface',
                    'lanes', 'lanes:forward', 'lanes:backward',
                    'cycleway', 'cycleway:both', 'cycleway:left', 'cycleway:right',
                    'bicycle', 'bicycle:conditional',
                    'sidewalk', 'sidewalk:left', 'sidewalk:right', 'foot',
                    'psv', 'bus', 'bus:lanes:forward', 'bus:lanes:backward',
                    'vehicle:lanes:backward', 'vehicle:lanes:forward',
                    'footway',
            ]
export_path = 'C:/DATA/CLOUD STORAGE/polybox/Research/SNMan/SNMan Shared/qgis_previews/'
CRS = 'epsg:2056'      # CH1905+ projected CRS

ox.utils.config(useful_tags_way=osm_tags)

custom_filter = [
    (
        f'["highway"]["area"!~"yes"]["access"!~"private"]'
        f'["highway"!~"abandoned|bridleway|bus_guideway|construction|corridor|elevator|'
        f'escalator|planned|platform|proposed|raceway|'
        f'steps"]'
        f'["service"!~"alley|driveway|emergency_access|parking|parking_aisle|private"]'
        f'["access"!~"no"]'
    ),
    (
        f'["bicycle"]["bicycle"!~"no"]["area"!~"yes"]'
    ),
    (
        f'["bicycle:conditional"]["area"!~"yes"]'
    ),
    (
        f'["bus"="yes"]["area"!~"yes"]'
    ),
    (
        f'["psv"="yes"]["area"!~"yes"]'
    ),
    (
        f'["highway"="footway"]["area"!~"yes"]'
    )
]

print('Get data from OSM server')
street_graph = ox.graph_from_place(
    #'Zurich, Zurich, Switzerland',
    #'Wipkingen, Zurich, Switzerland',
    #'Industriequartier, Zurich, Switzerland',
    'Seebach, Zurich, Switzerland',
    #'Altstetten, Zurich, Switzerland',
    #'Altstadt, Zurich, Switzerland',
    #'Oberstrass, Zurich, Switzerland',
    custom_filter=custom_filter,
    simplify=True,
    simplify_strict=False,
    retain_all=True,
    one_edge_per_direction=False,
)

if 1:
    print('Convert CRS of street graph to 2056')
    snman.convert_crs_of_street_graph(street_graph, CRS)
    nodes = copy.copy(street_graph.nodes)

if 1:
    print('Split through edges in intersections')
    # must be run a few times for including buffers of newly added nodes
    # TODO: Resolve shapely deprecation warnings
    intersections = None
    for i in range(6):
        intersections, a = snman.split_through_edges_in_intersections(street_graph, INTERSECTION_TOLERANCE)

    print('Update precalculated attributes')
    snman.update_precalculated_attributes(street_graph)

    print('Add connections between components in intersections')
    snman.connect_components_in_intersections(street_graph, intersections)

    print('Save intersection geometries into a file')
    snman.export_gdf(intersections, export_path + 'intersections.gpkg', columns=['geometry'])

if 1:
    print('Save raw street graph')
    snman.export_streetgraph(street_graph, export_path + 'raw_edges.gpkg', export_path + 'raw_nodes.gpkg')

if 1:
    print('Consolidate intersections')
    street_graph = ox.simplification.consolidate_intersections(
        street_graph, tolerance=INTERSECTION_TOLERANCE, rebuild_graph=True, dead_ends=True, reconnect_edges=True
    )

if 1:
    print('Generate lanes')
    snman.generate_lanes(street_graph)

if 1:
    print('Normalize edge directions, enforce direction from lower to higher node id')
    snman.normalize_edge_directions(street_graph)

#snman.export_streetgraph(street_graph, export_path + 'edges_early.gpkg', export_path + 'nodes_early.gpkg')

if 1:
    print('Convert into an undirected graph')
    street_graph = ox.utils_graph.get_undirected(street_graph)

if 1:
    print('Identify hierarchy')
    snman.add_hierarchy(street_graph)

if 1:
    print('Merge parallel and consecutive edges, repeat a few times')
    for i in range(5):
        snman.merge_parallel_edges(street_graph)
        snman.merge_consecutive_edges(street_graph)
        pass

if 1:
    print('Add lane stats')
    snman.generate_lane_stats(street_graph)

if 0:
    #TODO: Improve performance with geodataframe operations
    print('Add public transport')
    pt_network = snman.import_shp_to_gdf("C:/DATA/CLOUD STORAGE/polybox/Research/SNMan/SNMan Shared/stadt_zuerich_open_data/Linien_des_offentlichen_Verkehrs_-OGD/ZVV_LINIEN_GEN_L.shp")
    snman.match_pt(street_graph, pt_network)

if 1:
    print('Add lane stats')
    snman.generate_lane_stats(street_graph)

if 1:
    print('Update OSM tags')
    snman.update_osm_tags(street_graph)

if 1:
    print('Set given lanes')
    snman.set_given_lanes(street_graph)

if 1:
    print('Create directed graph of given lanes')
    given_lanes_graph = snman.create_given_lanes_graph(street_graph)

if 1:
    print('Export network without lanes')
    snman.export_streetgraph(street_graph, export_path + 'edges.gpkg', export_path + 'nodes.gpkg')

if 1:
    print('Export network with lanes')
    #TODO: Fix problems with saving lanes as a GeoPackage
    snman.export_streetgraph_with_lanes(street_graph, 'ln_desc', export_path + 'edges_lanes.shp')

if 1:
    print('Export network with given lanes')
    snman.export_streetgraph_with_lanes(street_graph, 'given_lanes', export_path + 'edges_given_lanes.shp')

if 0:
    print('Export given lanes')
    snman.export_streetgraph(given_lanes_graph, export_path + 'given_lanes.gpkg', export_path + 'given_lanes_nodes.gpkg')

if 1:
    print('Export OSM XML')
    snman.export_osm_xml(street_graph, export_path + 'new_network.osm',{
        'lanes', 'lanes:forward', 'lanes:backward', 'lanes:both_ways',
        'cycleway', 'cycleway:lane', 'cycleway:left', 'cycleway:left:lane', 'cycleway:right', 'cycleway:right:lane',
        'bus:lanes:backward', 'bus:lanes:forward', 'vehicle:lanes:backward', 'vehicle:lanes:forward',
        'maxspeed'
    })


print('Done!')
