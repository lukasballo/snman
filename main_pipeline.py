# import osmnx as ox
# import matplotlib.pyplot as plt
# import fiona
import geopandas as gpd
import copy
import shapely as shp
import networkx as nx
import osmnx as ox
import snman

print('Starting...')

# Constants
INTERSECTION_TOLERANCE = 10

inputs_path = 'C:/DATA/CLOUD STORAGE/polybox/Research/SNMan/SNMan Shared/inputs/'
export_path = 'C:/DATA/CLOUD STORAGE/polybox/Research/SNMan/SNMan Shared/qgis_previews/'

ox.utils.config(useful_tags_way=snman.constants.OSM_TAGS)

print('Get data from OSM server')
G = snman.osmnx.graph_from_place(
    'Seebach, Zurich, Switzerland',
    custom_filter=snman.constants.OSM_FILTER,
    simplify=True,
    simplify_strict=False,
    retain_all=True,
    one_edge_per_direction=False,
)

print('Convert CRS of street graph to 2056')
snman.convert_crs_of_street_graph(G, snman.constants.CRS)
nodes = copy.copy(G.nodes)

print('Load regions')
regions = snman.load_regions(inputs_path + 'regions/regions.shp', default_region=True, street_graph=G)

print('Load manual intersections')
given_intersections_gdf = snman.load_intersections(
    inputs_path + 'intersection_polygons/intersection_polygons.shp',
    inputs_path + 'intersection_points/intersection_points.shp'
)

print('Detect intersections')
intersections_gdf = snman.simplification.merge_nodes_geometric(G, 10, given_intersections_gdf=given_intersections_gdf)

print('Save intersection geometries into a file')
snman.export_gdf(intersections_gdf, export_path + 'intersections.gpkg', columns=['geometry'])
snman.export_gdf(intersections_gdf, export_path + 'intersections_polygons.gpkg', columns=['geometry'])
snman.export_gdf(gpd.GeoDataFrame(intersections_gdf, geometry='point_geometry'), export_path + 'intersections_points.gpkg', columns=['point_geometry'])

if 0:
    # must be run a few times for including buffers of newly added nodes
    for i in range(6):
        print('Split through edges in intersections')
        intersections = snman.split_through_edges_in_intersections(G, INTERSECTION_TOLERANCE, regions=regions)

        print('Update precalculated attributes')
        snman.update_precalculated_attributes(G)

        print('Add connections between components in intersections')
        snman.connect_components_in_intersections(G, intersections)



if 1:
    print('Save raw street graph')
    snman.export_streetgraph(G, export_path + 'raw_edges.gpkg', export_path + 'raw_nodes.gpkg')

if 1:
    print('Consolidate intersections')
    G = snman.simplification.consolidate_intersections(
        G, intersections_gdf,
        reconnect_edges=True
    )

if 1:
    print('Generate lanes')
    snman.generate_lanes(G)

if 1:
    print('Normalize edge directions, enforce direction from lower to higher node id')
    snman.normalize_edge_directions(G)

if 1:
    print('Convert into an undirected graph')
    G = ox.utils_graph.get_undirected(G)

if 1:
    print('Identify hierarchy')
    snman.add_hierarchy(G)

if 1:
    print('Merge parallel and consecutive edges, repeat a few times')
    for i in range(5):
        snman.merge_parallel_edges(G)
        snman.merge_consecutive_edges(G)
        pass

if 1:
    print('Add lane stats')
    snman.generate_lane_stats(G)

if 0:
    #TODO: Improve performance with geodataframe operations
    print('Add public transport')
    pt_network = snman.import_shp_to_gdf("C:/DATA/CLOUD STORAGE/polybox/Research/SNMan/SNMan Shared/stadt_zuerich_open_data/Linien_des_offentlichen_Verkehrs_-OGD/ZVV_LINIEN_GEN_L.shp")
    snman.match_pt(G, pt_network)

if 0:
    print('Update OSM tags')
    snman.update_osm_tags(G)

if 1:
    print('Set given lanes')
    snman.set_given_lanes(G)

if 1:
    print('Create directed graph of given lanes')
    G_minimal_graph_input = snman.create_given_lanes_graph(G)

if 1:
    print('Export network without lanes')
    snman.export_streetgraph(G, export_path + 'edges.gpkg', export_path + 'nodes.gpkg')

if 1:
    print('Export network with lanes')
    #TODO: Fix problems with saving lanes as a GeoPackage
    snman.export_streetgraph_with_lanes(G, 'ln_desc', export_path + 'edges_lanes.shp')

if 1:
    print('Export network with given lanes')
    snman.export_streetgraph_with_lanes(G, 'given_lanes', export_path + 'edges_given_lanes.shp')

if 1:
    print('Export given lanes')
    snman.export_streetgraph(G_minimal_graph_input, export_path + 'given_lanes.gpkg', export_path + 'given_lanes_nodes.gpkg')

if 0:
    print('Export OSM XML')
    snman.export_osm_xml(G, export_path + 'new_network.osm',{
        'lanes', 'lanes:forward', 'lanes:backward', 'lanes:both_ways',
        'cycleway', 'cycleway:lane', 'cycleway:left', 'cycleway:left:lane', 'cycleway:right', 'cycleway:right:lane',
        'bus:lanes:backward', 'bus:lanes:forward', 'vehicle:lanes:backward', 'vehicle:lanes:forward',
        'maxspeed'
    })

print('Link elimination')
G_minimal_graph_output = snman.link_elimination(G_minimal_graph_input)

print('Export minimal graph - output')
snman.export_streetgraph(G_minimal_graph_output, export_path + 'minimal_graph_out_edges.gpkg', export_path + 'minimal_graph_out_nodes.gpkg')


print('Done!')
