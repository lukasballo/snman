import geopandas as gpd
import copy
import shapely as shp
import networkx as nx
import snman.osmnx_customized as oxc
import snman

print('Starting...')

# Constants
INTERSECTION_TOLERANCE = 10
inputs_path = 'C:/DATA/CLOUD STORAGE/polybox/Research/SNMan/SNMan Shared/inputs/'
export_path = 'C:/DATA/CLOUD STORAGE/polybox/Research/SNMan/SNMan Shared/qgis_previews/'
oxc.settings.useful_tags_way = snman.constants.OSM_TAGS

# =====================================================================================
# LOAD DATA
# =====================================================================================

print('Load perimeters')
perimeters = snman.load_perimeters(inputs_path + 'perimeters/perimeters.shp')

print('Get data from OSM server')
# At this step, simplification means only removing degree=2 edges
G = oxc.graph_from_polygon(
    perimeters.loc['zollikerberg']['geometry'],
    custom_filter=snman.constants.OSM_FILTER,
    simplify=True,
    simplify_strict=False,
    retain_all=True,
    one_edge_per_direction=False
)

print('Prepare graph')
snman.prepare_graph(G)

print('Convert CRS of street graph to 2056')
snman.convert_crs_of_street_graph(G, snman.constants.CRS)
nodes = copy.copy(G.nodes)

print('Load regions')
regions = snman.load_regions(inputs_path + 'regions/regions.shp', default_tolerance=10, street_graph=G)

print('Load manual intersections')
given_intersections_gdf = snman.load_intersections(
    inputs_path + 'intersection_polygons/intersection_polygons.shp',
    inputs_path + 'intersection_points/intersection_points.shp'
)

# =====================================================================================
# CONSOLIDATE INTERSECTIONS
# =====================================================================================

print('Detect intersections')
intersections_gdf = snman.simplification.merge_nodes_geometric(
    G, INTERSECTION_TOLERANCE,
    given_intersections_gdf=given_intersections_gdf,
    regions=regions
)

print('Save intersection geometries into a file')
snman.export_gdf(intersections_gdf, export_path + 'intersections_polygons.gpkg', columns=['geometry'])

if 1:
    # must be run a few times for including buffers of newly added nodes
    for i in range(3):
        print('Split through edges in intersections')
        intersections = snman.split_through_edges_in_intersections(G, intersections_gdf)

        print('Add layers to nodes')
        snman.graph_tools._add_layers_to_nodes(G)

        print('Update precalculated attributes')
        snman.update_precalculated_attributes(G)

        print('Detect intersections')
        intersections_gdf = snman.simplification.merge_nodes_geometric(
            G, INTERSECTION_TOLERANCE,
            given_intersections_gdf=given_intersections_gdf,
            regions=regions
        )

        print('Add connections between components in intersections')
        snman.connect_components_in_intersections(G, intersections_gdf, separate_layers=True)

    print('Save intersection geometries into a file')
    snman.export_gdf(intersections_gdf, export_path + 'intersections_polygons.gpkg', columns=['geometry'])

if 1:
    print('Save raw street graph')
    snman.export_streetgraph(G, export_path + 'raw_edges.gpkg', export_path + 'raw_nodes.gpkg')

if 1:
    print('Consolidate intersections')
    G = snman.simplification.consolidate_intersections(
        G, intersections_gdf,
        reconnect_edges=True
    )

# =====================================================================================
# ENRICH AND ADJUST GRAPH
# =====================================================================================

if 1:
    print('Generate lanes')
    snman.generate_lanes(G)

if 1:
    print('Normalize edge directions, enforce direction from lower to higher node id')
    snman.normalize_edge_directions(G)

if 1:
    print('Convert into an undirected graph')
    G = oxc.utils_graph.get_undirected(G)

if 1:
    print('Identify hierarchy')
    snman.add_hierarchy(G)

# =====================================================================================
# CONSOLIDATE PARALLEL AND CONSECUTIVE EDGES
# =====================================================================================

if 1:
    print('Merge parallel and consecutive edges, repeat a few times')
    for i in range(5):
        snman.merge_parallel_edges(G)
        snman.merge_consecutive_edges(G)
        pass

if 1:
    print('Simplify link geometries')
    for id, edge in G.edges.items():
        edge['geometry'] = edge['geometry'].simplify(25, preserve_topology=False)

if 1:
    print('Add lane stats')
    snman.generate_lane_stats(G)

# =====================================================================================
# ENRICH
# =====================================================================================

if 0:
    #TODO: Improve performance
    print('Add public transport')
    pt_network = snman.import_shp_to_gdf("C:/DATA/CLOUD STORAGE/polybox/Research/SNMan/SNMan Shared/stadt_zuerich_open_data/Linien_des_offentlichen_Verkehrs_-OGD/ZVV_LINIEN_GEN_L.shp")
    snman.match_pt(G, pt_network)

if 1:
    print('Update OSM tags')
    snman.update_osm_tags(G)

if 1:
    print('Add elevation')
    spn = oxc.stats.count_streets_per_node(G, nodes=G.nodes)
    nx.set_node_attributes(G, values=spn, name="street_count")
    G = oxc.elevation.add_node_elevations_raster(G, inputs_path + 'ch_dhm_25/2056/ch_dhm_2056.tif', cpus=1)
    G = oxc.elevation.add_edge_grades(G, add_absolute=False)

# =====================================================================================
# GIVEN LANES
# =====================================================================================

if 1:
    print('Set given lanes')
    snman.set_given_lanes(G, bidirectional_for_dead_ends=False)

if 1:
    print('Create directed graph of given lanes')
    G_minimal_graph_input = snman.create_given_lanes_graph(G, hierarchies_to_remove={snman.hierarchy.HIGHWAY})

# =====================================================================================
# VARIA
# =====================================================================================

if 1:
    print('Keep only the largest connected component')
    snman.graph_tools.add_connected_component_ids(G)
    G = snman.graph_tools.keep_only_the_largest_connected_component(G)

# =====================================================================================
# EXPORT
# =====================================================================================

if 1:
    print('Export network with given lanes')
    snman.export_streetgraph_with_lanes(G, 'given_lanes', export_path + 'edges_given_lanes.gpkg')

if 1:
    print('Export given lanes')
    snman.export_streetgraph(G_minimal_graph_input, export_path + 'given_lanes.gpkg', export_path + 'given_lanes_nodes.gpkg')

# =====================================================================================
# REBUILD AND EXPORT
# =====================================================================================

if 1:
    print('Link elimination')
    G_minimal_graph_output = snman.link_elimination(G_minimal_graph_input)

    print('Export minimal graph - output')
    snman.export_streetgraph(G_minimal_graph_output, export_path + 'minimal_graph_out_edges.gpkg', export_path + 'minimal_graph_out_nodes.gpkg')

    print('Rebuild lanes according to the OWTOP graph')
    snman.owtop.rebuild_lanes_from_owtop_graph(G, G_minimal_graph_output, hierarchies_to_protect={snman.hierarchy.HIGHWAY})

    print('Export rebuilt network with lanes')
    snman.export_streetgraph_with_lanes(G, 'ln_desc_after', export_path + 'edges_lanes_after.shp')


# =====================================================================================
# EXPORT
# =====================================================================================



#TODO-FAILURE: Remove inplace CRS conversion of G
if 0:
    print('Export OSM XML')
    snman.export_osm_xml(G, export_path + 'new_network.osm',{
        'highway', 'lanes', 'lanes:forward', 'lanes:backward', 'lanes:both_ways',
        'cycleway', 'cycleway:lane', 'cycleway:left', 'cycleway:left:lane', 'cycleway:right', 'cycleway:right:lane',
        'bus:lanes:backward', 'bus:lanes:forward', 'vehicle:lanes:backward', 'vehicle:lanes:forward',
        'maxspeed', 'oneway',
        '_connected_component'
    }, uv_tags=True, tag_all_nodes=True)

if 1:
    print('Export network without lanes')
    snman.export_streetgraph(G, export_path + 'edges.gpkg', export_path + 'nodes.gpkg',
        edge_columns=snman.constants.EXPORT_EDGE_COLUMNS
    )

if 1:
    print('Export network without lanes')
    snman.export_streetgraph(G, export_path + 'edges_all_attributes.gpkg', export_path + 'nodes_all_attributes.gpkg')

if 1:
    print('Export network with lanes')
    #TODO: Fix problems with saving lanes as a GeoPackage
    snman.export_streetgraph_with_lanes(G, 'ln_desc', export_path + 'edges_lanes.shp')



print('Done!')
