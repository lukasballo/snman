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
crs = 'epsg:2056'      # CH1905+ projected CRS

ox.utils.config(useful_tags_way=osm_tags)

# TODO: Include highway=living_street and service roads
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
    'Zurich, Zurich, Switzerland',
    #'Seebach, Zurich, Switzerland',
    #'Altstetten, Zurich, Switzerland',
    #'Altstadt, Zurich, Switzerland',
    #'Oberstrass, Zurich, Switzerland',
    custom_filter=custom_filter,
    simplify=True,
    simplify_strict=False,
    retain_all=True,
    one_edge_per_direction=False,
)

print('Convert CRS of street graph to 2056')
snman.convert_crs_of_street_graph(street_graph, crs)
nodes = copy.copy(street_graph.nodes)

print('Normalize edge directions, enforce direction from lower to higher node id')
snman.normalize_edge_directions(street_graph)

for i in range(3):
    intersections_geom = snman.osmnx.simplification._merge_nodes_geometric(
        street_graph,
        INTERSECTION_TOLERANCE
        )
    edges = street_graph.edges(data=True, keys=True)
    edges_geom = gpd.GeoSeries(
        [edge[3].get('geometry') for edge in edges],
        crs=crs
        )

    intersections = gpd.GeoDataFrame({'geometry': intersections_geom},geometry='geometry')
    intersections['intersection_geometry'] = intersections_geom
    intersections['intersection_centroid'] = intersections.centroid

    edges = gpd.GeoDataFrame({'nx': edges, 'geometry': edges_geom}, geometry='geometry')
    edges['edge_geometry'] = edges['geometry']
    # create new column with endpoints of each edge,
    # we will need them to detect if an edge start/ends in an intersection buffer
    edges['edge_endpoints'] = edges.apply(
        lambda x: shp.ops.MultiPoint([
            x['geometry'].coords[0],
            x['geometry'].coords[-1]
            ])
            if isinstance(x.get('geometry'), shp.ops.LineString) and not x.get('geometry').is_empty
            else None,
        axis=1
        )

    a = gpd.sjoin(intersections, edges, how="inner", predicate="intersects", lsuffix='i', rsuffix='e')
    a['edge_endpoint_is_in_intersection'] = a.apply(
        lambda x: x['intersection_geometry'].intersects(x['edge_endpoints']),
        axis=1
    )
    a = a[a['edge_endpoint_is_in_intersection']==False]
    a['split_point'] = a.apply(
        lambda x: shp.ops.nearest_points(x['intersection_centroid'], x['edge_geometry'])[1],
        axis=1
        )
    a.apply(
        lambda x: snman.graph_tools._split_edge(street_graph, x['nx'], x['split_point']),
        axis=1
        )


street_count = snman.osmnx.stats.count_streets_per_node(street_graph)
nx.set_node_attributes(street_graph, street_count, name="street_count")


# Save the intersection polygons into a file
snman.export_gdf(intersections, export_path + 'intersections.gpkg', columns=['geometry'])




if 1:
    print('Consolidate intersections')
    street_graph = ox.simplification.consolidate_intersections(
        street_graph, tolerance=INTERSECTION_TOLERANCE, rebuild_graph=True, dead_ends=True, reconnect_edges=True
    )

print('Generate lanes')
snman.generate_lanes(street_graph)

print('Normalize edge directions, enforce direction from lower to higher node id')
snman.normalize_edge_directions(street_graph)

snman.export_streetgraph(street_graph, export_path + 'edges_early.gpkg', export_path + 'nodes_early.gpkg')


print('Convert into an undirected graph')
street_graph = ox.utils_graph.get_undirected(street_graph)




print('Identify hierarchy')
snman.add_hierarchy(street_graph)

snman.export_streetgraph(street_graph, export_path + 'raw_edges.gpkg', export_path + 'raw_nodes.gpkg')

if 1:
    print('Merge parallel and consecutive edges, repeat a few times')
    for i in range(5):
        snman.merge_parallel_edges(street_graph)
        snman.merge_consecutive_edges(street_graph)
        pass


print('Add lane stats')
snman.generate_lane_stats(street_graph)

"""
print('Fixing one-sided intersections')

for node_id in list(street_graph.nodes):
    node = [node_id, street_graph.nodes(data=True)[node_id]]
    #edge = [34,45,0,street_graph.get_edge_data(34,45,0)]
    edges = street_graph.edges(data=True, keys=True)
    for edge in copy.deepcopy(edges):
        # Exclude all edges from/to the examined node
        if node[0] in edge[0:1]:
            continue
        node_geom = shp.ops.Point(node[1].get('x'), node[1].get('y'))
        edge_geom = edge[3].get('geometry', shp.ops.LineString())
        distance = node_geom.distance(edge_geom)
        if distance < 25:
            snman.graph_tools._split_edge(street_graph, edge, node)

for i in range(5):
    snman.merge_parallel_edges(street_graph)
    snman.merge_consecutive_edges(street_graph)
    pass

"""


if 0:
    print('Add public transport')
    pt_network = snman.import_shp_to_gdf("C:/DATA/CLOUD STORAGE/polybox/Research/SNMan/SNMan Shared/stadt_zuerich_open_data/Linien_des_offentlichen_Verkehrs_-OGD/ZVV_LINIEN_GEN_L.shp")
    snman.match_pt(street_graph, pt_network)

print('Add lane stats')
snman.generate_lane_stats(street_graph)

print('Update OSM tags')
snman.update_osm_tags(street_graph)

print('Set given lanes')
snman.set_given_lanes(street_graph)

print('Create directed graph of given lanes')
given_lanes_graph = snman.create_given_lanes_graph(street_graph)

print('Export network without lanes')
snman.export_streetgraph(street_graph, export_path + 'edges.gpkg', export_path + 'nodes.gpkg')

print('Export network with lanes')
#TODO: Fix problems with saving lanes as a GeoPackage
snman.export_streetgraph_with_lanes(street_graph, 'ln_desc', export_path + 'edges_lanes.shp')

#print('Export network with given lanes')
#snman.export_streetgraph_with_lanes(street_graph, 'given_lanes', export_path + 'edges_given_lanes.shp')

#print('Export given lanes')
#snman.export_streetgraph(given_lanes_graph, export_path + 'given_lanes.gpkg', export_path + 'given_lanes_nodes.gpkg')

print('Export OSM XML')
snman.export_osm_xml(street_graph, export_path + 'new_network.osm',{
    'lanes', 'lanes:forward', 'lanes:backward', 'lanes:both_ways',
    'cycleway', 'cycleway:lane', 'cycleway:left', 'cycleway:left:lane', 'cycleway:right', 'cycleway:right:lane',
    'maxspeed'
})

print('Done!')
