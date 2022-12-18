from .constants import *
from . import lanes, io, constants
from shapely.ops import substring
import shapely
import networkx as nx
import copy
import osmnx
import geopandas as gpd
import pandas as pd
import itertools
from . import osmnx_customized as oxc


def _edge_hash(edge):
    return '-'.join(map(str, edge[0:3]))


def _remove_edge_from_list(edges, edge_to_remove, dead_ends=True):
    edges_cleaned = []
    for idx, candidate in enumerate(edges):
        if (_edge_hash(candidate) != _edge_hash(edge_to_remove)
                and not (not dead_ends and candidate[3].get('dead_end'))):
            edges_cleaned.append(candidate)
    return edges_cleaned


def _get_neighbors(graph, edge, dead_ends=True):
    adjacent_nodes = edge[0:2]
    u_neighbors = list(graph.edges(nbunch=adjacent_nodes[0], data=True, keys=True))
    v_neighbors = list(graph.edges(nbunch=adjacent_nodes[1], data=True, keys=True))
    # Remove this edge from the neighbors
    u_neighbors = _remove_edge_from_list(u_neighbors, edge, dead_ends=dead_ends)
    v_neighbors = _remove_edge_from_list(v_neighbors, edge, dead_ends=dead_ends)
    return [
        _unique_edges(u_neighbors),
        _unique_edges(v_neighbors),
        _unique_edges(u_neighbors + v_neighbors)
    ]


def _unique_edges(edges):
    unique_edges = []
    hashes = set()
    for edge in edges:
        this_hash = _edge_hash(edge)
        this_hash in hashes or unique_edges.append(edge)
        hashes.add(this_hash)
    return unique_edges


def normalize_edge_directions(street_graph):
    edges = list(street_graph.edges(data=True, keys=True))
    for edge in edges:
        if edge[0] > edge[1]:
            _reverse_edge(street_graph, edge)
            pass


def _reverse_edge(street_graph, edge, reverse_topology=True):

    """
    Flip the edge direction

    Parameters
    ----------
    street_graph : nx.MultiDiGraph

    edge : Tuple

    reverse_topology : Boolean
        Also flip the start and end node. Automatically false if the street_graph is undirected

    """

    u = edge[0]
    v = edge[1]
    key = edge[2]
    data = edge[3]

    if not street_graph.has_edge(u, v, key):
        return edge

    # Don't reverse topology if the graph is not directed
    reverse_topology = reverse_topology and nx.is_directed(street_graph)

    # Remove the old edge
    if reverse_topology:
        street_graph.remove_edge(u, v, key)

    # Reverse the lanes
    if data.get('ln_desc') is not None:
        data['ln_desc'] = lanes._reverse_lanes(data['ln_desc'])

    # Reverse the geometry
    if data.get('geometry') and data.get('geometry') != shapely.ops.LineString():
        #print(data)
        #print(data['geometry'].wkt)
        data['geometry'] = substring(data['geometry'], 1, 0, normalized=True)

    # Flip the reversed flag
    data[constants.KEY_REVERSED] = not data.get(constants.KEY_REVERSED, False)

    # Add the new edge
    if reverse_topology:
        key = street_graph.add_edge(v, u, **data)
    else:
        nx.set_edge_attributes(street_graph, {(u, v, key): data})

    # Return the resulting edge
    if reverse_topology:
        return v, u, key, data
    else:
        return edge


def _split_edge(street_graph, u, v, key, split_point):

    # Don't continue if the edge does not exist
    if not street_graph.has_edge(u,v,key):
        #print('edge does not exist:', u, v, key)
        return False

    edge_data = street_graph.get_edge_data(u,v,key)

    # Assign a new node id
    split_node_id = max(street_graph.nodes) + 1

    # duplicate the existing edge data as a basis for the data of the two new edges
    edge1_data = copy.deepcopy(edge_data)
    edge2_data = copy.deepcopy(edge_data)

    # split geometry
    line = edge_data.get('geometry', False)
    if line != False:
        node_point = split_point
        # snap the split point onto the line
        split_point = shapely.ops.nearest_points(node_point, line)[1]
        # Just using the point to split the line does not work because it is not exactly on the line.
        # So we create a small circle around it
        split_circle = split_point.buffer(1)
        new_lines = shapely.ops.split(line, split_circle)

        # Stop if the splitting was not successful
        if len(new_lines.geoms) != 3:
            return False

        # Extend the geometries to the split node
        edge1_data['geometry'] = shapely.ops.LineString(
            list(new_lines.geoms[0].coords) + list(node_point.coords)
        )
        edge1_data['_split'] = 1
        edge2_data['geometry'] = shapely.ops.LineString(
            list(node_point.coords) + list(new_lines.geoms[2].coords)
        )
        edge2_data['_split'] = 2

    # Split topology
    street_graph.add_node(split_node_id, x=split_point.x, y=split_point.y, _split_node = True)
    street_graph.remove_edge(u, v, key)
    key1 = street_graph.add_edge(u, split_node_id, **edge1_data)
    key2 = street_graph.add_edge(split_node_id, v, **edge2_data)

    # Reverse the new edge data if the new node has changed their topological direction
    if not nx.is_directed(street_graph):
        if u > split_node_id:
            _reverse_edge(street_graph, (u, split_node_id, key1, edge1_data))
        if split_node_id > v:
            _reverse_edge(street_graph, (v, split_node_id, key2, edge2_data))


def split_through_edges_in_intersections(G, intersections):

    intersections = intersections.copy()
    intersections['ix_geometry'] = intersections['geometry']
    intersections['ix_centroid'] = intersections.centroid

    edges = oxc.graph_to_gdfs(G, nodes=False)
    edges['e_geometry'] = edges['geometry']
    #edges['u'] = edges.index.get_level_values('u')
    #edges['v'] = edges.index.get_level_values('v')
    #edges['key'] = edges.index.get_level_values('key')

    # create new column with endpoints of each edge,
    # we will need them to detect if an edge start/ends in an intersection buffer
    edges['e_endpoints'] = edges.apply(
        lambda x: shapely.ops.MultiPoint([
            x['geometry'].coords[0],
            x['geometry'].coords[-1]
        ])
        if isinstance(x.get('geometry'), shapely.ops.LineString) and not x.get('geometry').is_empty
        else None,
        axis=1
    )

    # build pairs of intersections and intersecting edges, keep index from edges
    a = gpd.sjoin(edges, intersections, how="inner", predicate="intersects", lsuffix='e', rsuffix='i')
    # stop here if there are no edge/intersection pairs
    if len(a) == 0:
        return
    # add a new column telling us whether the edge starts/ends within the intersection
    a['edge_endpoint_in_intersection'] = a.apply(
        lambda x: x['ix_geometry'].intersects(x['e_endpoints']),
        axis=1
    )
    # filter for only those intersection/edge pairs where the edge does not start/end in the intersection
    a = a[a['edge_endpoint_in_intersection'] == False]

    # find out where should the edge be split
    a['split_point'] = a.apply(
        lambda x: shapely.ops.nearest_points(x['ix_centroid'], x['e_geometry'])[1],
        axis=1
    )

    a.apply(
        lambda row: _split_edge(G, *row.name, row.split_point),
        axis=1
    )


def connect_components_in_intersections(street_graph, intersections, separate_layers=True):
    """
    Creates connections between weakly connected components so that they can be merged into a single intersection
    In this process, we use the following restrictions:
    - ignoring dead-end nodes (so that we don't connect adjacent dead-ends)
    - connecting only components that are on the same layer, e.g. on the same level of a multi-level intersection

    Parameters
    ----------
    street_graph : nx.MultiDiGraph

    intersections : gpd.GeoDataFrame
        Output from the function split_through_edges_in_intersections()
    """

    # get the nodes as a geodataframe
    node_points = osmnx.graph_to_gdfs(street_graph, edges=False)[["geometry", "street_count", "highway"]]
    # eliminate dead ends from the process to keep them disconnected
    node_points = node_points.query('street_count != 1')
    # assign each node to an intersection (polygon)
    gdf = gpd.sjoin(node_points, intersections, how="left", predicate="within")
    # clean up the columns of the resulting geodataframe (cluster=id of the intersection geometry)
    gdf = gdf.drop(columns="geometry").rename(columns={"index_right": "cluster"})

    groups = gdf.groupby("cluster")
    for cluster_label, nodes_subset in groups:
        if len(nodes_subset) > 1:
            # identify all (weakly connected) components in cluster
            # wccs is a list of sets containing id's of the nodes in each weakly connected component
            wccs = list(nx.weakly_connected_components(street_graph.subgraph(nodes_subset.index)))
            # skip this node if there are not at least two components (one component = no need for further connections)
            if len(wccs) <= 1:
                continue
            nodes = []
            layers = {}
            for wcc in wccs:
                # select a node that will be used to make the connector edge
                first_node = min(wcc)
                nodes.append(first_node)
                # find out the layers of each wcc
                this_wcc_layers = [a.get('layer',0) for u,v,k,a in street_graph.edges(wcc, keys=True, data=True)]
                if this_wcc_layers != []:
                    layer = max(set(this_wcc_layers), key = this_wcc_layers.count)
                else:
                    layer = None
                layers[first_node] = layer
            # iterate over all combinations of the selected nodes to create a complete from-to mesh
            for combination in itertools.combinations(nodes,2):
                u = min(combination)
                v = max(combination)
                layer_u = layers[u]
                layer_v = layers[v]
                # skip this connector if the layer sets don't match
                if separate_layers and layer_u != layer_v:
                    continue
                all_nodes = dict(street_graph.nodes(data=True))
                node_u = all_nodes.get(u)
                node_v = all_nodes.get(v)
                geom = ''
                if node_u and node_v:
                    geom = shapely.ops.LineString((
                        shapely.ops.Point(node_u.get('x'), node_u.get('y')),
                        shapely.ops.Point(node_v.get('x'), node_v.get('y'))
                    ))
                street_graph.add_edge(
                    u,v,geometry=geom,
                    _components_connector=True, layer=max(layers.items())[1], _layers=str(layers), osmid='0'
                )

def update_precalculated_attributes(street_graph):
    street_count = osmnx.stats.count_streets_per_node(street_graph)
    nx.set_node_attributes(street_graph, street_count, name="street_count")