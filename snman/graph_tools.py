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


def create_OSM_graph(G):
    pass


def create_directed_graph(G):
    pass


def add_connected_component_ids(G):
    """
    For directed graphs: Adds IDs of weakly and strongly connected components to all edges
    For undirected graphs: Adds IDs of connected components to all edges

    Parameters
    ----------
    G : nx.MultiDiGraph
    """


    if G.is_directed():

        wccs = nx.weakly_connected_components(G)
        for wcc, nodes in enumerate(wccs):
            subgraph = nx.subgraph(G, nodes)
            nx.set_edge_attributes(subgraph, wcc, '_weakly_connected_component')

        sccs = nx.strongly_connected_components(G)
        for scc, nodes in enumerate(sccs):
            subgraph = nx.subgraph(G, nodes)
            nx.set_edge_attributes(subgraph, wcc, '_strongly_connected_component')

    else:

        ccs = nx.connected_components(G)
        for wcc, nodes in enumerate(ccs):
            subgraph = nx.subgraph(G, nodes)
            nx.set_edge_attributes(subgraph, wcc, '_connected_component')


def keep_only_the_largest_connected_component(G, weak:bool=False) -> nx.MultiGraph:
    """
    Remove all nodes and edges that are disconnected from the largest connected component.
    For directed graphs, strong connectedness will be considered, unless weak=True
        :param weak: use weakly connected component in case of a directed graph
        :returns: copy of subgraph representing the largest connected component
    """

    if G.is_directed():
        if weak:
            nodes = max(nx.weakly_connected_components(G), key=len)
        else:
            nodes = max(nx.strongly_connected_components(G), key=len)
    else:
        nodes = max(nx.connected_components(G), key=len)

    return G.subgraph(nodes).copy()


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


def normalize_edge_directions(G):
    edges = list(G.edges(data=True, keys=True))
    for edge in edges:
        if edge[0] > edge[1]:
            _reverse_edge(G, edge)
            pass


def _reverse_edge(G, edge, reverse_topology=True):

    """
    Flip the edge direction

    Parameters
    ----------
    G : nx.MultiDiGraph

    edge : Tuple

    reverse_topology : Boolean
        Also flip the start and end node. Automatically false if the G is undirected

    """

    u = edge[0]
    v = edge[1]
    key = edge[2]
    data = edge[3]

    if not G.has_edge(u, v, key):
        return edge

    # Don't reverse topology if the graph is not directed
    reverse_topology = reverse_topology and nx.is_directed(G)

    # Remove the old edge
    if reverse_topology:
        G.remove_edge(u, v, key)

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
        key = G.add_edge(v, u, **data)
    else:
        nx.set_edge_attributes(G, {(u, v, key): data})

    # Return the resulting edge
    if reverse_topology:
        return v, u, key, data
    else:
        return edge


def _split_edge(G, u, v, key, split_point):

    # Don't continue if the edge does not exist
    if not G.has_edge(u,v,key):
        #print('edge does not exist:', u, v, key)
        return False

    edge_data = G.get_edge_data(u,v,key)

    # Assign a new node id
    split_node_id = max(G.nodes) + 1

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
    G.add_node(split_node_id, x=split_point.x, y=split_point.y, _split_node = True)
    G.remove_edge(u, v, key)
    key1 = G.add_edge(u, split_node_id, **edge1_data)
    key2 = G.add_edge(split_node_id, v, **edge2_data)

    # Reverse the new edge data if the new node has changed their topological direction
    if not nx.is_directed(G):
        if u > split_node_id:
            _reverse_edge(G, (u, split_node_id, key1, edge1_data))
        if split_node_id > v:
            _reverse_edge(G, (v, split_node_id, key2, edge2_data))


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

    # stop here if there are no instances to process
    if len(a) == 0:
        return

    # find out where should the edge be split
    a['split_point'] = a.apply(
        lambda x: shapely.ops.nearest_points(x['ix_centroid'], x['e_geometry'])[1],
        axis=1
    )

    a.apply(
        lambda row: _split_edge(G, *row.name, row.split_point),
        axis=1
    )


def _is_motorized(edge):
    return edge.get('highway') not in ACTIVE_HIGHWAY_VALUES


def _add_layers_to_nodes(G):
    node_layers = {id: {'motorized': set(), 'active': set()} for id, node in G.nodes.items()}
    for edge_id, edge in G.edges.items():
        node_ids = [edge_id[0], edge_id[1]]
        mode = 'motorized' if _is_motorized(edge) else 'active'
        for node_id in node_ids:
            node_layers[node_id][mode].add(edge.get('layer',0))
    nx.set_node_attributes(G, node_layers, 'layers')


def _are_node_layers_compatible(layers1, layers2):

    if None in [layers1, layers2]:
        return False

    # filter out empty sets (i.e. modes with no layers)
    nodes = [
        list(
            filter(
                lambda x: x != set(),
                [node['motorized'], node['active']]
            )
        ) for node in [layers1, layers2]]

    # stop here if one of the nodes has no layers for any mode
    if 0 in [len(node) for node in nodes]:
        return False

    # keep only the layers of the highest mode of each node
    nodes = [node[0] for node in nodes]

    # check if the layers for the highest mode of each node overlap
    return not nodes[0].isdisjoint(nodes[1])


def connect_components_in_intersections(G, intersections, separate_layers=True):
    """
    Creates connections between weakly connected components so that they can be merged into a single intersection
    In this process, we use the following restrictions:
    - ignoring dead-end nodes (so that we don't connect adjacent dead-ends)
    - connecting only components that are on the same layer, e.g. on the same level of a multi-level intersection

    Parameters
    ----------
    G : nx.MultiDiGraph

    intersections : gpd.GeoDataFrame
        Output from the function split_through_edges_in_intersections()
    """

    # get the nodes as a geodataframe
    node_points = osmnx.graph_to_gdfs(G, edges=False)[["geometry", "street_count", "highway"]]
    # eliminate dead ends from the process to keep them as they are
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
            wccs = list(nx.weakly_connected_components(G.subgraph(nodes_subset.index)))
            # skip this node if there are not at least two components (one component = no need for further connections)
            if len(wccs) <= 1:
                continue
            nodes = []
            layers = {}
            # iterate over all combinations of the selected nodes (one from each wcc) to create a complete from-to mesh
            for combination in itertools.combinations(range(len(wccs)),2):

                a = combination[0]
                b = combination[1]

                # create gdfs of nodes from each wcc
                a = oxc.graph_to_gdfs(G.subgraph(wccs[a]), edges=False)
                b = oxc.graph_to_gdfs(G.subgraph(wccs[b]), edges=False)

                # join with nearest points
                joined = gpd.sjoin_nearest(a, b, distance_col="distance")
                joined['node_a'] = joined.index
                joined = joined.rename(columns={"index_right": "node_b"})

                # get the closest pair
                closest_pair = joined.sort_values(by='distance', ascending=True).to_dict('records')[0]
                a = closest_pair['node_a']
                b = closest_pair['node_b']

                a_data = G.nodes[a]
                b_data = G.nodes[b]

                # skip this connector if the layer sets don't match
                if (separate_layers
                    and not _are_node_layers_compatible(a_data.get('layers'), b_data.get('layers'))
                ):
                    continue

                geom = shapely.ops.LineString((
                    shapely.ops.Point(a_data.get('x'), a_data.get('y')),
                    shapely.ops.Point(b_data.get('x'), b_data.get('y'))
                ))

                G.add_edge(
                    a, b, geometry=geom,
                    osmid=0,
                    _components_connector=True
                )

def update_precalculated_attributes(G):
    street_count = osmnx.stats.count_streets_per_node(G)
    nx.set_node_attributes(G, street_count, name="street_count")