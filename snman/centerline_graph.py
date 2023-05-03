from . import osmnx_customized as oxc
from . import any_graph, lane_config
from .constants import *
import geopandas as gpd
import networkx as nx
import shapely
import pyproj


def surrogate_missing_edge_geometries(G):
    for uvk, data in G.edges.items():
        if 'geometry' not in data:
            data['geometry'] = shapely.LineString([
                shapely.Point(get_node_point(G, uvk[0])),
                shapely.Point(get_node_point(G, uvk[1]))
            ])


def get_node_point(G, node):
    data = G.nodes[node]
    return shapely.Point(data['x'], data['y'])


def add_connected_component_ids(Gc):
    """
    For directed graphs: Adds IDs of weakly ('_weakly_connected_component')
    and strongly ('_strongly_connected_component') connected components to all edges
    For undirected graphs: Adds IDs of connected components ('_connected_component') to all edges

    Parameters
    ----------
    Gc : nx.MultiGraph
    Returns
    -------
    None
    """

    if Gc.is_directed():

        wccs = nx.weakly_connected_components(Gc)
        for wcc, nodes in enumerate(wccs):
            subgraph = nx.subgraph(Gc, nodes)
            nx.set_edge_attributes(subgraph, wcc, '_weakly_connected_component')

        sccs = nx.strongly_connected_components(Gc)
        for scc, nodes in enumerate(sccs):
            subgraph = nx.subgraph(Gc, nodes)
            nx.set_edge_attributes(subgraph, wcc, '_strongly_connected_component')

    else:

        ccs = nx.connected_components(Gc)
        for wcc, nodes in enumerate(ccs):
            subgraph = nx.subgraph(Gc, nodes)
            nx.set_edge_attributes(subgraph, wcc, '_connected_component')


def keep_only_the_largest_connected_component(Gc, weak=False):
    """
    Remove all nodes and edges that are disconnected from the largest connected component.
    For directed graphs, strong connectedness will be considered, unless weak=True

    Parameters
    ----------
    Gc : nx.MultiGraph
        street graph
    weak : bool
        use weakly connected component in case of a directed graph

    Returns
    -------
    H : copy of subgraph representing the largest connected component
    """

    if Gc.is_directed():
        if weak:
            nodes = max(nx.weakly_connected_components(Gc), key=len)
        else:
            nodes = max(nx.strongly_connected_components(Gc), key=len)
    else:
        nodes = max(nx.connected_components(Gc), key=len)

    return Gc.subgraph(nodes).copy()


def update_precalculated_attributes(G):
    """
    Update edge attributes in a street graph that have been pre-calculated by osmnx
    when the graph was created

    Parameters
    ----------
    G : nx.MultiGraph or nx.MultiDiGraph
        street graph

    Returns
    -------
    None
    """
    street_count = oxc.stats.count_streets_per_node(G)
    nx.set_node_attributes(G, street_count, name="street_count")


def convert_crs_of_centerline_graph(Gc, to_crs):
    """
    Convert the coordinate reference system of the geometries in a graph.
    The graph must follow the convention of networkx:
        * edges have a 'geometry' attribute with a shapely geometry
        * nodes have 'x' and 'y' attributes holding the coordinates as floats
        * the graph has a 'crs' attribute

    Parameters
    ----------
    Gc : nx.MultiGraph
        centerline graph
    to_crs : int
        target crs

    Returns
    -------
    None
    """

    # Initialize the CRS transformer
    from_crs = pyproj.CRS(Gc.graph['crs'])
    to_crs = pyproj.CRS(to_crs)
    project = pyproj.Transformer.from_crs(from_crs, to_crs, always_xy=True).transform

    # Update the street_graph's metadata
    Gc.graph["crs"] = to_crs

    # Transform the geometry of all edges
    for edge in Gc.edges(data=True, keys=True):
        if "geometry" in edge[3]:
            edge[3]["geometry"] = shapely.ops.transform(project, edge[3]["geometry"])

    # Transform the geometry of all nodes
    for id, data in Gc.nodes.items():
        geom = shapely.geometry.Point(data.get('x'), data.get('y'))
        geom = shapely.ops.transform(project, geom)
        data['x'] = geom.x
        data['y'] = geom.y


def centerline_graph_to_lane_graph(Gc, mode, lanes_attribute=KEY_LANES_DESCRIPTION):

    # initialize and copy graph attributes
    Gl = nx.MultiDiGraph()
    Gl.graph = Gc.graph

    for uvk, data in Gc.edges.items():
        u, v, k = uvk

        lanes_list = data.get(lanes_attribute)
        for lane in lanes_list:
            lp = lanes._lane_properties(lane)
            length = data['length']
            cost = _calculate_lane_cost(lane, length, mode)
            only_active_modes = lp.modes.issubset(ACTIVE_MODES)
            only_active_modes_length = length * only_active_modes

            common_attributes = {
                'length': length,
                'cost': cost,
                'only_active_modes': only_active_modes,
                'only_active_modes_length': only_active_modes_length
            }

            if mode in lp.modes:
                if lp.direction in [DIRECTION_FORWARD, DIRECTION_BOTH]:
                    Gl.add_edge(u, v, lane=lane, **common_attributes)
                if lp.direction in [DIRECTION_BACKWARD, DIRECTION_BOTH]:
                    Gl.add_edge(v, u, lane=lanes.reverse_lane(lane), **common_attributes)

    # take over the node attributes from the street graph
    nx.set_node_attributes(Gl, dict(Gc.nodes))

    return Gl


def _remove_edge_from_list(edges, edge_to_remove, dead_ends=True):
    """
    Remove edge from a list using its unique string id

    Parameters
    ----------
    edges : list
    edge_to_remove : list
    dead_ends : bool
        include dead ends in the resulting list

    Returns
    -------
    list
    """

    edges_cleaned = []
    for idx, candidate in enumerate(edges):
        if (candidate[0:3] != edge_to_remove[0:3]
                and not (not dead_ends and candidate[3].get('dead_end'))):
            edges_cleaned.append(candidate)
    return edges_cleaned


def _get_neighbors(G, edge, dead_ends=True):
    """
    Return the neighbor edges of a given edge. The result is a list of three lists:
        * neighbors at the side of u node
        * neighbors at the side of v node
        * neighbors from both sides

    Parameters
    ----------
    G : nx.MultiGraph
        street graph
    edge : tuple
        the complete tuple of an edge
    dead_ends : bool
        include dead ends int the results

    Returns
    -------
    list
        a list of edges
    """

    adjacent_nodes = edge[0:2]
    u_neighbors = list(G.edges(nbunch=adjacent_nodes[0], data=True, keys=True))
    v_neighbors = list(G.edges(nbunch=adjacent_nodes[1], data=True, keys=True))
    # Remove this edge from the neighbors
    u_neighbors = _remove_edge_from_list(u_neighbors, edge, dead_ends=dead_ends)
    v_neighbors = _remove_edge_from_list(v_neighbors, edge, dead_ends=dead_ends)
    return [
        _unique_edges(u_neighbors),
        _unique_edges(v_neighbors),
        _unique_edges(u_neighbors + v_neighbors)
    ]


def _unique_edges(edges):
    """
    Remove duplicates from a list of edges

    Parameters
    ----------
    edges : list

    Returns
    -------
    list
    """

    unique_edges = []
    ids = set()
    for edge in edges:
        this_id = edge[0:3]
        this_id in ids or unique_edges.append(edge)
        ids.add(this_id)
    return unique_edges
