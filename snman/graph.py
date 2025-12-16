import networkx as nx
import numpy as np
from . import osmnx_customized as oxc
from . import utils


def match_points_to_nodes(points_gdf, G):
    """
    Match points in a GeoDataFrame to nearest nodes in the graph.

    Parameters
    ----------
    points_gdf : gpd.GeoDataFrame
        GeoDataFrame with point geometries
    G : nx.MultiDiGraph
        Street graph

    Returns
    -------
    gpd.GeoDataFrame
        GeoDataFrame with added 'node' and 'point_to_node_dist' columns
    """
    nodes = oxc.nearest_nodes(
        G,
        list(map(lambda geom: geom.x, points_gdf.geometry)),
        list(map(lambda geom: geom.y, points_gdf.geometry)),
        return_dist=True
    )
    points_gdf['node'] = nodes[0]
    points_gdf['point_to_node_dist'] = nodes[1]
    return points_gdf


def match_point_to_node(point, G):
    """
    Match a single point to the nearest node in the graph.

    Parameters
    ----------
    point : shapely.Point
        Point geometry
    G : nx.MultiDiGraph
        Street graph

    Returns
    -------
    int
        Nearest node ID
    """
    nodes = oxc.nearest_nodes(
        G,
        [point.x],
        [point.y],
        return_dist=False
    )
    return nodes[0]


def weak_neighbors(G, node):
    """
    Returns neighbors of a node considering both incoming and outcoming edges

    Parameters
    ----------
    G : nx.DiGraph or nx.MultiDiGraph
    node : int

    Returns
    -------
    set
    """

    edges = list(G.in_edges(node)) + list(G.out_edges(node))
    return set(utils.flatten_list(edges)).difference({node})


def cost_increase_by_edge_removal(G, u, v, k, weight):
    """
    Returns the increase of cost between u and v if the (u,v,k) edge is removed.

    Parameters
    ----------
    G : nx.MultiDiGraph
        graph, this function works only for a MultiDiGraph
    u : int
    v : int
    k : int
    weight : str
        which attribute should be used as weight

    Returns
    -------
    float
        Cost increase (or np.nan if edge doesn't exist, or np.inf if path becomes disconnected)
    """

    H = G.edge_subgraph(
        filter(
            lambda candidate_uvk:
                not set(candidate_uvk[0:2]).isdisjoint({u, v}),
            G.edges
        )
    )

    cost_before_removal = nx.shortest_path_length(
        H,
        u, v, weight=weight
    )

    if not H.has_edge(u, v, k):
        return np.nan

    try:
        cost_after_removal = nx.shortest_path_length(
            H.edge_subgraph(set(G.edges).difference({(u, v, k)})),
            u, v, weight=weight
        )
    except (nx.NetworkXNoPath, nx.NodeNotFound):
        cost_after_removal = np.inf

    return cost_after_removal - cost_before_removal


def keep_only_the_largest_connected_component(G, weak=False):
    """
    Remove all nodes and edges that are disconnected from the largest connected component.
    For directed graphs, strong connectedness will be considered, unless weak=True

    Parameters
    ----------
    G : nx.MultiGraph
        street graph
    weak : bool
        use weakly connected component in case of a directed graph

    Returns
    -------
    H : copy of subgraph representing the largest connected component
    """

    if len(G.nodes) == 0:
        return G.subgraph([])

    if G.is_directed():
        if weak:
            nodes = max(nx.weakly_connected_components(G), key=len)
        else:
            nodes = max(nx.strongly_connected_components(G), key=len)
    else:
        nodes = max(nx.connected_components(G), key=len)

    return G.subgraph(nodes).copy()


def safe_degree(G, node):
    """
    For directed graphs, returns the sum of in_degree and out_degree.
    For undirected graphs, return the degree.

    Parameters
    ----------
    G : nx.Graph, nx.DiGraph, nx.MultiGraph, nx.MultiDiGraph
    node : int

    Returns
    -------
    int
    """

    if G.is_directed():
        return G.degree(node)
    else:
        return G.in_degree(node) + G.out_degree(node)


def safe_remove_edge(G, u, v, k, remove_dangling_nodes=True):
    """
    Removes an edge from the graph:
    - without throwing error if the edge does not exist
    - optional: removing any dangling nodes

    Parameters
    ----------
    G : nx.Graph, nx.DiGraph, nx.MultiGraph, nx.MultiDiGraph
    u : int
    v : int
    k : int

    Returns
    -------
    None
    """

    if not G.has_edge(u, v, k):
        return

    G.remove_edge(u, v, k)

    if remove_dangling_nodes:
        for node in [u, v]:
            if safe_degree(G, node) == 0:
                G.remove_node(node)


def plot_scc(G):
    """
    Plot the graph with nodes colored by strongly connected component.

    Parameters
    ----------
    G : nx.MultiDiGraph
        Street graph

    Returns
    -------
    None
    """
    # Calculate strongly connected components
    components = list(nx.strongly_connected_components(G))

    # Assign a different color to each component
    colors = {}
    for i, component in enumerate(components):
        for node in component:
            colors[node] = i

    # Draw the graph with nodes colored by component
    oxc.plot_graph(G, node_color=[colors[node] for node in G.nodes()])


def remove_isolated_nodes(G):
    """
    Removes all nodes that are not connected to any other nodes

    Parameters
    ----------
    G: nx.Graph

    Returns
    -------
    None
    """

    G.remove_nodes_from(list(nx.isolates(G)))


def apply_function_to_each_edge(G, function):
    """
    Applies the given function to each edge

    Parameters
    ----------
    G: nx.Graph
    function: function

    Returns
    -------
    None
    """

    for uvk, data in G.edges.items():
        function(G, uvk)


def apply_function_to_each_node(G, function):
    """
    Applies the given function to each node

    Parameters
    ----------
    G: nx.Graph
    function: function

    Returns
    -------
    None
    """

    for n, data in G.nodes.items():
        function(G, n)


class SNManGenericGraph:

    def apply_to_edges(self, fn):
        """
        Apply a function to all edges.

        Parameters
        ----------
        fn : function
            Function with signature fn(u, v, k, data) to apply to each edge

        Returns
        -------
        None
        """
        pass


class SNManMultiDiGraph(nx.MultiDiGraph, SNManGenericGraph):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
