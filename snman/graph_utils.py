import networkx as nx


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
