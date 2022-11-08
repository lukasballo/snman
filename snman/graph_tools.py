from . import lanes
from .constants import *
from shapely.ops import substring
import networkx as nx


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
            edge[3][KEY_REVERSED] = True
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

    reverse_topology = reverse_topology and nx.is_directed(street_graph)

    # Remove the old edge
    if reverse_topology:
        street_graph.remove_edge(u, v, key)

    # Reverse the lanes
    data['ln_desc'] = lanes._reverse_lanes(data['ln_desc'])

    # Reverse the geometry
    if 'geometry' in data:
        data['geometry'] = substring(data['geometry'], 1, 0, normalized=True)

    # Add the new edge
    if reverse_topology:
        key = street_graph.add_edge(v, u, **data)

    # Return the resulting edge
    if reverse_topology:
        return v, u, key, data
    else:
        return edge
