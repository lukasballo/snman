from . import lanes
from shapely.ops import substring

def _edge_hash(edge):
    return '-'.join(map(str, edge[0:3]))


def _remove_edge_from_list(edges, edge_to_remove, dead_ends=True):
    edges_cleaned = []
    for idx, candidate in enumerate(edges):
        if (_edge_hash(candidate) != _edge_hash(edge_to_remove)
                and not (dead_ends == False and candidate[3]['dead_end'] == True)):
            edges_cleaned.append(candidate)
    return edges_cleaned


def _get_neighbors(graph, edge, dead_ends=True):
    adjacent_nodes = edge[0:2]
    u_neighbors =\
        list(graph.in_edges(nbunch=adjacent_nodes[0], data=True, keys=True))\
        + list(graph.out_edges(nbunch=adjacent_nodes[0], data=True, keys=True))
    v_neighbors =\
        list(graph.in_edges(nbunch=adjacent_nodes[1], data=True, keys=True))\
        + list(graph.out_edges(nbunch=adjacent_nodes[1], data=True, keys=True))
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


def _reverse_edge(street_graph, edge):

    u = edge[0]
    v = edge[1]
    key = edge[2]
    data = edge[3]

    # Remove the old edge
    street_graph.remove_edge(u, v, key)

    # Reverse the lanes
    data['ln_desc'] = lanes._reverse_lanes(data['ln_desc'])

    # Reverse the geometry
    data['geometry'] = substring(data['geometry'], 1, 0, normalized=True)

    # Add the new edge
    key = street_graph.add_edge(v, u, **data)

    # Return a tuple with a proper edge format
    return v, u, key, data
