from . import lanes
from .constants import *
from shapely.ops import substring
import shapely
import networkx as nx
import copy


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
    if data.get('ln_desc'):
        data['ln_desc'] = lanes._reverse_lanes(data['ln_desc'])

    # Reverse the geometry
    if data.get('geometry') and data.get('geometry') != shapely.ops.LineString():
        data['geometry'] = substring(data['geometry'], 1, 0, normalized=True)

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


def _split_edge(street_graph, edge, split_point):

    u = edge[0]
    v = edge[1]
    key = edge[2]
    edge_data = edge[3]

    # Don't continue if the edge does not exist
    if not street_graph.has_edge(u,v,key):
        return False

    split_node_id = len(street_graph.nodes)

    edge1_data = copy.copy(edge_data)
    edge2_data = copy.copy(edge_data)

    # Split geometry
    line = edge_data.get('geometry', False)
    if line != False:
        node_point = split_point
        split_point = shapely.ops.nearest_points(node_point, line)[1]
        # Just using the point to split the line does not work because it is not exactly on the line.
        # So we create a small circle around it
        split_circle = split_point.buffer(1)
        new_lines = shapely.ops.split(line, split_circle)

        # Stop if the splitting was not successful
        if len(new_lines) != 3:
            return False

        # Extend the geometries to the split node
        edge1_data['geometry'] = shapely.ops.LineString(
            list(new_lines[0].coords) + list(node_point.coords)
        )
        edge2_data['geometry'] = shapely.ops.LineString(
            list(node_point.coords) + list(new_lines[2].coords)
        )

    # Split topology
    street_graph.add_node(split_node_id, x=split_point.x, y=split_point.y)
    street_graph.remove_edge(u, v, key)
    key1 = street_graph.add_edge(u, split_node_id, **edge1_data)
    key2 = street_graph.add_edge(split_node_id, v, **edge2_data)

    # Reverse the new edge data if the new node has changed their topological direction
    if u > split_node_id:
        _reverse_edge(street_graph, (u, split_node_id, key1, edge1_data))
    if split_node_id > v:
        _reverse_edge(street_graph, (v, split_node_id, key2, edge2_data))


"""
def _split_edge(street_graph, edge, split_node):

    u = edge[0]
    v = edge[1]
    key = edge[2]
    edge_data = edge[3]
    split_node_data = split_node[1]

    edge1_data = copy.copy(edge_data)
    edge2_data = copy.copy(edge_data)

    # Split geometry
    line = edge_data.get('geometry', False)
    if line != False:
        node_point = shapely.ops.Point(split_node_data.get('x'), split_node_data.get('y'))
        split_point = shapely.ops.nearest_points(node_point, line)[1]
        # Just using the point to split the line does not work because it is not exactly on the line.
        # So we create a small circle around it
        split_circle = split_point.buffer(1)
        new_lines = shapely.ops.split(line, split_circle)

        # Stop if the splitting was not successful
        if len(new_lines) != 3:
            return False

        # Extend the geometries to the split node
        edge1_data['geometry'] = shapely.ops.LineString(
            list(new_lines[0].coords) + list(node_point.coords)
        )
        edge2_data['geometry'] = shapely.ops.LineString(
            list(node_point.coords) + list(new_lines[2].coords)
        )

    # Split topology
    street_graph.remove_edge(u, v, key)
    key1 = street_graph.add_edge(u, split_node[0], **edge1_data)
    key2 = street_graph.add_edge(split_node[0], v, **edge2_data)

    # Reverse the new edge data if the new node has changed their topological direction
    if u > split_node[0]:
        _reverse_edge(street_graph, [u, split_node[0], key1, edge1_data], reverse_topology=False)
    if split_node[0] > v:
        _reverse_edge(street_graph, [v, split_node[0], key2, edge2_data], reverse_topology=False)
"""

