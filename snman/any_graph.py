import copy
import networkx as nx
import shapely as shp
from .constants import *
from . import lane_config
import geopandas as gpd


def split_edge(G, u, v, key, split_points):

    if not G.has_edge(u, v, key):
        return

    edge_data = G.get_edge_data(u, v, key)
    edge_linestring = edge_data.get('geometry', False)

    if edge_linestring is False:
        return

    # snap split points to the edge linestring
    split_points = shp.ops.nearest_points(edge_linestring, split_points)[0]
    # make a small buffer around the points to deal with numeric errors
    split_circles = shp.MultiPoint(split_points).buffer(0.1)
    # split the edge linestring using the buffer circles
    segments = shp.ops.split(edge_linestring, split_circles)

    # generate a list of node points
    node_points = list(map(lambda segment: segment.centroid, list(segments.geoms)[1::2]))
    node_points = [shp.Point(edge_linestring.coords[0])] + node_points + [shp.Point(edge_linestring.coords[-1])]

    # generate a list of new node ids
    start_node_id = max(G.nodes) + 1
    node_ids = list(range(start_node_id, start_node_id + len(node_points)-2))
    node_ids = [u] + node_ids + [v]

    # generate a list of edge linestrings
    edge_linestrings = list(segments.geoms)[0::2]

    # zip
    nodes = list(zip(node_ids, node_points))

    if not len(nodes) == len(split_points) + 2 == len(edge_linestrings) + 1:
        print('not len(nodes) == len(split_points) + 2 == len(edge_linestrings) + 1')
        print(nodes)
        print(split_points)
        print(edge_linestrings)
        return

    new_edges = []

    for i, node in enumerate(nodes):
        # add node into the graph
        if i != 0 and i != len(node_points)-1:
            G.add_node(node[0], x=node[1].x, y=node[1].y, _split_node=True)
        # add edge into the graph
        if i != 0:
            previous_node = nodes[i-1]
            new_u = previous_node[0]
            new_v = node[0]
            new_edge_data = copy.deepcopy(edge_data)
            # extend the geometry to the exact node points
            new_edge_data['geometry'] = shp.LineString(
                [previous_node[1]] +
                list(edge_linestrings[i-1].coords) +
                [node[1]]
            )
            new_key = G.add_edge(new_u, new_v, **new_edge_data)
            new_edges.append((new_u, new_v, new_key, new_edge_data))

            # in undirected graphs, the edge topology might have reversed implicitly
            # in such cases, we need to reverse everything else in the edge as well
            if not nx.is_directed(G) and new_u > new_v:
                reverse_edge(G, new_u, new_v, new_key, reverse_topology=False)

    # remove the original edge
    G.remove_edge(u, v, key)

    return new_edges


def split_edge_old(G, u, v, key, split_point):
    """
    Split the edge at a given point and create two child edges, together with a new node in between

    Parameters
    ----------
    G : nx.MultiGraph
        street graph
    u : int
        edge to be split - u node
    v : int
        edge to be split - v node
    key : int
        edge to be split - key
    split_point : shapely.geometry.Point

    Returns
    -------
    None
    """

    # Don't continue if the edge does not exist
    if not G.has_edge(u,v,key):
        #print('edge does not exist:', u, v, key)
        return False

    edge_data = G.get_edge_data(u, v, key)

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
        split_point = shp.ops.nearest_points(node_point, line)[1]
        # Just using the point to split the line does not work because it is not exactly on the line.
        # So we create a small circle around it
        split_circle = split_point.buffer(1)
        new_lines = shp.ops.split(line, split_circle)

        # Stop if the splitting was not successful
        if len(new_lines.geoms) != 3:
            return False

        # Extend the geometries to the split node
        edge1_data['geometry'] = shp.ops.LineString(
            list(new_lines.geoms[0].coords) + list(node_point.coords)
        )
        edge1_data['_split'] = 1
        edge1_data['_split_point'] = split_point.wkt
        edge2_data['geometry'] = shp.ops.LineString(
            list(node_point.coords) + list(new_lines.geoms[2].coords)
        )
        edge2_data['_split'] = 2
        edge1_data['_split_point'] = split_point.wkt

    # Split topology
    G.add_node(split_node_id, x=split_point.x, y=split_point.y, _split_node=True)
    G.remove_edge(u, v, key)
    key1 = G.add_edge(u, split_node_id, **edge1_data)
    key2 = G.add_edge(split_node_id, v, **edge2_data)

    # Reverse the new edge data if the new node has changed their topological direction
    if not nx.is_directed(G):
        if u > split_node_id:
            reverse_edge(G, (u, split_node_id, key1, edge1_data))
        if split_node_id > v:
            reverse_edge(G, (v, split_node_id, key2, edge2_data))

    edge1 = (u, split_node_id, key1, edge1_data)
    edge2 = (split_node_id, v, key2, edge2_data)

    # return a tuple of the two resulting edges
    return edge1, edge2


def reverse_edge(G, u, v, key, reverse_topology=True):
    """
    Flip the edge direction, including lanes and geometry

    Parameters
    ----------
    G : nx.MultiDiGraph
        street graph
    edge : tuple
        the edge to be reversed
    reverse_topology : boolean
        flip the start and end node, automatically false if the G is undirected, be careful when using this
        as is may corrupt the graph by creating inconsistencies between the topological direction and the geometry
    """

    if not G.has_edge(u, v, key):
        return

    data = G.get_edge_data(u, v, key)

    # don't reverse topology if the graph is not directed
    reverse_topology = reverse_topology and nx.is_directed(G)

    # remove the old edge
    if reverse_topology:
        G.remove_edge(u, v, key)

    # reverse lanes
    if data.get('ln_desc') is not None:
        data['ln_desc'] = lane_config.reverse_lanes(data['ln_desc'])

    # reverse sensors
    sensors_forward = data.get('sensors_forward', [])
    sensors_backward = data.get('sensors_backward', [])
    data['sensors_forward'] = sensors_backward
    data['sensors_backward'] = sensors_forward

    # reverse geometry
    if data.get('geometry') and data.get('geometry') != shp.ops.LineString():
        data['geometry'] = shp.ops.substring(data['geometry'], 1, 0, normalized=True)

    # flip the reversed flag
    data[KEY_REVERSED] = not data.get(KEY_REVERSED, False)

    # add the new edge
    if reverse_topology:
        key = G.add_edge(v, u, **data)
    else:
        nx.set_edge_attributes(G, {(u, v, key): data})

    # return the resulting edge
    if reverse_topology:
        return v, u, key, data
    else:
        return u, v, key, data
