from . import osmnx_customized as oxc
from . import graph, space_allocation, _errors, geometry_tools
from .constants import *
import geopandas as gpd
import networkx as nx
import pyproj
import copy
import shapely
import numpy as np


def get_subgraph_with_invalid_geometries(G):
    filtered_edges = filter(
        lambda e: e[3].get('geometry') is not None and e[3].get('geometry').is_valid > 2,
        G.edges(keys=True, data=True)
    )
    return nx.Graph(filtered_edges)


def get_subgraph_with_empty_geometries(G):
    filtered_edges = filter(
        lambda e: e[3].get('geometry') is None or e[3].get('geometry').is_empty > 2,
        G.edges(keys=True, data=True)
    )
    return nx.Graph(filtered_edges)


def organize_edge_directions(G, method='lower_to_higher_node_id', key_lanes_description=KEY_LANES_DESCRIPTION):
    """
    Ensure that all edges in a street graph have the direction from the lower to the higher node id.
    Edges with a different direction will be reversed, including the lanes and their geometry.

    Parameters
    ----------
    G : nx.MultiDiGraph
        OSM graph
    method : str
        - lower_to_higher_node_id: each edge will be digitized from the lower to the higher node id
        - by_osm_convention: one-way edges will be digitized according to their lane direction
        - by_top_order_lanes: each edges will be digitized such that more >=50% of one-way lanes with the highest order
        point forward

    Returns
    -------
    None
    """

    edges = list(G.edges(data=True, keys=True))
    for edge in edges:

        if method == 'lower_to_higher_node_id':
            if edge[0] > edge[1]:
                reverse_edge(G, *edge[0:3])

        elif method == 'by_osm_convention':
            if space_allocation.is_backward_oneway_street(edge[3].get(key_lanes_description)):
                reverse_edge(G, *edge[0:3])

        elif method == 'by_top_order_lanes':
            if space_allocation.is_backward_by_top_order_lanes(edge[3].get(key_lanes_description)):
                reverse_edge(G, *edge[0:3])

        else:
            raise _errors.OptionNotImplemented('Method ' + str(method) + ' is not valid')


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


def add_connected_component_ids(G):
    """
    For directed graphs: Adds IDs of weakly ('_weakly_connected_component')
    and strongly ('_strongly_connected_component') connected components to all edges
    For undirected graphs: Adds IDs of connected components ('_connected_component') to all edges

    Parameters
    ----------
    G : nx.MultiGraph
    Returns
    -------
    None
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


def convert_crs(G, to_crs):
    """
    Convert the coordinate reference system of the geometries in a graph.
    The graph must follow the convention of networkx:
        * edges have a 'geometry' attribute with a shapely geometry
        * nodes have 'x' and 'y' attributes holding the coordinates as floats
        * the graph has a 'crs' attribute

    Parameters
    ----------
    G : nx.MultiGraph
        centerline graph
    to_crs : int
        target crs

    Returns
    -------
    None
    """

    # Initialize the CRS transformer
    from_crs = pyproj.CRS(G.graph['crs'])
    to_crs = pyproj.CRS(to_crs)
    project = pyproj.Transformer.from_crs(from_crs, to_crs, always_xy=True).transform

    # Update the street_graph's metadata
    G.graph["crs"] = to_crs

    # Transform the geometry of all edges
    for edge in G.edges(data=True, keys=True):
        if "geometry" in edge[3]:
            edge[3]["geometry"] = shapely.ops.transform(project, edge[3]["geometry"])

    # Transform the geometry of all nodes
    for id, data in G.nodes.items():
        geom = shapely.geometry.Point(data.get('x'), data.get('y'))
        geom = shapely.ops.transform(project, geom)
        data['x'] = geom.x
        data['y'] = geom.y


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


def split_edge(G, u, v, key, split_points):

    if not G.has_edge(u, v, key):
        return

    edge_data = G.get_edge_data(u, v, key)
    edge_linestring = edge_data.get('geometry', False)

    if edge_linestring is False:
        return

    # snap split points to the edge linestring
    split_points = shapely.ops.nearest_points(edge_linestring, split_points)[0]
    # make a small buffer around the points to deal with numeric errors
    split_circles = shapely.MultiPoint(split_points).buffer(0.5)
    # split the edge linestring using the buffer circles
    segments = shapely.ops.split(edge_linestring, split_circles)

    # generate a list of new node points (where the linestring has been split)
    node_points = list(map(lambda segment: segment.centroid, list(segments.geoms)[1::2]))
    # add first and last node of the linestring to the new node points
    node_points = [shapely.Point(edge_linestring.coords[0])] + node_points + [shapely.Point(edge_linestring.coords[-1])]

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
        print(u, v, key)
        print(edge_data.get('geometry', shapely.Point()).wkt)
        print(shapely.MultiPoint(split_points).wkt)
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
            new_edge_data['geometry'] = shapely.LineString(
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


def reverse_edge(
        G, u, v, key,
        reverse_topology=True,
        lane_description_keys={KEY_LANES_DESCRIPTION, KEY_LANES_DESCRIPTION_AFTER, KEY_GIVEN_LANES_DESCRIPTION}
):
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
    lane_description_keys : set
        which lane description keys should be reversed
    """

    if not G.has_edge(u, v, key):
        return

    data = G.get_edge_data(u, v, key)

    # don't reverse topology if the graph is not directed
    reverse_topology = reverse_topology and nx.is_directed(G)

    # reverse lanes
    for lane_description_key in lane_description_keys:
        allocation = data.get(lane_description_key)
        if allocation is not None:
            allocation.reverse_allocation()

    # TODO: auto identify attributes to be reversed based on parts of their name, e.g, '>'/'<', or 'forward'/'backward'
    # reverse sensors
    sensors_forward = data.get('sensors_forward', [])
    sensors_backward = data.get('sensors_backward', [])
    data['sensors_forward'] = sensors_backward
    data['sensors_backward'] = sensors_forward

    # reverse geometry
    if data.get('geometry') and data.get('geometry') != shapely.ops.LineString():
        try:
            data['geometry'] = geometry_tools.reverse_linestring(data['geometry'])
        except AttributeError:
            pass

    # flip the reversed flag
    data[KEY_REVERSED] = not data.get(KEY_REVERSED, False)

    # remove the old edge
    if reverse_topology:
        G.remove_edge(u, v, key)

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


def delete_edges_without_lanes(G, lane_description_key=KEY_LANES_DESCRIPTION):
    edges_without_lanes = dict(filter(lambda x: len(x[1]) == 0, nx.get_edge_attributes(G, lane_description_key).items()))
    G.remove_edges_from(edges_without_lanes)


def filter_lanes_by_modes(G, modes, lane_description_key=KEY_LANES_DESCRIPTION, delete_empty_edges=True, **kwargs):

    for uvk, data in G.edges.items():
        lanes = data.get(lane_description_key, [])
        data[lane_description_key] = space_allocation.filter_lanes_by_modes(lanes, modes, **kwargs)

    if delete_empty_edges:
        delete_edges_without_lanes(G, lane_description_key=lane_description_key)
        graph.remove_isolated_nodes(G)

    return G


def filter_by_hierarchy(G, hierarchy_levels):
    edges = dict(filter(lambda x: x[1] in hierarchy_levels, nx.get_edge_attributes(G, 'hierarchy').items()))
    return G.edge_subgraph(edges).copy()


def calculate_edge_cost(G, u, v, k, direction, mode, lanes_description=KEY_LANES_DESCRIPTION, include_tentative=False):

    cost_list = []
    uvk = (u, v, k)
    data = G.edges[uvk]
    for lane in data[lanes_description]:
        length = data['length']
        slope = data['grade']
        cost = space_allocation._calculate_lane_cost(
            lane, length, slope, mode, direction, include_tentative=include_tentative
        )
        cost_list.append(cost)

    return min(cost_list) if len(cost_list) > 0 else np.inf


def add_edge_costs(G, lanes_description=KEY_LANES_DESCRIPTION):

    for uvk, data in G.edges.items():
        for mode in MODES:
            for direction in ['<', '>']:
                cost = calculate_edge_cost(G, *uvk, direction, mode, lanes_description=lanes_description)
                key = '_'.join(['cost', lanes_description, mode, direction])
                data[key] = cost


def add_pseudo_cycling_lanes(G, lanes_description=KEY_LANES_DESCRIPTION):

    for uvk, data in G.edges.items():
        lanes = data[lanes_description]
        lane_types = [lane.lanetype for lane in lanes]
        has_m_lanes = LANETYPE_MOTORIZED in lane_types

        # ignore if there are no lanes for private cars
        if not has_m_lanes:
            continue

        for direction in [DIRECTION_FORWARD, DIRECTION_BACKWARD]:
            cycling_cost = calculate_edge_cost(G, *uvk, direction, MODE_CYCLING, lanes_description=lanes_description)

            if cycling_cost == np.Inf:
                lanes.append(LANETYPE_CYCLING_PSEUDO + direction)


def clone(G, edges=True):

    H = nx.MultiDiGraph()
    H.graph = G.graph
    H.nodes = G.nodes

    if edges:
        for uvk, data in G.edges.items():
            H.add_edge(*uvk, **copy.deepcopy(data))

    return H


def separate_edges_for_lane_directions(G, lanes_key=KEY_LANES_DESCRIPTION):

    H = clone(G, edges=False)

    for uvk, data in G.edges.items():

        u, v, k = uvk
        lanes = data[lanes_key]
        lanes_forward = []
        lanes_backward = []

        for lane in lanes:
            lp = space_allocation._lane_properties(lane)
            if lp.direction == DIRECTION_FORWARD:
                lanes_forward.append(lane)
            elif lp.direction == DIRECTION_BACKWARD:
                lanes_backward.append(lane)
            # convert bidirectional lanes into two separate oneway lanes
            elif lp.direction == DIRECTION_BOTH:
                lp_forward = space_allocation._lane_properties(lane)
                lp_forward.direction = DIRECTION_FORWARD
                lanes_forward.append(str(lp_forward))
                lp_backward = space_allocation._lane_properties(lane)
                lp_backward.direction = DIRECTION_BACKWARD
                lanes_backward.append(str(lp_backward))

        new_data = copy.deepcopy(data)
        del new_data[lanes_key]
        del new_data['sensors_forward'], new_data['sensors_backward']

        if len(lanes_forward) > 0:
            H.add_edge(
                u, v, **new_data,
                **{lanes_key: lanes_forward},
                sensors_forward=data['sensors_forward'],
                sensors_backward=[]
            )

        if len(lanes_backward) > 0:
            k = H.add_edge(
                u, v, **new_data,
                **{lanes_key: lanes_backward},
                sensors_forward=[],
                sensors_backward=data['sensors_backward']
            )
            reverse_edge(H, u, v, k)

    organize_edge_directions(H, method='by_osm_convention')
    space_allocation.update_osm_tags(H, lanes_description_key=lanes_key)
    add_edge_costs(H, lanes_description=lanes_key)

    return H


class StreetGraph(graph.SNManMultiDiGraph):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
