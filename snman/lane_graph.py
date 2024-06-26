import copy

from . import osmnx_customized as oxc
from . import space_allocation, geometry_tools, graph, street_graph
from .constants import *
import networkx as nx
import numpy as np


def create_lane_graph(G, lanes_attribute=KEY_LANES_DESCRIPTION):
    """
    Creates a new lane graph, derived from the street graph
    Parameters
    ----------
    G : street_graph.StreetGraph
    lanes_attribute : str
        which attribute should be used for the lane description

    Returns
    -------

    """

    G = copy.deepcopy(G)

    # initialize lane graph and copy attributes from street graph
    L = LaneGraph()
    L.graph = G.graph

    for uvk, data in G.edges.items():
        u, v, k = uvk

        lanes_list = data.get(lanes_attribute, space_allocation.SpaceAllocation([]))

        for i, lane in enumerate(lanes_list):

            lane = copy.deepcopy(lane)

            if lane == '':
                continue

            reverse = lane.direction in {DIRECTION_BACKWARD, DIRECTION_BACKWARD_OPTIONAL}

            if reverse:
                lane.reverse_direction()

            attributes = {}
            length = data['length']
            slope = data.get('grade', 0)
            lane_id = '-'.join([str(u), str(v), str(k), str(i)])
            attributes['length'] = length

            attributes['u_G'] = u
            attributes['v_G'] = v
            attributes['k_G'] = k
            attributes['lane_id_within_street'] = i

            attributes['primary_mode'] = lane.get_primary_mode()
            #attributes['lane_id'] = lane_id
            attributes['lanetype'] = lane.lanetype
            #attributes['direction'] = lane.direction
            attributes['width'] = lane.width
            attributes['osm_highway'] = data.get('highway')
            attributes['maxspeed'] = data.get('maxspeed')

            #attributes['status'] = lane.status
            #attributes['fixed'] = lane.status == STATUS_FIXED
            #attributes['mandatory_lane'] = lane.status == STATUS_ONE_DIRECTION_MANDATORY
            #attributes['coupled_with_opposite_direction'] = lane.direction == DIRECTION_BOTH
            if not reverse:
                costs = {}
                for mode in MODES:
                    cost = space_allocation._calculate_lane_cost(lane, length, slope, mode)
                    costs['cost_' + mode] = cost
                L.add_edge(
                    u, v, lane_id, **{**attributes, **costs}, lane=lane, backward=0, instance=1,
                    geometry=data.get('geometry')
                )
            if reverse:
                costs = {}
                for mode in MODES:
                    cost = space_allocation._calculate_lane_cost(lane, length, -slope, mode)
                    costs['cost_' + mode] = cost
                L.add_edge(
                    v, u, lane_id, **{**attributes, **costs}, lane=lane, backward=1, instance=1,
                    geometry=geometry_tools.reverse_linestring(data.get('geometry'))
                )
            if lane.direction in [DIRECTION_BOTH, DIRECTION_TBD]:
                costs = {}
                for mode in MODES:
                    cost = space_allocation._calculate_lane_cost(lane, length, slope, mode)
                    costs['cost_' + mode] = cost
                L.add_edge(
                    u, v, lane_id, **{**attributes, **costs},
                    lane=lane, backward=0, instance=1,
                    geometry=data.get('geometry')
                )
                costs = {}
                for mode in MODES:
                    cost = space_allocation._calculate_lane_cost(lane, length, -slope, mode)
                    costs['cost_' + mode] = cost
                opposite_lane = copy.copy(lane)
                L.add_edge(
                    v, u, lane_id, **{**attributes, **costs},
                    lane=opposite_lane, backward=1, instance=2,
                    geometry=geometry_tools.reverse_linestring(data.get('geometry'))
                )

    # take over the node attributes from the street graph
    nx.set_node_attributes(L, dict(G.nodes))

    return L



def calculate_stats(L, mode):
    """
    Calculates a set of standardized measures for a lane graph:

    * **usable_N_nodes**: number of nodes that are accessible for the given mode
    * **usable_N_edges**: number of edges that are accessible for the given mode
    * **convex_hull_km2**: an approximate area of the network using a convex hull of all accessible nodes
    * **usable_lane_km**: total length of lanes that are accessible for the given mode,
        bidirectional lanes count only once,
        pseudo-lanes with width=0 don't count
    * **usable_surface_km2**: total surface of lanes that are accessible for the given mode
    * **as_primary_lane_km**: total length of lanes where the given mode is primary
        (see constants.LANE_TYPES, the primary mode is the first one in the mode list of each lane type),
        bidirectional lanes count only once,
        pseudo-lanes with width=0 don't count
    * **as_primary_mode_lane_surface_km2**: total surface of lanes where the given mode is primary
    * **avg_betweenness_centrality_norm**: average normalized betweenness centrality of edges,
        the edge weights consider the scaling of cycling distanced by comfort
    * **avg_shortest_path_vod_km**: average length of the shortest path,
        for cycling, the length is scaled using the comfort factors (cycling infrastructure will lead to lower values),
        this means that you cannot compare the absolute average shortest path across modes, only the change before/after
    * **avg_shortest_path_km**: average length of the shortest path, without scaling for cycling

    Parameters
    ----------
    L : nx.MultiDiGraph
        lane graph
    mode : str
        for which mode should the stats be calculated, see constants.MODES

    Returns
    -------
    dict

    """

    L = copy.deepcopy(L)

    # calculate the approximate area as a convex hull of all nodes
    points_gpd = oxc.graph_to_gdfs(L, edges=False)
    area_km2 = points_gpd.geometry.unary_union.convex_hull.area / pow(1000, 2)

    # set betweenness centrality
    nx.set_edge_attributes(L, nx.edge_betweenness_centrality(L, normalized=True, weight='cost_'+mode), 'bc')

    return {
        '_L': copy.deepcopy(L),
        '_mode': mode,
        'usable_N_nodes': len(L.nodes),
        'usable_N_edges':
            round(
                sum(
                    [
                        1
                        for uvk, e
                        in L.edges.items()
                        if e['instance'] == 1
                    ]
                ),
                3
            ),
        'convex_hull_km2': area_km2,
        'usable_lane_km':
            round(
                sum(
                    [
                        e['length']
                        for uvk, e
                        in L.edges.items()
                        if e['lane'].width > 0 and e['instance'] == 1
                    ]
                ) / 1000,
                3
            ),
        'usable_lane_surface_km2':
            round(
                sum(
                    [
                        e['length'] * e['lane'].width
                        for uvk, e
                        in L.edges.items()
                        if e['instance'] == 1
                    ]
                ) / pow(1000, 2),
                3
            ),
        'as_primary_mode_lane_km':
            round(
                sum(
                    [
                        e['length'] * (e['primary_mode'] == mode)
                        for uvk, e
                        in L.edges.items()
                        if e['lane'].width > 0 and e['instance'] == 1
                    ]
                ) / 1000,
                3
            ),
        'as_primary_mode_lane_surface_km2':
            round(
                sum(
                    [
                        e['length'] * e['lane'].width * (e['primary_mode'] == mode)
                        for uvk, e
                        in L.edges.items()
                        if e['instance'] == 1
                    ]
                ) / pow(1000, 2),
                3
            ),
        'as_primary_mode_lane_surface_km2':
            round(
                sum(
                    [
                        e['length'] * e['lane'].width * (e['primary_mode'] == mode)
                        for uvk, e
                        in L.edges.items()
                        if e['instance'] == 1
                    ]
                ) / pow(1000, 2),
                3
            ),
        'avg_betweenness_centrality_norm':
            round(
                sum(
                    [
                        e['bc']
                        for uvk, e
                        in L.edges.items()
                        if e['instance'] == 1
                    ]
                ) /
                sum(
                    [
                        1
                        for uvk, e
                        in L.edges.items()
                        if e['instance'] == 1
                    ]
                ),
                5
            ),
        'avg_shortest_path_vod_km':
            round(
                nx.average_shortest_path_length(L, 'cost_' + mode) / 1000,
                3
            ),
        'avg_shortest_path_km':
            round(
                nx.average_shortest_path_length(L, 'length') / 1000,
                3
            ),
    }


def get_street_lanes(L, u_G, v_G, k_G, direction=None, only_first_instance=False):
    """
    Returns a dictionary of edges from the lane graph that correspond to a particular street in the street graph

    Parameters
    ----------
    L: nx.MultiDiGraph
        lane graph
    u_G: int
        u in street graph
    v_G: int
        v in street graph
    k_G: int
        k in street graph
    direction: str
        returns only the lanes in forward or backward direction if specified,
        returns all lanes if None
    only_first_instance:
        of true, returns only the first instance of every bidirectional lane
    Returns
    -------
    dict
    """

    forward = dict(L.get_edge_data(u_G, v_G, default={}))
    forward = {(u_G, v_G, k): value for k, value in forward.items()}
    backward = dict(L.get_edge_data(v_G, u_G, default={}))
    backward = {(v_G, u_G, k): value for k, value in backward.items()}

    if direction is None:
        # concat forward and backward lanes
        lanes = forward
        lanes.update(backward)
    elif direction == DIRECTION_FORWARD:
        lanes = forward
    elif direction == DIRECTION_BACKWARD:
        lanes = backward

    # filter for those matching the k_G key which means they belong to the correct (parallel) street
    lanes = {uvk: data for uvk, data in lanes.items() if data['k_G'] == k_G}
    if only_first_instance:
        lanes = {uvk: data for uvk, data in lanes.items() if data['instance'] == 1}
    return lanes


def calculate_street_width(L, u_G, v_G, k_G):
    """
    Calculates the width of a street in a lane graph

    Parameters
    ----------
    L: nx.MultiDiGraph
        lane graph
    u_G: int
        u in street graph
    v_G: int
        v in street graph
    k_G:int
        k in street graph

    Returns
    -------
    int
    """

    lanes = get_street_lanes(L, u_G, v_G, k_G)
    widths = map(lambda lane: lane['lane'].width * (lane['instance'] == 1), lanes.values())
    total_width = sum(widths)
    return total_width


def get_horizontal_position_of_lane(L, u, v, k):

    data = L.edges[(u, v, k)]

    lanes = get_street_lanes(L, data['u_G'], data['v_G'], data['k_G'])
    lanes = {uvk[2]: lane for uvk, lane in lanes.items()}
    lanes = {k: lane for k, lane in sorted(lanes.items())}
    widths = [lane['lane'].width for lane in lanes.values()]
    widths = np.cumsum(widths) - np.array(widths)/2 - sum(widths) / 2
    widths = dict(zip(lanes.keys(), widths))

    return widths[k]


def get_modes_of_street(L, u_G, v_G, k_G):
    M = L.subgraph([u_G, v_G])
    M = M.edge_subgraph([uvk for uvk, data in M.edges.items() if data['k_G'] == k_G])
    modes_M = [data['lane'].get_modes() for uvk, data in M.edges.items()]
    modes_M = set.union(*modes_M) if len(modes_M) > 0 else set()
    return modes_M


def get_lanes_by_mode(L, u_G, v_G, k_G, mode):
    M = L.subgraph([u_G, v_G])
    M = M.edge_subgraph(
        [uvk for uvk, data in M.edges.items() if (data['k_G'] == k_G and mode in data['lane'].get_modes())]
    )
    return dict(M.edges.items())


def get_lanes_by_filter(L, u_G, v_G, k_G, filter=lambda x: True, direction=DIRECTION_BOTH):

    M = L.subgraph([u_G, v_G])

    if direction == DIRECTION_FORWARD:
        M = M.edge_subgraph(
            [uvk for uvk, data in M.edges.items() if uvk[0:2] == (u_G, v_G)]
        )
    if direction == DIRECTION_BACKWARD:
        # exclude cyclical edges so that they are not double counted
        M = M.edge_subgraph(
            [uvk for uvk, data in M.edges.items() if uvk[0:2] == (v_G, u_G) and v_G != u_G]
        )

    M = M.edge_subgraph(
        [uvk for uvk, data in M.edges.items() if (data['k_G'] == k_G and filter(data))]
    )

    return dict(M.edges.items())


def get_dependent_parking_lanes(L, u, v, k):
    # u, v, k = (2253,2248,'2253-2248-0-3')
    data = L.edges[(u, v, k)]
    u_G = data['u_G']
    v_G = data['v_G']
    k_G = data['k_G']

    # reduce the graph to the necessary size
    M = L.subgraph([u_G, v_G])

    modes_M = get_modes_of_street(M, u_G, v_G, k_G)

    # create another subgraph with the chosen lane removed
    N = L.edge_subgraph(
        [uvk for uvk, data in M.edges.items() if uvk != (u, v, k)]
    )

    modes_N = get_modes_of_street(N, u_G, v_G, k_G)

    if (
            MODE_CAR_PARKING in modes_M and MODE_PRIVATE_CARS in modes_M and
            MODE_CAR_PARKING in modes_N and MODE_PRIVATE_CARS not in modes_N
    ):
        return get_lanes_by_mode(M, u_G, v_G, k_G, MODE_CAR_PARKING)

    else:
        return dict()


def merge_lanes_and_equalize_widths(L, u_G, v_G, k_G, filter=lambda lane: True):
    """
    Gathers all lanes that satisfy the filter.
    Then deletes all duplicate lanes in the same direction and equalizes the widths of the remaining lanes.

    Example usage: Keep only one cycling lane in each direction and make the resulting lanes equally wide.

    Parameters
    ----------
    L: lane_graph.LaneGraph
    u_G: int
    v_G: int
    k_G: int
    filter: function
    """
    lanes_by_direction = {}
    total_width = 0
    for direction in [DIRECTION_BACKWARD, DIRECTION_FORWARD]:
        lanes = get_lanes_by_filter(
            L, u_G, v_G, k_G,
            filter=filter,
            direction=direction
        )
        if len(lanes) > 0:
            lanes_by_direction[direction] = lanes
            total_width += sum([data['lane'].width for uvk, data in lanes.items()])

    for direction, lanes in lanes_by_direction.items():
        new_width = total_width / len(lanes_by_direction)
        lanes[list(lanes.keys())[0]]['lane'].width = new_width
        lanes[list(lanes.keys())[0]]['width'] = new_width
        for uvk in list(lanes.keys())[1:]:
            L.remove_edge(*uvk)
            print('removed', uvk)


class LaneGraph(graph.SNManMultiDiGraph):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

