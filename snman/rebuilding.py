import copy, math
import networkx as nx
import numpy as np
import geopandas as gpd
from . import space_allocation, hierarchy, street_graph, graph, io, merge_edges, lane_graph, access_graph
from .constants import *
from . import osmnx_customized as oxc


def multi_set_needed_node_access(G, source_lanes_attribute=KEY_LANES_DESCRIPTION):
    """
    Assigns a set of attributes to every node defining whether it must remain accessible to each mode.

    Following attributes will be added to each node:
     - needs_access_by_foot : bool
     - needs_access_by_private_cars : bool
     - needs_access_by_transit : bool
     - needs_access_by_cycling : bool

    This is a built-in default function that can be replaced by a custom one in the rebuilding process.

    Parameters
    ----------
    G
    source_lanes_attribute

    Returns
    -------

    """
    for mode in MODES:
        H = copy.deepcopy(G)
        H = street_graph.filter_lanes_by_modes(H, {mode}, source_lanes_attribute)
        for i, data in H.nodes.items():
            G.nodes[i]['needs_access_by_' + mode] = True


def multi_set_given_lanes(
        G,
        source_lanes_attribute=KEY_LANES_DESCRIPTION,
        target_lanes_attribute=KEY_GIVEN_LANES_DESCRIPTION,
        public_transit_mode='mandatory_like_existing',
        parking_mode='mandatory_like_existing',
        motorized_traffic_on_all_streets=False,
        motorized_traffic_road_hierarchies=(hierarchy.LOCAL_ROAD, hierarchy.MAIN_ROAD, hierarchy.HIGHWAY),
        cycling_infra_road_hierarchies=(hierarchy.LOCAL_ROAD, hierarchy.MAIN_ROAD),
        hierarchies_to_fix=(hierarchy.PATHWAY, hierarchy.SERVICE)
):

    """
    Builds a wishlist of lanes for every street. The 'optional' lanes can be removed in course of the rebuilding.
    This is a built-in default function that can be replaced by a custom one in the rebuilding process.

    Parameters
    ----------
    G: nx.MultiDiGraph
        street graph
    source_lanes_attribute: str
        key of the existing lanes
    target_lanes_attribute: str
        key of the rebuilt lanes
    public_transit_mode: str
        - 'mandatory_like_existing': keep all public transit lanes as they are today
        - 'none': no public transit lanes
        - 'add_dedicated': every street will get dedicated public transit lanes if there are any routes running
    parking_mode: str
        - 'mandatory_like_existing': keep all parking as it is
        - 'optional_like_existing': reduce parking if space is needed for modes with higher priority
        - 'none': no parking
    motorized_traffic_on_all_streets: bool
        if True, at least one lane for private cars will be enforced on every street
    motorized_traffic_road_hierarchies: tuple
        which hierarchies should have infrastructure for private cars
    cycling_infra_road_hierarchies: tuple
        which hierarchies should have infrastructure for cycling
    hierarchies_to_fix: tuple
        which hierarchies should be left as they are, i.e., their existing lanes are also their 'given' lanes

    Returns
    -------
    None

    """

    for uvk, data in G.edges.items():

        source_lanes = data[source_lanes_attribute]
        target_lanes = space_allocation.SpaceAllocation([])

        # for streets with "hierarchy to fix" keep everything as it is
        if data.get('hierarchy') in hierarchies_to_fix:
            target_lanes = copy.deepcopy(source_lanes)

        else:
            # -- Add dedicated public transit lanes --

            if 1:
                # all_dedicated: add optional pt lane to every street with pt
                if public_transit_mode == 'all_dedicated':
                    if data.get('pt_backward'):
                        target_lanes.append(
                            space_allocation.Lane(LANETYPE_DEDICATED_PT, DIRECTION_BACKWARD, status=STATUS_OPTIONAL)
                        )
                    if data.get('pt_forward'):
                        target_lanes.append(
                            space_allocation.Lane(LANETYPE_DEDICATED_PT, DIRECTION_FORWARD, status=STATUS_OPTIONAL)
                        )
                # mandatory_like_existing: all existing pt lanes remain as they are
                elif public_transit_mode == 'mandatory_like_existing':
                    existing_pt_lanes = space_allocation.filter_lanes_by_modes(
                        source_lanes, {MODE_TRANSIT}, operator='exact'
                    )
                    target_lanes += existing_pt_lanes
                # none: no pt lanes
                elif public_transit_mode == 'none':
                    pass
                else:
                    print('pt_mode unknown:', public_transit_mode)
                    return

            # -- Add car lanes --
            if 1:
                # save target lanes so far back to the street graph
                data[target_lanes_attribute] = target_lanes
                car_cost_forward = street_graph.calculate_edge_cost(
                    G, *uvk, DIRECTION_FORWARD, MODE_PRIVATE_CARS, lanes_description=target_lanes_attribute
                )
                car_cost_backward = street_graph.calculate_edge_cost(
                    G, *uvk, DIRECTION_BACKWARD, MODE_PRIVATE_CARS, lanes_description=target_lanes_attribute
                )

                # motorized_traffic_on_all_streets=True: ensure every street is accessible to car traffic
                # if it is not already
                if motorized_traffic_on_all_streets:
                    if data.get('hierarchy') in motorized_traffic_road_hierarchies:
                        if math.inf in [car_cost_backward, car_cost_forward]:
                            target_lanes.append(
                                space_allocation.Lane(LANETYPE_MOTORIZED, DIRECTION_TBD, status=STATUS_FIXED)
                            )
                # motorized_traffic_on_all_streets=False: add optional lane to every street
                else:
                    if data.get('hierarchy') in motorized_traffic_road_hierarchies:
                        if math.inf in [car_cost_backward, car_cost_forward]:
                            target_lanes.append(
                                space_allocation.Lane(LANETYPE_MOTORIZED, DIRECTION_TBD, status=STATUS_OPTIONAL)
                            )

            # -- Add cycling lanes --

            if 1:
                if data.get('hierarchy') in cycling_infra_road_hierarchies:
                    if data.get('hierarchy') == hierarchy.MAIN_ROAD:
                        # twice the standard width in each direction on main roads
                        factor = 2
                    else:
                        # standard width on other roads
                        factor = 2

                    target_lanes += [
                        space_allocation.Lane(
                            LANETYPE_CYCLING_LANE, DIRECTION_BACKWARD,
                            status=STATUS_OPTIONAL, width=1.3
                        )
                        for i in range(factor)
                    ]

                    target_lanes += [
                        space_allocation.Lane(
                            LANETYPE_CYCLING_LANE, DIRECTION_FORWARD,
                            status=STATUS_OPTIONAL, width=1.3
                        )
                        for i in range(factor)
                    ]

            # -- Add parking --
            if 1:
                if parking_mode == 'none':
                    pass
                elif parking_mode in ['mandatory_like_existing', 'optional_like_existing', 'existing_by_need']:
                    if parking_mode == 'mandatory_like_existing':
                        status = STATUS_FIXED
                    elif parking_mode == 'optional_like_existing':
                        status = STATUS_OPTIONAL
                    elif parking_mode == 'existing_by_need':
                        status = STATUS_BY_NEED
                    existing_parking_lanes = space_allocation.filter_lanes_by_modes(
                        source_lanes, {MODE_CAR_PARKING}, operator='exact'
                    )
                    for lane in existing_parking_lanes:
                        lane = copy.copy(lane)
                        lane.status = status
                        target_lanes.append(lane)
                elif parking_mode == 'everywhere_by_need':
                    target_lanes.append(
                        space_allocation.Lane(LANETYPE_PARKING_PARALLEL, DIRECTION_FORWARD, status=STATUS_BY_NEED)
                    )
                else:
                    print('parking mode not implemented:', parking_mode)
                    return

            # -- Add non-traffic space --

            #target_lanes.append(
            #    space_allocation.Lane(LANETYPE_NON_TRAFFIC, DIRECTION_FORWARD, status=STATUS_OPTIONAL)
            #)

        data[target_lanes_attribute] = target_lanes


def get_effective_subgraph(L, weight, node_inclusion):
    """
    Returns a subgraph respecting the weights and relevant nodes.
    The resulting subgraph excludes all edges with data[weight] == inf
    and all nodes where data.get(node_inclusion) == False

    Parameters
    ----------
    L
    weight
    node_inclusion

    Returns
    -------

    """
    # keep only the nodes to be included
    M = L.subgraph(
        [i for i, data in L.nodes.items() if data.get(node_inclusion)]
    )
    # keep only the edges without infinite weights
    M = M.edge_subgraph(
        [uvk for uvk, data in M.edges.items() if data.get(weight) != np.inf]
    )
    return M


def is_strongly_connected_plus(L, weight, node_inclusion, exclude_edges=()):
    """
    Checks if the lane graph would still fulfill the 'strongly connected' requirement after removing a set of edges.
    A helper function for the built-in, as well as any custom helper function.

    Parameters
    ----------
    G
    weight
    node_inclusion
    exclude_edges

    Returns
    -------

    """

    M = get_effective_subgraph(L, weight, node_inclusion)
    M.edge_subgraph(
        [uvk for uvk, data in M.edges.items() if uvk not in exclude_edges]
    )

    return nx.is_strongly_connected(M)


def check_if_edges_disconnect_graph(L, edges):
    """
    Checks if removing a set of edges increases
    the number of strongly connected components.

    Parameters
    ----------
    L
    edges

    Returns
    -------
    bool
    """

    M = L.edge_subgraph(
        [uvk for uvk, data in L.edges.items() if uvk not in edges]
    )

    scc_before = nx.number_strongly_connected_components(L)
    scc_after = nx.number_strongly_connected_components(M)

    return scc_after > scc_before


def check_if_edge_breaks_transit(G, L, u, v, k):
    """
    Returns True if the removal of an edge in the lane graph would break the public transit network.
    Otherwise, returns False, including the case if the public transit network is already broken on this street.

    Parameters
    ----------
    G : street_graph.StreetGraph
    L : lane_graph.LaneGraph
    u : int
        u in lane graph
    v : int
        v in lane graph
    k : str
        key in lane graph

    Returns
    -------
    bool

    """

    L_data = L.edges[(u, v, k)]
    u_G = L_data['u_G']
    v_G = L_data['v_G']
    k_G = L_data['k_G']
    G_data = G.edges[(u_G, v_G, k_G)]

    pt_forward = G_data.get('pt_forward', False)
    pt_backward = G_data.get('pt_backward', False)

    keys = set(L.edges.keys())
    keys.remove((u, v, k))
    M = L.edge_subgraph(keys)

    breaks_transit = False
    for direction in [DIRECTION_FORWARD, DIRECTION_BACKWARD]:

        direction_word = 'forward' if direction == DIRECTION_FORWARD else 'backward'

        has_pt = G_data.get(f'pt_{direction_word}', False)
        if not has_pt:
            continue

        lanes_before = lane_graph.get_lanes_by_filter(
            L, u_G, v_G, k_G,
            filter=lambda x: MODE_TRANSIT in x['lane'].get_modes(),
            direction=direction)

        lanes_after = lane_graph.get_lanes_by_filter(
            M, u_G, v_G, k_G,
            filter=lambda x: MODE_TRANSIT in x['lane'].get_modes(),
            direction=direction)

        breaks_transit = breaks_transit or (len(lanes_before) > 0) != (len(lanes_after) > 0)

    return breaks_transit


def _remove_car_lanes(
        L, L_existing,
        G, width_attribute,
        A,
        verbose=False
):
    """
    a helper for multi_rebuild(), takes care of the car lanes removal
    """

    iteration_step = 1
    uv_changed = None
    while True:

        if verbose:
            print('iteration', iteration_step)

        # calculate betweenness centrality
        #print(uv_changed)
        if iteration_step == 1 or (uv_changed is not None):
            #print('calculate bc')

            # create a subgraph that contains only car lanes, this is needed to reduce computing time
            # for betweenness centrality
            M = L.edge_subgraph([
                uvk for uvk, data in L.edges.items()
                if data['lanetype'] in {LANETYPE_MOTORIZED, LANETYPE_HIGHWAY}
            ])

            #print(len(L.nodes), len(L.edges))
            #print(len(M.nodes), len(M.edges))

            weight = 'cost_' + MODE_PRIVATE_CARS
            n_nodes = len(M.nodes)
            k = int(n_nodes / 10) if n_nodes > 300 else n_nodes
            bc = nx.edge_betweenness_centrality(M, k, weight=weight, seed=9)
            nx.set_edge_attributes(L, bc, 'bc_' + MODE_PRIVATE_CARS)

        #print('calculate_excess_width')

        # calculate excess width (how much wider are the new lanes than the old ones)
        for uvk, data in G.edges.items():
            width_before = data[width_attribute]
            width_after = lane_graph.calculate_street_width(L, *uvk)
            data['_after_width_total_m'] = width_after
            data['_after_excess_width_m'] = width_after - width_before

        #print('list of unfixed')

        # create a list of unfixed car edges, these are removal candidates
        removal_candidates_car = list(filter(
            lambda edge:
            (edge[1]['lane'].status != STATUS_FIXED or edge[1]['lane'].direction == DIRECTION_TBD)
            and edge[1]['lanetype'] in {LANETYPE_MOTORIZED, LANETYPE_HIGHWAY},
            # and G.edges[(edge[1]['u_G'],edge[1]['v_G'],edge[1]['k_G'])]['_after_excess_width_m'] > 0,
            L.edges.items()
        ))

        # stop here if no removal candidates exist
        if len(removal_candidates_car) == 0:
            break

        #print('sort')

        # sort by excess width
        removal_candidates_car = sorted(
            removal_candidates_car,
            key=lambda x: (
                G.edges[(x[1]['u_G'], x[1]['v_G'], x[1]['k_G'])]['_after_excess_width_m'],
                1 - x[1]['bc_private_cars']
            )
        )

        removal_candidate_car = removal_candidates_car[-1]
        remove_edge_uvk = removal_candidate_car[0]
        remove_edge_data = removal_candidate_car[1]
        remove_edge_street_data = G.edges[(remove_edge_data['u_G'], remove_edge_data['v_G'], remove_edge_data['k_G'])]

        opposite_edge_uvk = (remove_edge_uvk[1], remove_edge_uvk[0], remove_edge_uvk[2])

        P = L.edge_subgraph(
            [uvk for uvk, data in L.edges.items() if data.get('cost_car_parking') != np.inf]
        )
        P = P.subgraph(
            [node for node, degree in dict(G.degree()).items() if degree > 0]
        )

        for i, data in L.nodes.items():
            data['_cars_and_parking'] = data.get('needs_access_by_private_cars', False) or i in P.nodes

        # would the removal of this edge disconnect the graph?
        M = get_effective_subgraph(L, 'cost_private_cars', 'needs_access_by_private_cars')
        edges_disconnect_graph = check_if_edges_disconnect_graph(M, [remove_edge_uvk])

        # would the removal of this edge violate the public transit requirement?
        edge_breaks_transit = check_if_edge_breaks_transit(G, L, *remove_edge_uvk)

        #print('check')

        # is the last direction of a mandatory lane with direction tbd?
        is_last_direction_of_mandatory_lane = \
            remove_edge_data['lane'].status == STATUS_FIXED and remove_edge_data['lane'].direction == DIRECTION_TBD \
            and not L.has_edge(*opposite_edge_uvk)

        # is the last car lane of a street with parking?
        dependent_parking_lanes = lane_graph.get_dependent_parking_lanes(L, *remove_edge_uvk)

        if (
                dependent_parking_lanes == {}
                or access_graph.effect_of_parking_removal_on_underassignment(A, dependent_parking_lanes.keys()) == 0
        ):
            removes_needed_parking = False
        else:
            removes_needed_parking = True
            L.edges[remove_edge_uvk]['lane'].status = STATUS_FIXED

        #print('remove')

        if (
                not edges_disconnect_graph
                and not edge_breaks_transit
                and not is_last_direction_of_mandatory_lane
                and not removes_needed_parking
        ):

            # execute the removal
            if verbose:
                print('removed', remove_edge_uvk)
            L.remove_edge(*remove_edge_uvk)
            uv_changed = remove_edge_uvk[0:2]

            if dependent_parking_lanes != {}:
                L.remove_edges_from(dependent_parking_lanes.keys())
                A.remove_nodes_from(dependent_parking_lanes.keys())

        else:
            # mark as fixed
            if verbose:
                print(remove_edge_uvk, edges_disconnect_graph, is_last_direction_of_mandatory_lane)
            L.edges[remove_edge_uvk]['lane'].status = STATUS_FIXED
            if L.edges[remove_edge_uvk]['lane'].direction == DIRECTION_TBD:
                L.edges[remove_edge_uvk]['lane'].direction = DIRECTION_FORWARD
            L.edges[remove_edge_uvk]['twin_factor'] = 1

            L.edges[remove_edge_uvk]['_disconnect_graph'] = edges_disconnect_graph
            L.edges[remove_edge_uvk]['_breaks_transit'] = edge_breaks_transit
            L.edges[remove_edge_uvk]['_last_dir_of_mandatory'] = is_last_direction_of_mandatory_lane
            L.edges[remove_edge_uvk]['_removes_needed_parking'] = removes_needed_parking

            if verbose:
                print('fixed', remove_edge_uvk)
            uv_changed = None

        iteration_step += 1


def _narrow_down_car_lanes(
        L, L_existing,
        G, width_attribute,
        verbose=False
):
    """
    a helper for multi_rebuild(), replaces two full car lanes with a bidirectional one on local roads
    """

    for uvk, data in G.edges.items():

        # apply only to local roads and dead ends
        if data['hierarchy'] not in (hierarchy.LOCAL_ROAD, hierarchy.DEAD_END):
            continue

        u, v, k = uvk
        lanes = lane_graph.get_street_lanes(L, *uvk)
        motorized_lanes = {uvk: data for uvk, data in lanes.items() if data['lanetype'] == LANETYPE_MOTORIZED}
        has_forward = False
        has_backward = False
        instance = 1
        for lane_u, lane_v, lane_k in motorized_lanes.keys():
            if lane_u == u and lane_v == v:
                has_forward = True
            if lane_u == v and lane_v == u:
                has_backward = True
        if has_forward and has_backward and len(motorized_lanes) == 2:
            lane_id = list(motorized_lanes.keys())[0][2]
            template_data = L.edges[list(motorized_lanes.keys())[0]]
            del template_data['instance']
            template_data['direction'] = DIRECTION_BOTH
            template_data['width'] = LANE_TYPES[LANETYPE_MOTORIZED + DIRECTION_BOTH]['width']
            template_data['lane'] = space_allocation.Lane(LANETYPE_MOTORIZED, DIRECTION_BOTH)
            template_data['twin_factor'] = 0.5
            template_data['lane_id'] = lane_id
            for motorized_lane_uvk, motorized_lane_data in motorized_lanes.items():
                L.remove_edge(*motorized_lane_uvk)
                L.add_edge(
                    *motorized_lane_uvk[0:2], lane_id,
                    **template_data, instance=instance
                )
                instance += 1


def _remove_parking(
        L, L_existing,
        G, width_attribute,
        A,
        verbose=False
):
    """
    a helper for multi_rebuild(), takes care of the parking removal
    """

    i = 1
    while True:

        if verbose:
            print('iteration', i)
        i += 1

        # calculate excess width (how much wider are the new lanes than the old ones)
        for uvk, data in G.edges.items():
            width_before = data[width_attribute]
            width_after = lane_graph.calculate_street_width(L, *uvk)
            data['_after_width_total_m'] = width_after
            data['_after_excess_width_m'] = width_after - width_before

        # create a list of unfixed edges, these are removal candidates
        removal_candidates_parking = list(filter(
            lambda edge:
            edge[1]['lane'].status != STATUS_FIXED
            and edge[1]['lanetype'] in {LANETYPE_PARKING_PARALLEL}
            and G.edges[(edge[1]['u_G'], edge[1]['v_G'], edge[1]['k_G'])]['_after_excess_width_m'] > 0,
            L.edges.items()
        ))

        # stop here if no removal candidates exist
        if len(removal_candidates_parking) == 0:
            break

        # sort by excess width
        removal_candidates_parking = sorted(
            removal_candidates_parking,
            key=lambda x: G.edges[(x[1]['u_G'], x[1]['v_G'], x[1]['k_G'])]['_after_excess_width_m'],
        )

        removal_candidate_parking = removal_candidates_parking[-1]
        remove_edge_uvk = removal_candidate_parking[0]

        if access_graph.effect_of_parking_removal_on_underassignment(A, [remove_edge_uvk]) == 0:
            L.remove_edge(*remove_edge_uvk)
            A.remove_node(remove_edge_uvk)
            access_graph.gravity_model(A)

        else:
            L.edges[remove_edge_uvk]['lane'].status = STATUS_FIXED


def _remove_cycling_lanes(
        L, L_existing,
        G, width_attribute,
        verbose=False
):
    """
    a helper for multi_rebuild(), takes care of cycling lanes removal
    """
    iteration_step = 1
    # a list to keep track of uv pairs where any edges in the lange graph have changed
    uv_changed = []

    while True:

        if verbose:
            print('iteration', iteration_step)

        #print('calculate excess width')
        # calculate excess width (how much wider are the new lanes than the old ones)
        for uvk, data in G.edges.items():
            width_before = data[width_attribute]
            width_after = lane_graph.calculate_street_width(L, *uvk)
            data['_after_width_total_m'] = width_after
            data['_after_excess_width_m'] = width_after - width_before


        #print('calculate cost inrease on removal')
        # calculate cost increase on removal for each link
        # at first, do it for all edges but then only for those uv pairs where something changed
        for uvk, data in L.edges.items():
            if (iteration_step == 1 or uvk[0:2] == uv_changed or uvk[1::-1] == uv_changed)\
                    and data['lane'].lanetype == LANETYPE_CYCLING_LANE:
                data['_cost_increase_by_removal'] = graph.cost_increase_by_edge_removal(L, *uvk, 'cost_cycling')

        # create a list of unfixed edges, these are removal candidates
        removal_candidates_cycling = list(filter(
            lambda edge:
                # is unfixed
                edge[1]['lane'].status != STATUS_FIXED
                # is a cycling lane
                and edge[1]['lane'].lanetype == LANETYPE_CYCLING_LANE
                # is part of a street that has excess width
                and G.edges[(edge[1]['u_G'], edge[1]['v_G'], edge[1]['k_G'])]['_after_excess_width_m'] > 0,
            L.edges.items()
        ))

        # stop here if no removal candidates exist
        if len(removal_candidates_cycling) == 0:
            break

        #print('sort')
        # sort the removal candidates
        removal_candidates_cycling = sorted(
            removal_candidates_cycling,
            key=lambda x: (
                # those on streets with the largest excess width
                G.edges[(x[1]['u_G'], x[1]['v_G'], x[1]['k_G'])]['_after_excess_width_m'],
                # within these, those with the lowest importance
                x[1]['_cost_increase_by_removal'] * -1
            )
        )

        removal_candidate_cycling = removal_candidates_cycling[-1]
        remove_edge_uvk = removal_candidate_cycling[0]

        #print('strongly connected?')
        M = get_effective_subgraph(L, 'cost_cycling', 'needs_access_by_cycling')
        edges_disconnect_graph = check_if_edges_disconnect_graph(M, [remove_edge_uvk])

        #print('remove')
        if not edges_disconnect_graph:
            L.remove_edge(*remove_edge_uvk)
            uv_changed = remove_edge_uvk[0:2]
            if verbose:
                print('removed', remove_edge_uvk)
        else:
            L.edges[remove_edge_uvk]['lane'].status = STATUS_FIXED
            if verbose:
                print('fixed', remove_edge_uvk)

        iteration_step += 1


def _merge_transit_with_car_lanes(
        L, L_existing,
        G, width_attribute,
        verbose=False
):
    """
    a helper for multi_rebuild(), takes care of cycling lanes removal
    """

    i = 1
    while True:

        if verbose:
            print('iteration', i)
        i += 1

        # calculate excess width (how much wider are the new lanes than the old ones)
        for uvk, data in G.edges.items():
            width_before = data[width_attribute]
            width_after = lane_graph.calculate_street_width(L, *uvk)
            data['_after_width_total_m'] = width_after
            data['_after_excess_width_m'] = width_after - width_before

            # create a list of unfixed edges, these are removal candidates
        removal_candidates_pt = list(filter(
            lambda edge:
            edge[1]['lane'].status != STATUS_FIXED
            and edge[1]['lanetype'] in {LANETYPE_DEDICATED_PT}
            and G.edges[(edge[1]['u_G'], edge[1]['v_G'], edge[1]['k_G'])]['_after_excess_width_m'] > 0,
            L.edges.items()
        ))

        # stop here if no removal candidates exist
        if len(removal_candidates_pt) == 0:
            break

            # sort by excess width
        removal_candidates_pt = sorted(
            removal_candidates_pt,
            key=lambda x: G.edges[(x[1]['u_G'], x[1]['v_G'], x[1]['k_G'])]['_after_excess_width_m'],
        )

        removal_candidate_pt = removal_candidates_pt[-1]
        remove_edge_uvk = removal_candidate_pt[0]
        remove_edge_k_G = removal_candidate_pt[1]['k_G']

        # look for a parallel car lane
        parallel_car_lanes = lane_graph.get_street_lanes(L, *remove_edge_uvk[0:2], remove_edge_k_G,
                                                         direction=DIRECTION_FORWARD)
        parallel_car_lanes = filter(
            lambda edge:
            edge[1]['k_G'] == remove_edge_k_G
            and edge[1]['lanetype'] in {LANETYPE_MOTORIZED, LANETYPE_HIGHWAY},
            parallel_car_lanes.items()
        )
        # fix a parallel car lane (if there is such)
        parallel_car_lanes = list(sorted(
            parallel_car_lanes,
            key=lambda x: (x[1]['lane'] == STATUS_FIXED) * 1
        ))
        if len(parallel_car_lanes) > 0:
            parallel_car_lane = parallel_car_lanes[-1]
            parallel_car_lane_uvk = parallel_car_lane[0]
            L.edges[parallel_car_lane_uvk]['lane'].status = STATUS_FIXED
            L.remove_edge(*remove_edge_uvk)
            if verbose:
                print('removed', remove_edge_uvk)
        else:
            L.edges[remove_edge_uvk]['lane'].status = STATUS_FIXED
            if verbose:
                print('fixed', remove_edge_uvk)

    return L


def _adjust_width_of_lanes(
        L, L_existing,
        G, width_attribute,
        lanetypes=None,
        verbose=False
):
    """
    Extend the width of cycling lanes to fill out the remaining space

    Parameters
    ----------
    L: lane_graph.LaneGraph
    L_existing: lane_graph.LaneGraph
    G: street_graph.StreetGraph
    width_attribute: str
    verbose: bool
    """

    for uvk, data in G.edges.items():

        width_before = data[width_attribute]
        width_after = lane_graph.calculate_street_width(L, *uvk)
        excess_width = width_after - width_before
        data['_after_width_total_m'] = width_after
        data['_after_excess_width_m'] = excess_width

        if excess_width < 0:
            lanes = lane_graph.get_street_lanes(L, *uvk)
            # filter
            if lanetypes:
                lanes = {uvk: data for uvk, data in lanes.items() if data['lane'].lanetype in lanetypes}
            total_width = sum([data['lane'].width for uvk, data in lanes.items()])
            for lane_uvk, lane_data in lanes.items():
                if verbose:
                    print(uvk, -excess_width)
                lane_data['lane'].width = lane_data['lane'].width - excess_width * lane_data['lane'].width / total_width
                lane_data['width'] = lane_data['lane'].width
                lane_data['lane'].status = STATUS_FIXED


def multi_rebuild(
        L, L_existing,
        G, width_attribute,
        verbose=False
):
    """
    Default rebuilding function based on a heuristic of removing links from a lane graph.
    The changes are made inplace on the lane graph L.

    Parameters
    ----------
    L: nx.MultiLaneGraph
        given lane graph with a wishlist of edges, some of them
        will be removed throughout the process
    L_existing: nx.MultiLaneGraph
        existing lane graph for reference, not needed here but may be useful
        for alternative variations of this function
    G: nx.MultiLaneGraph
        street graph
    width_attribute: str
        the attribute key to find the width of each street in the street graph
    verbose: bool

    Returns
    -------
    None
    """

    if verbose:
        print('---- removing car lanes ------')
    _remove_car_lanes(L, L_existing, G, width_attribute, verbose)

    if verbose:
        print('---- narrow down car lanes ------')
    _narrow_down_car_lanes(L, L_existing, G, width_attribute, verbose)

    if verbose:
        print('---- removing parking ------')
    _remove_parking(L, L_existing, G, width_attribute, verbose)

    if verbose:
        print('---- removing cycling lanes ------')
    _remove_cycling_lanes(L, L_existing, G, width_attribute, verbose)

    if verbose:
        print('---- merging transit lanes with car lanes ------')
    _merge_transit_with_car_lanes(L, L_existing, G, width_attribute, verbose)

    if verbose:
        print('---- adjusting non-traffic spaces ------')
    _adjust_non_traffic(L, L_existing, G, width_attribute, verbose)

    return L


def rebuild_streets_based_on_lane_graph(
        G,
        L,
        target_lane_key=KEY_LANES_DESCRIPTION_AFTER,
        hierarchies_to_protect=[]
):
    for G_uvk, G_data in G.edges.items():

        # skip this edge if its hierarchy is protected
        if G_data['hierarchy'] in hierarchies_to_protect:
            continue

        lanes_description = []
        for direction in [DIRECTION_BACKWARD, DIRECTION_FORWARD]:
            lanes = lane_graph.get_street_lanes(L, *G_uvk, direction=direction)
            for L_uvk, L_data in lanes.items():
                lane = (L_data['lane'])
                # for bidirectional lanes, work only with the first instance
                if lane.direction in {DIRECTION_BOTH, DIRECTION_BOTH_OPTIONAL}:
                    # check if the opposite edge is present
                    if L.has_edge(L_uvk[1], L_uvk[0], L_uvk[2]):
                        # if opposite edge is present, add a bidirectional lane
                        if L_data['instance'] == 1:
                            lane.direction = DIRECTION_BOTH
                            lane.width = L_data['width']
                            lanes_description.append(lane)
                        else:
                            # do nothing with the second instance
                            pass
                    else:
                        # if not, add a single-direction lane
                        lane.direction = DIRECTION_FORWARD
                        lane.width = L_data['width'] * 2/3
                        lanes_description.append(lane)

                # for all other lanes
                else:
                    lane.direction = direction
                    lane.width = L_data['width']
                    lanes_description.append(lane)
        G_data[target_lane_key] = lanes_description


def multi_rebuild_regions(
        G,
        rebuilding_regions_gdf,
        width_attribute=KEY_LANES_DESCRIPTION + '_width_total_m',
        rebuilding_function=multi_rebuild,
        given_lanes_function=multi_set_given_lanes,
        public_transit_mode='mandatory_like_existing',
        parking_mode='mandatory_like_existing',
        needed_node_access_function=multi_set_needed_node_access,
        add_fix_hierarchies={hierarchy.PATHWAY, hierarchy.OTHER_HIERARCHY},
        existing_lanes_attribute=KEY_LANES_DESCRIPTION,
        given_lanes_attribute=KEY_GIVEN_LANES_DESCRIPTION,
        target_lanes_attribute=KEY_LANES_DESCRIPTION_AFTER,
        export_L=None, export_H=None,
        export_when='before_and_after',
        verbose=False,
):
    """
    Process each rebuilding region. By default, the redesign process is defined by the built-in functions
    multi_set_given_lanes() and multi_rebuild().

    However, they can be overridden by custom functions to generate alternative design scenarios.

    Parameters
    ----------
    G: nx.GeoDataFrame
        street graph
    rebuilding_regions_gdf: gpd.GeoDataFrame
    width_attribute: str
        the attribute where the available width of each street in the street graph is stored
    rebuilding_function: function
        The function that will be used to rebuild the network, i.e. perform the optimization.
        Per default, a built-in heuristic will be used based on link elimination from the lane graph.
        You can pass a custom function to perform this task differently or compare different approaches.
        Make sure that your function takes the same arguments as multi_rebuild().
    given_lanes_function: function
        the function that will be used to define the given lanes, i.e. a wishlist of lanes for each edge
    public_transit_mode: str
        specifies the design of public transit infrastructure, see multi_set_given_lanes()
    parking_mode: str
        specifies the design of parking infrastructure, see multi_set_given_lanes()
    needed_node_access_function
    existing_lanes_attribute: str
        the attribute with existing lanes
    given_lanes_attribute: str
        the attribute with 'given' lanes which is a wishlist of lanes created by the given_lanes_function
    target_lanes_attribute
        the attribute where the resulting lanes will be stored
    verbose: bool
        for debugging

    Returns
    -------
    None
    """

    # initialize the target lanes attribute as a copy of the given lanes
    nx.set_edge_attributes(G, nx.get_edge_attributes(G, existing_lanes_attribute), target_lanes_attribute)
    # ensure consistent edge directions
    street_graph.organize_edge_directions(G)

    for i, rebuilding_region in rebuilding_regions_gdf.iterrows():

        print('rebuilding region', i)

        if len(rebuilding_region['hierarchies_to_include']) > 0:
            hierarchies_to_include = rebuilding_region['hierarchies_to_include']
        else:
            hierarchies_to_include = hierarchy.HIERARCHIES

        print('include', hierarchies_to_include)

        hierarchies_to_fix = (
            hierarchy.HIERARCHIES.difference(hierarchies_to_include)
            .union(rebuilding_region['hierarchies_to_fix'])
            .union(add_fix_hierarchies)
        )

        print('fix', hierarchies_to_fix)

        # make a graph cutout based on the region geometry and skip this region if the resulting subgraph is empty
        H = oxc.truncate.truncate_graph_polygon(G, rebuilding_region.geometry, quadrat_width=30, retain_all=True)
        if len(H.edges) == 0:
            continue

        print('after cutout')
        print(H.nodes)

        # keep only hierarchies to include
        H = street_graph.filter_by_hierarchy(H, hierarchies_to_include)

        print('after hierarchies')
        print(H.nodes)

        # set given lanes and required access for nodes according to network rules
        given_lanes_function(
            H,
            hierarchies_to_fix=hierarchies_to_fix,
            motorized_traffic_on_all_streets=rebuilding_region['keep_all_streets'],
            public_transit_mode=public_transit_mode,
            parking_mode=parking_mode
        )

        print('after given lanes')
        print(H.nodes)

        needed_node_access_function(H)

        print('after needed nodes')
        print(H.nodes)

        # simplify the graph by removing intermediate nodes
        merge_edges.reset_intermediate_nodes(H)
        merge_edges.merge_consecutive_edges(H, distinction_attributes={KEY_LANES_DESCRIPTION_AFTER})

        print('after merging')
        print(H.nodes)

        # make lane graph and ensure it is strongly connected
        L = lane_graph.create_lane_graph(H, KEY_GIVEN_LANES_DESCRIPTION)
        #L = graph.keep_only_the_largest_connected_component(L)

        print('L')
        print(L.nodes)

        # export the lane graphs (before rebuilding) for debugging purposes
        if export_when in ['before_and_after', 'before']:
            if export_L:
                io.export_street_graph(L, export_L[0], export_L[1])
            if export_H:
                io.export_street_graph(H, export_H[0], export_H[1])

        # execute the multi rebuilding function
        L = rebuilding_function(L, None, H, width_attribute, verbose=verbose)

        # use the resulting lane graph (with edges that have not been removed) to rebuild the street graph
        rebuild_streets_based_on_lane_graph(
            H,
            L,
            hierarchies_to_protect=hierarchies_to_fix
        )

        # reconstruct the original street graph with intermediary nodes
        merge_edges.reconstruct_consecutive_edges(H)
        street_graph.organize_edge_directions(H)

        # export the lane graphs (after rebuilding) for debugging purposes
        if export_when in ['before_and_after', 'after']:
            if export_L:
                io.export_street_graph(L, *export_L)
            if export_H:
                io.export_street_graph(H, *export_H)

        # write rebuilt lanes from the subgraph into the main graph
        nx.set_edge_attributes(G, nx.get_edge_attributes(H, KEY_LANES_DESCRIPTION_AFTER), KEY_LANES_DESCRIPTION_AFTER)
