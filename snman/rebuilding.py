import copy, math
import networkx as nx
import numpy as np
import geopandas as gpd
from numpy.f2py.auxfuncs import throw_error

from . import space_allocation, hierarchy, street_graph, graph, io, merge_edges, lane_graph, access_graph, _errors, utils
from .constants import *
from . import osmnx_customized as oxc


def multi_set_needed_node_access(
        G,
        source_lanes_attribute=KEY_LANES_DESCRIPTION,
        method='maintain_access_to_nodes',
        modes=MODES
):
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
    for mode in modes:

        H = copy.deepcopy(G)
        if method == 'maintain_access_to_nodes':
            H = street_graph.filter_lanes_by_modes(H, {mode}, lane_description_key=source_lanes_attribute)

        elif method == 'access_to_parking':
            street_graph.filter_lanes_by_function(
                H,
                lambda lane: lane.status == STATUS_FIXED, lane_description_key=source_lanes_attribute
            )

        else:
            pass

        for i, data in H.nodes.items():
            G.nodes[i]['needs_access_by_' + mode] = True


def multi_set_given_lanes(
        G,
        source_lanes_attribute=KEY_LANES_DESCRIPTION,
        target_lanes_attribute=KEY_GIVEN_LANES_DESCRIPTION,
        forced_given_lanes_attribute=KEY_FORCED_GIVEN_LANES_DESCRIPTION,
        public_transit_mode='mandatory_like_existing',
        parking_mode='mandatory_like_existing',
        bicycle_parking_mode='mandatory_like_existing',
        cycling_infrastructure_mode='optional',
        cycling_lane_width=0.9,
        motorized_traffic_on_all_streets=False,
        motorized_traffic_road_hierarchies=(hierarchy.LOCAL_ROAD, hierarchy.MAIN_ROAD, hierarchy.HIGHWAY, hierarchy.SERVICE),
        motorized_traffic_lane_mode='separate_lanes',
        hierarchies_with_cycling_lanes=(hierarchy.LOCAL_ROAD, hierarchy.MAIN_ROAD),
        hierarchies_to_fix=(hierarchy.PATHWAY),
        non_traffic_mode='mandatory_like_existing'
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
    non_traffic_mode: str
        - 'mandatory_like_existing': keel all non-traffic spaces as they are

    Returns
    -------
    None

    """

    for uvk, data in G.edges.items():

        edge_hierarchy = data.get('hierarchy')
        normal_lane_width = normal_lane_width_by_hierarchy.get(
            edge_hierarchy,
            normal_lane_width_by_hierarchy[hierarchy.LOCAL_ROAD]
        )

        source_lanes = data[source_lanes_attribute]
        target_lanes = space_allocation.SpaceAllocation([])

        # for streets with "hierarchy to fix" keep everything as it is
        if data.get('hierarchy') in hierarchies_to_fix:
            target_lanes = copy.deepcopy(source_lanes)

        # transfer forced given lanes if they are defined
        elif type(data.get(forced_given_lanes_attribute)) == space_allocation.SpaceAllocation:
            target_lanes = data.get(forced_given_lanes_attribute)

        else:
            # -- Add dedicated public transit lanes --

            if 1:
                # all_dedicated: add optional pt lane to every street with pt
                if public_transit_mode == 'all_dedicated':
                    if data.get('pt_backward'):
                        target_lanes.append(
                            space_allocation.Lane(LANETYPE_DEDICATED_PT, DIRECTION_BACKWARD, status=STATUS_OPTIONAL, width=normal_lane_width)
                        )
                    if data.get('pt_forward'):
                        target_lanes.append(
                            space_allocation.Lane(LANETYPE_DEDICATED_PT, DIRECTION_FORWARD, status=STATUS_OPTIONAL, width=normal_lane_width)
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

                if motorized_traffic_lane_mode == 'separate_lanes':
                    use_direction = DIRECTION_TBD
                elif motorized_traffic_lane_mode == 'bidirectional_lanes':
                    use_direction = DIRECTION_BOTH
                else:
                    raise _errors.OptionNotImplemented(
                        f"Motorized traffic lane mode {motorized_traffic_lane_mode} is not implemented."
                    )

                # motorized_traffic_on_all_streets=True: ensure every street is accessible to car traffic
                # if it is not already
                if motorized_traffic_on_all_streets:
                    if data.get('hierarchy') in motorized_traffic_road_hierarchies:
                        if math.inf in [car_cost_backward, car_cost_forward]:
                            if use_direction != DIRECTION_BOTH:
                                target_lanes.append(
                                    space_allocation.Lane(LANETYPE_MOTORIZED, use_direction, status=STATUS_FIXED, width=normal_lane_width)
                                )
                            else:
                                target_lanes.append(
                                    space_allocation.Lane(LANETYPE_MOTORIZED, use_direction, status=STATUS_ONE_DIRECTION_MANDATORY)
                                )
                # motorized_traffic_on_all_streets=False: add optional lane to every street
                else:
                    if data.get('hierarchy') in motorized_traffic_road_hierarchies:
                        if math.inf in [car_cost_backward, car_cost_forward]:
                            if use_direction != DIRECTION_BOTH:
                                target_lanes.append(
                                    space_allocation.Lane(LANETYPE_MOTORIZED, use_direction, status=STATUS_OPTIONAL, width=normal_lane_width)
                                )
                            else:
                                target_lanes.append(
                                    space_allocation.Lane(LANETYPE_MOTORIZED, use_direction, status=STATUS_OPTIONAL)
                                )

            # -- Add cycling lanes --

            # in this mode, all mixed cycling and pedestrian paths (X lanes) will be kept
            if cycling_infrastructure_mode=='x_mandatory_like_existing':
                existing_x_lanes = space_allocation.filter_lanes_by_lanetypes(
                    source_lanes, {LANETYPE_FOOT_CYCLING_MIXED},
                )
                target_lanes += existing_x_lanes

            if 1:
                if data.get('hierarchy') in hierarchies_with_cycling_lanes:
                    if data.get('hierarchy') == hierarchy.MAIN_ROAD:
                        # twice the standard width in each direction on main roads
                        factor = 1
                    else:
                        # standard width on other roads
                        factor = 1

                    # Bike lane widths are set to a small value to avoid removing substandard infrastructure
                    # TODO: Make the bike lane widths user-defined
                    target_lanes += [
                        space_allocation.Lane(
                            LANETYPE_CYCLING_LANE, DIRECTION_BACKWARD,
                            status=STATUS_OPTIONAL, width=cycling_lane_width
                        )
                        for i in range(factor)
                    ]

                    target_lanes += [
                        space_allocation.Lane(
                            LANETYPE_CYCLING_LANE, DIRECTION_FORWARD,
                            status=STATUS_OPTIONAL, width=cycling_lane_width
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
                        lane = copy.deepcopy(lane)
                        lane.status = status
                        target_lanes.append(lane)
                elif parking_mode == 'everywhere_by_need':
                    target_lanes.append(
                        space_allocation.Lane(LANETYPE_PARKING_PARALLEL, DIRECTION_FORWARD, status=STATUS_BY_NEED)
                    )
                    #print(space_allocation.Lane(LANETYPE_PARKING_PARALLEL, DIRECTION_FORWARD, status=STATUS_BY_NEED))
                    #print(target_lanes)
                else:
                    print('parking mode not implemented:', parking_mode)
                    return

            # -- Add bicycle parking --

            if bicycle_parking_mode == 'mandatory_like_existing':
                existing_bicycle_parking_lanes = space_allocation.filter_lanes_by_lanetypes(
                    source_lanes, {LANETYPE_BICYCLE_PARKING}
                )
                for lane in existing_bicycle_parking_lanes:
                    lane = copy.deepcopy(lane)
                    lane.status = STATUS_FIXED
                    target_lanes.append(lane)

            # -- Add non-traffic space --

            if non_traffic_mode == 'mandatory_like_existing':
                existing_non_traffic_lanes = space_allocation.filter_lanes_by_lanetypes(
                    source_lanes, {LANETYPE_NON_TRAFFIC}
                )
                for lane in existing_non_traffic_lanes:
                    lane = copy.deepcopy(lane)
                    lane.status = STATUS_FIXED
                    target_lanes.append(lane)


            # reorder the lanes to match a convention
            seed_side = space_allocation.assign_seed_side(data['geometry'])
            target_lanes = space_allocation._reorder_lanes_on_edge(target_lanes, seed_side=seed_side)

        data[target_lanes_attribute] = target_lanes
        #print(target_lanes)


def get_effective_subgraph(L, weight=None, node_inclusion=None):
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

    M = copy.deepcopy(L)

    if node_inclusion:
        # keep only the nodes to be included
        M = M.subgraph(
            [i for i, data in M.nodes.items() if data.get(node_inclusion)]
        )

    if weight:
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


def get_number_of_scc_with_inclusion_attribute(L, inclusion_attribute):

    L = copy.deepcopy(L)

    scc = list(nx.strongly_connected_components(L))

    # Create a mapping from node to SCC ID
    node_to_scc_id = {}
    for idx, component in enumerate(scc):
        for node in component:
            node_to_scc_id[node] = idx

    # Add SCC ID as a node attribute
    nx.set_node_attributes(L, node_to_scc_id, name='scc_id')

    scc_ids = set()
    for i, data in L.nodes.items():
        if data[inclusion_attribute] == True:
            scc_ids.add(data['scc_id'])

    #print(scc_ids)
    return len(scc_ids)

def check_if_edges_disconnect_graph(L, edges, inclusion_attribute=None):
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

    if inclusion_attribute:
        scc_before = get_number_of_scc_with_inclusion_attribute(L, inclusion_attribute)
        scc_after = get_number_of_scc_with_inclusion_attribute(M, inclusion_attribute)

    else:
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


def check_if_has_dependent_parking_lanes(L, u, v, k):
    """
    Returns true if the lane has at least one dependent parking lane,
    i.e., being the last car lane with a parallel parking lane.
    In that case, removing it would break access to the parking lane.

    Parameters
    ----------
    L: lane_graph.LaneGraph
    u: int
    v: int
    k: str

    Returns
    -------
    bool

    """

    return len(lane_graph.get_dependent_parking_lanes(L, u, v, k)) > 0


def _remove_car_lanes(
        L, L_existing,
        G, width_attribute,
        A,
        mode='lowest_bc',
        stop_at_step=np.inf,
        verbose=False
):
    """
    a helper for multi_rebuild(), takes care of the car lanes removal

    Parameters
    ----------
    L: lane_graph.LaneGraph
    L_existing: lane_graph.LaneGraph
    G: street_graph.StreetGraph
    width_attribute: str
    A: access_graph.AccessGraph
    mode: str
    stop_at_step: int
    verbose: bool
    """

    iteration_step = 1
    uv_changed = None
    while iteration_step <= stop_at_step:

        if verbose:
            print('iteration', iteration_step)

        # calculate betweenness centrality
        if iteration_step == 1 or (uv_changed is not None):

            # create a subgraph that contains only car lanes, this is needed to reduce computing time
            # for betweenness centrality
            M = L.edge_subgraph([
                uvk for uvk, data in L.edges.items()
                if data['lanetype'] in {LANETYPE_MOTORIZED, LANETYPE_HIGHWAY}
            ])

            n_nodes = len(M.nodes)
            k = int(n_nodes / 10) if n_nodes > 300 else n_nodes
            bc = nx.edge_betweenness_centrality(M, k, weight='cost_' + MODE_PRIVATE_CARS, seed=9)
            nx.set_edge_attributes(L, bc, 'bc_' + MODE_PRIVATE_CARS)


        # calculate excess width (how much wider are the new lanes than the old ones)
        for uvk, data in G.edges.items():
            width_before = data[width_attribute]
            width_after = lane_graph.calculate_street_width(L, *uvk)
            data['_after_width_total_m'] = width_after
            data['_after_excess_width_m'] = width_after - width_before

        # create a list of unfixed car edges, these are removal candidates
        removal_candidates_car = list(filter(
            lambda edge:
                (edge[1]['lane'].status != STATUS_FIXED or edge[1]['lane'].direction == DIRECTION_TBD)
                and edge[1]['lanetype'] in {LANETYPE_MOTORIZED, LANETYPE_HIGHWAY},
                #and G.edges[(edge[1]['u_G'], edge[1]['v_G'], edge[1]['k_G'])]['_after_excess_width_m'] > 0,
                L.edges.items()
        ))

        # stop here if no removal candidates exist
        if len(removal_candidates_car) == 0:
            break

        # sort removal candidates

        if mode == 'lowest_bc':
            removal_candidates_car = sorted(
                removal_candidates_car,
                key=lambda x: (
                    #G.edges[(x[1]['u_G'], x[1]['v_G'], x[1]['k_G'])]['_after_excess_width_m'] > 0,
                    x[1]['bc_private_cars']
                )
            )

        elif mode == 'highest_bc':
            removal_candidates_car = sorted(
                removal_candidates_car,
                key=lambda x: (
                    #G.edges[(x[1]['u_G'], x[1]['v_G'], x[1]['k_G'])]['_after_excess_width_m'] > 0,
                    -x[1]['bc_private_cars']
                )
            )

        else:
            raise NotImplementedError(f'Mode {mode} is not implemented')


        removal_candidate_car = removal_candidates_car[0]
        remove_edge_uvk = removal_candidate_car[0]
        remove_edge_data = removal_candidate_car[1]

        opposite_edge_uvk = (remove_edge_uvk[1], remove_edge_uvk[0], remove_edge_uvk[2])

        # make a subgraph containing only parking lanes
        P = L.edge_subgraph(
            [uvk for uvk, data in L.edges.items() if MODE_CAR_PARKING in data['lane'].get_modes()]
        )
        P = P.subgraph(
            [node for node, degree in dict(G.degree()).items() if degree > 0]
        )

        # set a flag whether each node must maintain access
        # (either due to required access or because of touching a parking lane)
        for i, data in L.nodes.items():
            data['_cars_and_parking'] = data.get('needs_access_by_private_cars', False) or i in P.nodes

        # checks
        # would the removal of this edge disconnect the graph?
        M = get_effective_subgraph(L, weight='cost_private_cars')
        edges_disconnect_graph = check_if_edges_disconnect_graph(
            M, [remove_edge_uvk], 'needs_access_by_private_cars'
        )

        # would the removal of this edge violate the public transit requirement?
        edge_breaks_transit = check_if_edge_breaks_transit(G, L, *remove_edge_uvk)

        # is the last direction of a mandatory lane with direction tbd?
        is_last_direction_of_mandatory_lane = \
            remove_edge_data['lane'].status == STATUS_FIXED and remove_edge_data['lane'].direction == DIRECTION_TBD \
            and not L.has_edge(*opposite_edge_uvk)

        # has dependent parking lanes?
        has_dependent_parking_lanes = check_if_has_dependent_parking_lanes(L, *remove_edge_uvk)

        if (
            edges_disconnect_graph
            or edge_breaks_transit
            or is_last_direction_of_mandatory_lane
            or has_dependent_parking_lanes
        ):

            # mark as fixed
            L.edges[remove_edge_uvk]['lane'].status = STATUS_FIXED

            # if direction tbd
            if L.edges[remove_edge_uvk]['lane'].direction == DIRECTION_TBD:
                # convert to forward and make twin factor 1
                L.edges[remove_edge_uvk]['lane'].direction = DIRECTION_FORWARD
                L.edges[remove_edge_uvk]['twin_factor'] = 1
                L.edges[remove_edge_uvk]['instance'] = 1

                # if there is an opposite instance that was already converted into a forward lane, assign a new key
                # (having two forward lanes with the same key causes trouble downstream)
                if opposite_edge_uvk in L.edges and L.edges[opposite_edge_uvk]['lane'].direction == DIRECTION_FORWARD:
                    remove_edge_uvk_old = copy.deepcopy(remove_edge_uvk)
                    # tuple -> list, so that we can do item assignment
                    remove_edge_uvk = list(remove_edge_uvk)
                    remove_edge_uvk[2] = remove_edge_uvk[2] + 'a'
                    # to change the key, we need to create a new edge and delete the old one
                    L.add_edge(*remove_edge_uvk, **L.edges[remove_edge_uvk_old])
                    print(remove_edge_uvk)
                    L.remove_edge(*remove_edge_uvk_old)

            # write the check results into the edge attributes for debugging
            L.edges[remove_edge_uvk]['_disconnects_graph'] = str(edges_disconnect_graph)
            L.edges[remove_edge_uvk]['_breaks_transit'] = str(edge_breaks_transit)
            L.edges[remove_edge_uvk]['_last_dir_of_mandatory'] = str(is_last_direction_of_mandatory_lane)
            L.edges[remove_edge_uvk]['_has_dependent_parking'] = str(has_dependent_parking_lanes)

            if verbose:
                print('fixed', remove_edge_uvk)
            uv_changed = None

        else:

            # remove
            if verbose:
                print('removed', remove_edge_uvk)

            # if direction is both and opposite instance exists, make the opposite instance forward and instance no. 1
            if (
                    L.edges[remove_edge_uvk]['lane'].direction == DIRECTION_BOTH
                    and L.has_edge(*opposite_edge_uvk)
            ):
                L.edges[opposite_edge_uvk]['lane'].direction = DIRECTION_FORWARD
                L.edges[opposite_edge_uvk]['lane'].set_default_width()
                L.edges[opposite_edge_uvk]['instance'] = 1
                # if the original lane has one direction mandatory, make the remaining lane fixed
                if L.edges[remove_edge_uvk]['lane'].status in [STATUS_FIXED, STATUS_ONE_DIRECTION_MANDATORY]:
                    L.edges[opposite_edge_uvk]['lane'].direction = DIRECTION_FORWARD
                    L.edges[opposite_edge_uvk]['lane'].set_default_width()
                    L.edges[opposite_edge_uvk]['lane'].status = STATUS_FIXED
                    L.edges[opposite_edge_uvk]['_last_dir_of_mandatory_from_bidir'] = str(True)

            L.remove_edge(*remove_edge_uvk)
            uv_changed = remove_edge_uvk[0:2]

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
        motorized_lanes = {uvk: data for uvk, data in lanes.items() if data['lane'].lanetype == LANETYPE_MOTORIZED}
        has_forward = False
        has_backward = False
        instance = 1
        for lane_u, lane_v, lane_k in motorized_lanes.keys():
            if lane_u == u and lane_v == v:
                has_forward = True
            if lane_u == v and lane_v == u:
                has_backward = True
        if has_forward and has_backward and len(motorized_lanes) == 2:
            # take the k of the first lane as the future id
            lane_id = list(motorized_lanes.keys())[0][2]
            # take the first lane as a template
            template_data = L.edges[list(motorized_lanes.keys())[0]]
            del template_data['instance']
            template_data['lane'] = space_allocation.Lane(LANETYPE_MOTORIZED, DIRECTION_BOTH)
            template_data['width'] = template_data['lane'].width
            template_data['direction'] = template_data['lane'].direction
            template_data['twin_factor'] = 0.5
            template_data['lane_id'] = lane_id
            for motorized_lane_uvk, motorized_lane_data in motorized_lanes.items():
                L.remove_edge(*motorized_lane_uvk)
                L.add_edge(
                    *motorized_lane_uvk[0:2], lane_id,
                    **copy.deepcopy(template_data), instance=instance
                )
                instance += 1


def _remove_parking(
        L, L_existing,
        G, width_attribute,
        A,
        gravity_iterations=15,
        lanetype=LANETYPE_PARKING_PARALLEL,
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
                and edge[1]['lanetype'] in {lanetype},
                #and G.edges[(edge[1]['u_G'], edge[1]['v_G'], edge[1]['k_G'])]['_after_excess_width_m'] > 0,
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

        if access_graph.effect_of_parking_removal_on_underassignment(
                A, [remove_edge_uvk], iterations=gravity_iterations
        ) == 0:
            L.remove_edge(*remove_edge_uvk)
            A.remove_node(remove_edge_uvk)
            access_graph.gravity_model(A, iterations=gravity_iterations)

            if verbose:
                print('removed', remove_edge_uvk)

        else:
            L.edges[remove_edge_uvk]['lane'].status = STATUS_FIXED

            if verbose:
                print('fixed', remove_edge_uvk)


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
        only_too_narrow=False,
        verbose=False
):
    """
    Extend the width of cycling lanes to fill out the remaining space
    TODO: Merge with space_allocation.adjust_widths_to_street()

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

        # check if criteria is met
        if only_too_narrow is True and excess_width < 0:
            pass
        elif only_too_narrow is False and excess_width != 0:
            pass
        else:
            continue

        lanes = lane_graph.get_street_lanes(L, *uvk)
        # filter
        if lanetypes:
            lanes = {uvk: data for uvk, data in lanes.items() if data['lane'].lanetype in lanetypes}
        total_width = sum([data['lane'].width for uvk, data in lanes.items()])

        for lane_uvk, lane_data in lanes.items():
            if verbose:
                print(uvk, -excess_width)
            if total_width == 0:
                lane_data['lane'].width = lane_data['lane'].width - excess_width / len(lanes)
            else:
                lane_data['lane'].width = lane_data['lane'].width - excess_width * lane_data['lane'].width / total_width
            lane_data['width'] = lane_data['lane'].width
            lane_data['lane'].status = STATUS_FIXED


"""
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
"""


def multi_rebuild_region(
        G,
        polygon,
        car_access_needs=None,
        bicycle_access_needs=None,
        hierarchies_to_include=hierarchy.HIERARCHIES,
        hierarchies_to_fix=(),
        add_fix_hierarchies=(hierarchy.PATHWAY, hierarchy.OTHER_HIERARCHY),
        motorized_traffic_on_all_streets=False,
        given_lanes_function=multi_set_given_lanes,
        remove_car_lanes_mode=None,
        motorized_traffic_lane_mode='separate_lanes',
        public_transit_mode='mandatory_like_existing',
        cycling_infrastructure_mode='optional',
        parking_mode='everywhere_by_need',
        car_parking_radius=200,
        bicycle_parking_radius=200,
        gravity_iterations=15,
        needed_node_access_function=multi_set_needed_node_access,
        given_lanes_attribute=KEY_GIVEN_LANES_DESCRIPTION,
        target_lanes_attribute=KEY_LANES_DESCRIPTION_AFTER,
        forced_given_lanes_attribute=KEY_FORCED_GIVEN_LANES_DESCRIPTION,
        fill_out_with_non_traffic=True,
        # export_L=None, export_H=None,
        # export_when='before_and_after',
        verbose=False,
        save_steps_path=None,
        save_steps_scaling_factor=1,
        **kwargs
):

    print('truncating street graph')
    H = oxc.truncate.truncate_graph_polygon(G, polygon, quadrat_width=100, retain_all=True)
    H = street_graph.filter_by_hierarchy(H, hierarchies_to_include)
    H = copy.deepcopy(H)
    print('done truncating')

    print('GIVEN LANES')

    # set given lanes and required access for nodes according to network rules
    given_lanes_function(
        H,
        source_lanes_attribute=target_lanes_attribute,
        target_lanes_attribute=given_lanes_attribute,
        forced_given_lanes_attribute=forced_given_lanes_attribute,
        hierarchies_to_fix=tuple(set(hierarchies_to_fix) | set(add_fix_hierarchies)),
        motorized_traffic_on_all_streets=motorized_traffic_on_all_streets,
        motorized_traffic_lane_mode=motorized_traffic_lane_mode,
        public_transit_mode=public_transit_mode,
        parking_mode=parking_mode,
        cycling_infrastructure_mode=cycling_infrastructure_mode,
        **kwargs
    )

    # for each node, set whether it must be accessible by each mode
    needed_node_access_function(H)

    # for i, data in H.nodes.items():
    #    print(i, 'needs_access_by_private_cars' in data.keys())

    L = lane_graph.create_lane_graph(H)
    if save_steps_path:
        io.export_HLA(save_steps_path, '0', L=L, scaling_factor=save_steps_scaling_factor)

    """
    if not hierarchy.HIERARCHIES <= hierarchies_to_include:
        # simplify the graph by removing intermediate nodes
        merge_edges.reset_intermediate_nodes(H)
        merge_edges.merge_consecutive_edges(H, distinction_attributes={given_lanes_attribute}, min_width=True)
    """

    L = lane_graph.create_lane_graph(H, lanes_attribute=given_lanes_attribute)

    if len(L.edges) == 0:
        print('no edges')
        return

    A = access_graph.create_access_graph(
        L, gpd.clip(car_access_needs, polygon), radius=car_parking_radius
    )
    B = access_graph.create_access_graph(
        L, gpd.clip(bicycle_access_needs, polygon), radius=bicycle_parking_radius,
        lanetype=LANETYPE_BICYCLE_PARKING, vehicle_length=BICYCLE_WIDTH_IN_RACK
    )
    access_graph.gravity_model(A, iterations=gravity_iterations)
    access_graph.gravity_model(B, iterations=gravity_iterations)

    # Given lanes
    if save_steps_path:
        io.export_HLA(save_steps_path, '1', H=H, L=L, A=A, B=B, scaling_factor=save_steps_scaling_factor)

    # Remove lanes in steps
    print('REMOVE CAR PARKING')
    _remove_parking(
        L, None, H, 'width', A, gravity_iterations=gravity_iterations, verbose=True
    )
    print('REMOVE BICYCLE PARKING')
    _remove_parking(
        L, None, H, 'width', B, gravity_iterations=gravity_iterations, verbose=True,
        lanetype=LANETYPE_BICYCLE_PARKING
    )

    # determine which nodes must be accessible by car, e.g., due to parking lanes
    M = copy.deepcopy(L)
    for uvk, data in L.edges.items():
        if not (
                data['lane'].lanetype == LANETYPE_PARKING_PARALLEL
                or (data['lane'].status in [STATUS_FIXED, STATUS_ONE_DIRECTION_MANDATORY] and data['lane'].lanetype == LANETYPE_MOTORIZED)
        ):
            M.remove_edge(*uvk)

    graph.remove_isolated_nodes(M)

    for i, data in L.nodes.items():
        data['needs_access_by_private_cars'] = M.has_node(i)

    if save_steps_path:
        io.export_HLA(save_steps_path, '2', L=L, A=A, B=B, scaling_factor=save_steps_scaling_factor)

    print('REMOVE CAR LANES')
    _remove_car_lanes(L, None, H, 'width', A, remove_car_lanes_mode, verbose=True)
    if save_steps_path:
        io.export_HLA(save_steps_path, '3', H=H, L=L, A=A, B=B, scaling_factor=save_steps_scaling_factor)
    #if motorized_traffic_lane_mode == 'bidirectional_lanes':
    #    print('NARROW DOWN CAR LANES')
    #    _remove_car_lanes(L, None, H, 'width', A, remove_car_lanes_mode, verbose=True)
    #    if save_steps_path:
    #        io.export_HLA(save_steps_path, '3a', H=H, L=L, A=A, scaling_factor=save_steps_scaling_factor)
    print('REMOVE CYCLING LANES')
    _remove_cycling_lanes(L, None, H, 'width', verbose=True)
    if save_steps_path:
        io.export_HLA(save_steps_path, '4', H=H, L=L, A=A, B=B, scaling_factor=save_steps_scaling_factor)
    print('MERGE TRANSIT AND CAR LANES')
    _merge_transit_with_car_lanes(L, None, H, 'width', verbose=True)
    if save_steps_path:
        io.export_HLA(save_steps_path, '5', H=H, L=L, A=A, B=B, scaling_factor=save_steps_scaling_factor)
    print('ADJUST WIDTH OF LANES')
    # at first, try to fill out with non-traffic
    if fill_out_with_non_traffic:
        _adjust_width_of_lanes(
            L, None, H, 'width', lanetypes={LANETYPE_NON_TRAFFIC},
            only_too_narrow=True, verbose=True
        )
    # second, try wo wide cycling lanes
    _adjust_width_of_lanes(
        L, None, H, 'width', lanetypes={LANETYPE_CYCLING_LANE},
        only_too_narrow=True, verbose=True
    )
    # if not possible, widen all other lanes proportionally
    _adjust_width_of_lanes(L, None, H, 'width', verbose=True)
    if save_steps_path:
        io.export_HLA(save_steps_path, '6', L=L, A=A, B=B, scaling_factor=save_steps_scaling_factor)

    for uvk, data in H.edges.items():
        lane_graph.merge_lanes_and_equalize_widths(
            L, *uvk,
            filter=lambda lane: lane['lane'].lanetype == LANETYPE_CYCLING_LANE
        )

    street_graph.lane_graph_to_street_graph(
        H, L,
        KEY_LANES_DESCRIPTION_AFTER
    )

    if save_steps_path:
        io.export_HLA(save_steps_path, '7', H=H, L=L, A=A, B=B, scaling_factor=save_steps_scaling_factor)

    """
    if not hierarchy.HIERARCHIES <= hierarchies_to_include:
        # reconstruct the original street graph with intermediary nodes
        merge_edges.reconstruct_consecutive_edges(H)
        street_graph.organize_edge_directions(H, key_lanes_description=KEY_LANES_DESCRIPTION_AFTER)
    """

    if save_steps_path:
        io.export_HLA(save_steps_path, '8', H=H, L=L, A=A, B=B, scaling_factor=save_steps_scaling_factor)

    I = copy.deepcopy(H)

    # remove edges that must remain fixed
    for uvk, data in H.edges.items():
        if data['hierarchy'] in hierarchies_to_fix:
            I.remove_edge(*uvk)

    # write rebuilt lanes from the subgraph into the main graph
    nx.set_edge_attributes(
        G,
        nx.get_edge_attributes(I, target_lanes_attribute),
        target_lanes_attribute
    )

    return H, L, A


def multi_rebuild_regions(
        G,
        rebuilding_regions_gdf,
        rebuilding_function=multi_rebuild_region,
        source_lanes_attribute=KEY_LANES_DESCRIPTION,
        target_lanes_attribute=KEY_LANES_DESCRIPTION_AFTER,
        **kwargs
):
    """
    Rebuilds the street space allocation by applying the given rebuilding regions, with their own settings.

    Parameters
    ----------
    G: street_graph.StreetGraph
    rebuilding_regions_gdf: gpd.GeoDataFrame
    rebuilding_function: function
        Here you can provide an own rebuilding function, as a plugin.
        Otherwise, the default one, based on betweenness centrality will be used.
    source_lanes_attribute: str
        In which attribute are the status quo lanes stored
    target_lanes_attribute: str
        Into which attribute should we write the rebuilt lanes
    kwargs: **kwargs
        additional arguments for the **multi_rebuild_region()** function
    """

    # initialize the target lanes attribute as a copy of the given lanes
    nx.set_edge_attributes(
        G,
        copy.deepcopy(nx.get_edge_attributes(G, source_lanes_attribute)),
        target_lanes_attribute
    )

    # H, L, A = multi_rebuild_region(G, rebuilding_regions_gdf['geometry'][0], access_needs)

    rebuilding_regions_gdf = rebuilding_regions_gdf.sort_values(by='order')

    HLAs = rebuilding_regions_gdf.apply(
        lambda region: rebuilding_function(
            G,
            region['geometry'],
            hierarchies_to_include=region['hierarchies_to_include'] if len(
                region['hierarchies_to_include']) > 0 else hierarchy.HIERARCHIES,
            hierarchies_to_fix=region['hierarchies_to_fix'],
            motorized_traffic_on_all_streets=region['motorized_traffic_on_all_streets'],
            parking_mode=region['parking_mode'],
            car_parking_radius=region['parking_radius'],
            remove_car_lanes_mode=region['remove_car_lanes_mode'],
            motorized_traffic_lane_mode=region['motorized_traffic_lane_mode'],
            hierarchies_with_cycling_lanes=region['hierarchies_with_cycling_lanes'],
            cycling_infrastructure_mode=region['cycling_infrastructure_mode'],
            **kwargs
        ),
        axis=1
    )

    return HLAs