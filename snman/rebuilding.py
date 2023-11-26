import copy, math
import networkx as nx
import geopandas as gpd
from . import utils, distribution, space_allocation, hierarchy, street_graph, graph, io, merge_edges, lane_graph
from .constants import *
from . import osmnx_customized as oxc


def link_elimination(L, verbose=False):
    """
    **Deprecated, use multi_rebuild()**

    Generating a network fo one-way streets. A greedy algorithm that sequentially removes links from the graph
    until no link can be removed without losing strong connectivity.

    In literature, the problem is sometimes referred to as One-Way Traffic Organization problem (OWTOP).

    The solution in this function is inspired by a cycling network design algorithm in

    Steinacker, C., D.-M. Storch, M. Timme and M. Schröder (2022)
    Demand-driven design of bicycle infrastructure networks for improved urban bikeability,
    Nature Computational Science, DOI: https://doi.org/10.1038/s43588-022-00318-w.


    Parameters
    ----------
    L: nx.DiGraph
        lane graph of "given" car lanes based on the high-level design rules applied in distribution.set_given_lanes
    verbose : bool
        print internal details during the process

    Returns
    -------
    L : nx.DiGraph
        lane graph of car lanes after rebuilding
    """

    # Get the giant weakly connected component (remove any unconnected parts)
    gcc = sorted(nx.weakly_connected_components(L), key=len, reverse=True)[0]
    L = L.subgraph(gcc).copy()

    if verbose:
        print('Initialized graph has ', len(L.nodes), ' nodes and ', len(L.edges), ' edges')

    if not nx.is_strongly_connected(L):
        print('Initialized graph is not strongly connected')
        return

    def opposite_direction_exists(L, u, v, k):
        return L.has_edge(v, u, k)

    # Remove edges
    i = 0
    while True:
        i += 1
        if verbose and i % 10 == 0:
            print('Iteration ', i)

        # calculate betweenness centrality
        bc = nx.edge_betweenness_centrality(L, weight='cost_' + MODE_PRIVATE_CARS)
        nx.set_edge_attributes(L, bc, 'bc')
        removal_candidates = [edge for edge in list(L.edges.items()) if edge[1].get('fixed', False) == False]

        # stop here if no removal candidates exist
        if len(removal_candidates) == 0:
            break

        # find the edge with the lowest bc
        removal_candidates = sorted(removal_candidates, key=lambda x: x[1]['bc'])
        removal_candidate = list(removal_candidates)[0]
        remove_edge_id = removal_candidate[0]

        # remove the edge but add it back if the conditions get violated
        L.remove_edge(*remove_edge_id)
        if \
                not nx.is_strongly_connected(L) \
                        or (removal_candidate[1].get('mandatory_lane') and not opposite_direction_exists(L,
                                                                                                         *remove_edge_id)):
            attributes = copy.deepcopy(removal_candidate[1])
            attributes['fixed'] = True
            L.add_edge(*remove_edge_id, **attributes)
        else:
            pass

    return L


def rebuild_regions(
        G,
        rebuilding_regions_gdf,
        rebuilding_function=link_elimination,
        source_lanes_attribute=KEY_LANES_DESCRIPTION,
        target_lanes_attribute=KEY_LANES_DESCRIPTION_AFTER,
        verbose=False,
        export_L=None,
        export_H=None
):
    """
    ** Deprecated. Use multi_rebuild_regions() **

    Iterates through each one of the *rebuilding_regions* and applies the given *rebuilding_function*
    to the lane within its area and road hierarchies.

    Parameters
    ----------
    G : nx.MultiDiGraph
        street graph
    rebuilding_regions_gdf : gpd.GeoDataFrame
        regions to be rebuilt
    rebuilding_function : function
        a function to be applied to the lane graph within each region,
        see *link_elimination()* as example for the required inputs and outputs
    source_lanes_attribute : str
        key of the source lane description, where should the lanes be loaded from
    target_lanes_attribute : str
        key of the target lane description, where should the resulting lanes be saved
    verbose : bool
        show internal process information if true
    export_L : list
        a list of two strings, paths to save the edges and nodes of the intermediary lane graphs, use for debugging
    export_H : list
        a list of two strings, paths to save the edges and nodes of the intermediary street graphs, use for debugging

    Returns
    -------
    None
    """

    # initialize lanes after rebuild
    nx.set_edge_attributes(G, nx.get_edge_attributes(G, source_lanes_attribute), target_lanes_attribute)
    # ensure consistent edge directions
    street_graph.organize_edge_directions(G)

    for i, rebuilding_region in rebuilding_regions_gdf.iterrows():

        print('Rebuilding region', i, rebuilding_region['description'])

        hierarchies_to_fix = hierarchy.HIERARCHIES.difference(rebuilding_region['hierarchies_to_include']).union(
            rebuilding_region['hierarchies_to_fix'])

        polygon = rebuilding_region['geometry']

        # make a graph cutout based on the region geometry and skip this region if the resulting subgraph is empty
        H = oxc.truncate.truncate_graph_polygon(G, polygon, quadrat_width=100, retain_all=True)
        if len(H.edges) == 0:
            continue

        # keep only hierarchies to include
        if len(rebuilding_region['hierarchies_to_include']) > 0:
            H = street_graph.filter_by_hierarchy(H, rebuilding_region['hierarchies_to_include'])

        # set given lanes according to network rules
        distribution.set_given_lanes(
            H,
            maintain_motorized_access_on_street_level=rebuilding_region['keep_all_streets'],
            hierarchies_to_fix=rebuilding_region['hierarchies_to_fix']
        )

        # keep only the car lanes
        H = street_graph.filter_lanes_by_modes(H, {MODE_PRIVATE_CARS}, lane_description_key=KEY_GIVEN_LANES_DESCRIPTION)

        # simplify the graph by removing intermediate nodes
        merge_edges.reset_intermediate_nodes(H)
        merge_edges.merge_consecutive_edges(H, distinction_attributes={KEY_LANES_DESCRIPTION_AFTER})

        # make lane graph and ensure it is strongly connected
        L = lane_graph.create_lane_graph(H, KEY_GIVEN_LANES_DESCRIPTION)
        L = graph.keep_only_the_largest_connected_component(L)

        # export the lane graphs (before rebuilding) for debugging purposes
        if export_L:
            io.export_street_graph(L, export_L[0], export_L[1])
        if export_H:
            io.export_street_graph(H, export_H[0], export_H[1])

        # execute the rebuilding function
        L = rebuilding_function(L, verbose=verbose)

        # export the lane graphs (after rebuilding) for debugging purposes
        if export_L:
            io.export_street_graph(L, export_L[0], export_L[1])
        if export_H:
            io.export_street_graph(H, export_H[0], export_H[1])

        # use the resulting lane graph to rebuild the street graph
        rebuild_lanes_from_owtop_graph(
            H,
            L,
            hierarchies_to_protect=rebuilding_region['hierarchies_to_fix'],
            source_lanes_attribute=KEY_LANES_DESCRIPTION_AFTER,
            target_lanes_attribute=KEY_LANES_DESCRIPTION_AFTER
        )

        # reconstruct the original street graph before removing the intermediary nodes
        merge_edges.reconstruct_consecutive_edges(H)
        street_graph.organize_edge_directions(H)

        # export the street subgraph (after rebuilding) for debugging purposes
        if export_H:
            io.export_street_graph(H, export_H[0], export_H[1])

        # write rebuilt lanes from the subgraph into the main graph
        nx.set_edge_attributes(G, nx.get_edge_attributes(H, KEY_LANES_DESCRIPTION_AFTER), KEY_LANES_DESCRIPTION_AFTER)


def rebuild_lanes_from_owtop_graph(
        G,
        O,
        hierarchies_to_protect=[],
        source_lanes_attribute=KEY_LANES_DESCRIPTION,
        target_lanes_attribute=KEY_LANES_DESCRIPTION_AFTER
):
    """
    **Deprecated, use multi_rebuild()**

    Update lanes in the street graph to match the topology in the owtop graph

    Parameters
    ----------
    G : nx.MultiGraph
        street graph
    O : nx.DiGraph
        owtop graph
    hierarchies_to_protect : list
        which street hierarchies should not be changed
    source_lanes_attribute : str
        attribute holding the lanes that should be used as input
    target_lanes_attribute : str
        attribute holding the lanes that should be used as output

    Returns
    -------
    None
    """

    n_car_lanes = {}
    for id, data in G.edges.items():
        u, v, k = id
        n_car_lanes[(u, v)] = O.has_edge(u, v) * 1
        n_car_lanes[(v, u)] = O.has_edge(v, u) * 1

    # iterate over all streets
    for id, data in G.edges.items():
        u, v, k = id
        lanes_before = data[source_lanes_attribute]
        lanes_after = lanes_before.copy()

        # don't touch this street if its hierarchy is protected
        if data['hierarchy'] in hierarchies_to_protect:
            continue

        # iterate over all lanes
        for i, l in enumerate(lanes_before):

            # M- lanes
            if l == LANETYPE_MOTORIZED + DIRECTION_BOTH:

                # if there are connections in both directions
                if n_car_lanes[(u, v)] >= 1 and n_car_lanes[(v, u)] >= 1:
                    # keep the lane as it is
                    n_car_lanes[(u, v)] -= 1
                    n_car_lanes[(v, u)] -= 1

                # if there is only a forward connection
                elif n_car_lanes[(u, v)] >= 1 and n_car_lanes[(v, u)] == 0:
                    # turn it into [L<,M>]
                    lanes_after[i] = [
                        LANETYPE_CYCLING_LANE + DIRECTION_BACKWARD,
                        LANETYPE_MOTORIZED + DIRECTION_FORWARD,
                    ]
                    n_car_lanes[(u, v)] -= 1

                # if there is only a backward connection
                elif n_car_lanes[(u, v)] == 0 and n_car_lanes[(v, u)] >= 1:
                    # convert it into [M<,L>]
                    lanes_after[i] = [
                        LANETYPE_MOTORIZED + DIRECTION_BACKWARD,
                        LANETYPE_CYCLING_LANE + DIRECTION_FORWARD,
                    ]
                    n_car_lanes[(v, u)] -= 1

                # if there is no connection
                else:
                    # convert it into [P<,P>,P>]
                    lanes_after[i] = [
                        LANETYPE_CYCLING_TRACK + DIRECTION_BACKWARD,
                        LANETYPE_CYCLING_TRACK + DIRECTION_FORWARD,
                        LANETYPE_CYCLING_TRACK + DIRECTION_FORWARD,
                    ]

            # M< and M> lanes
            if l in [
                LANETYPE_MOTORIZED + DIRECTION_BACKWARD,
                LANETYPE_MOTORIZED + DIRECTION_FORWARD
            ]:

                # if forward connection exists (process forward first)
                if n_car_lanes[(v, u)] >= 1:
                    # convert into M<
                    lanes_after[i] = LANETYPE_MOTORIZED + DIRECTION_BACKWARD
                    n_car_lanes[(v, u)] -= 1

                # if backward connection exists
                elif n_car_lanes[(u, v)] >= 1:
                    # convert into M>
                    lanes_after[i] = LANETYPE_MOTORIZED + DIRECTION_FORWARD
                    n_car_lanes[(u, v)] -= 1

                else:
                    # convert it into [P<,P>]
                    lanes_after[i] = [
                        LANETYPE_CYCLING_TRACK + DIRECTION_BACKWARD,
                        LANETYPE_CYCLING_TRACK + DIRECTION_FORWARD,
                    ]

        # try to convert single-direction lanes to bidirectional ones
        for i, l in enumerate(lanes_before):

            # try to convert M< to M-
            if l == LANETYPE_MOTORIZED + DIRECTION_BACKWARD:
                if n_car_lanes[(v, u)] >= 1:
                    lanes_after[i] = str(
                        space_allocation._lane_properties(LANETYPE_MOTORIZED + DIRECTION_BOTH + str(
                            space_allocation._lane_properties(l).width
                        ))
                    )
                    n_car_lanes[(v, u)] -= 1

            # try to convert M> to M-
            if l == LANETYPE_MOTORIZED + DIRECTION_FORWARD:
                if n_car_lanes[(u, v)] >= 1:
                    lanes_after[i] = str(
                        space_allocation._lane_properties(LANETYPE_MOTORIZED + DIRECTION_BOTH + str(
                            space_allocation._lane_properties(l).width
                        ))
                    )
                    n_car_lanes[(u, v)] -= 1

        # add the remaining access directions as additional bidirectional lanes
        n_bidirectional = min([n_car_lanes[(v, u)], n_car_lanes[(u, v)]])
        lanes_after.extend([LANETYPE_MOTORIZED + DIRECTION_BOTH] * n_bidirectional)
        n_car_lanes[(v, u)] -= n_bidirectional
        n_car_lanes[(u, v)] -= n_bidirectional

        # ...and single-direction lanes
        lanes_after.extend([LANETYPE_MOTORIZED + DIRECTION_BACKWARD] * n_car_lanes[(v, u)])
        n_car_lanes[(v, u)] = 0
        lanes_after.extend([LANETYPE_MOTORIZED + DIRECTION_FORWARD] * n_car_lanes[(u, v)])
        n_car_lanes[(u, v)] = 0

        lanes_after = list(utils.flatten_list(lanes_after))
        data[target_lanes_attribute] = lanes_after


def multi_set_needed_node_access(G, source_lanes_attribute=KEY_LANES_DESCRIPTION):
    """
    Assigns a set of flags to every node defining whether it must remain accessible to each mode.
    This is a built-in default function that can be replaced by a custom one in the rebuilding process.

    Parameters
    ----------
    G
    source_lanes_attribute

    Returns
    -------

    """
    for mode in MODES:
        H = street_graph.filter_lanes_by_modes(G, {mode}, source_lanes_attribute)
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
        hierarchies_to_fix=(hierarchy.PATHWAY,)
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
        target_lanes = []

        lane_stats = space_allocation._lane_stats(source_lanes)

        # for streets with "hierarchy to fix" keep everything as it is
        if data.get('hierarchy') in hierarchies_to_fix:
            target_lanes = source_lanes

        else:

            # -- Add dedicated public transit lanes --

            # all_dedicated: add optional pt lane to every street with pt
            if public_transit_mode == 'all_dedicated':
                if data.get('pt_backward'):
                    target_lanes += [LANETYPE_DEDICATED_PT + DIRECTION_BACKWARD_OPTIONAL]
                if data.get('pt_forward'):
                    target_lanes += [LANETYPE_DEDICATED_PT + DIRECTION_FORWARD_OPTIONAL]
            # mandatory_like_existing: all existing pt lanes remain as they are
            elif public_transit_mode == 'mandatory_like_existing':
                existing_pt_lanes = space_allocation.filter_lanes_by_modes(source_lanes, {MODE_TRANSIT}, exact=True)
                target_lanes += existing_pt_lanes
            # none: no pt lanes
            elif public_transit_mode == 'none':
                pass
            else:
                print('pt_mode unknown:', public_transit_mode)
                return

            # -- Ensure connectivity for public transit --

            # save target lanes so far back to the street graph
            data[target_lanes_attribute] = target_lanes
            # on every edge with public transit, add a motorized traffic lane if no connection is possible
            pt_cost_forward = street_graph.calculate_edge_cost(
                G, *uvk, DIRECTION_FORWARD, MODE_TRANSIT, lanes_description=target_lanes_attribute
            )
            if data.get('pt_forward') and pt_cost_forward == math.inf:
                target_lanes += [LANETYPE_MOTORIZED + DIRECTION_FORWARD]
            pt_cost_backward = street_graph.calculate_edge_cost(
                G, *uvk, DIRECTION_BACKWARD, MODE_TRANSIT, lanes_description=target_lanes_attribute
            )
            if data.get('pt_backward') and pt_cost_backward == math.inf:
                target_lanes += [LANETYPE_MOTORIZED + DIRECTION_BACKWARD]

            # -- Add car lanes --

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
                        target_lanes += [LANETYPE_MOTORIZED + DIRECTION_TBD]
            # motorized_traffic_on_all_streets=False: add optional lane to every street
            else:
                if data.get('hierarchy') in motorized_traffic_road_hierarchies:
                    if math.inf in [car_cost_backward, car_cost_forward]:
                        target_lanes += [LANETYPE_MOTORIZED + DIRECTION_TBD_OPTIONAL]

            # -- Add cycling lanes --


            if data.get('hierarchy') in cycling_infra_road_hierarchies:
                if data.get('hierarchy') == hierarchy.MAIN_ROAD:
                    # twice the standard width in each direction on main roads
                    factor = 2
                else:
                    # standard width on other roads
                    factor = 1
                target_lanes += [LANETYPE_CYCLING_LANE + DIRECTION_BACKWARD_OPTIONAL] * factor
                target_lanes += [LANETYPE_CYCLING_LANE + DIRECTION_FORWARD_OPTIONAL] * factor

            # -- Add parking --

            if parking_mode == 'none':
                pass
            else:
                if parking_mode == 'mandatory_like_existing':
                    direction = DIRECTION_BOTH
                elif parking_mode == 'optional_like_existing':
                    direction = DIRECTION_BOTH_OPTIONAL
                else:
                    print('parking mode not implemented:', parking_mode)
                    return
                existing_parking_lanes = space_allocation.filter_lanes_by_modes(
                    source_lanes, {MODE_CAR_PARKING}, exact=True
                )
                for lane in existing_parking_lanes:
                    lp = space_allocation._lane_properties(lane)
                    lp.direction = direction
                    target_lanes += [str(lp)]

            # -- Add non-traffic space --

            target_lanes += [LANETYPE_NON_TRAFFIC + DIRECTION_BACKWARD, LANETYPE_NON_TRAFFIC + DIRECTION_FORWARD]

        data[target_lanes_attribute] = target_lanes


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

    # read from L and delete from M
    M = copy.deepcopy(L)

    # exclude edges based on their weight
    for uvk, data in L.edges.items():
        weight_value = data.get(weight)
        if weight_value == math.inf:
            M.remove_edge(uvk[0], uvk[1], key=uvk[2])

    # io.export_street_graph(M, export_path + 'L_edges.gpkg', export_path + 'L_nodes.gpkg')

    N = copy.deepcopy(M)

    # exclude nodes based on a given attribute
    for i, data in N.nodes.items():
        degree = N.degree(i)
        to_be_included = data.get(node_inclusion, False)
        if degree == 0 and to_be_included == False:
            M.remove_node(i)

    # exclude edges based on an optional list
    for uvk in exclude_edges:
        M.remove_edge(*uvk)

    return nx.is_strongly_connected(M)


def _remove_car_lanes(
        L, L_existing,
        G, width_attribute,
        verbose
):
    """
    a helper for multi_rebuild(), takes care of the car lanes removal
    """

    i = 1
    while True:

        if verbose:
            print('iteration', i)
        i += 1

        # calculate betweenness centrality
        bc = nx.edge_betweenness_centrality(L, weight='cost_' + MODE_PRIVATE_CARS)
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
            edge[1].get('fixed', False) == False
            and edge[1]['lanetype'] in {LANETYPE_MOTORIZED, LANETYPE_HIGHWAY},
            # and G.edges[(edge[1]['u_G'],edge[1]['v_G'],edge[1]['k_G'])]['_after_excess_width_m'] > 0,
            L.edges.items()
        ))

        # stop here if no removal candidates exist
        if len(removal_candidates_car) == 0:
            break

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

        remove_edge_uvks_to_test = [remove_edge_uvk]
        opposite_edge_uvk = (remove_edge_uvk[1], remove_edge_uvk[0], remove_edge_uvk[2])

        # add connected edges that would need to be removed as well
        if remove_edge_data['coupled_with_opposite_direction']:
            remove_edge_uvks_to_test.append(opposite_edge_uvk)

        # is still strongly connected?
        is_strongly_connected = is_strongly_connected_plus(
            L, 'cost_private_cars', 'needs_access_by_private_cars',
            exclude_edges=remove_edge_uvks_to_test
        )

        # is the last direction of a mandatory lane with direction tbd?
        is_last_direction_of_mandatory_lane = \
            remove_edge_data['mandatory_lane'] == True \
            and not L.has_edge(*opposite_edge_uvk)

        if is_strongly_connected and not is_last_direction_of_mandatory_lane:

            L.remove_edge(*remove_edge_uvk)
            if verbose:
                print('removed', remove_edge_uvk)

        else:
            if verbose:
                print(remove_edge_uvk, is_strongly_connected, is_last_direction_of_mandatory_lane)
            L.edges[remove_edge_uvk]['fixed'] = True
            L.edges[remove_edge_uvk]['twin_factor'] = 1
            if verbose:
                print('fixed', remove_edge_uvk)


def _remove_parking(
        L, L_existing,
        G, width_attribute,
        verbose
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
            edge[1].get('fixed', False) == False
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
        opposite_edge_uvk = (remove_edge_uvk[1], remove_edge_uvk[0], remove_edge_uvk[2])

        L.remove_edge(*remove_edge_uvk)
        L.remove_edge(*opposite_edge_uvk)


def _remove_cycling_lanes(
        L, L_existing,
        G, width_attribute,
        verbose
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
        removal_candidates_cycling = list(filter(
            lambda edge:
            edge[1].get('fixed', False) == False
            and edge[1]['lanetype'] in {LANETYPE_CYCLING_LANE}
            and G.edges[(edge[1]['u_G'], edge[1]['v_G'], edge[1]['k_G'])]['_after_excess_width_m'] > 0,
            L.edges.items()
        ))

        # stop here if no removal candidates exist
        if len(removal_candidates_cycling) == 0:
            break

        # sort by excess width
        removal_candidates_cycling = sorted(
            removal_candidates_cycling,
            key=lambda x: G.edges[(x[1]['u_G'], x[1]['v_G'], x[1]['k_G'])]['_after_excess_width_m'],
        )

        removal_candidate_cycling = removal_candidates_cycling[-1]
        remove_edge_uvk = removal_candidate_cycling[0]

        # is still strongly connected?
        sc = is_strongly_connected_plus(
            L, 'cost_cycling', 'needs_access_by_cycling',
            exclude_edges={remove_edge_uvk}
        )

        if sc:
            L.remove_edge(*remove_edge_uvk)
            if verbose:
                print('removed', remove_edge_uvk)
        else:
            L.edges[remove_edge_uvk]['fixed'] = True
            if verbose:
                print('fixed', remove_edge_uvk)


def _merge_transit_with_car_lanes(
        L, L_existing,
        G, width_attribute,
        verbose
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
            edge[1].get('fixed', False) == False
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
            key=lambda x: x[1].get('fixed', False) * 1
        ))
        if len(parallel_car_lanes) > 0:
            parallel_car_lane = parallel_car_lanes[-1]
            parallel_car_lane_uvk = parallel_car_lane[0]
            L.edges[parallel_car_lane_uvk]['fixed'] = True
            L.remove_edge(*remove_edge_uvk)
            if verbose:
                print('removed', remove_edge_uvk)
        else:
            L.edges[remove_edge_uvk]['fixed'] = True
            if verbose:
                print('fixed', remove_edge_uvk)

    return L


def _adjust_non_traffic(
        L, L_existing,
        G, width_attribute,
        verbose
):
    for uvk, data in G.edges.items():

        width_before = data[width_attribute]
        width_after = lane_graph.calculate_street_width(L, *uvk)
        excess_width = width_after - width_before
        data['_after_width_total_m'] = width_after
        data['_after_excess_width_m'] = excess_width

        if excess_width < 0:
            lanes = lane_graph.get_street_lanes(L, *uvk)
            # filter
            lanes = {uvk: data for uvk, data in lanes.items() if data['lanetype'] == LANETYPE_NON_TRAFFIC}
            for lane_uvk, lane_data in lanes.items():
                if verbose:
                    print(uvk, -excess_width)
                lane_data['width'] = -excess_width / len(lanes)


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
                lp = space_allocation._lane_properties(L_data['lane'])
                # for bidirectional lanes, work only with the first instance
                if lp.direction in {DIRECTION_BOTH, DIRECTION_BOTH_OPTIONAL}:
                    if L_data['instance'] == 1:
                        lp.direction = DIRECTION_BOTH
                        lp.width = L_data['width']
                        lanes_description += [str(lp)]
                # for all other lanes
                else:
                    lp.direction = direction
                    lp.width = L_data['width']
                    lanes_description += [str(lp)]
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
        export_when=None,
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
        H = oxc.truncate.truncate_graph_polygon(G, rebuilding_region.geometry, quadrat_width=100, retain_all=True)
        if len(H.edges) == 0:
            continue

        # keep only hierarchies to include
        H = street_graph.filter_by_hierarchy(H, hierarchies_to_include)

        # set given lanes and required access for nodes according to network rules
        given_lanes_function(
            H,
            hierarchies_to_fix=hierarchies_to_fix,
            motorized_traffic_on_all_streets=rebuilding_region['keep_all_streets'],
            public_transit_mode=public_transit_mode,
            parking_mode=parking_mode
        )
        needed_node_access_function(H)

        # simplify the graph by removing intermediate nodes
        merge_edges.reset_intermediate_nodes(H)
        merge_edges.merge_consecutive_edges(H, distinction_attributes={KEY_LANES_DESCRIPTION_AFTER})

        # make lane graph and ensure it is strongly connected
        L = lane_graph.create_lane_graph(H, KEY_GIVEN_LANES_DESCRIPTION)
        L = graph.keep_only_the_largest_connected_component(L)

        # export the lane graphs (before rebuilding) for debugging purposes
        if export_when in [None, 'before']:
            if export_L:
                io.export_street_graph(L, export_L[0], export_L[1])
            if export_H:
                io.export_street_graph(H, export_H[0], export_H[1])

        # execute the multi rebuilding function
        rebuilding_function(L, None, H, width_attribute, verbose=verbose)

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
        if export_when in [None, 'after']:
            if export_L:
                io.export_street_graph(L, *export_L)
            if export_H:
                io.export_street_graph(H, *export_H)

        # write rebuilt lanes from the subgraph into the main graph
        nx.set_edge_attributes(G, nx.get_edge_attributes(H, KEY_LANES_DESCRIPTION_AFTER), KEY_LANES_DESCRIPTION_AFTER)
