import copy

import networkx as nx
from . import utils, distribution, space_allocation, hierarchy, street_graph, graph_utils, io
from .constants import *
from . import osmnx_customized as oxc


def rebuild_regions(
        G,
        rebuilding_regions_gdf,
        source_lanes_attribute=KEY_LANES_DESCRIPTION,
        target_lanes_attribute=KEY_LANES_DESCRIPTION_AFTER,
        verbose=False,
        export_L=False
):

    # initialize lanes after rebuild
    nx.set_edge_attributes(G, nx.get_edge_attributes(G, source_lanes_attribute), target_lanes_attribute)

    for i, rebuilding_region in rebuilding_regions_gdf.iterrows():

        print('Rebuilding region', i, rebuilding_region['description'])

        hierarchies_to_fix = hierarchy.HIERARCHIES.difference(rebuilding_region['hierarchies_to_include']).union(
            rebuilding_region['hierarchies_to_fix'])

        polygon = rebuilding_region['geometry']
        # make subgraph
        H = oxc.truncate.truncate_graph_polygon(G, polygon, quadrat_width=100, retain_all=True)
        # keep only hierarchies toi include
        if len(rebuilding_region['hierarchies_to_include']) > 0:
            H = street_graph.filter_by_hierarchy(H, rebuilding_region['hierarchies_to_include'])
        # set given lanes according to network rules
        distribution.set_given_lanes(
            H,
            maintain_motorized_access_on_street_level=rebuilding_region['keep_all_streets'],
            hierarchies_to_fix=rebuilding_region['hierarchies_to_fix']
        )
        # get only the car lanes
        H = street_graph.filter_lanes_by_modes(H, {MODE_PRIVATE_CARS})
        # make lane graph
        L = street_graph.to_lane_graph(H, KEY_GIVEN_LANES_DESCRIPTION)
        # make sure that the graph is strongly connected
        L = graph_utils.keep_only_the_largest_connected_component(L)

        if export_L:
            io.export_street_graph(L, export_L[0], export_L[1])

        # link elimination
        L = link_elimination(L, verbose=verbose)

        if export_L:
            io.export_street_graph(L, export_L[0], export_L[1])

        # rebuild lanes in subgraph based on new lane graph
        rebuild_lanes_from_owtop_graph(
            H,
            L,
            hierarchies_to_protect=rebuilding_region['hierarchies_to_fix'],
            source_lanes_attribute=KEY_LANES_DESCRIPTION_AFTER,
            target_lanes_attribute=KEY_LANES_DESCRIPTION_AFTER
        )
        # write rebuilt lanes from subgraph into the main graph
        nx.set_edge_attributes(G, nx.get_edge_attributes(H, KEY_LANES_DESCRIPTION_AFTER), KEY_LANES_DESCRIPTION_AFTER)


def link_elimination(L, verbose=False):
    """
    Generating a network fo one-way streets. A greedy algorithm that sequentially removes links from the graph
    until no link can be removed without losing strong connectivity.

    The problem is referred to in the literature as One-Way Traffic Organization problem (OWTOP).

    Parameters
    ----------
    L: nx.DiGraph
        lane graph
    verbose : bool
        print internal details during the process

    Returns
    -------
    L : nx.DiGraph
        a copy of the graph after link elimination
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
        bc = nx.edge_betweenness_centrality(L, weight='cost_'+MODE_PRIVATE_CARS)
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
            not nx.is_strongly_connected(L)\
            or (removal_candidate[1].get('mandatory_lane') and not opposite_direction_exists(L, *remove_edge_id)):
            attributes = copy.deepcopy(removal_candidate[1])
            attributes['fixed'] = True
            L.add_edge(*remove_edge_id, **attributes)
        else:
            pass


    return L


def rebuild_lanes_from_owtop_graph(
        G,
        O,
        hierarchies_to_protect=[],
        source_lanes_attribute=KEY_LANES_DESCRIPTION,
        target_lanes_attribute=KEY_LANES_DESCRIPTION_AFTER
):
    """
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
