import networkx as nx
import matplotlib.pyplot as plt
from . import constants, hierarchy, utils


def link_elimination(O, keep_all_streets=True):

    # Get the giant weakly connected component (remove any unconnected parts)
    gcc = sorted(nx.weakly_connected_components(O), key=len, reverse=True)[0]
    O = O.subgraph(gcc).copy()

    # Add complementary edges
    for i, data in O.edges.items():
        if data.get('fixed', False) == False:
            O.add_edge(i[1], i[0], fixed=False)

    print('Initialized graph has ', len(O.nodes), ' and ', len(O.edges), ' edges')

    if not nx.is_strongly_connected(O):
        print('Initialized graph is not strongly connected')
        return

    def opposite_direction_exists(O, *edge_id):
        return O.has_edge(edge_id[1], edge_id[0])

    # Remove edges
    i = 0
    while True:
        i+=1
        print('Iteration ', i)
        # Calculate betweenness centrality
        bc = nx.edge_betweenness_centrality(O)
        nx.set_edge_attributes(O, bc, 'bc')
        edges = [edge for edge in list(O.edges.items()) if edge[1].get('fixed', False) == False]

        # Finish the loop if no edge candidates exist
        if len(edges) == 0:
            break
        edges = sorted(edges, key=lambda x: x[1]['bc'])
        edge = list(edges)[0]
        edge_id = edge[0]

        # Check if the new graph is strongly connected
        H = O.copy()
        H.remove_edge(*edge_id)
        if nx.is_strongly_connected(H) and (opposite_direction_exists(O, *edge_id) or not keep_all_streets):
            O.remove_edge(*edge_id)
        else:
            nx.set_edge_attributes(O, {edge_id: {'fixed': True}})

    return O


def rebuild_lanes_from_owtop_graph(G, O, hierarchies_to_protect=[]):
    """
    Rebuild lanes in the street graph according to the topology of the OWTOP graph
    :param G: Street Graph
    :param O: OWTOP graph
    :return:
    """

    n_car_lanes = {}
    for id, data in G.edges.items():
        u, v, k = id
        n_car_lanes[(u, v)] = O.has_edge(u, v) * 1
        n_car_lanes[(v, u)] = O.has_edge(v, u) * 1

    for id, data in G.edges.items():
        u, v, k = id
        lanes_before = data['ln_desc']
        lanes_after = lanes_before.copy()

        if data['hierarchy'] not in hierarchies_to_protect:

            for i, l in enumerate(lanes_before):

                # M- lanes
                if l == constants.LANETYPE_MOTORIZED + constants.DIRECTION_BOTH:

                    if n_car_lanes[(u, v)] >= 1 and n_car_lanes[(v, u)] >= 1:
                        # keep the lane as it is
                        n_car_lanes[(u, v)] -= 1
                        n_car_lanes[(v, u)] -= 1

                    elif n_car_lanes[(u, v)] >= 1 and n_car_lanes[(v, u)] == 0:
                        # turn it into [L<,M>]
                        lanes_after[i] = [
                            constants.LANETYPE_CYCLING_LANE + constants.DIRECTION_BACKWARD,
                            constants.LANETYPE_MOTORIZED + constants.DIRECTION_FORWARD,
                        ]
                        n_car_lanes[(u, v)] -= 1

                    elif n_car_lanes[(u, v)] == 0 and n_car_lanes[(v, u)] >= 1:
                        # convert it into [M<,L>]
                        lanes_after[i] = [
                            constants.LANETYPE_MOTORIZED + constants.DIRECTION_BACKWARD,
                            constants.LANETYPE_CYCLING_LANE + constants.DIRECTION_FORWARD,
                        ]
                        n_car_lanes[(v, u)] -= 1

                    else:
                        # convert it into [P<,P>]
                        lanes_after[i] = [
                            constants.LANETYPE_CYCLING_TRACK + constants.DIRECTION_BACKWARD,
                            constants.LANETYPE_CYCLING_TRACK + constants.DIRECTION_FORWARD,
                        ]

                # M< and M> lanes
                if l in [
                    constants.LANETYPE_MOTORIZED + constants.DIRECTION_BACKWARD,
                    constants.LANETYPE_MOTORIZED + constants.DIRECTION_FORWARD
                ]:

                    if n_car_lanes[(v, u)] >= 1:
                        # convert into M<
                        lanes_after[i] = constants.LANETYPE_MOTORIZED + constants.DIRECTION_BACKWARD
                        n_car_lanes[(v, u)] -= 1

                    elif n_car_lanes[(u, v)] >= 1:
                        # convert into M>
                        lanes_after[i] = constants.LANETYPE_MOTORIZED + constants.DIRECTION_FORWARD
                        n_car_lanes[(u, v)] -= 1

                    else:
                        # convert it into [P<,P>]
                        lanes_after[i] = [
                            constants.LANETYPE_CYCLING_TRACK + constants.DIRECTION_BACKWARD,
                            constants.LANETYPE_CYCLING_TRACK + constants.DIRECTION_FORWARD,
                        ]

        data['ln_desc_after'] = list(utils.flatten(lanes_after))


