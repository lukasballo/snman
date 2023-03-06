import networkx as nx
import matplotlib.pyplot as plt
from . import constants, hierarchy, utils, distribution
from . import osmnx_customized as oxc


def rebuild_regions(G, rebuilding_regions_gdf, initialize_ln_desc_after=True, **kwargs):

    if initialize_ln_desc_after:
        nx.set_edge_attributes(G, nx.get_edge_attributes(G, 'ln_desc'), 'ln_desc_after')

    for idx, data in rebuilding_regions_gdf[rebuilding_regions_gdf['active'] == True].iterrows():
        _rebuild_region(G, data['geometry'], data['hierarchies_to_include'], data['hierarchies_to_fix'], **kwargs)


def _rebuild_region(G, polygon, hierarchies_to_include, hierarchies_to_fix, **kwargs):

    # create a subgraph with only those edges that should be reorganized
    H = oxc.truncate.truncate_graph_polygon(G, polygon, quadrat_width=100, retain_all=True)
    if hierarchies_to_include is not []:
        filtered_edges = dict(filter(lambda key_value: key_value[1]['hierarchy']
            not in hierarchies_to_include, H.edges.items()))
        H.remove_edges_from(filtered_edges.keys())

    # initialize the input for link elimination
    H_minimal_graph_input = distribution.create_given_lanes_graph(H, hierarchies_to_fix=hierarchies_to_fix)
    #snman.export_streetgraph(H_minimal_graph_input, export_path + 'given_lanes.gpkg', export_path + 'given_lanes_nodes.gpkg')

    # run the link elimination
    H_minimal_graph_output = link_elimination(H_minimal_graph_input, **kwargs)
    #snman.export_streetgraph(H_minimal_graph_output, export_path + 'minimal_graph_out_edges.gpkg', export_path + 'minimal_graph_out_nodes.gpkg')

    # apply the link elimination output to the subgraph graph
    rebuild_lanes_from_owtop_graph(H, H_minimal_graph_output, hierarchies_to_protect=hierarchies_to_fix)

    # write the reorganized lanes from subgraph H into the main graph G
    nx.set_edge_attributes(G, nx.get_edge_attributes(H,'ln_desc_after'), 'ln_desc_after')

def link_elimination(O, keep_all_streets=True, verbose=False):

    # Get the giant weakly connected component (remove any unconnected parts)
    gcc = sorted(nx.weakly_connected_components(O), key=len, reverse=True)[0]
    O = O.subgraph(gcc).copy()

    # Add complementary edges
    for i, data in O.edges.items():
        if data.get('fixed', False) == False:
            O.add_edge(i[1], i[0], fixed=False)

    if verbose:
        print('Initialized graph has ', len(O.nodes), ' nodes and ', len(O.edges), ' edges')

    if not nx.is_strongly_connected(O):
        print('Initialized graph is not strongly connected')
        return

    def opposite_direction_exists(O, *edge_id):
        return O.has_edge(edge_id[1], edge_id[0])

    # Remove edges
    i = 0
    while True:
        i+=1
        if verbose:
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


