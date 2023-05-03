import networkx as nx
from . import constants, lane_config, hierarchy
from .constants import *


def set_given_lanes(
        G,
        source_lanes_attribute=constants.KEY_LANES_DESCRIPTION,
        target_lanes_attribute=constants.KEY_GIVEN_LANES_DESCRIPTION,
        bidirectional_for_dead_ends=True
):
    """
    Sets which lanes are given due to external policy definitions
    e.g. dedicated lanes for public transport, bidirectional lanes for cars, etc.

    Parameters
    ----------
    G : nx.MultiGraph
        street graph
    bidirectional_for_dead_ends : bool
        automatically enforce bidirectional streets for all dead ends (speeds up the link elimination process but
        can lead to errors if the dead end detection is buggy)

    Returns
    -------
    None
    """

    # TODO: Add support for dedicated transit lanes

    for id, data in G.edges.items():

        source_lanes = data[source_lanes_attribute]
        target_lanes = []

        lane_stats = lane_config._lane_stats(source_lanes)

        # for roads with public transport...
        if data.get('pt_tram') or data.get('pt_bus'):
            # ...keep forward lane if there is no dedicated forward/both ways lane
            if lane_stats.n_lanes_dedicated_pt_forward + lane_stats.n_lanes_dedicated_pt_both_ways == 0:
                target_lanes += [LANETYPE_MOTORIZED + DIRECTION_FORWARD]
            # ...keep backward lane if there is no dedicated backward/both ways lane
            if lane_stats.n_lanes_dedicated_pt_backward + lane_stats.n_lanes_dedicated_pt_both_ways == 0:
                target_lanes += [LANETYPE_MOTORIZED + DIRECTION_BACKWARD]

        # for normal roads, keep one single-direction lane
        elif data.get('hierarchy') in [hierarchy.MAIN_ROAD, hierarchy.LOCAL_ROAD, hierarchy.HIGHWAY]:
            target_lanes += [LANETYPE_MOTORIZED + DIRECTION_TBD]

        # for dead ends, create a single bi-directional lane
        elif data.get('hierarchy') == hierarchy.DEAD_END:
            if bidirectional_for_dead_ends:
                target_lanes += [LANETYPE_MOTORIZED + DIRECTION_BOTH]
            else:
                target_lanes += [LANETYPE_MOTORIZED + DIRECTION_TBD]

        data[target_lanes_attribute] = target_lanes


def create_given_lanes_graph(
        G,
        hierarchies_to_remove=[],
        hierarchies_to_fix=[],
        source_lanes_attribute=constants.KEY_LANES_DESCRIPTION
    ):
    """
    Returns a directed graph of given (mandatory) lanes. Lanes with changeable direction are marked with an attribute

    Parameters
    ----------
    G : nx.MultiGraph
        street graph
    hierarchies_to_remove : list
        streets of which hierarchies should be removed from the graph (i.e. not considered in the further process)
    hierarchies_to_fix : list
        streets of which hierarchies should be automatically converted into fixed edges (i.e. considered but not
        changed in the further process)
    source_lanes_attribute : str
        attribute holding the lanes that should be used as a starting point for edges from the hierarchies to fix

    Returns
    -------
    H : nx.DiGraph
        a graph of given lanes to be used in the rebuilding process
    """
    H = nx.DiGraph()
    H.graph['crs'] = G.graph['crs']
    H.add_nodes_from(G.nodes.items())

    for id, data in G.edges.items():
        u = id[0]
        v = id[1]

        if data.get('hierarchy') in hierarchies_to_fix:
            # use the existing lanes
            lanes_list = data.get(source_lanes_attribute, [])
        else:
            # use the "given" necessary lanes
            lanes_list = data.get(constants.KEY_GIVEN_LANES_DESCRIPTION, [])

        for lane in lanes_list:
            lane_properties = lane_config._lane_properties(lane)

            if data.get('hierarchy') in hierarchies_to_remove:
                continue

            if lane_properties.direction in [DIRECTION_FORWARD, DIRECTION_BOTH]:
                H.add_edge(u, v, fixed=True)

            if lane_properties.direction in [DIRECTION_BACKWARD, DIRECTION_BOTH]:
                H.add_edge(v, u, fixed=True)

            if lane_properties.direction in [DIRECTION_TBD]:
                H.add_edge(u, v, fixed=False)

    return H
