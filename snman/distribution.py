import networkx as nx
from . import constants, lanes, hierarchy


def set_given_lanes(G, bidirectional_for_dead_ends=True):
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
        data[constants.KEY_GIVEN_LANES_DESCRIPTION] = []

        # For roads with public transport, keep both directions
        if data.get('pt_tram') or data.get('pt_bus'):
            data[constants.KEY_GIVEN_LANES_DESCRIPTION] += [
                lanes.LANETYPE_MOTORIZED + lanes.DIRECTION_BACKWARD,
                lanes.LANETYPE_MOTORIZED + lanes.DIRECTION_FORWARD
            ]

        # For normal roads, keep one single-direction lane
        elif data.get('hierarchy') in [hierarchy.MAIN_ROAD, hierarchy.LOCAL_ROAD, hierarchy.HIGHWAY]:
            data[constants.KEY_GIVEN_LANES_DESCRIPTION] += [lanes.LANETYPE_MOTORIZED + lanes.DIRECTION_TBD]

        # For dead ends, create a single bi-directional lane
        elif data.get('hierarchy') == hierarchy.DEAD_END:
            if bidirectional_for_dead_ends:
                data[constants.KEY_GIVEN_LANES_DESCRIPTION] += [lanes.LANETYPE_MOTORIZED + lanes.DIRECTION_BOTH]
            else:
                data[constants.KEY_GIVEN_LANES_DESCRIPTION] += [lanes.LANETYPE_MOTORIZED + lanes.DIRECTION_TBD]


def create_given_lanes_graph(G, hierarchies_to_remove=[], hierarchies_to_fix=[]):
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
            lanes_list = data.get(constants.KEY_LANES_DESCRIPTION, [])
        else:
            # use the "given" necessary lanes
            lanes_list = data.get(constants.KEY_GIVEN_LANES_DESCRIPTION, [])

        for lane in lanes_list:
            lane_properties = lanes._lane_properties(lane)

            if data.get('hierarchy') in hierarchies_to_remove:
                continue

            if lane_properties.direction in [lanes.DIRECTION_FORWARD, lanes.DIRECTION_BOTH]:
                H.add_edge(u, v, fixed=True)

            if lane_properties.direction in [lanes.DIRECTION_BACKWARD, lanes.DIRECTION_BOTH]:
                H.add_edge(v, u, fixed=True)

            if lane_properties.direction in [lanes.DIRECTION_TBD]:
                H.add_edge(u, v, fixed=False)

    return H
