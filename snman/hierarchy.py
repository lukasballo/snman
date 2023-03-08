# Street Hierarchy levels
HIGHWAY = '0_highway'
MAIN_ROAD = '1_main_road'
LOCAL_ROAD = '2_local_road'
DEAD_END = '3_dead_end'
PATHWAY = '4_path'
OTHER_HIERARCHY = '9_other'


def add_hierarchy(G, iterations=20):
    """
    Label all streets with hierarchy levels
    TODO: speed-up the dead-end detection algorithm and remove user-defined number of iterations

    Parameters
    ----------
    G : nx.MultiGraph
        street graph
    iterations : int
        number of iterations for dead end detection

    Returns
    -------
    None
    """

    _identify_dead_ends(G, iterations)

    for edge in G.edges(data=True, keys=True):
        _add_edge_hierarchy(edge)


def _add_edge_hierarchy(edge):
    """
    Label one edge with its hierarchy level

    Parameters
    ----------
    edge : tuple
        the complete edge tuple

    Returns
    -------
    None
    """

    edge_data = edge[3]

    edge_data['hierarchy'] = OTHER_HIERARCHY

    if edge_data.get('highway') in {'primary', 'primary_link', 'secondary', 'secondary_link', 'tertiary', 'tertiary_link'}:
        edge_data['hierarchy'] = MAIN_ROAD

    if edge_data.get('highway') in {'residential', 'living_street', 'unclassified', 'service'}:
        edge_data['hierarchy'] = LOCAL_ROAD

    if edge_data.get('highway') in {'path', 'footway', 'cycleway'}:
        edge_data['hierarchy'] = PATHWAY

    if edge_data.get('highway') in {'construction', 'track'}:
        edge_data['hierarchy'] = OTHER_HIERARCHY

    if edge_data.get('dead_end'):
        edge_data['hierarchy'] = DEAD_END

    if edge_data.get('highway') in {'motorway', 'motorway_link', 'trunk', 'trunk_link'}:
        edge_data['hierarchy'] = HIGHWAY


def _identify_dead_ends(G, iterations):
    """
    Parameters
    ----------
    G : nx.MultiGraph
        street graph
    iterations : int

    Returns
    -------
    None
    """

    # TODO: Auto detect how many iterations are necessary to get through the entire graph
    # TODO: Change the algorithm, it seems to be inefficient and delivers wring results in some cases
    for i in range(iterations):

        for node in G.nodes():

            # Get a list of all adjacent edges that are:
            # - not yet labeled as dead ends
            # - accessible for cars
            adjacent_non_dead_end_edges = []
            for edge in list(G.edges(nbunch=node, data=True, keys=True)):
                if (not edge[3].get('dead_end')) and edge[3].get('highway')\
                        not in ['path', 'footway', 'track']:
                    adjacent_non_dead_end_edges.append(edge)

            # Only 1 adjacent edge -> mark this edge as dead end
            if len(adjacent_non_dead_end_edges) == 1:
                edge = adjacent_non_dead_end_edges[0]
                edge[3]['dead_end'] = True
