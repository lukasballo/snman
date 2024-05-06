# Street Hierarchy levels
HIGHWAY = '0_highway'
MAIN_ROAD = '1_main_road'
LOCAL_ROAD = '2_local_road'
DEAD_END = '3_dead_end'
PATHWAY = '4_path'
SERVICE = '5_service'
OTHER_HIERARCHY = '9_other'

HIERARCHIES = {HIGHWAY, MAIN_ROAD, LOCAL_ROAD, DEAD_END, PATHWAY, OTHER_HIERARCHY, SERVICE}

HIGHWAY_OSM = {'motorway', 'motorway_link', 'trunk', 'trunk_link'}

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

    if edge_data.get('highway') in {'primary', 'primary_link', 'secondary', 'secondary_link', 'tertiary', 'tertiary_link'}:
        edge_data['hierarchy'] = MAIN_ROAD
    elif edge_data.get('highway') in {'residential', 'living_street', 'unclassified'}:
        edge_data['hierarchy'] = LOCAL_ROAD
    elif edge_data.get('highway') in {'path', 'footway', 'cycleway'}:
        edge_data['hierarchy'] = PATHWAY
    elif edge_data.get('highway') in {'construction', 'track'}:
        edge_data['hierarchy'] = OTHER_HIERARCHY
    elif edge_data.get('dead_end'):
        edge_data['hierarchy'] = DEAD_END
    elif edge_data.get('highway') in HIGHWAY_OSM:
        edge_data['hierarchy'] = HIGHWAY
    elif edge_data.get('highway') == 'service':
        edge_data['hierarchy'] = SERVICE
    else:
        edge_data['hierarchy'] = OTHER_HIERARCHY

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

    # TODO: does not work properly!
    return

    for i in range(iterations):

        for node in G.nodes():

            # Get a list of all adjacent edges that are:
            # - not yet labeled as dead ends
            # - accessible for cars
            adjacent_non_dead_end_edges = []
            for edge in list(G.in_edges(nbunch=node, data=True, keys=True)) + list(G.out_edges(nbunch=node, data=True, keys=True)):
                if (not edge[3].get('dead_end')) and edge[3].get('highway')\
                        not in ['path', 'footway', 'track']:
                    adjacent_non_dead_end_edges.append(edge)

            # Only 1 adjacent edge -> mark this edge as dead end
            if len(adjacent_non_dead_end_edges) == 1:
                edge = adjacent_non_dead_end_edges[0]
                edge[3]['dead_end'] = True
