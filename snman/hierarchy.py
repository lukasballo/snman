from . import graph_tools


def add_hierarchy(street_graph):

    _identify_dead_ends(street_graph)

    for edge in street_graph.edges(data=True, keys=True):
        _add_edge_hierarchy(street_graph, edge)


def _add_edge_hierarchy(street_graph, edge):

    edge_data = edge[3]

    edge_data['hierarchy'] = 'other'

    if edge_data.get('dead_end'):
        edge_data['hierarchy'] = 'dead_end'

    if edge_data.get('highway') in {'motorway', 'motorway_link', 'trunk', 'trunk_link'}:
        edge_data['hierarchy'] = 'highway'

    if edge_data.get('highway') in {'primary', 'primary_link', 'secondary', 'secondary_link', 'tertiary', 'tertiary_link'}:
        edge_data['hierarchy'] = 'main'

    if edge_data.get('highway') in {'residential', 'living_street', 'unclassified'}:
        edge_data['hierarchy'] = 'local'

    if edge_data.get('highway') in {'path', 'footway', 'cycleway'}:
        edge_data['hierarchy'] = 'path'

    if edge_data.get('highway') in {'construction', 'track'}:
        edge_data['hierarchy'] = 'other'


def _identify_dead_ends(street_graph):
    pass


def _identify_dead_ends(graph):

    # TODO: Auto detect how many iterations are necessary to get through the entire graph
    for i in range(100):

        for node in graph.nodes():

            # Get a list of all adjacent edges that are not dead ends
            adjacent_non_dead_end_edges = [
                edge for edge in list(graph.edges(nbunch=node, data=True, keys=True))
                if not edge[3].get('dead_end')
            ]

            # Only 1 adjacent edge -> mark this edge as dead end
            if len(adjacent_non_dead_end_edges) == 1:
                edge = adjacent_non_dead_end_edges[0]
                edge[3]['dead_end'] = True


"""
def _dead_ends_process_edge(graph, edge):

    edge_data = edge[3]

    # Mark as processed
    edge_data['processed'] = True

    # Get neighbor edges on both sides, exclude edges already identified as dead ends
    neighbors = graph_tools._get_neighbors(graph, edge, dead_ends=False)
    edge_data['n_neigh_u'] = len(neighbors[0])
    edge_data['n_neigh_v'] = len(neighbors[1])

    # Initialize
    next_edges = []

    # No neighbors on one side means a dead end
    edge_data['dead_end'] = min((len(neighbors[0])), (len(neighbors[1]))) == 0

    # Find neighbor edges that should be processed next
    for neighbor in neighbors[2]:
        neighbor_data = neighbor[3]
        if not neighbor_data.get('processed'):
            next_edges.append(neighbor)

    # return next edges to be processed
    return next_edges
"""
