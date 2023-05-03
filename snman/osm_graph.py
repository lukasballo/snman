from . import any_graph, lane_config, centerline_graph


def prepare_graph(Gm):
    """
    A set of operations needed to make the street graph created by osmnx ready for snman operations

    Parameters
    ----------
    Gm : nx.MultiDiGraph
        OSM Graph

    Returns
    -------
    None
    """

    # ensure consistent data types
    for id, edge in Gm.edges.items():

        maxspeed = edge.get('maxspeed', '')
        edge['maxspeed'] = int(maxspeed) if maxspeed.isdigit() else -1

        layer = edge.get('layer', '')
        # isdigit only supports positive numbers, so we need to remove any '-' first
        edge['layer'] = int(layer) if layer.lstrip('-').isdigit() else 0

    centerline_graph.surrogate_missing_edge_geometries(Gm)


def organize_edge_directions(Gm, method='lower_to_higher_node_id'):
    """
    Ensure that all edges in a street graph have the direction from the lower to the higher node id.
    Edges with a different direction will be reversed, including the lanes and their geometry.

    Parameters
    ----------
    Gm : nx.MultiDiGraph
        OSM graph
    method : str
        - lower_to_higher_node_id: each edge will be digitized from the lower to the higher node id
        - by_osm_convention: one-way nodes will be digitized according to their lane direction

    Returns
    -------
    None
    """

    edges = list(Gm.edges(data=True, keys=True))
    for edge in edges:
        if method == 'lower_to_higher_node_id':
            if edge[0] > edge[1]:
                any_graph.reverse_edge(Gm, *edge[0:3])
        elif method == 'by_osm_convention':
            if lane_config._is_backward_oneway_street(edge[3].get('ln_desc')):
                any_graph.reverse_edge(Gm, *edge[0:3])
        else:
            raise 'Reorganization method not implemented: ' + str(method)
