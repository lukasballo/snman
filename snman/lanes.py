import math
import networkx as nx
from . import config


def generate_lanes(street_graph):
    """
    Reverse-engineer the lanes of each street edge and store them as a list in the attribute 'ln_desc'
    Naming convention:
    c = cycling,
    m = motorized traffic,
    > = forward,
    < = backward,
    - = both directions

    Parameters
    ----------
    street_graph : nx.MultiGraph
        contains the street network

    Returns
    -------
    None
    """
    for edge in street_graph.edges(data=True, keys=True):
        edge_data = edge[3]
        edge_data['ln_desc'] = _generate_lanes_for_edge(edge_data)


def _generate_lanes_for_edge(edge):
    """
    Reverse-engineer the lanes for one edge

    Parameters
    ----------
    edge : tuple
        (u, v, key, data)

    Returns
    -------
    lane_list : list
        a list of lanes, following the convention described under generate_lanes
    """
    lanes_list = []

    n_motorized_lanes = int(edge.get('lanes', -1))
    n_motorized_lanes_forward = int(edge.get('lanes:forward', -1))
    n_motorized_lanes_backward = int(edge.get('lanes:backward', -1))

    if edge['highway'] == 'footway' or edge['highway'] == 'path':
        lanes_list.append('c-')
    else:
        if edge['oneway'] == 1:
            if edge.get('cycleway:left') == 'lane' or edge.get('cycleway:both') == 'lane':
                lanes_list.append('c>')
            if n_motorized_lanes >= 1:
                lanes_list.extend(['m>'] * n_motorized_lanes)
            else:
                lanes_list.extend(['m>'])
            if edge.get('cycleway:right') == 'lane' or edge.get('cycleway:both') == 'lane' or edge.get('cycleway') == 'lane':
                lanes_list.append('c>')
        else:
            if edge.get('cycleway:left') == 'lane' or edge.get('cycleway:both') or edge.get('cycleway') == 'lane':
                lanes_list.append('c<')
            if n_motorized_lanes >= 2:
                if n_motorized_lanes_backward >= 0 and n_motorized_lanes_forward >= 0:
                    lanes_list.extend(['m<'] * n_motorized_lanes_backward)
                    lanes_list.extend(['m>'] * n_motorized_lanes_forward)
                else:
                    lanes_list.extend(['m<'] * math.floor(n_motorized_lanes / 2))
                    lanes_list.extend(['m>'] * math.ceil(n_motorized_lanes / 2))
            else:
                lanes_list.append('m-')
            if edge.get('cycleway:right') == 'lane' or edge.get('cycleway:both') or edge.get('cycleway') == 'lane':
                lanes_list.append('c>')

    #return ' | '.join(lanes_list)
    return lanes_list


def _reverse_lanes(lanes):
    """
    Reverse the order and direction of all lanes

    Parameters
    ----------
    lanes : list
        a list of lanes, following the convention described under generate_lanes

    Returns
    -------
    reversed_lanes : list
        lanes, with reversed order and directions
    """
    reversed_lanes = lanes
    reversed_lanes = [lane.replace('>', '>>') for lane in reversed_lanes]
    reversed_lanes = [lane.replace('<', '<<') for lane in reversed_lanes]
    reversed_lanes = [lane.replace('>>', '<') for lane in reversed_lanes]
    reversed_lanes = [lane.replace('<<', '>') for lane in reversed_lanes]
    reversed_lanes.reverse()
    return reversed_lanes


def generate_lane_stats(street_graph):
    """
    Add lane statistics to all edges for the street graph

    Params
    ------
    street_graph : nx.MultiGraph

    Returns
    -------
    None
    """
    for edge in street_graph.edges(data=True, keys=True):
        edge_data = edge[3]
        _generate_lane_stats_for_edge(edge_data)


def _generate_lane_stats_for_edge(edge):
    lanes_list = edge.get('ln_desc', '')
    n_lanes_cycling = lanes_list.count('c>') + lanes_list.count('c<') + lanes_list.count('c-') * 2
    n_lanes_motorized = lanes_list.count('m>') + lanes_list.count('m<') + lanes_list.count('m-') * 1.5
    width_cycling_m = n_lanes_cycling * config.lane_width_cycling_m
    width_motorized_m = n_lanes_motorized * config.lane_width_motorized_m
    width_total_m = width_cycling_m + width_motorized_m
    try:
        proportion_cycling = width_cycling_m / width_total_m
    except ZeroDivisionError:
        proportion_cycling = None

    edge['n_ln_cyc'] = n_lanes_cycling
    edge['n_ln_mot'] = n_lanes_motorized
    edge['w_cyc_m'] = width_cycling_m
    edge['w_mot_m'] = width_motorized_m
    edge['w_tot_m'] = width_total_m
    edge['prop_cyc'] = proportion_cycling

