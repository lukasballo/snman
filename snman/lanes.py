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

    if edge.get('highway') == 'footway' or edge.get('highway') == 'path':
        lanes_list.append('c-')
    else:
        if edge.get('oneway') == 1:
            if edge.get('cycleway:left') == 'lane' or edge.get('cycleway:both') == 'lane':
                lanes_list.append('c>')
            if n_motorized_lanes >= 1:
                lanes_list.extend(['m>'] * n_motorized_lanes)
            else:
                lanes_list.extend(['m>'])
            if edge.get('cycleway:right') == 'lane' or edge.get('cycleway:both') == 'lane' or edge.get('cycleway') == 'lane':
                lanes_list.append('c>')
        else:
            if edge.get('cycleway:left') == 'lane' or edge.get('cycleway:both') == 'lane' or edge.get('cycleway') == 'lane':
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
            if edge.get('cycleway:right') == 'lane' or edge.get('cycleway:both') == 'lane' or edge.get('cycleway') == 'lane':
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
    lanes = edge.get('ln_desc', [])

    width_cycling = 0
    width_motorized = 0
    width_total = 0

    for lane in lanes:
        lane_properties = _get_lane_properties(lane)
        if lane_properties['type'] == 'm':
            width_motorized += lane_properties['width']
        if lane_properties['type'] == 'c':
            width_cycling += lane_properties['width']
        width_total += lane_properties['width']

    try:
        proportion_cycling = width_cycling / width_total
    except ZeroDivisionError:
        proportion_cycling = None

    edge['w_cyc_m'] = width_cycling
    edge['w_mot_m'] = width_motorized
    edge['w_tot_m'] = width_total
    edge['prop_cyc'] = proportion_cycling

def _get_lane_properties(lane_description):
    return {
        'width': config.default_lane_widths_m.get(lane_description, 0),
        'type': lane_description[0:-1],
        'direction': lane_description[-1]
    }


