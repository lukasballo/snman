import math
import networkx as nx
from .constants import *


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
        lanes_list.append(LANETYPE_CYCLING_PATH + DIRECTION_BOTH)
    else:
        if edge.get('oneway') == 1:
            if edge.get('cycleway:left') == 'lane' or edge.get('cycleway:both') == 'lane':
                lanes_list.append(LANETYPE_CYCLING_LANE + DIRECTION_FORWARD)
            if n_motorized_lanes >= 1:
                lanes_list.extend([LANETYPE_MOTORIZED + DIRECTION_FORWARD] * n_motorized_lanes)
            else:
                lanes_list.extend([LANETYPE_MOTORIZED + DIRECTION_FORWARD])
            if edge.get('cycleway:right') == 'lane' or edge.get('cycleway:both') == 'lane' or edge.get('cycleway') == 'lane':
                lanes_list.append(LANETYPE_CYCLING_LANE + DIRECTION_FORWARD)
        else:
            if edge.get('cycleway:left') == 'lane' or edge.get('cycleway:both') == 'lane' or edge.get('cycleway') == 'lane':
                lanes_list.append(LANETYPE_CYCLING_LANE + DIRECTION_BACKWARD)
            if n_motorized_lanes >= 2:
                if n_motorized_lanes_backward >= 0 and n_motorized_lanes_forward >= 0:
                    lanes_list.extend([LANETYPE_MOTORIZED + DIRECTION_BACKWARD] * n_motorized_lanes_backward)
                    lanes_list.extend([LANETYPE_MOTORIZED + DIRECTION_FORWARD] * n_motorized_lanes_forward)
                else:
                    lanes_list.extend([LANETYPE_MOTORIZED + DIRECTION_BACKWARD] * math.floor(n_motorized_lanes / 2))
                    lanes_list.extend([LANETYPE_MOTORIZED + DIRECTION_FORWARD] * math.ceil(n_motorized_lanes / 2))
            else:
                lanes_list.append(LANETYPE_MOTORIZED + DIRECTION_BOTH)
            if edge.get('cycleway:right') == 'lane' or edge.get('cycleway:both') == 'lane' or edge.get('cycleway') == 'lane':
                lanes_list.append(LANETYPE_CYCLING_LANE + DIRECTION_FORWARD)

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
    # We use >> and << as temporary symbols during the process
    reversed_lanes = [lane.replace(DIRECTION_FORWARD, '>>') for lane in reversed_lanes]
    reversed_lanes = [lane.replace(DIRECTION_BACKWARD, '<<') for lane in reversed_lanes]
    reversed_lanes = [lane.replace('>>', DIRECTION_BACKWARD) for lane in reversed_lanes]
    reversed_lanes = [lane.replace('<<', DIRECTION_FORWARD) for lane in reversed_lanes]
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
    lanes = edge.get(LANES_DESCRIPTION_KEY, [])

    width_cycling = 0
    width_motorized = 0
    width_total = 0

    for lane in lanes:
        lane_properties = _lane_properties(lane)
        if lane_properties.lanetype == LANETYPE_MOTORIZED:
            width_motorized += lane_properties.width
        if lane_properties.lanetype == LANETYPE_CYCLING_LANE:
            width_cycling += lane_properties.width
        width_total += lane_properties.width

    try:
        proportion_cycling = width_cycling / width_total
    except ZeroDivisionError:
        proportion_cycling = None

    edge['w_cyc_m'] = width_cycling
    edge['w_mot_m'] = width_motorized
    edge['w_tot_m'] = width_total
    edge['prop_cyc'] = proportion_cycling

class _lane_properties:

    width = None
    lanetype = None
    direction = None
    motorized = None
    private_cars = None
    dedicated_cycling = None
    dedicated_cycling_lane = None
    dedicated_cycling_path = None

    def __init__(self, lane_description):

        self.width = DEFAULT_LANE_WIDTHS.get(lane_description, 0)
        self.lanetype = lane_description[0:-1]
        self.direction = lane_description[-1]
        self.motorized = lane_description[0:-1] in [LANETYPE_MOTORIZED, LANETYPE_DEDICATED_PT]
        self.private_cars = lane_description[0:-1] == LANETYPE_MOTORIZED
        self.dedicated_cycling = lane_description[0:-1] in [LANETYPE_CYCLING_PATH, LANETYPE_CYCLING_LANE]
        self.dedicated_cycling_lane = lane_description[0:-1] == LANETYPE_CYCLING_LANE
        self.dedicated_cycling_path = lane_description[0:-1] == LANETYPE_CYCLING_PATH

class _lane_stats:

    n_lanes_motorized_forward = 0
    n_lanes_motorized_backward = 0
    n_lanes_motorized_both_ways = 0
    n_lanes_motorized_direction_tbd = 0

    n_lanes_dedicated_cycling_forward = 0
    n_lanes_dedicated_cycling_backward = 0
    n_lanes_dedicated_cycling_both_ways = 0
    n_lanes_dedicated_cycling_direction_tbd = 0

    def __init__(self,lanes_description):
        for lane in lanes_description:
            lane_properties = _lane_properties(lane)
            direction = lane_properties.direction

            # Motorized Lanes
            if lane_properties.motorized:
                if direction == DIRECTION_FORWARD:
                    self.n_lanes_motorized_forward += 1
                elif direction == DIRECTION_BACKWARD:
                    self.n_lanes_motorized_backward += 1
                elif direction == DIRECTION_BOTH:
                    self.n_lanes_motorized_both_ways += 1
                elif direction == DIRECTION_TBD:
                    self.n_lanes_motorized_direction_tbd += 1

            # Cycling
            if lane_properties.dedicated_cycling:
                if direction == DIRECTION_FORWARD:
                    self.n_lanes_dedicated_cycling_forward += 1
                elif direction == DIRECTION_BACKWARD:
                    self.n_lanes_dedicated_cycling_backward += 1
                elif direction == DIRECTION_BOTH:
                    self.n_lanes_dedicated_cycling_both_ways += 1
                elif direction == DIRECTION_TBD:
                    self.n_lanes_dedicated_cycling_direction_tbd += 1


        self.n_lanes_motorized = \
            self.n_lanes_motorized_forward + self.n_lanes_motorized_backward\
            + self.n_lanes_motorized_both_ways + self.n_lanes_motorized_direction_tbd

def update_osm_tags(street_graph, lanes_description_key=LANES_DESCRIPTION_KEY):
    for edge in street_graph.edges(data=True, keys=True):
        _update_osm_tags_for_edge(edge, lanes_description_key)

def _update_osm_tags_for_edge(edge, lanes_description_key):

    data = edge[3]
    lane_stats = _lane_stats(data.get(lanes_description_key, []))

    # Clean the tags before updating
    data['lanes'] = None
    data['lanes:forward'] = None
    data['lanes:backward'] = None
    data['lanes:both_ways'] = None

    # Motorized lanes
    if lane_stats.n_lanes_motorized:
        data['lanes'] = lane_stats.n_lanes_motorized
    if lane_stats.n_lanes_motorized_forward:
        data['lanes:forward'] = lane_stats.n_lanes_motorized_forward
    if lane_stats.n_lanes_motorized_backward:
        data['lanes:backward'] = lane_stats.n_lanes_motorized_backward
    if lane_stats.n_lanes_motorized_both_ways > 0:
        data['lanes:both_ways'] = lane_stats.n_lanes_motorized_both_ways

    # Clean the tags before updating
    data['cycleway'] = None
    data['cycleway:lane'] = None
    data['cycleway:right'] = None
    data['cycleway:right:lane'] = None
    data['cycleway:left'] = None
    data['cycleway:left:lane'] = None

    # Cycling Infrastructure
    # TODO: For now, we only assume advisory lanes
    if lane_stats.n_lanes_dedicated_cycling_both_ways > 0 \
        or (lane_stats.n_lanes_dedicated_cycling_forward > 0 and lane_stats.n_lanes_dedicated_cycling_backward > 0):
        data['cycleway'] = 'lane'
        data['cycleway:lane'] = 'advisory'
    elif lane_stats.n_lanes_dedicated_cycling_forward > 0:
        data['cycleway:right'] = 'lane'
        data['cycleway:right:lane'] = 'advisory'
    elif lane_stats.n_lanes_dedicated_cycling_backward > 0:
        data['cycleway:left'] = 'lane'
        data['cycleway:left:lane'] = 'advisory'
