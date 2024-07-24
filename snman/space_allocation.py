from .constants import *
import math
import networkx as nx
from . import utils, hierarchy
import numpy as np
import copy
import itertools as it


def generate_lanes(G, attr=KEY_LANES_DESCRIPTION):
    """
    Reverse-engineer the lanes of each street edge and store them as a list in an attribute

    Parameters
    ----------
    G : nx.MultiDiGraph
        street graph
    attr : str
        in which attribute should the lanes be stored

    Returns
    -------
    None
    """

    for edge in G.edges(data=True, keys=True):
        edge_data = edge[3]
        edge_data[attr] = _generate_lanes_for_edge(edge_data)


def _generate_lanes_for_edge(edge):
    # TODO: add pedestrian/cycling paths with cycling=designated as separate cycling and pedestrian paths
    """
    Reverse-engineer the lanes for one edge

    Resulting format: {lane type} + {lane direction}
        * lane type: see LANETYPE_* under .constants
        * lane direction: see DIRECTION_* under .constants
        * example: M> (lane for motorized traffic, forward)

    Parameters
    ----------
    edge : dict
        the data dictionary of an edge

    Returns
    -------
    lane_list : list
        a list of lanes, following the format described above
    """

    # PART 1: INITIALIZE VARIABLES

    # left/right lanes: cycling lanes that are not included in the osm lanes tag
    left_lanes_list = []
    forward_lanes_list = []
    both_dir_lanes_list = []
    backward_lanes_list = []
    right_lanes_list = []

    # Reverse forward/backward if the edge has been reversed in the conversion into undirected graph
    if edge.get(KEY_REVERSED, False) == True:
        _DIRECTION_FORWARD = DIRECTION_BACKWARD
        _DIRECTION_BACKWARD = DIRECTION_FORWARD
    else:
        _DIRECTION_FORWARD = DIRECTION_FORWARD
        _DIRECTION_BACKWARD = DIRECTION_BACKWARD

    # motorized lanes: replace with highway lanes if the road is a highway (motorway in osm terminology)
    if edge.get('highway') in hierarchy.HIGHWAY_OSM:
        _LANETYPE_MOTORIZED = LANETYPE_HIGHWAY
    else:
        _LANETYPE_MOTORIZED = LANETYPE_MOTORIZED

    # is this street oneway?
    oneway = edge.get('oneway', False) or edge.get('junction', False) == 'roundabout'

    # get the general lane count, fill in with default value according to the exact case
    if oneway:
        # 1 for all oneway streets
        n_lanes = int(float(edge.get('lanes', 1)))
    elif edge.get('maxspeed', -1) >= 50 or edge.get('highway') in {'primary', 'secondary'}:
        # 2 for speeds of at least 50 kmh or secondary/primary roads
        n_lanes = int(float(edge.get('lanes', 2)))
    else:
        # 1 for slower roads
        n_lanes = int(float(edge.get('lanes', 1)))

    # initialize detailed lane counts
    n_lanes_forward = utils.safe_int(edge.get('lanes:forward', 0), fallback_value=0)
    n_lanes_backward = utils.safe_int(edge.get('lanes:backward', 0), fallback_value=0)
    n_lanes_both = 0

    n_lanes_motorized = 0
    n_lanes_motorized_forward = 0
    n_lanes_motorized_backward = 0
    n_lanes_motorized_both = 0

    n_lanes_dedicated_pt = 0
    n_lanes_dedicated_pt_forward = 0
    n_lanes_dedicated_pt_backward = 0
    n_lanes_dedicated_pt_both = 0

    # PART 2: RECONSTRUCT LANE COUNTS

    # if forward/backward lanes are defined, make sure that the total lane count is consistent
    if n_lanes_forward or n_lanes_backward:
        n_lanes = max([n_lanes_forward + n_lanes_backward, n_lanes])

    # forward and backward lanes explicitly defined and consistent
    if n_lanes == n_lanes_forward + n_lanes_backward:
        pass
    # two-way street with more than 1 lane but no explicit lane counts forward/backward are defined
    elif n_lanes > 1 and not oneway:
        n_lanes_forward = math.floor(n_lanes / 2)
        n_lanes_backward = math.ceil(n_lanes / 2)
    # two-way street with exactly one lane
    elif n_lanes == 1 and not oneway:
        n_lanes_both = 1
    # oneway street with n_lanes defined
    elif n_lanes > 0 and oneway:
        n_lanes_forward = n_lanes

    # If the edge is dedicated for public transport
    if (edge.get('highway') == 'service' or edge.get('access') == 'no') and (
            edge.get('psv') == 'yes' or edge.get('bus') == 'yes'):
        if oneway:
            n_lanes_dedicated_pt = max([n_lanes, 1])
            n_lanes_dedicated_pt_forward = max([n_lanes_forward, 1])
        else:
            n_lanes_dedicated_pt = max([n_lanes, 2])
            n_lanes_dedicated_pt_forward = max([n_lanes_forward, 1])
            n_lanes_dedicated_pt_backward = max([n_lanes_backward, 1])
    else:
        n_lanes_motorized = n_lanes
        n_lanes_motorized_forward = n_lanes_forward
        n_lanes_motorized_backward = n_lanes_backward
        n_lanes_motorized_both = n_lanes_both

    # PART 3: GENERATE LANES

    # Cycleway [with walking]
    if edge.get('highway') in CYCLING_HIGHWAY_VALUES:
        # walking allowed without segregation
        if edge.get('foot') in {'yes', 'designated'} and edge.get('segregated') != 'yes':
            # oneway for cyclists
            if edge.get('oneway') in {'yes', 1} or edge.get('oneway:bicycle') in {'yes', 1}:
                both_dir_lanes_list.append(Lane(LANETYPE_FOOT_CYCLING_MIXED, _DIRECTION_FORWARD))
            # both ways for cyclists
            else:
                both_dir_lanes_list.append(Lane(LANETYPE_FOOT_CYCLING_MIXED, DIRECTION_BOTH))
        # walking allowed with segregation
        elif edge.get('foot') in {'yes', 'designated'} and edge.get('segregated') == 'yes':
            # oneway for cyclists
            if edge.get('oneway') in {'yes', 1} or edge.get('oneway:bicycle') in {'yes', 1}:
                both_dir_lanes_list.append(Lane(LANETYPE_FOOT, DIRECTION_BOTH))
                forward_lanes_list.append(Lane(LANETYPE_CYCLING_TRACK, _DIRECTION_FORWARD))
            # both ways for cyclists
            else:
                both_dir_lanes_list.append(Lane(LANETYPE_FOOT, DIRECTION_BOTH))
                both_dir_lanes_list.append(Lane(LANETYPE_CYCLING_TRACK, DIRECTION_BOTH))
        # walking not allowed
        else:
            # oneway for cyclists
            if edge.get('oneway') in {'yes', 1}:
                forward_lanes_list.append(Lane(LANETYPE_CYCLING_TRACK, _DIRECTION_FORWARD))
            # both ways for cyclists
            else:
                both_dir_lanes_list.append(Lane(LANETYPE_CYCLING_TRACK, DIRECTION_BOTH))

    # Walkway [with cycling]
    elif edge.get('highway') in PEDESTRIAN_HIGHWAY_VALUES:
        # cycling allowed without segregation
        if (edge.get('bicycle') in {'yes', 'designated'} or edge.get('bicycle:conditional')) \
                and edge.get('segregated') != 'yes':
            # oneway for cyclists
            if edge.get('oneway') in {'yes', 1} or edge.get('oneway:bicycle') in {'yes', 1}:
                both_dir_lanes_list.append(Lane(LANETYPE_FOOT_CYCLING_MIXED, _DIRECTION_FORWARD))
            # both ways for cyclists
            else:
                both_dir_lanes_list.append(Lane(LANETYPE_FOOT_CYCLING_MIXED, DIRECTION_BOTH))
        # cycling allowed with segregation
        elif (edge.get('bicycle') in {'yes', 'designated'} or edge.get('bicycle:conditional')) \
                and edge.get('segregated') == 'yes':
            # oneway for cyclists
            if edge.get('oneway') in {'yes', 1} or edge.get('oneway:bicycle') in {'yes', 1}:
                both_dir_lanes_list.append(Lane(LANETYPE_FOOT + DIRECTION_BOTH))
                forward_lanes_list.append(Lane(LANETYPE_CYCLING_TRACK + _DIRECTION_FORWARD))
            # both ways for cyclists
            else:
                both_dir_lanes_list.append(Lane(LANETYPE_FOOT, DIRECTION_BOTH))
                both_dir_lanes_list.append(Lane(LANETYPE_CYCLING_TRACK, DIRECTION_BOTH))
        # cycling not allowed
        else:
            both_dir_lanes_list.append(Lane(LANETYPE_FOOT, DIRECTION_BOTH))

    # Normal road
    elif edge.get('highway') not in {'platform'}:
        # Add sidewalk left (inactive to prevent double sidewalks)
        # if edge.get('sidewalk') in {'left', 'both'}:
        #    left_lanes_list.extend([LANETYPE_FOOT + DIRECTION_BOTH])
        # Add cycling lane left
        if edge.get('cycleway:left') == 'lane' \
                or edge.get('cycleway:both') == 'lane' \
                or edge.get('cycleway') == 'lane':
            left_lanes_list.append(Lane(LANETYPE_CYCLING_LANE, _DIRECTION_BACKWARD))

        # Add cycling lane right
        if edge.get('cycleway:right') == 'lane' \
                or edge.get('cycleway:both') == 'lane' \
                or edge.get('cycleway') == 'lane':
            right_lanes_list.append(Lane(LANETYPE_CYCLING_LANE, _DIRECTION_FORWARD))

        # Add cycling track left
        if edge.get('cycleway:left') == 'track' \
                or edge.get('cycleway:both') == 'track' \
                or edge.get('cycleway') == 'track':
            left_lanes_list.append(Lane(LANETYPE_CYCLING_TRACK, _DIRECTION_BACKWARD))

        # Add cycling track right
        if edge.get('cycleway:right') == 'track' \
                or edge.get('cycleway:both') == 'track' \
                or edge.get('cycleway') == 'track':
            right_lanes_list.append(Lane(LANETYPE_CYCLING_TRACK, _DIRECTION_FORWARD))

        # Add cycling allowed in opposite direction
        if edge.get('cycleway') == 'opposite':
            left_lanes_list.append(Lane(LANETYPE_CYCLING_PSEUDO, _DIRECTION_BACKWARD))

        # add parking lanes
        if edge.get('parking:left') == 'lane' or edge.get('parking:both') == 'lane':
            left_lanes_list.append(Lane(LANETYPE_PARKING_PARALLEL, DIRECTION_BOTH))
        if edge.get('parking:right') == 'lane' or edge.get('parking:both') == 'lane':
            right_lanes_list.append(Lane(LANETYPE_PARKING_PARALLEL, DIRECTION_BOTH))

        # Add sidewalk right
        # if edge.get('sidewalk') in {'right', 'both'}:
        #    right_lanes_list.extend([LANETYPE_FOOT + DIRECTION_BOTH])

        backward_lanes_list.extend(
            [Lane(_LANETYPE_MOTORIZED, _DIRECTION_BACKWARD) for i in range(n_lanes_motorized_backward)]
        )
        backward_lanes_list.extend(
            [Lane(LANETYPE_DEDICATED_PT, _DIRECTION_BACKWARD) for i in range(n_lanes_dedicated_pt_backward)]
        )

        start_both_dir_lanes = len(backward_lanes_list)
        both_dir_lanes_list.extend(
            [Lane(LANETYPE_DEDICATED_PT, DIRECTION_BOTH) for i in range(n_lanes_dedicated_pt_both)]
        )
        both_dir_lanes_list.extend(
            [Lane(_LANETYPE_MOTORIZED, DIRECTION_BOTH) for i in range(n_lanes_motorized_both)]
        )

        start_forward_lanes = start_both_dir_lanes + len(both_dir_lanes_list)
        forward_lanes_list.extend(
            [Lane(LANETYPE_DEDICATED_PT, _DIRECTION_FORWARD) for i in range(n_lanes_dedicated_pt_forward)]
        )
        forward_lanes_list.extend(
            [Lane(_LANETYPE_MOTORIZED, _DIRECTION_FORWARD) for i in range(n_lanes_motorized_forward)]
        )

    # Everything else
    else:
        pass

    # PART 4: APPLY DEDICATED LANES BASED ON DETAILED LANE DESCRIPTIONS

    # apply bus lanes
    # first, try to use the "bus:lanes:*" tag. if not possible, use the "vehicle:lanes:*" tag
    # but never use both as this could lead to conflicts
    if edge.get('bus:lanes'):
        osm_bus_lanes = edge.get('bus:lanes', '').split('|')
        for i, lane in enumerate(osm_bus_lanes):
            if lane == 'designated' and i < len(backward_lanes_list):
                backward_lanes_list[i] = Lane(LANETYPE_DEDICATED_PT, _DIRECTION_BACKWARD)
            elif lane == 'designated' and i < start_forward_lanes:
                both_dir_lanes_list[i-start_both_dir_lanes] = Lane(LANETYPE_DEDICATED_PT, DIRECTION_BOTH)
            elif lane == 'designated' and i < start_forward_lanes + len(forward_lanes_list):
                forward_lanes_list[i-start_forward_lanes] = Lane(LANETYPE_DEDICATED_PT, _DIRECTION_FORWARD)

    elif edge.get('vehicle:lanes'):
        osm_bus_lanes = edge.get('vehicle:lanes', '').split('|')
        for i, lane in enumerate(osm_bus_lanes):
            if lane == 'no' and i < len(backward_lanes_list):
                backward_lanes_list[i] = Lane(LANETYPE_DEDICATED_PT, _DIRECTION_BACKWARD)
            elif lane == 'no' and i < start_forward_lanes:
                both_dir_lanes_list[i-start_both_dir_lanes] = Lane(LANETYPE_DEDICATED_PT, DIRECTION_BOTH)
            elif lane == 'no' and i < start_forward_lanes + len(forward_lanes_list):
                forward_lanes_list[i-start_forward_lanes] = Lane(LANETYPE_DEDICATED_PT, _DIRECTION_FORWARD)

    if edge.get('bus:lanes:forward'):
        osm_bus_lanes_forward = edge.get('bus:lanes:forward', '').split('|')
        for i, lane in enumerate(osm_bus_lanes_forward):
            if lane == 'designated' and i < len(forward_lanes_list):
                forward_lanes_list[i] = Lane(LANETYPE_DEDICATED_PT, _DIRECTION_FORWARD)

    elif edge.get('vehicle:lanes:forward'):
        osm_vehicle_lanes_forward = edge.get('vehicle:lanes:forward', '').split('|')
        for i, lane in enumerate(osm_vehicle_lanes_forward):
            if lane == 'no' and i < len(forward_lanes_list):
                forward_lanes_list[i] = Lane(LANETYPE_DEDICATED_PT, _DIRECTION_FORWARD)

    elif edge.get('busway:right') or edge.get('busway:both') or edge.get('busway'):
        utils.set_last_or_append(forward_lanes_list, Lane(LANETYPE_DEDICATED_PT, _DIRECTION_FORWARD))

    # do the same in backward direction
    if edge.get('bus:lanes:backward'):
        osm_bus_lanes_backward = edge.get('bus:lanes:backward', '').split('|')
        # reverse the lane order because we are  still working in the forward direction
        osm_bus_lanes_backward.reverse()
        for i, lane in enumerate(osm_bus_lanes_backward):
            if lane == 'designated' and i < len(backward_lanes_list):
                backward_lanes_list[i] = Lane(LANETYPE_DEDICATED_PT, _DIRECTION_BACKWARD)

    elif edge.get('vehicle:lanes:backward'):
        osm_vehicle_lanes_backward = edge.get('vehicle:lanes:backward', '').split('|')
        # reverse the lane order because we are still working in the forward direction
        osm_vehicle_lanes_backward.reverse()
        for i, lane in enumerate(osm_vehicle_lanes_backward):
            if lane == 'no' and i < len(backward_lanes_list):
                backward_lanes_list[i] = Lane(LANETYPE_DEDICATED_PT, _DIRECTION_BACKWARD)

    elif edge.get('busway:left') or edge.get('busway:both') or edge.get('busway'):
        utils.set_last_or_append(backward_lanes_list, Lane(LANETYPE_DEDICATED_PT, _DIRECTION_BACKWARD))

    # PART 5: RETURN

    # merge all pars of the road section
    lanes_list = left_lanes_list + backward_lanes_list + both_dir_lanes_list + forward_lanes_list + right_lanes_list
    return SpaceAllocation(lanes_list)


def generate_lane_stats(G, lanes_attribute=KEY_LANES_DESCRIPTION):
    """
    Add lane statistics to all edges for the street graph

    Parameters
    ------
    G : nx.MultiGraph
        street graph
    lanes_attribute : str
        which attribute describing the lanes should be used
        (e.g., lanes in status quo or lanes after rebuilding)

    Returns
    -------
    None
    """

    for edge in G.edges(data=True, keys=True):
        edge_data = edge[3]
        _generate_lane_stats_for_edge(edge_data, lanes_attribute)


def _generate_lane_stats_for_edge(edge, lanes_attribute=KEY_LANES_DESCRIPTION):
    # TODO: Generate stats for both status quo and after rebuilding
    """
    Add lane statistics to one edge. Following attributes are added:
        * width_cycling_m -> total width for cycling in meters
        * width_motorized_m -> total width for motorized traffic in meters
        * width_total_m -> total width of all lanes
        * n_lanes_motorized -> number of lanes for motorized traffic

    Parameters
    ----------
    edge : dict
        the data dictionary of an edge
    lanes_attribute : str
        which attribute describing the lanes should be used
        (e.g., lanes in status quo or lanes after rebuilding)

    Returns
    -------
    None
    """

    lanes = edge.get(lanes_attribute, [])
    lane_stats = _lane_stats(lanes)

    width_cycling = 0
    width_motorized = 0
    width_total = 0

    for lane in lanes:
        if lane.lanetype == LANETYPE_MOTORIZED:
            width_motorized += lane.width
        if lane.lanetype == LANETYPE_CYCLING_LANE:
            width_cycling += lane.width
        width_total += lane.width

    try:
        proportion_cycling = width_cycling / width_total
    except ZeroDivisionError:
        proportion_cycling = None

    # basic stats
    edge[lanes_attribute + '_width_cycling_m'] = width_cycling
    edge[lanes_attribute + '_width_motorized_m'] = width_motorized
    edge[lanes_attribute + '_width_total_m'] = width_total
    edge[lanes_attribute + '_n_lanes_motorized'] = lane_stats.n_lanes_motorized

    # description of best cycling option in each direction
    for user_dir_name, user_dir_description in {'forward': DIRECTION_FORWARD, 'backward': DIRECTION_BACKWARD}.items():
        for lane_direction in [user_dir_description, DIRECTION_BOTH]:
            for lanetype in CYCLING_QUALITY_HIERARCHY:
                lane_description = lanetype + lane_direction
                if lane_description in lanes:
                    edge[lanes_attribute + '_cycling_' + user_dir_name] = lane_description
                    break


    def get_standard_width(self):
        return LANE_TYPES[self.lanetype + self.direction]['width']

    def has_atypical_width(self):
        if self.width == self.get_standard_width():
            return False
        else:
            return True

    def __str__(self):
        if self.has_atypical_width():
            return self.lanetype + self.direction + str(self.width)
        else:
            return self.lanetype + self.direction


class _lane_stats:
    """
    A class for a standardized set of statistics over all lanes
    """

    n_lanes_motorized = 0
    n_lanes_motorized_forward = 0
    n_lanes_motorized_backward = 0
    n_lanes_motorized_both_ways = 0
    n_lanes_motorized_direction_tbd = 0

    n_lanes_private_cars = 0
    n_lanes_private_cars_forward = 0
    n_lanes_private_cars_backward = 0
    n_lanes_private_cars_both_ways = 0
    n_lanes_private_cars_direction_tbd = 0

    n_lanes_dedicated_pt = 0
    n_lanes_dedicated_pt_forward = 0
    n_lanes_dedicated_pt_backward = 0
    n_lanes_dedicated_pt_both_ways = 0
    n_lanes_dedicated_pt_direction_tbd = 0

    n_lanes_dedicated_cycling_lanes_forward = 0
    n_lanes_dedicated_cycling_lanes_backward = 0
    n_lanes_dedicated_cycling_lanes_both_ways = 0
    n_lanes_dedicated_cycling_lanes_direction_tbd = 0

    n_lanes_dedicated_cycling_tracks_forward = 0
    n_lanes_dedicated_cycling_tracks_backward = 0
    n_lanes_dedicated_cycling_tracks_both_ways = 0
    n_lanes_dedicated_cycling_tracks_direction_tbd = 0

    modes = ()

    def __init__(self, lanes):
        """
        Decodes a list of lanes into a set of statistics

        Parameters
        ----------
        lanes_description : SpaceAllocation
            a list of lanes following the format described in _generate_lanes_for_edge
        """

        self.modes = set()
        for lane in lanes:

            if lane == '':
                continue

            direction = lane.direction

            self.modes.update(lane.get_modes())

            # Motorized Lanes
            if lane.has_large_vehicles():
                if direction == DIRECTION_FORWARD:
                    self.n_lanes_motorized_forward += 1
                elif direction == DIRECTION_BACKWARD:
                    self.n_lanes_motorized_backward += 1
                elif direction == DIRECTION_BOTH:
                    self.n_lanes_motorized_both_ways += 1
                elif direction == DIRECTION_TBD:
                    self.n_lanes_motorized_direction_tbd += 1

            # Private cars
            if lane.has_mode(MODE_PRIVATE_CARS):
                if direction == DIRECTION_FORWARD:
                    self.n_lanes_private_cars_forward += 1
                elif direction == DIRECTION_BACKWARD:
                    self.n_lanes_private_cars_backward += 1
                elif direction == DIRECTION_BOTH:
                    self.n_lanes_private_cars_both_ways += 1
                elif direction == DIRECTION_TBD:
                    self.n_lanes_private_cars_direction_tbd += 1

            # PT lanes
            if lane.lanetype == LANETYPE_DEDICATED_PT:
                if direction == DIRECTION_FORWARD:
                    self.n_lanes_dedicated_pt_forward += 1
                if direction == DIRECTION_BACKWARD:
                    self.n_lanes_dedicated_pt_backward += 1
                if direction == DIRECTION_BOTH:
                    self.n_lanes_dedicated_pt_both_ways += 1
                if direction == DIRECTION_TBD:
                    self.n_lanes_dedicated_pt_direction_tbd += 1

            # Cycling lanes
            if lane.lanetype == LANETYPE_CYCLING_LANE:
                if direction == DIRECTION_FORWARD:
                    self.n_lanes_dedicated_cycling_lanes_forward += 1
                elif direction == DIRECTION_BACKWARD:
                    self.n_lanes_dedicated_cycling_lanes_backward += 1
                elif direction == DIRECTION_BOTH:
                    self.n_lanes_dedicated_cycling_lanes_both_ways += 1
                elif direction == DIRECTION_TBD:
                    self.n_lanes_dedicated_cycling_lanes_direction_tbd += 1

            # Cycling paths
            if lane.lanetype in {LANETYPE_CYCLING_TRACK, LANETYPE_FOOT_CYCLING_MIXED}:
                if direction == DIRECTION_FORWARD:
                    self.n_lanes_dedicated_cycling_tracks_forward += 1
                elif direction == DIRECTION_BACKWARD:
                    self.n_lanes_dedicated_cycling_tracks_backward += 1
                elif direction == DIRECTION_BOTH:
                    self.n_lanes_dedicated_cycling_tracks_both_ways += 1
                elif direction == DIRECTION_TBD:
                    self.n_lanes_dedicated_cycling_tracks_direction_tbd += 1

        self.n_lanes_motorized = \
            self.n_lanes_motorized_forward + self.n_lanes_motorized_backward \
            + self.n_lanes_motorized_both_ways + self.n_lanes_motorized_direction_tbd

        self.n_lanes_private_cars = \
            self.n_lanes_private_cars_forward + self.n_lanes_private_cars_backward \
            + self.n_lanes_private_cars_both_ways + self.n_lanes_private_cars_direction_tbd

        self.n_lanes_dedicated_pt = \
            self.n_lanes_dedicated_pt_forward + self.n_lanes_dedicated_pt_backward \
            + self.n_lanes_dedicated_pt_both_ways + self.n_lanes_dedicated_pt_direction_tbd


def _get_opposite_direction(direction):
    if direction == DIRECTION_FORWARD:
        return DIRECTION_BACKWARD
    elif direction == DIRECTION_BACKWARD:
        return DIRECTION_FORWARD
    else:
        return direction


def update_osm_tags(G, lanes_description_key=KEY_LANES_DESCRIPTION):
    """
    Update the osm tags of all edges to match their current lanes. This is necessary after the simplification when
    multiple edges are merged into single edge

    Parameters
    ----------
    G : nx.MultiGraph
        street graph
    lanes_description_key : str
        which attribute should be used as a source of lane data

    Returns
    -------
    None
    """
    for edge in G.edges(data=True, keys=True):
        _update_osm_tags_for_edge(edge, lanes_description_key)


def _update_osm_tags_for_edge(edge, lanes_description_key):
    """
    Update OSM tags of one edge to match its current lanes

    Parameters
    ----------
    edge
    lanes_description_key

    Returns
    -------
    None
    """

    data = edge[3]
    lane_stats = _lane_stats(data.get(lanes_description_key, []))

    # Clean the tags before updating
    data['lanes'] = None
    data['lanes:forward'] = None
    data['lanes:backward'] = None
    data['lanes:both_ways'] = None
    data['oneway'] = None

    # Update the highway tag
    if lane_stats.n_lanes_private_cars > 0:
        # leave the highway tag as it is
        pass
    elif lane_stats.n_lanes_dedicated_pt > 0:
        # set to service road
        data['highway'] = 'service'
    else:
        # set to path
        data['highway'] = 'path'

    # Clean the tags before updating
    data['bus:lanes:backward'] = None
    data['bus:lanes:forward'] = None
    data['vehicle:lanes:backward'] = None
    data['vehicle:lanes:forward'] = None

    # Motorized lanes
    if lane_stats.n_lanes_motorized > 0:

        data['lanes'] = lane_stats.n_lanes_motorized
        if lane_stats.n_lanes_motorized_forward:
            data['lanes:forward'] = lane_stats.n_lanes_motorized_forward
        if lane_stats.n_lanes_motorized_backward:
            data['lanes:backward'] = lane_stats.n_lanes_motorized_backward
        if lane_stats.n_lanes_motorized_both_ways > 0:
            data['lanes:both_ways'] = lane_stats.n_lanes_motorized_both_ways
        if (
            lane_stats.n_lanes_motorized_forward > 0
            and lane_stats.n_lanes_motorized_backward + lane_stats.n_lanes_motorized_both_ways == 0
        ):
            data['oneway'] = 'yes'
        elif (
            lane_stats.n_lanes_motorized_backward > 0
            and lane_stats.n_lanes_motorized_forward + lane_stats.n_lanes_motorized_both_ways == 0
        ):
            data['oneway'] = '-1'
        else:
            data['oneway'] = 'no'

        # PT lanes
        if lane_stats.n_lanes_dedicated_pt_both_ways > 0 or lane_stats.n_lanes_dedicated_pt_backward > 0:
            data['bus:lanes:backward'] = '|'.join(
                ['designated'] * lane_stats.n_lanes_dedicated_pt_both_ways +
                ['designated'] * lane_stats.n_lanes_dedicated_pt_backward +
                ['permissive'] * lane_stats.n_lanes_private_cars_both_ways +
                ['permissive'] * lane_stats.n_lanes_private_cars_backward
            )
            data['vehicle:lanes:backward'] = '|'.join(
                ['no'] * lane_stats.n_lanes_dedicated_pt_both_ways +
                ['no'] * lane_stats.n_lanes_dedicated_pt_backward +
                ['yes'] * lane_stats.n_lanes_private_cars_both_ways +
                ['yes'] * lane_stats.n_lanes_private_cars_backward
            )

        if lane_stats.n_lanes_dedicated_pt_both_ways > 0 or lane_stats.n_lanes_dedicated_pt_forward > 0:
            data['bus:lanes:forward'] = '|'.join(
                ['designated'] * lane_stats.n_lanes_dedicated_pt_both_ways +
                ['designated'] * lane_stats.n_lanes_dedicated_pt_forward +
                ['permissive'] * lane_stats.n_lanes_private_cars_both_ways +
                ['permissive'] * lane_stats.n_lanes_private_cars_forward
            )
            data['vehicle:lanes:forward'] = '|'.join(
                ['no'] * lane_stats.n_lanes_dedicated_pt_both_ways +
                ['no'] * lane_stats.n_lanes_dedicated_pt_forward +
                ['yes'] * lane_stats.n_lanes_private_cars_both_ways +
                ['yes'] * lane_stats.n_lanes_private_cars_forward
            )

    # Clean the tags before updating
    data['cycleway'] = None
    data['cycleway:lane'] = None
    data['cycleway:right'] = None
    data['cycleway:right:lane'] = None
    data['cycleway:left'] = None
    data['cycleway:left:lane'] = None

    # Cycling lanes
    if lane_stats.n_lanes_dedicated_cycling_lanes_both_ways > 0 \
            or (
            lane_stats.n_lanes_dedicated_cycling_lanes_forward > 0 and lane_stats.n_lanes_dedicated_cycling_lanes_backward > 0):
        # both directions
        data['cycleway'] = 'lane'
        data['cycleway:lane'] = 'advisory'
    elif lane_stats.n_lanes_dedicated_cycling_lanes_forward > 0:
        # only forward
        data['cycleway:right'] = 'lane'
        data['cycleway:right:lane'] = 'advisory'
    elif lane_stats.n_lanes_dedicated_cycling_lanes_backward > 0:
        # only backward
        data['cycleway:left'] = 'lane'
        data['cycleway:left:lane'] = 'advisory'

    # Cycling tracks
    if lane_stats.n_lanes_dedicated_cycling_tracks_both_ways > 0 \
            or (
            lane_stats.n_lanes_dedicated_cycling_tracks_forward > 0 and lane_stats.n_lanes_dedicated_cycling_tracks_backward > 0):
        # both directions
        data['cycleway'] = 'track'
    elif lane_stats.n_lanes_dedicated_cycling_tracks_forward > 0:
        # only forward
        data['cycleway:right'] = 'track'
    elif lane_stats.n_lanes_dedicated_cycling_tracks_backward > 0:
        # only backward
        data['cycleway:left'] = 'track'

    maxspeed = data.get('maxspeed', -1)
    if maxspeed == -1 and 'maxspeed' in data:
        del data['maxspeed']


def is_backward_oneway_street(lanes):
    """
    Returns true if the given lane configuration represents a one-way street that is digitized in opposite direction.
    This is needed for preparing the graph for osm export where one-way streets are typically digitized
    in the forward direction.

    Parameters
    ----------
    lanes: list
        a list of lanes following the format described in _generate_lanes_for_edge

    Returns
    -------
    bool
    """

    ls = _lane_stats(lanes)

    if ls.n_lanes_motorized > 0:
        if (
            ls.n_lanes_motorized_forward == 0
            and ls.n_lanes_motorized_both_ways == 0
            and ls.n_lanes_motorized_direction_tbd == 0
            and ls.n_lanes_motorized_backward > 0
        ):
            return True

    elif (
            ls.n_lanes_dedicated_cycling_lanes_forward + ls.n_lanes_dedicated_cycling_tracks_forward == 0
            and ls.n_lanes_dedicated_cycling_lanes_both_ways + ls.n_lanes_dedicated_cycling_tracks_both_ways == 0
            and ls.n_lanes_dedicated_cycling_lanes_direction_tbd + ls.n_lanes_dedicated_cycling_tracks_direction_tbd == 0
            and ls.n_lanes_dedicated_cycling_lanes_backward + ls.n_lanes_dedicated_cycling_tracks_backward > 0
    ):
        return True

    return False


def is_backward_by_top_order_lanes(lanes):
    """
    This is the basis of a custom direction organization algorithm:
    A lane is considered as backward if more than 50% of the top-order lanes point backward.
    Lanes without an explicit direction are ignored.

    Parameters
    ----------
    lanes : SpaceAllocation

    Returns
    -------
    bool
    """

    lane_order_list = []
    for lane in lanes:
        if lane.direction in (DIRECTION_FORWARD, DIRECTION_BACKWARD):
            lane_order_list.append(lane.get_order())

    if len(lane_order_list) == 0:
        return False

    top_order = min(lane_order_list)

    balance = 0
    for lane in lanes:
        if lane.get_order() != top_order:
            continue
        if lane.direction == DIRECTION_FORWARD:
            balance += 1
        else:
            balance -= 1

    return balance < 0


def add_pseudo_contraflow_cycling(G, lane_attribute=KEY_LANES_DESCRIPTION):

    for uvk, data in G.edges.items():

        lanes = data[lane_attribute]

        # remove any existing S lanes
        for i, lane in enumerate(copy.copy(lanes)):
            if lane.lanetype == LANETYPE_CYCLING_PSEUDO:
                lanes.remove(i)

        # add contraflow S lanes if there are M lanes but no cycling infra
        has_m_lanes = False
        cycling_forward = False
        for i, lane in enumerate(copy.copy(lanes)):
            if lane.get_primary_mode() == MODE_PRIVATE_CARS:
                has_m_lanes = True




def reorder_lanes(G, lanes_attribute=KEY_LANES_DESCRIPTION, **kwargs):
    """
    Reorder lanes in the street graph: Puts the lanes on each street into a consistent order and consolidates cycling
    infrastructure.

    Parameters
    ----------
    G
    lanes_attribute

    Returns
    -------

    """

    for uvk, data in G.edges.items():
        data[lanes_attribute] = _reorder_lanes_on_edge(data[lanes_attribute], **kwargs)


def _reorder_lanes_on_edge(lanes, how='standard'):
    """
    Reorder the lane list and consolidate cycling infrastructure of an edge

    Parameters
    ----------
    lanes : SpaceAllocation

    Returns
    -------
    SpaceAllocation
    """

    # prepare a structure of street parts
    sorted_lanes = {}
    for mode in MODES:
        sorted_lanes[mode] = {}
        for direction in DIRECTIONS:
            sorted_lanes[mode][direction] = []

    # sort by primary mode and direction
    for i, l in enumerate(lanes):
        sorted_lanes[l.get_primary_mode()][l.direction].append(l)

    for mode in MODES:
        for direction in DIRECTIONS:
            sorted_lanes[mode][direction] = list(sorted(sorted_lanes[mode][direction]))

    # calculate stats
    cycling_lanes = list(utils.flatten_list(sorted_lanes[MODE_CYCLING].values()))
    width_cycling_total = sum([l.width for l in cycling_lanes])
    width_mixed = sum([l.width for l in cycling_lanes if l.lanetype == LANETYPE_FOOT_CYCLING_MIXED])
    n_mixed_total = len([l for l in cycling_lanes if l.lanetype == LANETYPE_FOOT_CYCLING_MIXED])
    n_cars = len(list(utils.flatten_list(sorted_lanes[MODE_PRIVATE_CARS])))

    # consolidate parking lanes, merge all parking lanes into a single with adjusted width
    parking_lanes =(
            sorted_lanes[MODE_CAR_PARKING][DIRECTION_FORWARD]
            + sorted_lanes[MODE_CAR_PARKING][DIRECTION_BACKWARD]
    )
    width = sum([lane.width for lane in parking_lanes])
    if width > 0:
        consolidated_parking = Lane(LANETYPE_PARKING_PARALLEL, DIRECTION_FORWARD)
        consolidated_parking.width = width
        sorted_lanes[MODE_CAR_PARKING][DIRECTION_FORWARD] = [consolidated_parking]

    new_allocation = copy.copy(lanes)
    new_allocation[:] = []

    if how == 'bidirectional_cycling':

        ordered_lanes = [

            sorted_lanes[MODE_FOOT][DIRECTION_BOTH][0::2],

            sorted_lanes[MODE_NON_TRAFFIC][DIRECTION_BACKWARD],
            sorted_lanes[MODE_NON_TRAFFIC][DIRECTION_FORWARD],
            sorted_lanes[MODE_PRIVATE_CARS][DIRECTION_BACKWARD],
            sorted_lanes[MODE_TRANSIT][DIRECTION_BACKWARD],

            sorted_lanes[MODE_TRANSIT][DIRECTION_BOTH],
            sorted_lanes[MODE_PRIVATE_CARS][DIRECTION_BOTH],

            sorted_lanes[MODE_TRANSIT][DIRECTION_FORWARD],
            sorted_lanes[MODE_PRIVATE_CARS][DIRECTION_FORWARD],
            sorted_lanes[MODE_CAR_PARKING][DIRECTION_FORWARD],

            sorted_lanes[MODE_CYCLING][DIRECTION_BACKWARD],
            sorted_lanes[MODE_CYCLING][DIRECTION_FORWARD],
            sorted_lanes[MODE_CYCLING][DIRECTION_BOTH],
            sorted_lanes[MODE_FOOT][DIRECTION_BOTH][1::2]

        ]

    else:

        ordered_lanes = [

            sorted_lanes[MODE_FOOT][DIRECTION_BOTH][0::2],
            sorted_lanes[MODE_CYCLING][DIRECTION_BOTH][0::2],

            sorted_lanes[MODE_CYCLING][DIRECTION_BACKWARD],
            sorted_lanes[MODE_NON_TRAFFIC][DIRECTION_BACKWARD],
            sorted_lanes[MODE_NON_TRAFFIC][DIRECTION_FORWARD],
            sorted_lanes[MODE_PRIVATE_CARS][DIRECTION_BACKWARD],
            sorted_lanes[MODE_TRANSIT][DIRECTION_BACKWARD],

            sorted_lanes[MODE_TRANSIT][DIRECTION_BOTH],
            sorted_lanes[MODE_PRIVATE_CARS][DIRECTION_BOTH],

            sorted_lanes[MODE_TRANSIT][DIRECTION_FORWARD],
            sorted_lanes[MODE_PRIVATE_CARS][DIRECTION_FORWARD],
            sorted_lanes[MODE_CAR_PARKING][DIRECTION_FORWARD],
            sorted_lanes[MODE_CYCLING][DIRECTION_FORWARD],

            sorted_lanes[MODE_CYCLING][DIRECTION_BOTH][1::2],
            sorted_lanes[MODE_FOOT][DIRECTION_BOTH][1::2]

        ]

    new_allocation.extend(it.chain.from_iterable(ordered_lanes))

    return new_allocation


def normalize_cycling_lanes(G, lanes_key=KEY_LANES_DESCRIPTION):
    """
    Replaces all cycling paths with cycling lanes for simplicity in the masterplan.

    Parameters
    ----------
    G: nx.MultiDiGraph

    Returns
    -------
    none
    """

    for uvk, data in G.edges.items():
        lanes = data[lanes_key]
        for i, lane in enumerate(lanes):
            lane
            if lane.lanetype == LANETYPE_CYCLING_TRACK:
                lane.lanetype = LANETYPE_CYCLING_LANE
                lanes[i] = lane


def _calculate_lane_cost(lane, length, slope, mode, direction=DIRECTION_FORWARD, include_tentative=False):
    """
    Returns the cost of traversing this lane using the specified mode.
    The resulting cost is relative to other lanes with the same mode but is not comparable across modes.
    (Use a separate factor outside this function for proper mode choice)

    Parameters
    ----------
    lane: Lane
    length: float
        length in meters
    slope: float
        slope, e.g. 0.05 for 5% incline
    mode: str
        mode from constants.MODES
    direction: str
    include_tentative: bool
        if true, also lanes with tentative direction will be included

    Returns
    -------
    float
    """

    lane

    if include_tentative and lane.status != STATUS_FIXED:
        return np.Inf

    # if this lane can not carry the specified mode, return infinity
    if mode not in lane.get_modes():
        return np.Inf

    # if this lane is not accessible in the specified direction, return infinity
    if lane.direction not in {direction}.union(ALTERNATIVE_DIRECTIONS.get(direction, set())):
        return np.Inf

    # apply the cycling cost factor if the mode is cycling
    elif mode == MODE_CYCLING:
        # adjust the slope according to the direction
        if direction in {DIRECTION_FORWARD, DIRECTION_FORWARD_OPTIONAL}:
            pass
        elif direction in {DIRECTION_BACKWARD, DIRECTION_BACKWARD_OPTIONAL}:
            slope = -slope
        else:
            # cannot calculate cycling cost if direction is unknown
            slope = math.inf
        return length * (1 + lane.get_cycling_vod_factor() + CYCLING_SLOPE_VOD(slope))
    # otherwise, return just the length
    else:
        return length


def filter_lanes_by_modes(lanes, modes, operator='or'):
    """
    Returns a subset of lanes given filter criteria

    Parameters
    ----------
    lanes : SpaceAllocation
    modes : set
    operator : str
        - 'exact': modes must match the lane modes exactly
        - 'or': at least one of the modes must match a lane mode
        - 'and': all modes must match a lane mode

    Returns
    -------
    list

    """

    lanes_result = copy.copy(lanes)
    lanes_result[:] = []

    if operator == 'exact':
        lanes_result += [lane for lane in lanes if lane.get_modes() == modes]
    elif operator == 'or':
        lanes_result += [lane for lane in lanes if not lane.get_modes().isdisjoint(modes)]
    elif operator == 'and':
        lanes_result += [lane for lane in lanes if modes.issubset(lane.get_modes())]

    return lanes_result


def filter_lanes_by_function(lanes, filter_function):

    lanes_result = copy.copy(lanes)
    lanes_result[:] = []

    lanes_result += [lane for lane in lanes if filter_function(lane) is True]

    return lanes_result


def space_allocation_from_string(sa_string):
    lane_strings = sa_string.split(' | ')
    lanes = [
        Lane(s[0], s[1], status=s[2], width=utils.safe_float(s[3:])) if len(s) > 2
        else Lane(s[0], s[1]) if len(s) == 2 # fallback for old lane format
        else None
        for s in lane_strings
    ]
    return SpaceAllocation(lanes)


def lane_from_string(s):
    return Lane(s[0], s[1], status=s[2], width=utils.safe_float(s[3:]))


class Lane:
    def __init__(self, lanetype, direction, status=STATUS_FIXED, width=None):
        """
        Create a new lane

        Parameters
        ----------
        lanetype: str
            See `snman.constants.LANETYPE_*` for options.
        direction: str
            See `snman.constants.DIRECTION_*` for options.
        status: str
            Denotes whether the lane is fixed or still negotiable. See `snman.constants.STATUS_*` for options.
            This is used in the reallocation process.
        width: float
            Physical width of the street that can be allocated.
            If None, the default width from `snman.constants` will be used
        """

        self.lanetype = lanetype
        self.direction = direction
        self.status = status
        if width:
            self.width = width
        else:
            self.set_default_width()
        self.valid = True

    def __str__(self):
        return self.lanetype + self.direction + self.status + str(self.width)

    def __eq__(self, other):
        return str(self) == str(other)

    def __lt__(self, other):
        return str(self) < str(other)

    def set_default_width(self):
        self.width = LANE_TYPES[self.lanetype + self.direction]['width']

    def get_modes(self):
        return set(LANE_TYPES[self.lanetype + self.direction]['modes'])

    def get_modes_as_list(self):
        return LANE_TYPES[self.lanetype + self.direction]['modes']

    def has_mode(self, mode):
        return mode in self.get_modes()

    def get_primary_mode(self):
        return utils.get_nth_element_of_list(list(self.get_modes_as_list()), 0)

    def has_large_vehicles(self):
        return not {MODE_PRIVATE_CARS, MODE_TRANSIT}.isdisjoint(self.get_modes())

    def is_dedicated_for_cycling(self):
        has_cycling = MODE_CYCLING in self.get_modes()
        has_nothing_else_than_cycling_and_walking = self.get_modes().issubset({MODE_CYCLING, MODE_FOOT})
        return has_cycling and has_nothing_else_than_cycling_and_walking

    def is_cycling_infa(self):
        self.lanetype in CYCLING_INFRA

    def get_cycling_vod_factor(self):
        return LANE_TYPES[self.lanetype + self.direction]['cycling_vod']

    def get_order(self):
        return LANE_TYPES[self.lanetype + self.direction]['order']

    def reverse_direction(self):
        self.direction = _get_opposite_direction(self.direction)


class SpaceAllocation(list):
    def __init__(self, lanes, maxspeed=None, hierarchy=None):
        """
        A space allocation

        Parameters
        ----------
        lanes: list[Lane]
        maxspeed: int
        hierarchy: str
        """

        self.extend(lanes)
        # TODO: Move hierarchy and maxspeed from StreetGraoph Edge into SpaceAllocation
        self.maxspeed = maxspeed
        self.hierarchy = hierarchy

    def __str__(self):
        stringified_lanes = list(
            map(lambda item: str(item), self)
        )
        return ' | '.join(stringified_lanes)

    def __eq__(self, other):
        return (
            type(other) is SpaceAllocation
            and str(SpaceAllocation(sorted(self))) == str(SpaceAllocation(sorted(other)))
        )

    def reverse_allocation(self):
        """
        Reverse the order and direction of all lanes
        """
        self.reverse()
        [lane.reverse_direction() for lane in self]
