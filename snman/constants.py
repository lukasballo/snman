import osmnx as ox

ox.config(log_console=False, use_cache=True)

DEFAULT_INTERSECTION_TOLERANCE = 10
DEFAULT_SIMPLIFICATION_RADIUS = 35

DIRECTION_FORWARD = '>'
DIRECTION_BACKWARD = '<'
DIRECTION_BOTH = '-'
DIRECTION_TBD = '?'
DIRECTION_FORWARD_OPTIONAL = ')'
DIRECTION_BACKWARD_OPTIONAL = '('
DIRECTION_BOTH_OPTIONAL = '/'
DIRECTION_TBD_OPTIONAL = '*'
DIRECTIONS = {
    DIRECTION_FORWARD, DIRECTION_BACKWARD, DIRECTION_BOTH, DIRECTION_TBD,
    DIRECTION_FORWARD_OPTIONAL, DIRECTION_BACKWARD_OPTIONAL, DIRECTION_BOTH_OPTIONAL, DIRECTION_TBD_OPTIONAL
}
TENTATIVE_DIRECTIONS = {
    DIRECTION_TBD, DIRECTION_TBD_OPTIONAL, DIRECTION_BOTH_OPTIONAL,
    DIRECTION_FORWARD_OPTIONAL, DIRECTION_BACKWARD_OPTIONAL
}

ALTERNATIVE_DIRECTIONS = {
    DIRECTION_FORWARD: {
        DIRECTION_BOTH, DIRECTION_TBD, DIRECTION_BOTH_OPTIONAL, DIRECTION_FORWARD_OPTIONAL, DIRECTION_TBD_OPTIONAL
    },
    DIRECTION_BACKWARD: {
        DIRECTION_BOTH, DIRECTION_TBD, DIRECTION_BOTH_OPTIONAL, DIRECTION_BACKWARD_OPTIONAL, DIRECTION_TBD_OPTIONAL
    }
}

MODE_FOOT = 'foot'
MODE_CYCLING = 'cycling'
MODE_PRIVATE_CARS = 'private_cars'
MODE_TRANSIT = 'transit'
MODE_CAR_PARKING = 'car_parking'
MODE_NON_TRAFFIC = 'non_traffic'
MODES = {MODE_FOOT, MODE_CYCLING, MODE_PRIVATE_CARS, MODE_TRANSIT, MODE_CAR_PARKING, MODE_NON_TRAFFIC}
ACTIVE_MODES = {MODE_FOOT, MODE_CYCLING}
MOTORIZED_MODES = {MODE_PRIVATE_CARS, MODE_TRANSIT}

LANETYPE_MOTORIZED = 'M'            # A normal lane accessible to car, public transport, and cyclists
LANETYPE_HIGHWAY = 'H'              # A highway lane
LANETYPE_DEDICATED_PT = 'T'         # Only for public transport
LANETYPE_CYCLING_TRACK = 'P'        # Only for cyclists and separated from other modes
LANETYPE_CYCLING_LANE = 'L'         # Advisory cycling lane, in some cases also used by other traffic
LANETYPE_CYCLING_PSEUDO = 'S'       # Contraflow cycling in one-way streets without cycling infrastructure
LANETYPE_FOOT_CYCLING_MIXED = 'X'   # Mixed, for cyclists and pedestrians
LANETYPE_FOOT = 'F'                 # Pedestrians only
LANETYPE_PARKING_PARALLEL = 'R'     # On-street parking
LANETYPE_PARKING_PERPENDICULAR = 'N'    # On-street parking
LANETYPE_PARKING_DIAGONAL = 'D'     # On-street parking
LANETYPE_NON_TRAFFIC = 'Z'          # No traffic, e.g., greenery or community spaces

PARALLEL_PARKING_CAR_LENGTH = 5

STATUS_FIXED = '*'
STATUS_ONE_DIRECTION_MANDATORY = '%'
STATUS_OPTIONAL = '/'
STATUS_BY_NEED = '!'

# For assessing cycling quality, from best to worst
# TODO: To be replaced with ranking according to cycling comfort factors
CYCLING_QUALITY_HIERARCHY = [
    LANETYPE_CYCLING_TRACK,
    LANETYPE_FOOT_CYCLING_MIXED,
    LANETYPE_CYCLING_LANE,
    LANETYPE_MOTORIZED,
    LANETYPE_CYCLING_PSEUDO
]

CYCLING_INFRA = [
    LANETYPE_CYCLING_TRACK,
    LANETYPE_CYCLING_LANE,
    LANETYPE_FOOT_CYCLING_MIXED
]

# For sorting lanes in the crosssection
# TODO: To be replaced by sorting heuristics
LANETYPE_ORDER = ['H', 'T', 'M', 'L', 'P', 'S', 'F']
DIRECTION_ORDER = ['<', '-', '>']

KEY_LANES_DESCRIPTION = 'ln_desc'               # under which key is the existing lane configuration
KEY_LANES_DESCRIPTION_AFTER = 'ln_desc_after'   # under which key is the lane configuration after rebuilding
KEY_GIVEN_LANES_DESCRIPTION = 'given_lanes'
KEY_REVERSED = '_reversed'          # which key tells if the edge has been reversed

KEY_SENSORS_FORWARD = 'sensors_forward'
KEY_SENSORS_BACKWARD = 'sensors_backward'

LANE_TYPES = {

    LANETYPE_HIGHWAY + DIRECTION_FORWARD:
        {'width': 4.0, 'order': 0, 'cycling_vod': 0, 'modes': [MODE_PRIVATE_CARS, MODE_TRANSIT]},
    LANETYPE_HIGHWAY + DIRECTION_BACKWARD:
        {'width': 4.0, 'order': 0, 'cycling_vod': 0, 'modes': [MODE_PRIVATE_CARS, MODE_TRANSIT]},
    LANETYPE_HIGHWAY + DIRECTION_BOTH:
        {'width': 6.0, 'order': 0, 'cycling_vod': 0, 'modes': [MODE_PRIVATE_CARS, MODE_TRANSIT]},

    LANETYPE_MOTORIZED + DIRECTION_FORWARD:
        {'width': 3.0, 'order': 1, 'cycling_vod': 0,
         'modes': [MODE_PRIVATE_CARS, MODE_TRANSIT, MODE_CYCLING, MODE_FOOT]},
    LANETYPE_MOTORIZED + DIRECTION_BACKWARD:
        {'width': 3.0, 'order': 1, 'cycling_vod': 0,
         'modes': [MODE_PRIVATE_CARS, MODE_TRANSIT, MODE_CYCLING, MODE_FOOT]},
    LANETYPE_MOTORIZED + DIRECTION_BOTH:
        {'width': 4.5, 'order': 1, 'cycling_vod': 0,
         'modes': [MODE_PRIVATE_CARS, MODE_TRANSIT, MODE_CYCLING, MODE_FOOT]},

    # lane to be kept but with direction to be decided yet
    LANETYPE_MOTORIZED + DIRECTION_TBD:
        {'width': 3.0, 'order': 1, 'cycling_vod': 0,
         'modes': [MODE_PRIVATE_CARS, MODE_TRANSIT, MODE_CYCLING, MODE_FOOT]},

    LANETYPE_MOTORIZED + DIRECTION_FORWARD_OPTIONAL:
        {'width': 3.0, 'order': 1, 'cycling_vod': 0,
         'modes': [MODE_PRIVATE_CARS, MODE_TRANSIT, MODE_CYCLING, MODE_FOOT]},
    LANETYPE_MOTORIZED + DIRECTION_BACKWARD_OPTIONAL:
        {'width': 3.0, 'order': 1, 'cycling_vod': 0,
         'modes': [MODE_PRIVATE_CARS, MODE_TRANSIT, MODE_CYCLING, MODE_FOOT]},
    LANETYPE_MOTORIZED + DIRECTION_BOTH_OPTIONAL:
        {'width': 4.5, 'order': 1, 'cycling_vod': 0,
         'modes': [MODE_PRIVATE_CARS, MODE_TRANSIT, MODE_CYCLING, MODE_FOOT]},
    # optional lane with undecided direction
    LANETYPE_MOTORIZED + DIRECTION_TBD_OPTIONAL:
        {'width': 3.0, 'order': 1, 'cycling_vod': 0,
         'modes': [MODE_PRIVATE_CARS, MODE_TRANSIT, MODE_CYCLING, MODE_FOOT]},

    LANETYPE_DEDICATED_PT + DIRECTION_FORWARD:
        {'width': 3.0, 'order': 2, 'cycling_vod': 0, 'modes': [MODE_TRANSIT]},
    LANETYPE_DEDICATED_PT + DIRECTION_BACKWARD:
        {'width': 3.0, 'order': 2, 'cycling_vod': 0, 'modes': [MODE_TRANSIT]},
    LANETYPE_DEDICATED_PT + DIRECTION_BOTH:
        {'width': 4.5, 'order': 2, 'cycling_vod': 0, 'modes': [MODE_TRANSIT]},
    LANETYPE_DEDICATED_PT + DIRECTION_FORWARD_OPTIONAL:
        {'width': 3.0, 'order': 2, 'cycling_vod': 0, 'modes': [MODE_TRANSIT]},
    LANETYPE_DEDICATED_PT + DIRECTION_BACKWARD_OPTIONAL:
        {'width': 3.0, 'order': 2, 'cycling_vod': 0, 'modes': [MODE_TRANSIT]},

    LANETYPE_CYCLING_LANE + DIRECTION_FORWARD:
        {'width': 1.5, 'order': 3, 'cycling_vod': -0.51, 'modes': [MODE_CYCLING]},
    LANETYPE_CYCLING_LANE + DIRECTION_BACKWARD:
        {'width': 1.5, 'order': 3, 'cycling_vod': -0.51, 'modes': [MODE_CYCLING]},
    LANETYPE_CYCLING_LANE + DIRECTION_FORWARD_OPTIONAL:
        {'width': 1.5, 'order': 3, 'cycling_vod': -0.51, 'modes': [MODE_CYCLING]},
    LANETYPE_CYCLING_LANE + DIRECTION_BACKWARD_OPTIONAL:
        {'width': 1.5, 'order': 3, 'cycling_vod': -0.51, 'modes': [MODE_CYCLING]},
    LANETYPE_CYCLING_LANE + DIRECTION_BOTH:
        {'width': 2.0, 'order': 3, 'cycling_vod': -0.51, 'modes': [MODE_CYCLING]},

    LANETYPE_CYCLING_TRACK + DIRECTION_FORWARD:
        {'width': 1.5, 'order': 4, 'cycling_vod': -0.51, 'modes': [MODE_CYCLING]},
    LANETYPE_CYCLING_TRACK + DIRECTION_BACKWARD:
        {'width': 1.5, 'order': 4, 'cycling_vod': -0.51, 'modes': [MODE_CYCLING]},
    LANETYPE_CYCLING_TRACK + DIRECTION_BOTH:
        {'width': 2.5, 'order': 4, 'cycling_vod': -0.51, 'modes': [MODE_CYCLING]},

    LANETYPE_CYCLING_PSEUDO + DIRECTION_FORWARD:
        {'width': 0.0, 'order': 5, 'cycling_vod': 0, 'modes': [MODE_CYCLING]},
    LANETYPE_CYCLING_PSEUDO + DIRECTION_BACKWARD:
        {'width': 0.0, 'order': 5, 'cycling_vod': 0, 'modes': [MODE_CYCLING]},
    LANETYPE_CYCLING_PSEUDO + DIRECTION_BOTH:
        {'width': 0.0, 'order': 5, 'cycling_vod': 0, 'modes': [MODE_CYCLING]},

    LANETYPE_FOOT_CYCLING_MIXED + DIRECTION_FORWARD:
        {'width': 2.5, 'order': 6, 'cycling_vod': -0.51, 'modes': [MODE_CYCLING, MODE_FOOT]},
    LANETYPE_FOOT_CYCLING_MIXED + DIRECTION_BACKWARD:
        {'width': 2.5, 'order': 6, 'cycling_vod': -0.51, 'modes': [MODE_CYCLING, MODE_FOOT]},
    LANETYPE_FOOT_CYCLING_MIXED + DIRECTION_BOTH:
        {'width': 2.5, 'order': 6, 'cycling_vod': -0.51, 'modes': [MODE_CYCLING, MODE_FOOT]},

    LANETYPE_FOOT + DIRECTION_FORWARD:
        {'width': 1.8, 'order': 7, 'cycling_vod': 0, 'modes': [MODE_FOOT]},
    LANETYPE_FOOT + DIRECTION_BACKWARD:
        {'width': 1.8, 'order': 7, 'cycling_vod': 0, 'modes': [MODE_FOOT]},
    LANETYPE_FOOT + DIRECTION_BOTH:
        {'width': 1.8, 'order': 7, 'cycling_vod': 0, 'modes': [MODE_FOOT]},

    LANETYPE_PARKING_PARALLEL + DIRECTION_BOTH:
        {'width': 2, 'order': 8, 'cycling_vod': 0, 'modes': [MODE_CAR_PARKING]},
    LANETYPE_PARKING_PARALLEL + DIRECTION_BOTH_OPTIONAL:
        {'width': 2, 'order': 8, 'cycling_vod': 0, 'modes': [MODE_CAR_PARKING]},
    LANETYPE_PARKING_DIAGONAL + DIRECTION_BOTH:
        {'width': 4.5, 'order': 8, 'cycling_vod': 0, 'modes': [MODE_CAR_PARKING]},
    LANETYPE_PARKING_PERPENDICULAR + DIRECTION_BOTH:
        {'width': 6, 'order': 8, 'cycling_vod': 0, 'modes': [MODE_CAR_PARKING]},
    LANETYPE_PARKING_PARALLEL + DIRECTION_FORWARD:
        {'width': 2, 'order': 8, 'cycling_vod': 0, 'modes': [MODE_CAR_PARKING]},
    LANETYPE_PARKING_PARALLEL + DIRECTION_BACKWARD:
        {'width': 2, 'order': 8, 'cycling_vod': 0, 'modes': [MODE_CAR_PARKING]},

    LANETYPE_NON_TRAFFIC + DIRECTION_FORWARD:
        {'width': 0.0, 'order': 9, 'cycling_vod': 0, 'modes': [MODE_NON_TRAFFIC]},
    LANETYPE_NON_TRAFFIC + DIRECTION_BACKWARD:
        {'width': 0.0, 'order': 9, 'cycling_vod': 0, 'modes': [MODE_NON_TRAFFIC]},

}


def CYCLING_SLOPE_VOD(slope):
    """
    Returns the VoD parameter for a given slope. See this paper, page 8:

    Meister, A., M. Felder, B. Schmid and K.W. Axhausen (2023)
    Route choice modeling for cyclists on urban networks,
    *Transportation Research Part A: Policy and Practice*, **173**, 103723.

    Parameters
    ----------
    slope: float

    Returns
    -------
    int

    """
    slope = float(slope)
    if slope < 0.03:
        return 0
    if 0.02 <= slope < 0.06:
        return +0.55
    elif 0.06 <= slope < 0.10:
        return +3.11
    elif 0.10 <= slope:
        return +4.33


# Which highway=* values represent different infrastructures (primarily) for pedestrians and cyclists
PEDESTRIAN_HIGHWAY_VALUES = {'footway', 'path', 'track', 'pedestrian', 'steps'}
CYCLING_HIGHWAY_VALUES = {'cycleway'}
ACTIVE_HIGHWAY_VALUES = PEDESTRIAN_HIGHWAY_VALUES.union(CYCLING_HIGHWAY_VALUES)


OSM_HIGHWAY_VALUES = {
    'motorway':         {'level': None},
    'motorway_link':    {'level': None},
    'trunk':            {'level': None},
    'trunk_link':       {'level': None},
    'primary':          {'level': None},
    'primary_link':     {'level': None},
    'secondary':        {'level': None},
    'secondary_link':   {'level': None},
    'tertiary':         {'level': None},
    'tertiary_link':    {'level': None},
    'service':          {'level': None},
    'busway':           {'level': None},
    'unclassified':     {'level': None},
    'road':             {'level': None},
    'residential':      {'level': None},
    'living_street':    {'level': None},
    'track':            {'level': None},
    'path':             {'level': None},
    'cycleway':         {'level': None},
    'footway':          {'level': None},
    'pedestrian':       {'level': None},
    'steps':            {'level': None},
    'rest_area':        {'level': None},
}

level = 0
for val, attributes in OSM_HIGHWAY_VALUES.items():
    attributes['level'] = level
    level += 1


OSM_TAGS = {
    'bridge', 'tunnel', 'layer', 'oneway', 'oneway:bicycle', 'ref', 'name',
    'highway', 'maxspeed', 'service', 'access', 'area',
    'landuse', 'width', 'est_width', 'junction', 'surface',
    'lanes', 'lanes:forward', 'lanes:backward',
    'cycleway', 'cycleway:both', 'cycleway:left', 'cycleway:right',
    'bicycle', 'bicycle:conditional',
    'sidewalk', 'sidewalk:left', 'sidewalk:right', 'foot',
    'psv', 'bus', 'bus:lanes', 'bus:lanes:forward', 'bus:lanes:backward',
    'vehicle:lanes', 'vehicle:lanes:backward', 'vehicle:lanes:forward',
    'busway', 'busway:right', 'busway:left',
    'footway',
    'parking:left', 'parking:right', 'parking:both'
}

OSM_FILTER = [
    # regular roads
    (
        f'["highway"]["area"!~"yes"]["access"!~"private"]'
        f'["highway"!~"abandoned|bridleway|bus_guideway|corridor|elevator|'
        f'escalator|planned|platform|proposed|raceway|construction|footway|pedestrian|steps|path|service"]'
        f'["service"!~"alley|driveway|emergency_access|parking|parking_aisle|private"]'
        f'["access"!~"no"]'
    ),
    # bus roads marked as bus=yes
    (
        f'["highway"]["bus"="yes"]["area"!~"yes"]'
        f'["highway"!~"abandoned|bridleway|bus_guideway|corridor|elevator|'
        f'escalator|planned|platform|proposed|raceway|construction"]'
    ),
    # bus roads marked as psv=yes
    (
        f'["highway"]["psv"="yes"]["area"!~"yes"]'
        f'["highway"!~"abandoned|bridleway|bus_guideway|corridor|elevator|'
        f'escalator|planned|platform|proposed|raceway|construction"]'
    ),
    # cycling paths
    (
        f'["highway"]["bicycle"]["bicycle"!~"no"]["area"!~"yes"]'
        f'["highway"!~"abandoned|bridleway|bus_guideway|corridor|elevator|'
        f'escalator|planned|platform|proposed|raceway|construction"]'
    ),
    # conditional cycling paths
    (
        f'["highway"]["bicycle:conditional"]["area"!~"yes"]'
        f'["highway"!~"abandoned|bridleway|bus_guideway|corridor|elevator|'
        f'escalator|planned|platform|proposed|raceway|construction"]'
    ),
    # pedestrian paths incl. those mapped as area
    #(
    #    f'["highway"="footway|pedestrian"]'
    #    f'["highway"!~"abandoned|bridleway|bus_guideway|corridor|elevator|'
    #    f'escalator|planned|platform|proposed|raceway|construction"]'
    #)
]

DEFAULT_CRS = 'epsg:2056'

EXPORT_EDGE_COLUMNS = [
    'grade',
    'ln_desc',
    'width_total_m',
    'width_motorized_m',
    'length',
    'n_lanes_motorized',
    'cycling_forward',
    'cycling_backward',
    'maxspeed',
    'highway',
    'hierarchy',
    'layer',
    'adt_max_forward',
    'adt_max_backward',
    'adt_avg_forward',
    'adt_avg_backward'
]

EXPORT_OSM_TAGS = {
    'highway', 'maxspeed',
    'lanes', 'lanes:forward', 'lanes:backward', 'lanes:both_ways', 'oneway',
    'cycleway', 'cycleway:lane', 'cycleway:left', 'cycleway:left:lane', 'cycleway:right', 'cycleway:right:lane',
    'bus:lanes:backward', 'bus:lanes:forward', 'vehicle:lanes:backward', 'vehicle:lanes:forward',
    '_connected_component', 'sensors_forward', 'sensors_backward',
    'length'
}
