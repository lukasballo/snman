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
DIRECTION_TBD_OPTIONAL = '*'
DIRECTIONS = [
    DIRECTION_FORWARD, DIRECTION_BACKWARD, DIRECTION_BOTH,
    DIRECTION_TBD, DIRECTION_FORWARD_OPTIONAL, DIRECTION_BACKWARD_OPTIONAL, DIRECTION_TBD_OPTIONAL
]

MODE_FOOT = 'foot'
MODE_CYCLING = 'cycling'
MODE_PRIVATE_CARS = 'private_cars'
MODE_TRANSIT = 'transit'
MODES = {MODE_FOOT, MODE_CYCLING, MODE_PRIVATE_CARS, MODE_TRANSIT}
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

LANE_TYPES = {

    LANETYPE_HIGHWAY + DIRECTION_FORWARD:
        {'width': 4.0, 'order': 0, 'cycling_cost_factor': 1.0, 'modes': [MODE_PRIVATE_CARS, MODE_TRANSIT]},
    LANETYPE_HIGHWAY + DIRECTION_BACKWARD:
        {'width': 4.0, 'order': 0, 'cycling_cost_factor': 1.0, 'modes': [MODE_PRIVATE_CARS, MODE_TRANSIT]},
    LANETYPE_HIGHWAY + DIRECTION_BOTH:
        {'width': 6.0, 'order': 0, 'cycling_cost_factor': 1.0, 'modes': [MODE_PRIVATE_CARS, MODE_TRANSIT]},

    LANETYPE_MOTORIZED + DIRECTION_FORWARD:
        {'width': 3.0, 'order': 1, 'cycling_cost_factor': 1.0,
         'modes': [MODE_PRIVATE_CARS, MODE_TRANSIT, MODE_CYCLING, MODE_FOOT]},
    LANETYPE_MOTORIZED + DIRECTION_BACKWARD:
        {'width': 3.0, 'order': 1, 'cycling_cost_factor': 1.0,
         'modes': [MODE_PRIVATE_CARS, MODE_TRANSIT, MODE_CYCLING, MODE_FOOT]},
    LANETYPE_MOTORIZED + DIRECTION_BOTH:
        {'width': 4.5, 'order': 1, 'cycling_cost_factor': 1.0,
         'modes': [MODE_PRIVATE_CARS, MODE_TRANSIT, MODE_CYCLING, MODE_FOOT]},

    # lane to be kept but with direction to be decided yet
    LANETYPE_MOTORIZED + DIRECTION_TBD:
        {'width': 3.0, 'order': 1, 'cycling_cost_factor': 1.0,
         'modes': [MODE_PRIVATE_CARS, MODE_TRANSIT, MODE_CYCLING, MODE_FOOT]},

    LANETYPE_MOTORIZED + DIRECTION_FORWARD_OPTIONAL:
        {'width': 3.0, 'order': 1, 'cycling_cost_factor': 1.0,
         'modes': [MODE_PRIVATE_CARS, MODE_TRANSIT, MODE_CYCLING, MODE_FOOT]},
    LANETYPE_MOTORIZED + DIRECTION_BACKWARD_OPTIONAL:
        {'width': 3.0, 'order': 1, 'cycling_cost_factor': 1.0,
         'modes': [MODE_PRIVATE_CARS, MODE_TRANSIT, MODE_CYCLING, MODE_FOOT]},

    # optional lane with undecided direction
    LANETYPE_MOTORIZED + DIRECTION_TBD_OPTIONAL:
        {'width': 3.0, 'order': 1, 'cycling_cost_factor': 1.0,
         'modes': [MODE_PRIVATE_CARS, MODE_TRANSIT, MODE_CYCLING, MODE_FOOT]},

    LANETYPE_DEDICATED_PT + DIRECTION_FORWARD:
        {'width': 3.0, 'order': 2, 'cycling_cost_factor': 1.0, 'modes': [MODE_TRANSIT]},
    LANETYPE_DEDICATED_PT + DIRECTION_BACKWARD:
        {'width': 3.0, 'order': 2, 'cycling_cost_factor': 1.0, 'modes': [MODE_TRANSIT]},
    LANETYPE_DEDICATED_PT + DIRECTION_BOTH:
        {'width': 4.5, 'order': 2, 'cycling_cost_factor': 1.0, 'modes': [MODE_TRANSIT]},

    LANETYPE_CYCLING_LANE + DIRECTION_FORWARD:
        {'width': 1.5, 'order': 3, 'cycling_cost_factor': 0.5, 'modes': [MODE_CYCLING]},
    LANETYPE_CYCLING_LANE + DIRECTION_BACKWARD:
        {'width': 1.5, 'order': 3, 'cycling_cost_factor': 0.5, 'modes': [MODE_CYCLING]},
    LANETYPE_CYCLING_LANE + DIRECTION_BOTH:
        {'width': 2.0, 'order': 3, 'cycling_cost_factor': 0.5, 'modes': [MODE_CYCLING]},

    LANETYPE_CYCLING_TRACK + DIRECTION_FORWARD:
        {'width': 1.5, 'order': 4, 'cycling_cost_factor': 0.5, 'modes': [MODE_CYCLING]},
    LANETYPE_CYCLING_TRACK + DIRECTION_BACKWARD:
        {'width': 1.5, 'order': 4, 'cycling_cost_factor': 0.5, 'modes': [MODE_CYCLING]},
    LANETYPE_CYCLING_TRACK + DIRECTION_BOTH:
        {'width': 2.5, 'order': 4, 'cycling_cost_factor': 0.5, 'modes': [MODE_CYCLING]},

    LANETYPE_CYCLING_PSEUDO + DIRECTION_FORWARD:
        {'width': 0.0, 'order': 5, 'cycling_cost_factor': 1.0, 'modes': [MODE_CYCLING]},
    LANETYPE_CYCLING_PSEUDO + DIRECTION_BACKWARD:
        {'width': 0.0, 'order': 5, 'cycling_cost_factor': 1.0, 'modes': [MODE_CYCLING]},
    LANETYPE_CYCLING_PSEUDO + DIRECTION_BOTH:
        {'width': 0.0, 'order': 5, 'cycling_cost_factor': 1.0, 'modes': [MODE_CYCLING]},

    LANETYPE_FOOT_CYCLING_MIXED + DIRECTION_FORWARD:
        {'width': 2.5, 'order': 6, 'cycling_cost_factor': 0.5, 'modes': [MODE_CYCLING, MODE_FOOT]},
    LANETYPE_FOOT_CYCLING_MIXED + DIRECTION_BACKWARD:
        {'width': 2.5, 'order': 6, 'cycling_cost_factor': 0.5, 'modes': [MODE_CYCLING, MODE_FOOT]},
    LANETYPE_FOOT_CYCLING_MIXED + DIRECTION_BOTH:
        {'width': 2.5, 'order': 6, 'cycling_cost_factor': 0.5, 'modes': [MODE_CYCLING, MODE_FOOT]},

    LANETYPE_FOOT + DIRECTION_FORWARD:
        {'width': 1.8, 'order': 7, 'cycling_cost_factor': 1.0, 'modes': [MODE_FOOT]},
    LANETYPE_FOOT + DIRECTION_BACKWARD:
        {'width': 1.8, 'order': 7, 'cycling_cost_factor': 1.0, 'modes': [MODE_FOOT]},
    LANETYPE_FOOT + DIRECTION_BOTH:
        {'width': 1.8, 'order': 7, 'cycling_cost_factor': 1.0, 'modes': [MODE_FOOT]},
}

# Which highway=* values represent different infrastructures (primarily) for pedestrians and cyclists
PEDESTRIAN_HIGHWAY_VALUES = {'footway', 'path', 'track', 'pedestrian', 'steps'}
CYCLING_HIGHWAY_VALUES = {'cycleway'}
ACTIVE_HIGHWAY_VALUES = PEDESTRIAN_HIGHWAY_VALUES.union(CYCLING_HIGHWAY_VALUES)

OSM_TAGS = {
    'bridge', 'tunnel', 'layer', 'oneway', 'oneway:bicycle', 'ref', 'name',
    'highway', 'maxspeed', 'service', 'access', 'area',
    'landuse', 'width', 'est_width', 'junction', 'surface',
    'lanes', 'lanes:forward', 'lanes:backward',
    'cycleway', 'cycleway:both', 'cycleway:left', 'cycleway:right',
    'bicycle', 'bicycle:conditional',
    'sidewalk', 'sidewalk:left', 'sidewalk:right', 'foot',
    'psv', 'bus', 'bus:lanes', 'bus:lanes:forward', 'bus:lanes:backward',
    'vehicle:lanes:backward', 'vehicle:lanes:forward',
    'footway'
}

OSM_FILTER = [
    (
        f'["highway"]["area"!~"yes"]["access"!~"private"]'
        f'["highway"!~"abandoned|bridleway|bus_guideway|corridor|elevator|'
        f'escalator|planned|platform|proposed|raceway|construction|footway|pedestrian|steps|path"]'
        f'["service"!~"alley|driveway|emergency_access|parking|parking_aisle|private"]'
        f'["access"!~"no"]'
    ),
    (
        f'["highway"]["bicycle"]["bicycle"!~"no"]["area"!~"yes"]'
        f'["highway"!~"abandoned|bridleway|bus_guideway|corridor|elevator|'
        f'escalator|planned|platform|proposed|raceway|construction"]'
    ),
    (
        f'["highway"]["bicycle:conditional"]["area"!~"yes"]'
        f'["highway"!~"abandoned|bridleway|bus_guideway|corridor|elevator|'
        f'escalator|planned|platform|proposed|raceway|construction"]'
    ),
    (
        f'["highway"]["bus"="yes"]["area"!~"yes"]'
        f'["highway"!~"abandoned|bridleway|bus_guideway|corridor|elevator|'
        f'escalator|planned|platform|proposed|raceway|construction"]'
    ),
    (
        f'["highway"]["psv"="yes"]["area"!~"yes"]'
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
    '_connected_component'
}
