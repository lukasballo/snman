from snman import osmnx as ox

ox.config(log_console=False, use_cache=True)

DIRECTION_FORWARD = '>'
DIRECTION_BACKWARD = '<'
DIRECTION_BOTH = '-'
DIRECTION_TBD = '?'

LANETYPE_MOTORIZED = 'M'            # A normal lane accessible to car, public transport, and cyclists
LANETYPE_DEDICATED_PT = 'T'         # Only for public transport
LANETYPE_CYCLING_TRACK = 'P'        # Only for cyclists
LANETYPE_CYCLING_LANE = 'L'         # Advisory cycling lane, in some cases also used by other traffic
LANETYPE_CYCLING_PSEUDO = 'S'       # Allowing cyclists to travel despite missing separate infrastructure,
                                    # e.g. in the opposite direction of one-way streets
LANETYPE_FOOT_CYCLING_MIXED = 'X'   # Mixed, for cyclists and pedestrians
LANETYPE_FOOT = 'F'                 # Pedestrians only

LANETYPE_ORDER = ['T', 'M', 'L', 'P', 'S', 'F']
DIRECTION_ORDER = ['<', '-', '>']

KEY_LANES_DESCRIPTION = 'ln_desc'   # under which key of each edge is the existing lane configuration
KEY_GIVEN_LANES_DESCRIPTION = 'given_lanes'
KEY_REVERSED = '_reversed'          # which key tells if the edge has been reversed


# > = unidirectional lane with defined direction
# ? = unidirectional lane with direction yet to be defined
# - = bidirectional lane (e.g. local streets with light traffic or cycling paths without lanes)
DEFAULT_LANE_WIDTHS = {

    LANETYPE_MOTORIZED + DIRECTION_FORWARD: 3,
    LANETYPE_MOTORIZED + DIRECTION_BACKWARD: 3,
    LANETYPE_MOTORIZED + DIRECTION_TBD: 3,
    LANETYPE_MOTORIZED + DIRECTION_BOTH: 4.5,

    LANETYPE_DEDICATED_PT + DIRECTION_FORWARD: 3,
    LANETYPE_DEDICATED_PT + DIRECTION_BACKWARD: 3,
    LANETYPE_DEDICATED_PT + DIRECTION_TBD: 3,
    LANETYPE_DEDICATED_PT + DIRECTION_BOTH: 4.5,

    LANETYPE_CYCLING_LANE + DIRECTION_FORWARD: 1.3,
    LANETYPE_CYCLING_LANE + DIRECTION_BACKWARD: 1.3,
    LANETYPE_CYCLING_LANE + DIRECTION_TBD: 1.3,
    LANETYPE_CYCLING_LANE + DIRECTION_BOTH: 2.6,

    LANETYPE_CYCLING_TRACK + DIRECTION_FORWARD: 1.8,
    LANETYPE_CYCLING_TRACK + DIRECTION_BACKWARD: 1.8,
    LANETYPE_CYCLING_TRACK + DIRECTION_TBD: 1.8,
    LANETYPE_CYCLING_TRACK + DIRECTION_BOTH: 2,

    LANETYPE_CYCLING_PSEUDO + DIRECTION_FORWARD: 0,
    LANETYPE_CYCLING_PSEUDO + DIRECTION_BACKWARD: 0,
    LANETYPE_CYCLING_PSEUDO + DIRECTION_TBD: 0,
    LANETYPE_CYCLING_PSEUDO + DIRECTION_BOTH: 0,

    LANETYPE_FOOT_CYCLING_MIXED + DIRECTION_FORWARD: 2.5,
    LANETYPE_FOOT_CYCLING_MIXED + DIRECTION_BACKWARD: 2.5,
    LANETYPE_FOOT_CYCLING_MIXED + DIRECTION_TBD: 2.5,
    LANETYPE_FOOT_CYCLING_MIXED + DIRECTION_BOTH: 2.5,

    LANETYPE_FOOT + DIRECTION_FORWARD: 1.8,
    LANETYPE_FOOT + DIRECTION_BACKWARD: 1.8,
    LANETYPE_FOOT + DIRECTION_TBD: 1.8,
    LANETYPE_FOOT + DIRECTION_BOTH: 1.8,

}


OSM_TAGS = ['bridge', 'tunnel', 'layer', 'oneway', 'ref', 'name',
                    'highway', 'maxspeed', 'service', 'access', 'area',
                    'landuse', 'width', 'est_width', 'junction', 'surface',
                    'lanes', 'lanes:forward', 'lanes:backward',
                    'cycleway', 'cycleway:both', 'cycleway:left', 'cycleway:right',
                    'bicycle', 'bicycle:conditional',
                    'sidewalk', 'sidewalk:left', 'sidewalk:right', 'foot',
                    'psv', 'bus', 'bus:lanes:forward', 'bus:lanes:backward',
                    'vehicle:lanes:backward', 'vehicle:lanes:forward',
                    'footway',
            ]

OSM_FILTER = [
    (
        f'["highway"]["area"!~"yes"]["access"!~"private"]'
        f'["highway"!~"abandoned|bridleway|bus_guideway|construction|corridor|elevator|'
        f'escalator|planned|platform|proposed|raceway"]'
        f'["service"!~"alley|driveway|emergency_access|parking|parking_aisle|private"]'
        f'["access"!~"no"]'
    ),
    (
        f'["highway"]["bicycle"]["bicycle"!~"no"]["area"!~"yes"]'
    ),
    (
        f'["highway"]["bicycle:conditional"]["area"!~"yes"]'
    ),
    (
        f'["highway"]["bus"="yes"]["area"!~"yes"]'
    ),
    (
        f'["highway"]["psv"="yes"]["area"!~"yes"]'
    ),
    # pedestrian paths incl. those mapped as area
    (
        f'["highway"="footway|pedestrian"]'
    )
]

CRS = 'epsg:2056'