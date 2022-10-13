from snman import osmnx as ox

ox.config(log_console=False, use_cache=True)

# > = unidirectional lane with defined direction
# ? = unidirectional lane with direction yet to be defined
# - = bidirectional lane (e.g. local streets with light traffic or cycling paths without lanes)
default_lane_widths_m = {
    'm>': 3,
    'm<': 3,
    'm?': 3,
    'm-': 4.5,
    'c>': 1.3,
    'c<': 1.3,
    'c?': 1.3,
    'c-': 2.6
}

lane_descriptions = {
    'm' : 'motorized traffic',
    'c' : 'cycling'
}