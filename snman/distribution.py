import math
import networkx as nx
from . import config


def set_given_lanes(street_graph):
    """
    Sets which lanes are given due to external policy definitions
    e.g. dedicated lanes for public transport, bidirectional lanes for cars, etc.
    """

    #TODO: Make the hardcoded policy definition user-configurable

    for id, data in street_graph.edges.items():
        data['given_lanes'] = []

        if data.get('pt_tram') or data.get('pt_bus'):
            data['given_lanes'] += ['m<', 'm>']

        else:

            if data.get('hierarchy') in ['main', 'local']:
                data['given_lanes'] += ['m?']

            elif data.get('hierarchy') in ['dead_end']:
                data['given_lanes'] += ['m-']

        # In case of highways keep all lanes as they are
        if data.get('hierarchy') == 'highway':
            data['given_lanes'] = data.get('ln_desc')