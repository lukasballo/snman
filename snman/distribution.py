import math
import networkx as nx
from . import config, lanes


def set_given_lanes(street_graph):
    """
    Sets which lanes are given due to external policy definitions
    e.g. dedicated lanes for public transport, bidirectional lanes for cars, etc.
    """

    #TODO: Make the (currently hardcoded) policy definition user-configurable
    #TODO: Add support for dedicated transit lanes

    for id, data in street_graph.edges.items():
        data['given_lanes'] = []

        if data.get('pt_tram') or data.get('pt_bus'):
            data['given_lanes'] += ['m<', 'm>']

        else:

            if data.get('hierarchy') in ['1_main', '2_local']:
                data['given_lanes'] += ['m?']

            elif data.get('hierarchy') in ['3_dead_end']:
                data['given_lanes'] += ['m-']

        # In case of highways keep all lanes as they are
        if data.get('hierarchy') == '0_highway':
            data['given_lanes'] = data.get('ln_desc')

def create_given_lanes_graph(street_graph):
    """
    Returns a directed graph of given (mandatory) lanes. Lanes with changeable direction are marked with an attribute
    """
    given_lanes_graph = nx.DiGraph()
    given_lanes_graph.graph['crs'] = street_graph.graph['crs']
    given_lanes_graph.add_nodes_from(street_graph.nodes.items())

    for id, data in street_graph.edges.items():
        u = id[0]
        v = id[1]
        given_lanes = data.get('given_lanes',[])

        for lane in given_lanes:
            lane_properties = lanes._get_lane_properties(lane)

            if lane_properties['direction'] in ['>', '-']:
                given_lanes_graph.add_edge(u, v, fixed_direction=True)

            if lane_properties['direction'] in ['<', '-']:
                given_lanes_graph.add_edge(v, u, fixed_direction=True)

            if lane_properties['direction'] in ['?']:
                given_lanes_graph.add_edge(u, v, fixed_direction=False)

    return given_lanes_graph