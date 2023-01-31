from .constants import *
from . import graph_tools

def prepare_graph(G):

    # ensure consistent data types
    for id, edge in G.edges.items():

        maxspeed = edge.get('maxspeed', '')
        edge['maxspeed'] = int(maxspeed) if maxspeed.isdigit() else -1

        layer = edge.get('layer', '')
        edge['layer'] = int(layer) if layer.isdigit() else 0



    # add layers to graph
    #graph_tools._add_layers_to_nodes(G)