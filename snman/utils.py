from .constants import *
from . import graph_tools

def prepare_graph(G):

    # ensure consistent data types
    for id, edge in G.edges.items():
        maxspeed = edge.get('maxspeed')
        if maxspeed and maxspeed.isdigit():
            edge['maxspeed'] = int(maxspeed)
        else:
            edge['maxspeed'] = -1

    # add layers to graph
    #graph_tools._add_layers_to_nodes(G)