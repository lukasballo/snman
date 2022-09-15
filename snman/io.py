import networkx as nx
from . import osmnx as ox
from . import config

def export_to_shp(street_graph):
    """
    Exports street graph as a shape file

    Params
    ------
    street_graph : nx.MultiGraph
    """
    nodes, edges = ox.graph_to_gdfs(street_graph)
    edges.columns = [column[0:10] for column in edges.columns]
    #edges['ln_desc'] = edges['ln_desc'].apply(lambda ln_desc: print(type(ln_desc)))
    edges['ln_desc'] = edges['ln_desc'].apply(lambda ln_desc: ' | '.join(ln_desc))
    """
    for ln_desc in edges['ln_desc']:
        print(ln_desc)
        ' | '.join(ln_desc)
    """
    #edges['ln_desc'] = edges['ln_desc'].apply(lambda ln_desc: '*')
    edges.to_file(config.data_path + 'edges.shp', schema=None)
