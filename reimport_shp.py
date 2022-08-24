import osmnx_ebc as ox
import networkx as nx
import math
import config
import re
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import momepy

print("Reimport running...")

edges = gpd.read_file(config.data_path + 'edges.shp')
edges = edges.to_crs(2056)
edges['dead_end'] = False
edges['processed'] = False
edges['n_neigh_u'] = None
edges['n_neigh_v'] = None
G = momepy.gdf_to_nx(edges, approach="primal", directed=False)

edge_iterator = G.edges(data=True)
node_iterator = G.nodes

def get_edge_hash(edge):
    edge_data = edge[2]
    return str(edge_data['u']) + str(edge_data['v']) + str(edge_data['__key'])

def remove_edge_from_list(edges, edge_to_remove, dead_ends=True):
    edges_cleaned = []
    for idx, candidate in enumerate(edges):
        if(get_edge_hash(candidate) != get_edge_hash(edge_to_remove)
                and not(dead_ends == False and candidate[2]['dead_end'] == True)):
            edges_cleaned.append(candidate)

    return edges_cleaned

def get_neighbors(edge, dead_ends=True):
    adjacent_nodes = edge[0:2]
    u_neighbors = list(G.edges(nbunch=adjacent_nodes[0], data=True))
    v_neighbors = list(G.edges(nbunch=adjacent_nodes[1], data=True))
    u_neighbors = remove_edge_from_list(u_neighbors, edge, dead_ends=dead_ends)
    v_neighbors = remove_edge_from_list(v_neighbors, edge, dead_ends=dead_ends)
    return [u_neighbors, v_neighbors, u_neighbors + v_neighbors]

def process_edge(edge):
    #print('Processing Edge: ', edge)
    edge_data = edge[2]
    edge_data['processed'] = True
    neighbors = get_neighbors(edge, dead_ends=False)
    edge_data['n_neigh_u'] = len(neighbors[0])
    edge_data['n_neigh_v'] = len(neighbors[1])
    next_edges = []
    if min((len(neighbors[0])), (len(neighbors[1]))) == 0:
        edge_data['dead_end'] = True
    for neighbor in neighbors[2]:
        #print('nb: ', neighbor)
        neighbor_data = neighbor[2]
        if neighbor_data['processed'] == False:
            #process_edge(neighbor)
            next_edges.append(neighbor)
    # return next edges to be processed
    return next_edges

def process_edges(edges):
    next_edges = []
    for edge in edges:
        #print(edge)
        #print(type(edge[0]))
        #print(edge[0])
        #print('PROCESSING: ', edge)
        result = process_edge(edge[0])
        if (len(result) > 0):
            next_edges.append(result)
    return next_edges

next_edges = []

# Identify dead end nodes
for node in node_iterator:
    degree = G.degree(node)
    if(degree == 1):
        adjacent_edge = list(G.edges(nbunch=node, data=True))[0]
        result = process_edge(adjacent_edge)
        if (len(result) > 0):
            next_edges.append(result)

for step in range(0,100):
    next_edges = process_edges(next_edges)



"""
# Plot
positions = {n: [n[0], n[1]] for n in list(G.nodes)}
nx.draw(G, node_size=2)
plt.show()
"""

processed_nodes = momepy.nx_to_gdf(G)[0]
processed_nodes.to_file(config.data_path + 'processed_nodes.shp')

processed_edges = momepy.nx_to_gdf(G)[1]
processed_edges.to_file(config.data_path + 'processed_edges.shp')

print("Done")
