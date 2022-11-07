from snman import constants, graph_tools as lib
import geopandas as gpd
import momepy

print("Reimport running...")

edges = gpd.read_file(config.data_path + 'edges.shp')
edges = edges.to_crs(2056)
# Add columns
edges['dead_end'] = False
edges['processed'] = False
edges['n_neigh_u'] = None
edges['n_neigh_v'] = None
# Convert dataframe to a graph
G = momepy.gdf_to_nx(edges, approach="primal", directed=False)

edge_iterator = G.edges(data=True)
node_iterator = G.nodes

def process_edge(graph, edge):
    #print('Processing Edge: ', edge)
    edge_data = edge[2]
    edge_data['processed'] = True
    neighbors = lib.get_neighbors(graph, edge, dead_ends=False)
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

def process_edges(graph, edges):
    next_edges = []
    for edge in edges:
        result = process_edge(graph, edge[0])
        if (len(result) > 0):
            next_edges.append(result)
    return next_edges

next_edges = []

# Identify dead end nodes
for node in node_iterator:
    degree = G.degree(node)
    if(degree == 1):
        adjacent_edge = list(G.edges(nbunch=node, data=True))[0]
        result = process_edge(G, adjacent_edge)
        if (len(result) > 0):
            next_edges.append(result)

for step in range(0, 100):
    next_edges = process_edges(G, next_edges)



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
