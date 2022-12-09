import networkx as nx
import matplotlib.pyplot as plt
import snman

keep_all_streets = True

G = nx.MultiDiGraph()
G.add_nodes_from([1,2,3,4,5,6,7,8])

# Add edges with arbitrary directions
G.add_edge(2,1)
G.add_edge(2,3)
G.add_edge(3,4)
G.add_edge(1,4)
G.add_edge(3,5)

G.add_edge(5,6, fixed=True)
G.add_edge(6,7)
G.add_edge(7,8)
G.add_edge(5,8)

G.add_edge(2,7)


# Add complementary edges
for i, data in G.edges.items():
    G.add_edge(i[1], i[0], i[2])

def opposite_direction_exists(G, *edge_id):
    return G.has_edge(edge_id[1], edge_id[0])

# Try removing edges
while True:
    # Calculate betweenness centrality
    bc = nx.edge_betweenness_centrality(G)
    nx.set_edge_attributes(G, bc, 'bc')
    edges = [edge for edge in list(G.edges.items()) if edge[1].get('fixed', False) == False]

    # Finish the loop if no edge candidates exist
    if len(edges) == 0:
        break
    edges = sorted(edges, key=lambda x: x[1]['bc'])
    edge = list(edges)[0]
    edge_id = edge[0]

    # Check if the new graph is strongly connected
    H = G.copy()
    H.remove_edge(*edge_id)
    if nx.is_strongly_connected(H) and (opposite_direction_exists(G, *edge_id) or not keep_all_streets):
        G.remove_edge(*edge_id)
    else:
        nx.set_edge_attributes(G, {edge_id: {'fixed': True}})

#nx.draw_networkx(G)
#plt.show()
