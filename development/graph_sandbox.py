import networkx as nx
import matplotlib.pyplot as plt
import snman

keep_all_streets = True

G = nx.MultiDiGraph()
G.add_nodes_from([0,1,2,3,4,5])

pos = {
    0: (0,1),
    1: (1,1),
    2: (2,2),
    3: (2,0),
    4: (0,2),
    5: (1,2)
}

# Add edges with arbitrary directions
G.add_edge(0,1)
G.add_edge(1,2)
G.add_edge(2,3)
G.add_edge(1,3)

#G.add_edge(4,5)

wccs = nx.weakly_connected_components(G)
wcc_ids = {}
for wcc, nodes in enumerate(wccs):
    for node in nodes:
        wcc_ids[node] = wcc


nx.set_node_attributes(G, wcc, 'wcc')

# Plot
nx.draw_networkx(G, pos=pos)
nx.draw_networkx_labels(G, pos, font_size=7, labels=wcc_ids, verticalalignment='top')
plt.show()
