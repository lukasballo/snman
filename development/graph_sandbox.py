import networkx as nx
import matplotlib.pyplot as plt
import snman

keep_all_streets = True

G = nx.DiGraph()
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

G.add_edge(4,5)
G.add_edge(5,4)

wccs = nx.strongly_connected_components(G)
for wcc, nodes in enumerate(wccs):
    subgraph = nx.subgraph(G, nodes)
    nx.set_edge_attributes(subgraph, wcc, 'wcc')

nx.set_node_attributes(G, wcc, 'wcc')

edge_labels = nx.get_edge_attributes(G, "wcc")

# Plot
nx.draw_networkx(G, pos=pos)
#nx.draw_networkx_labels(G, pos, font_size=7, labels=wcc_ids, verticalalignment='top')
nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels)
plt.show()
