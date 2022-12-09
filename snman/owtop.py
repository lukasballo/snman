import networkx as nx
import matplotlib.pyplot as plt


def link_elimination(G, keep_all_streets=True):

    # Get the giant weakly connected component (remove any unconnected parts)
    gcc = sorted(nx.weakly_connected_components(G), key=len, reverse=True)[0]
    G = G.subgraph(gcc).copy()

    # Add complementary edges
    for i, data in G.edges.items():
        if data.get('fixed', False) == False:
            G.add_edge(i[1], i[0], fixed=False)

    print('Initialized graph has ',len(G.nodes),' and ',len(G.edges),' edges')

    if not nx.is_strongly_connected(G):
        print('Initialized graph is not strongly connected')
        return

    def opposite_direction_exists(G, *edge_id):
        return G.has_edge(edge_id[1], edge_id[0])

    # Remove edges
    i = 0
    while True:
        i+=1
        print('Iteration ', i)
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

    return G


