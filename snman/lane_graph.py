from . import osmnx_customized as oxc
import networkx as nx
import statistics as stats


def calculate_stats(L):

    points_gpd = oxc.graph_to_gdfs(L, edges=False)
    area_km2 = points_gpd.geometry.unary_union.convex_hull.area / pow(1000, 2)

    return {
        'N_nodes': len(L.nodes),
        'N_edges': len(L.edges),
        'N_strongly_connected_components': len(list(nx.strongly_connected_components(L))),
        'N_weakly_connected_components': len(list(nx.weakly_connected_components(L))),
        'lane_km':
            round(
                sum(nx.get_edge_attributes(L, 'length').values()) / 1000,
            3),
        'lane_km/km2':
            round(
                (sum(nx.get_edge_attributes(L, 'length').values()) / 1000) / area_km2,
            3),
        'active_modes_only_lane_km':
            round(
                sum(nx.get_edge_attributes(L, 'only_active_modes_length').values()) / 1000,
            3),
        'active_modes_only_lane_km/km2':
            round(
                (sum(nx.get_edge_attributes(L, 'only_active_modes_length').values()) / 1000) / area_km2,
            3),
        'avg_betweenness_centrality_norm':
            round(
                stats.mean(nx.betweenness_centrality(L, normalized=True).values()),
            5),
        'avg_shortest_path_length_km':
            round(
                nx.average_shortest_path_length(L, 'cost') / 1000,
            3),
    }