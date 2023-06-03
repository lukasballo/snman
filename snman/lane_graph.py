import copy

from . import osmnx_customized as oxc
import networkx as nx
import statistics as stats


def calculate_stats(L, mode):

    L = copy.deepcopy(L)

    # calculate the approximate area as a convex hull of all nodes
    points_gpd = oxc.graph_to_gdfs(L, edges=False)
    area_km2 = points_gpd.geometry.unary_union.convex_hull.area / pow(1000, 2)

    # set betweenness centrality
    nx.set_edge_attributes(L, nx.edge_betweenness_centrality(L, normalized=True, weight='cost_'+mode), 'bc')

    return {
        'usable_N_nodes': len(L.nodes),
        'usable_N_edges':
            round(
                sum([1 * e['twin_factor'] for uvk, e in L.edges.items()]),
                3
            ),
        'convex_hull_km2': area_km2,
        'usable_lane_km':
            round(
                sum([e['length'] * e['twin_factor'] for uvk, e in L.edges.items()]) / 1000,
                3
            ),
        'usable_lane_surface_km2':
            round(
                sum([e['length'] * e['width'] * e['twin_factor'] for uvk, e in L.edges.items()]) / pow(1000, 2),
                3
            ),
        'as_primary_mode_lane_km':
            round(
                sum([e['length'] * e['twin_factor'] * (e['primary_mode'] == mode) for uvk, e in L.edges.items()]) / 1000,
                3
            ),
        'as_primary_mode_lane_surface_km2':
            round(
                sum([e['length'] * e['width'] * e['twin_factor'] * (e['primary_mode'] == mode) for uvk, e in L.edges.items()]) / pow(1000, 2),
                3
            ),
        'avg_betweenness_centrality_norm':
            round(
                sum([e['bc'] * e['twin_factor'] for uvk, e in L.edges.items()]) / sum([e['twin_factor'] for uvk, e in L.edges.items()]),
                5
            ),
        'avg_shortest_path_length_km':
            round(
                nx.average_shortest_path_length(L, 'cost_' + mode) / 1000,
            3),
    }