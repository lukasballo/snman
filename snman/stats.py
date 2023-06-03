from .constants import *
from . import space_allocation, street_graph, graph_utils, lane_graph
from . import osmnx_customized as oxc
import copy
import networkx as nx
import pandas as pd


def street_sections(G, key_lanes_description=KEY_LANES_DESCRIPTION, weight='length'):

    G = copy.deepcopy(G)

    space_allocation.reorder_lanes(G, lanes_attribute=KEY_LANES_DESCRIPTION_AFTER)
    street_graph.organize_edge_directions(
        G,
        method='by_top_order_lanes',
        key_lanes_description=KEY_LANES_DESCRIPTION_AFTER
    )

    nodes, edges = oxc.graph_to_gdfs(G, nodes=True, edges=True)
    edges[key_lanes_description] = edges[key_lanes_description].apply(lambda x: '|'.join(x))
    result = edges[[key_lanes_description, weight]].groupby(key_lanes_description).sum(weight)
    result = result.sort_values(weight, ascending=False)

    return result


def street_sections_change(
        G,
        key_lanes_description=KEY_LANES_DESCRIPTION,
        key_lanes_description_after=KEY_LANES_DESCRIPTION_AFTER,
        weight='length'
):

    G = copy.deepcopy(G)


def network_metrics_for_all_measurement_regions(G, measurement_regions_gdf, plot_scc=False):

    metrics = measurement_regions_gdf.apply(
        lambda row: (
            row['name'],
            network_metrics(
                oxc.truncate.truncate_graph_polygon(G, row['geometry'], quadrat_width=100, retain_all=True),
                plot_scc=plot_scc
            )
        ),
        axis=1
    )

    return dict(list(metrics))


def network_metrics(G, plot_scc=False):

    if len(G.edges) == 0:
        return pd.DataFrame()

    configs = [
        {
            'label': 'CARS BEFORE',
            'lanes': KEY_LANES_DESCRIPTION,
            'mode': MODE_PRIVATE_CARS
        },
        {
            'label': 'CARS AFTER',
            'lanes': KEY_LANES_DESCRIPTION_AFTER,
            'mode': MODE_PRIVATE_CARS
        },
        {
            'label': 'CYCLING BEFORE',
            'lanes': KEY_LANES_DESCRIPTION,
            'mode': MODE_CYCLING
        },
        {
            'label': 'CYCLING AFTER',
            'lanes': KEY_LANES_DESCRIPTION_AFTER,
            'mode': MODE_CYCLING
        },
        {
            'label': 'TRANSIT BEFORE',
            'lanes': KEY_LANES_DESCRIPTION,
            'mode': MODE_TRANSIT
        },
        {
            'label': 'TRANSIT AFTER',
            'lanes': KEY_LANES_DESCRIPTION_AFTER,
            'mode': MODE_TRANSIT
        },
    ]

    df = pd.DataFrame(
        map(
            lambda config: _generate_metrics_for_one_config(
                G, config['label'], config['lanes'], config['mode'], plot_scc
            ),
            configs
        )
    )

    df = df.set_index(['label']).transpose()
    return df


def _generate_metrics_for_one_config(G, label, lanes_key, mode, plot_scc):
    L = street_graph.to_lane_graph(
        street_graph.filter_lanes_by_modes(G.copy(), {mode}, lane_description_key=lanes_key),
        lanes_attribute=lanes_key
    )
    L_lcc = graph_utils.keep_only_the_largest_connected_component(L)
    stats = lane_graph.calculate_stats(L_lcc, mode)
    stats['N_nodes_full'] = len(L.nodes)
    stats['N_edges_full'] = len(L.edges)
    stats['N_nodes_deleted'] = len(L.nodes) - len(L_lcc.nodes)
    stats['N_edges_deleted'] = len(L.edges) - len(L_lcc.edges)
    stats['label'] = label
    if plot_scc:
        print(label, 'N deleted =', stats['N_nodes_deleted'], 'E deleted =', stats['N_edges_deleted'])
        graph_utils.plot_scc(L)
    return stats
