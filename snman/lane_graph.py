import copy

from . import osmnx_customized as oxc
from . import space_allocation
from .constants import *
import networkx as nx


def create_lane_graph(G, lanes_attribute=KEY_LANES_DESCRIPTION):
    """
    Creates a new lane graph, derived from the street graph
    Parameters
    ----------
    G : nx.MultiDiGraph
        street graph
    lanes_attribute : str
        which attribute should be used for the lane description

    Returns
    -------

    """

    # initialize and copy graph attributes
    L = nx.MultiDiGraph()
    L.graph = G.graph

    for uvk, data in G.edges.items():
        u, v, k = uvk

        lanes_list = data.get(lanes_attribute)

        for i, lane in enumerate(lanes_list):

            lp = space_allocation._lane_properties(lane)
            reverse = lp.direction in {DIRECTION_BACKWARD, DIRECTION_BACKWARD_OPTIONAL}

            if reverse:
                lane = space_allocation.reverse_lane(lane)
                lp = space_allocation._lane_properties(lane)

            attributes = {}
            length = data['length']
            lane_id = '-'.join([str(u), str(v), str(k), str(i)])
            attributes['length'] = length

            for mode in MODES:
                cost = space_allocation._calculate_lane_cost(lane, length, mode)
                attributes['cost_' + mode] = cost

            attributes['primary_mode'] = lp.primary_mode
            attributes['lane_id'] = lane_id
            attributes['lanetype'] = lp.lanetype
            attributes['direction'] = lp.direction
            attributes['width'] = lp.width
            attributes['osm_highway'] = data.get('highway')
            attributes['maxspeed'] = data.get('maxspeed')

            attributes['fixed'] = lp.direction not in [
                DIRECTION_BACKWARD_OPTIONAL, DIRECTION_FORWARD_OPTIONAL, DIRECTION_TBD_OPTIONAL, DIRECTION_TBD
            ]
            attributes['mandatory_lane'] = lp.direction not in [
                DIRECTION_BACKWARD_OPTIONAL, DIRECTION_FORWARD_OPTIONAL, DIRECTION_TBD_OPTIONAL
            ]

            if not reverse:
                L.add_edge(u, v, lane_id, **attributes, lane=lane, twin_factor=1)
            if reverse:
                L.add_edge(v, u, lane_id, **attributes, lane=lane, twin_factor=1)
            if lp.direction in [
                DIRECTION_BOTH, DIRECTION_TBD, DIRECTION_TBD_OPTIONAL
            ]:
                L.add_edge(u, v, lane_id, **attributes, lane=lane, twin_factor=0.5)
                L.add_edge(v, u, lane_id, **attributes, lane=space_allocation.reverse_lane(lane), twin_factor=0.5)

    # take over the node attributes from the street graph
    nx.set_node_attributes(L, dict(G.nodes))

    return L



def calculate_stats(L, mode):
    """
    Calculates a set of standardized measures for a lane graph:

    * **usable_N_nodes**: number of nodes that are accessible for the given mode
    * **usable_N_edges**: number of edges that are accessible for the given mode
    * **convex_hull_km2**: an approximate area of the network using a convex hull of all accessible nodes
    * **usable_lane_km**: total length of lanes that are accessible for the given mode,
        bidirectional lanes count only once,
        pseudo-lanes with width=0 don't count
    * **usable_surface_km2**: total surface of lanes that are accessible for the given mode
    * **as_primary_lane_km**: total length of lanes where the given mode is primary
        (see constants.LANE_TYPES, the primary mode is the first one in the mode list of each lane type),
        bidirectional lanes count only once,
        pseudo-lanes with width=0 don't count
    * **as_primary_mode_lane_surface_km2**: total surface of lanes where the given mode is primary
    * **avg_betweenness_centrality_norm**: average normalized betweenness centrality of edges,
        the edge weights consider the scaling of cycling distanced by comfort
    * **avg_shortest_path_vod_km**: average length of the shortest path,
        for cycling, the length is scaled using the comfort factors (cycling infrastructure will lead to lower values),
        this means that you cannot compare the absolute average shortest path across modes, only the change before/after
    * **avg_shortest_path_km**: average length of the shortest path, without scaling for cycling

    Parameters
    ----------
    L : nx.MultiDiGraph
        lane graph
    mode : str
        for which mode should the stats be calculated, see constants.MODES

    Returns
    -------
    dict

    """

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
                sum(
                    [
                        1 * e['twin_factor']
                        for uvk, e
                        in L.edges.items()
                    ]
                ),
                3
            ),
        'convex_hull_km2': area_km2,
        'usable_lane_km':
            round(
                sum(
                    [
                        e['length'] * e['twin_factor']
                        for uvk, e
                        in L.edges.items()
                        if e['width'] > 0
                    ]
                ) / 1000,
                3
            ),
        'usable_lane_surface_km2':
            round(
                sum(
                    [
                        e['length'] * e['width'] * e['twin_factor']
                        for uvk, e
                        in L.edges.items()
                    ]
                ) / pow(1000, 2),
                3
            ),
        'as_primary_mode_lane_km':
            round(
                sum(
                    [
                        e['length'] * e['twin_factor'] * (e['primary_mode'] == mode)
                        for uvk, e
                        in L.edges.items()
                        if e['width'] > 0
                    ]
                ) / 1000,
                3
            ),
        'as_primary_mode_lane_surface_km2':
            round(
                sum(
                    [
                        e['length'] * e['width'] * e['twin_factor'] * (e['primary_mode'] == mode)
                        for uvk, e
                        in L.edges.items()
                    ]
                ) / pow(1000, 2),
                3
            ),
        'as_primary_mode_lane_surface_km2':
            round(
                sum(
                    [
                        e['length'] * e['width'] * e['twin_factor'] * (e['primary_mode'] == mode)
                        for uvk, e
                        in L.edges.items()
                    ]
                ) / pow(1000, 2),
                3
            ),
        'avg_betweenness_centrality_norm':
            round(
                sum(
                    [
                        e['bc'] * e['twin_factor']
                        for uvk, e
                        in L.edges.items()
                    ]
                ) /
                sum(
                    [
                        e['twin_factor']
                        for uvk, e
                        in L.edges.items()
                    ]
                ),
                5
            ),
        'avg_shortest_path_vod_km':
            round(
                nx.average_shortest_path_length(L, 'cost_' + mode) / 1000,
                3
            ),
        'avg_shortest_path_km':
            round(
                nx.average_shortest_path_length(L, 'length') / 1000,
                3
            ),
    }