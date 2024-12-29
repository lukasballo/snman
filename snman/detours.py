import networkx as nx
import pandas as pd
import geopandas as gpd
import shapely as shp
import itertools as it
from . import lane_graph, graph

def get_cheapest_edge_between_nodes(L, u, v, weight='cost_private_cars'):
    """
    Calculates the edge with the lowest cost between u and v.
    The is needed for assigning the correct lane to a path

    Parameters
    ----------
    L: lane_graph.LaneGraph
    u: int
    v: int
    weight: str

    Returns
    -------
    tuple
    """

    edges = L.get_edge_data(u, v)
    key = min(edges, key=lambda edge: edges[edge][weight])
    return (u, v, key), edges[key]


def shortest_path(L, u, v, weight='cost_private_cars', include_geometry=True, crs=2056):
    """
    Calculates the shortest path between u and v.

    Parameters
    ----------
    L: lane_graph.LaneGraph
    u: int
    v: int
    weight: str
    include_geometry: bool
        if True, shortest path will include geometry, otherwise only the value will be calculated
        (which is faster)
    crs: int

    Returns
    -------
    dict
    """

    if include_geometry:
        nodes = nx.dijkstra_path(L, u, v, weight)
        node_pairs = list(zip(nodes, nodes[1:]))

        # build path from node pairs
        path = map(lambda node_pair: get_cheapest_edge_between_nodes(L, *node_pair, weight=weight)[1], node_pairs)
        path = gpd.GeoDataFrame(path, geometry='geometry', crs=crs)

        return {
            'u': u,
            'v': v,
            'weight': weight,
            'total_cost': sum(path[weight]),
            'geometry': shp.ops.linemerge(path.unary_union)
        }

    else:
        total_cost = nx.dijkstra_path_length(L, u, v, weight)
        return {
            'u': u,
            'v': v,
            'weight': weight,
            'total_cost': total_cost
        }


def shortest_paths(L, uv_pairs, weight='cost_private_cars', include_geometry=True, crs=2056):
    """
    Calculates the shortest paths for pairs of u and v.

    Parameters
    ----------
    L: lane_graph.LaneGraph
    uv_pairs: list
        a list of (u,v) tuples
    weight: str
    crs: int

    Returns
    -------
    gpd.GeoDataFrame
    """

    paths = map(
        lambda uv_pair: shortest_path(L, *uv_pair, weight=weight, include_geometry=include_geometry),
        uv_pairs
    )
    if include_geometry:
        return gpd.GeoDataFrame(paths, geometry='geometry', crs=crs)
    else:
        return pd.DataFrame(paths)


def detour_metrics_for_origin(L, origin, destinations, weight='cost_private_cars', include_geometry=False, crs=2056):
    """
    Calculates detour metrics for a given origin.

    Parameters
    ----------
    L: lane_graph.LaneGraph
    origin: int
    destinations: list
    weight: str
    crs: int

    Returns
    -------

    """
    combinations = list(it.product(
        [origin],
        destinations
    ))
    #print(combinations)
    paths = shortest_paths(L, combinations, weight=weight, include_geometry=include_geometry, crs=crs)

    agg = {
        'mean_cost': ('total_cost', 'mean'),
        'median_cost': ('total_cost', 'median')
    }

    if include_geometry:
        agg['geometry'] = ('geometry', list)

    origin_metrics = paths.groupby('u').agg(**agg)

    if include_geometry:
        origin_metrics['geometry'] = origin_metrics['geometry'].apply(lambda geom: shp.MultiLineString(geom))

    return origin_metrics.iloc[0]
