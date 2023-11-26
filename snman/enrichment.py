import copy

import leuvenmapmatching.matcher.base
from leuvenmapmatching.matcher.distance import DistanceMatcher
from leuvenmapmatching.map.inmem import InMemMap
from leuvenmapmatching import visualization as mmviz
import shapely as shp
from statistics import mean
import pandas as pd
import geopandas as gpd
import networkx as nx
import warnings
from . import utils, space_allocation, street_graph, io
from .constants import *


def match_linestrings(
        G, source, column_configs,
        remove_short_overlaps=True, max_dist2=None,
        remove_sidetrips=True,
        lanes_key = KEY_LANES_DESCRIPTION,
        modes = None,
        _save_map=None,
        **distance_matcher_args
):
    """
    Make a spatial join of the graph edges and polylines in a GeoDataFrame.
    Copy selected attributes from the GeoFDataFrame to the graph edges.

    The spatial join is made using LeuvenMapMatching:
    https://github.com/wannesm/LeuvenMapMatching

    Parameters
    ----------
    G : nx.MultiGraph
        street graph, target
    source : gpd.GeoDataFrame
        data source
    column_configs : list
        a list of dictionaries, see example

    Examples
    --------

    >>> column_configs = [
    >>>     {'source_column': 'DTV_ALLE',   'target_column': 'adt_avg',         'agg': 'avg' },
    >>>     {'source_column': 'DTV_ALLE',   'target_column': 'adt_max',         'agg': 'max' },
    >>>     {'source_column': 'FROMNODENO', 'target_column': 'npvm_fromnodeno', 'agg': 'list'},
    >>>     {'source_column': 'TONODENO',   'target_column': 'npvm_tonodeno',   'agg': 'list'}
    >>> ]

    """

    # stop here if the source geodataframe is empty
    if len(source) == 0:
        return

    # create empty in-memory map for the mapmatching process
    map_con = InMemMap("source", use_latlon=False, use_rtree=True, index_edges=True, crs_xy=2056)

    # add nodes to the in-memory map
    # please note that lv works with lat, lon (reverse order)
    for id, data in G.nodes.items():
        map_con.add_node(id, (data['y'], data['x']))

    # keep only edges accessible to at least one of the provided modes
    if modes:
        H = street_graph.filter_lanes_by_modes(G, modes, lane_description_key=lanes_key)
    else:
        H = G

    if _save_map:
        I = copy.deepcopy(G)
        I.remove_edges_from(G.edges())

    # add edges to the in-memory map
    for uvk, data in H.edges.items():
        u = int(uvk[0])
        v = int(uvk[1])
        # add edges to the map
        map_con.add_edge(u, v)
        map_con.add_edge(v, u)

        if _save_map:
            I.add_edge(u, v)
            I.add_edge(v, u)

    if _save_map:
        io.export_street_graph(I, *_save_map)

    # create a matcher
    matcher = DistanceMatcher(map_con, **distance_matcher_args)

    # submit all source linestrings to the matcher
    source['node_pairs'] = source.apply(
        lambda x: _submit_linestring_to_matcher(
            matcher,
            x['geometry'],
            remove_short_overlaps,
            remove_sidetrips,
            max_dist2
        ), axis=1)

    # transfer the attributes as specified in "column_configs"
    # note that there may be multiple values that will be transferred to a single target edge,
    # so they need to be aggregated (see next step)
    for config in column_configs:

        edge_values = {}
        for index, edge in source.iterrows():
            value = edge[config['source_column']]
            node_pairs = edge['node_pairs']
            if len(node_pairs) < 1:
                continue
            for node_pair in node_pairs:
                u = node_pair[0]
                v = node_pair[1]
                if edge_values.get((u, v)) is None:
                    edge_values[(u, v)] = []
                edge_values[(u, v)].append(value)

        # aggregate the target data using one of the following functions
        for uvk, data in G.edges.items():
            if config['agg'] == 'avg':
                data[config['target_column'] + '_forward']  = mean(edge_values.get(uvk[:2], [0]))
                data[config['target_column'] + '_backward'] = mean(edge_values.get(uvk[:2][::-1], [0]))
            if config['agg'] == 'max':
                data[config['target_column'] + '_forward']  = max(edge_values.get(uvk[:2], [0]))
                data[config['target_column'] + '_backward'] = max(edge_values.get(uvk[:2][::-1], [0]))
            if config['agg'] == 'count':
                data[config['target_column'] + '_forward']  = len(edge_values.get(uvk[:2], []))
                data[config['target_column'] + '_backward'] = len(edge_values.get(uvk[:2][::-1], []))
            if config['agg'] == 'list':
                data[config['target_column'] + '_forward']  = str(edge_values.get(uvk[:2], []))
                data[config['target_column'] + '_backward'] = str(edge_values.get(uvk[:2][::-1], []))
            if config['agg'] == 'add_lanes':
                forward = utils.flatten_list(edge_values.get(uvk[:2], []))
                data[config['target_column']].extend(forward)
                backward = utils.flatten_list(edge_values.get(uvk[:2][::-1], []))
                backward = space_allocation.reverse_lanes(backward)
                data[config['target_column']] = backward + data[config['target_column']]
            if config['agg'] == 'replace_lanes':
                forward = list(utils.flatten_list(edge_values.get(uvk[:2], [])))
                backward = list(utils.flatten_list(edge_values.get(uvk[:2][::-1], [])))
                backward = space_allocation.reverse_lanes(backward)
                if len(forward+backward) > 0:
                    data[config['target_column']] = backward + forward
                    #print(data[config['target_column']])


def _submit_linestring_to_matcher(matcher, geom, remove_short_overlaps, remove_sidetrips, max_distance):
    """
    Matches one linestring, respecting the specialties of lvmapmatching,
    e.g., reversed coordinates

    Parameters
    ----------
    matcher : DistanceMatcher
    geom : shp.geometry.LineString

    Returns
    -------
    list
        matched node IDs

    """
    # convert any multilinestring to a linestring
    geom = utils.multilinestring_to_linestring(geom)
    path = geom.coords
    # coordinates must be reversed for lvmapmatching
    path = [coords[::-1] for coords in path]
    matcher.match(path)

    distances = matcher.path_all_distances()
    if len(distances) < 2:
        return []

    matched_nodes = matcher.path_pred_onlynodes

    # remove first/last node from target path if the overlap is short
    if remove_short_overlaps:
        for ij in [(0, 1), (-1, -2)]:

            if len(matched_nodes) < 2:
                return []
            # i: first/last node
            # j: second/second-last node
            i, j = ij
            # first/last point of the source (to be matched)
            source_i = shp.Point(geom.coords[i])
            # first/last matched node on the target network
            target_i = shp.Point(
                matcher.map.node_coordinates(
                    matched_nodes[i]
                )[::-1]
            )
            target_j = shp.Point(
                matcher.map.node_coordinates(
                    matched_nodes[j]
                )[::-1]
            )
            dist_source_i_target_i = source_i.distance(target_i)
            dist_source_i_target_j = source_i.distance(target_j)

            if dist_source_i_target_i > dist_source_i_target_j:
                # delete the point
                del matched_nodes[i]

    if max_distance:
        if max(distances) > max_distance:
            return []

    matched_node_pairs = [(matched_nodes[i], matched_nodes[i + 1]) for i in range(len(matched_nodes) - 1)]

    # see the german term 'stichfahrt' for sidetrips
    if remove_sidetrips and len(matched_node_pairs) >= 2:
        matched_node_pairs_filtered = []
        # iterate over all pairs except the last one
        for i, node_pair in enumerate(matched_node_pairs[:-2]):
            next_node_pair = matched_node_pairs[i+1]
            # check if the next pair is a reverse of this one which would make it a sidetrip
            if node_pair[::-1] != next_node_pair:
                matched_node_pairs_filtered.append(node_pair)
        # write the result back into the original list
        matched_node_pairs = matched_node_pairs_filtered

    return matched_node_pairs


def match_lane_edits(G, lane_edits, lanes_key=KEY_LANES_DESCRIPTION, **mapmatching_arguments):
    column_configs = [
        {'source_column': 'add_lanes', 'target_column': lanes_key, 'agg': 'add_lanes'}
    ]

    return match_linestrings(G, lane_edits, column_configs, lanes_key=lanes_key, **mapmatching_arguments)


def match_parking_spots(
        G, parking_spots,
        parking_space_length=7, parking_space_length_for_one_lane=24,
        remove_previous_parking=True
):

    # copy the parking spaces dataset and convert the points into zero-length linestrings
    # that can be matched onto the edges
    parking_spots = copy.deepcopy(parking_spots)
    parking_spots.geometry = parking_spots.geometry.apply(lambda x: shp.LineString([x, x]))

    parking_count_key = '_n_parking_spots'
    column_configs = [
        {'source_column': 'id1', 'target_column': parking_count_key, 'agg': 'count'}
    ]

    # create a copy of the street graph that only contains ground-level edges
    H = copy.deepcopy(G)
    for uvk, data in G.edges.items():
        if data.get('layer') != 0:
            H.remove_edge(*uvk)

    # match on the ground-level street graph
    match_linestrings(
        H, parking_spots, column_configs, remove_short_overlaps=False,
        modes=[MODE_PRIVATE_CARS, MODE_TRANSIT],
        max_dist=30, max_dist_init=30, max_lattice_width=5
    )

    # copy the matched attributes into the original street graph
    for _direction in ['_forward', '_backward']:
        nx.set_edge_attributes(
            G,
            nx.get_edge_attributes(H, parking_count_key + _direction),
            parking_count_key + _direction
        )

    for uvk, data in G.edges.items():
        # please note that the forward/backward distinction is meaningless and completely random in this case
        parking_count = data.get(parking_count_key + '_forward', 0) + data.get(parking_count_key + '_backward', 0)
        data['n_parking_spots'] = parking_count
        length = data['length']
        parking_space_density = utils.safe_division(parking_count, length)
        # estimating the number of parking lanes based on the number of parking spots that have been matched
        n_parking_lanes = round(parking_space_density / (1 / parking_space_length))
        if n_parking_lanes == 0 and parking_space_density > 1 / parking_space_length_for_one_lane:
            n_parking_lanes = 1

        if remove_previous_parking:
            data[KEY_LANES_DESCRIPTION] = space_allocation.filter_lanes_by_modes(
                data[KEY_LANES_DESCRIPTION],
                MODES.difference({MODE_CAR_PARKING})
            )

        data[KEY_LANES_DESCRIPTION].extend([LANETYPE_PARKING_PARALLEL + DIRECTION_BOTH] * n_parking_lanes)


def match_public_transit(G, pt_routes):
    """
    Match public transit routes onto the street graph using mapmatching.

    Parameters
    ----------
    G : nx.MultiDiGraph
        street graph
    routes : gpd.GeoDataFrame
        transit routes

    Returns
    -------

    """

    pt_routes = pt_routes.explode(ignore_index=True)

    column_configs = [
        {'source_column': 'LINIENNUMM', 'target_column': '_pt_routes', 'agg': 'list'},
        {'source_column': 'RICHTUNG', 'target_column': '_pt_directions', 'agg': 'list'},
    ]

    match_linestrings(
        G, pt_routes, column_configs, remove_short_overlaps=False,
        modes=(MODE_TRANSIT, MODE_PRIVATE_CARS),
        max_dist=200, max_dist_init=500, max_lattice_width=5
    )

    for uvk, data in G.edges.items():
        data['pt_forward'] = data['_pt_routes_forward'] != '[]'
        data['pt_backward'] = data['_pt_routes_backward'] != '[]'
        # for backward compatibility, to be removed later
        data['pt_bus'] = data['pt_forward'] or data['pt_backward']


def match_sensors(G, sensors_df):
    """
    Assign traffic sensors to edges in the street graph. Must be used with OSM graph.

    Parameters
    ----------
    G : nx.MultiDiGraph
        OSM graph
    sensors_df : pd.DataFrame
        sensors and their osm links

    Returns
    -------
    None
    """

    s = sensors_df.copy()
    s['id'] = s.index
    s = s.set_index(['u', 'v', 'osmid']).sort_index()

    for uvk, data in G.edges.items():
        u, v, k = uvk

        i_forward = (u, v, data['osmid'])
        if i_forward in s.index:
            data['sensors_forward'] = [s.loc[i_forward]['id']]
        else:
            data['sensors_forward'] = []

        i_backward = (v, u, data['osmid'])
        if i_backward in s.index:
            data['sensors_backward'] = [s.loc[i_backward]['id']]
        else:
            data['sensors_backward'] = []


def match_public_transit_by_buffers(G, pt_network):
    """
    Matches the public transport network onto the street network using simple buffers and spatial join.

    DEPRECATED: Use match_public_transit() instead.

    Parameters
    ----------
    G: nx.MultiGraph
        street graph
    pt_network: gpd.GeoDataFrame
        Linestrings of public transport routes, based on the open dataset of ZVV
        "Linien des öffentlichen Verkehrs (OGD) (kantonaler Datensatz)"
        https://data.stadt-zuerich.ch/dataset/ktzh_linien_des_oeffentlichen_verkehrs__ogd_,
        with following additional columns:
            * TYPE (Tram, Bus, Nightbus, Microbus)
            * ALIGNMENT (tunnel, <None>)
    """

    # Add buffers around the pt routes
    pt_network_buffers = pt_network.copy()
    pt_network_buffers.geometry = pt_network_buffers.geometry.buffer(15, resolution=16)

    for edge in G.edges(data=True, keys=True):

        edge_geometry = edge[3]['geometry']

        # Get all pt routes whose buffers intersect with this edge
        pt_routes = pt_network_buffers[
            pt_network_buffers.intersects(edge_geometry).tolist()
            and pt_network_buffers['ALIGNMENT'] != 'tunnel'
        ].copy()

        if len(pt_routes) == 0:
            continue

        # Calculate the overlapping length
        with warnings.catch_warnings():
            # we suppress warnings due to a known issue in the intersection function
            # see here: https://github.com/shapely/shapely/issues/1345
            warnings.simplefilter("ignore")
            pt_routes['intersection_length_prop'] = [
                edge_geometry.intersection(pt_route).length / edge_geometry.length if edge_geometry.length != 0 else 0
                for pt_route in pt_routes.geometry
            ]

        # Keep only those that overlap over a substantial part
        pt_routes = pt_routes[pt_routes['intersection_length_prop'] > 0.7]

        edge[3]['pt_routes'] = str(pt_routes.LINIENNUMM.tolist())
        # TODO: Distinguish directions of pt lines
        edge[3]['pt_tram'] = pt_routes['TYPE'].eq('Tram').any() * 1
        edge[3]['pt_bus'] = pt_routes['TYPE'].eq('Bus').any() * 1
        edge[3]['pt_mcrbus'] = pt_routes['TYPE'].eq('Microbus').any() * 1

