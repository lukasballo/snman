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

from numpy.f2py.auxfuncs import throw_error

from . import utils, space_allocation, street_graph, io, geometry_tools, _errors
from . import osmnx_customized as oxc
from .constants import *


def aggregate_points_to_nearest_nodes(G, points, attribute, aggregation, max_distance=100, fillna_value=0):
    """
    Matches points to nearest nodes and adds an attribute with an aggregated value from the points.

    Parameters
    ----------
    G: street_graph.StreetGraph
    points: gpd.GeoDataFrame
    attribute: str
        which attribute from the points geodataframe should be aggregated into the nodes
    aggregation: str
        see pandas.core.groupby.generic.DataFrameGroupBy for a complete list of aggregation functions,
        e.g., 'sum'
    max_distance: int
    fillna_value: any
        a value that will be written if the aggregation does not return anything
    """

    nodes = oxc.graph_to_gdfs(G, edges=False)
    nodes['id'] = nodes.index

    points = gpd.sjoin_nearest(
        points,
        nodes, how='inner', max_distance=max_distance, distance_col='dist')
    points.rename(columns={'index_right': 'node_id'}, inplace=True)
    points = points.groupby('node_id').agg({attribute: aggregation})

    nodes = nodes.merge(points, how='left', left_on='id', right_on='node_id')
    nodes.set_index('id', inplace=True)
    nodes[attribute].fillna(fillna_value, inplace=True)

    a = dict(nodes[attribute])
    attribute_name = attribute + '_' + aggregation
    nx.set_node_attributes(G, fillna_value, name=attribute_name)
    nx.set_node_attributes(G, a, name=attribute_name)


def match_linestrings(
        G, source, column_configs,
        remove_short_overlaps=True, max_dist2=None,
        remove_sidetrips=True,
        lanes_key = KEY_LANES_DESCRIPTION,
        modes = None,
        verbose = False,
        _save_map = None,
        max_dist=200,
        max_dist_init=400,
        max_lattice_width=20,
        **distance_matcher_args
):
    """
    Make a spatial join of the graph edges and polylines in a GeoDataFrame.
    Copy selected attributes from the GeoFDataFrame to the graph edges.

    The spatial join is made using LeuvenMapMatching:
    https://github.com/wannesm/LeuvenMapMatching

    # TODO: Clean up the aggregation functions

    Parameters
    ----------
    G : nx.MultiGraph
        street graph, target
    source : gpd.GeoDataFrame
        data source
    column_configs : list
        a list of dictionaries, following aggregation functions are supported:
        avg, max, list, add_lanes, replace_lanes

    Examples
    --------

    >>> column_configs = [
    >>>     {'source_column': 'DTV_ALLE',   'target_column': 'adt_avg',         'agg': 'avg' },
    >>>     {'source_column': 'DTV_ALLE',   'target_column': 'adt_max',         'agg': 'max' },
    >>>     {'source_column': 'FROMNODENO', 'target_column': 'npvm_fromnodeno', 'agg': 'list'},
    >>>     {'source_column': 'lanes',      'target_column': 'ln_desc',         'agg': 'replace_lanes'}
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
    H = copy.deepcopy(G)
    if modes:
        street_graph.filter_lanes_by_modes(H, modes, lane_description_key=lanes_key)

    if _save_map:
        I = copy.deepcopy(H)
        I.remove_edges_from(H.edges())

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
    matcher = DistanceMatcher(
        map_con,
        max_dist=max_dist,
        max_dist_init=max_dist_init,
        max_lattice_width=max_lattice_width,
        **distance_matcher_args
    )

    def submit_source(row):

        if verbose:
            print(f'submitting linestring {row.name}')

        return _submit_linestring_to_matcher(
            matcher,
            row['geometry'],
            remove_short_overlaps,
            remove_sidetrips,
            max_dist2
        )

    # submit all source linestrings to the matcher
    source['node_pairs'] = source.apply(submit_source, axis=1)

    #print(source['node_pairs'])

    # transfer the attributes as specified in "column_configs"
    # note that there may be multiple values that will be transferred to a single target edge,
    # so they need to be aggregated (see next step)
    for config in column_configs:

        edge_values = {}
        for index, edge in source.iterrows():
            # value from the source edge
            value = edge[config['source_column']]
            node_pairs = edge['node_pairs']
            if len(node_pairs) < 1:
                continue
            for node_pair in node_pairs:
                # target edge
                u = node_pair[0]
                v = node_pair[1]
                if edge_values.get((u, v)) is None:
                    edge_values[(u, v)] = []
                edge_values[(u, v)].append(
                    copy.deepcopy(value)
                )

        # aggregate the target data using one of the following functions
        for uvk, data in G.edges.items():

            if config['agg'] == 'avg':
                data[config['target_column'] + '_forward']  = mean(edge_values.get(uvk[:2], [0]))
                data[config['target_column'] + '_backward'] = mean(edge_values.get(uvk[:2][::-1], [0]))

            elif config['agg'] == 'max':
                data[config['target_column'] + '_forward']  = max(edge_values.get(uvk[:2], [0]))
                data[config['target_column'] + '_backward'] = max(edge_values.get(uvk[:2][::-1], [0]))

            elif config['agg'] == 'sum':
                data[config['target_column'] + '_forward'] = sum(edge_values.get(uvk[:2], [0]))
                data[config['target_column'] + '_backward'] = sum(edge_values.get(uvk[:2][::-1], [0]))

            elif config['agg'] == 'count':
                data[config['target_column'] + '_forward']  = len(edge_values.get(uvk[:2], []))
                data[config['target_column'] + '_backward'] = len(edge_values.get(uvk[:2][::-1], []))

            elif config['agg'] == 'list':
                data[config['target_column'] + '_forward']  = str(edge_values.get(uvk[:2], []))
                data[config['target_column'] + '_backward'] = str(edge_values.get(uvk[:2][::-1], []))

            elif config['agg'] == 'max_no_direction':
                forward  = max(edge_values.get(uvk[:2], [0]))
                backward = max(edge_values.get(uvk[:2][::-1], [0]))
                data[config['target_column']] = max([str(forward), str(backward)])

            elif config['agg'] == 'list_no_direction':
                forward  = edge_values.get(uvk[:2], [])
                backward = edge_values.get(uvk[:2][::-1], [])
                data[config['target_column']] = str(forward + backward)

            elif config['agg'] == 'has_match_no_direction':
                forward  = edge_values.get(uvk[:2], False)
                backward = edge_values.get(uvk[:2][::-1], False)
                data[config['target_column']] = not (forward is False and backward is False)

            elif config['agg'] == 'add_lanes':
                forward = space_allocation.SpaceAllocation(utils.flatten_list(edge_values.get(uvk[:2], [])))
                data[config['target_column']].extend(forward)
                backward = space_allocation.SpaceAllocation(utils.flatten_list(edge_values.get(uvk[:2][::-1], [])))
                backward.reverse_allocation()
                data[config['target_column']] = backward + data[config['target_column']]

            elif config['agg'] == 'replace_lanes':
                forward = space_allocation.SpaceAllocation(utils.flatten_list(edge_values.get(uvk[:2], [])))
                backward = space_allocation.SpaceAllocation(utils.flatten_list(edge_values.get(uvk[:2][::-1], [])))
                backward.reverse_allocation()
                if len(forward+backward) > 0:
                    data[config['target_column']] = space_allocation.SpaceAllocation(backward + forward)

            # delete the values that were used in this step. this avoids that they will be matched onto multiple edges
            if edge_values.get(uvk[:2], None) is not None:
                del edge_values[uvk[:2]]
            if edge_values.get(uvk[:2][::-1], None) is not None:
                del edge_values[uvk[:2][::-1]]


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

    #print('*')

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
    #print(f"matched node pairs: {matched_node_pairs}")

    # see the german term 'stichfahrt' for sidetrips
    if remove_sidetrips and len(matched_node_pairs) >= 2:
        matched_node_pairs_filtered = []
        # iterate over all pairs except the last one
        last_pair_was_removed = False
        for i, node_pair in enumerate(matched_node_pairs):
            next_node_pair = matched_node_pairs[i+1] if i+2 < len(matched_nodes) else None
            if last_pair_was_removed == True:
                last_pair_was_removed = False
                # (don't add the node pair)
            else:
                # check if the next pair is a reverse of this one which would make it a sidetrip
                if node_pair[::-1] == next_node_pair:
                    last_pair_was_removed=True
                    # (don't add the node pair)
                else:
                    matched_node_pairs_filtered.append(node_pair)
        # write the result back into the original list
        matched_node_pairs = matched_node_pairs_filtered

    #print(f"after sidetrips removal: {matched_node_pairs}")
    return matched_node_pairs


def match_lane_edits(
        G, manual_lanes,
        lanes_key=KEY_LANES_DESCRIPTION,
        source_column='lanes',
        **mapmatching_arguments
):
    """
    Replaces lanes on edges of the street graph based on a geodataframe of manually drawn lines with the new lanes.
    Each line in the geodataframe will be mapmathed onto the street graph automatically.

    Parameters
    ----------
    G: street_graph.StreetGraph
    manual_lanes: gpd.GeoDataFrame
        approximate lines with the necessary attribute containing the target lanes
    lanes_key: str
        which lanes key should be replaced
    mapmatching_arguments: **kwargs
        see leuven mapmatching
    """

    column_configs = [{'source_column': source_column, 'target_column': lanes_key, 'agg': 'replace_lanes'}]
    match_linestrings(
        G, manual_lanes, column_configs,
        **mapmatching_arguments
    )


def match_hierarchy_edits(
        G, manual_hierarchy,
        source_column='hierarchy',
        **mapmatching_arguments
):
    """
    Replaces hierarchies of the street graph based on a geodataframe of manually drawn lines with the new lanes.
    Each line in the geodataframe will be mapmathed onto the street graph automatically.

    Parameters
    ----------
    G: street_graph.StreetGraph
    manual_hierarchy: gpd.GeoDataFrame
        approximate lines with the necessary attribute containing the target hierarchy
    mapmatching_arguments: **kwargs
        see leuven mapmatching
    """

    # store the old hierarchy
    for uvk, data in G.edges.items():
        data['_hierarchy_old'] = data['hierarchy']

    # for all unmatched edges, the matcher will overwrite the hierarchy values with '0'
    column_configs = [{'source_column': source_column, 'target_column': 'hierarchy', 'agg': 'max_no_direction'}]
    match_linestrings(
        G, manual_hierarchy, column_configs,
        **mapmatching_arguments
    )

    # overwrite the '0' values with old hierarchy values
    for uvk, data in G.edges.items():
        if data['hierarchy'] == '0':
            data['hierarchy'] = data['_hierarchy_old']


def match_width(
        G, manual_hierarchy,
        source_column='width',
        **mapmatching_arguments
):
    """
    Replaces widths of the street graph based on a geodataframe of manually drawn lines with the new lanes.
    Each line in the geodataframe will be mapmathed onto the street graph automatically.

    Parameters
    ----------
    G: street_graph.StreetGraph
    manual_hierarchy: gpd.GeoDataFrame
        approximate lines with the necessary attribute containing the target width
    mapmatching_arguments: **kwargs
        see leuven mapmatching
    """

    # store the old hierarchy
    for uvk, data in G.edges.items():
        data['_width_old'] = data.get('width')

    # for all unmatched edges, the matcher will overwrite the hierarchy values with '0'
    column_configs = [{'source_column': source_column, 'target_column': 'width', 'agg': 'max_no_direction'}]
    match_linestrings(
        G, manual_hierarchy, column_configs,
        **mapmatching_arguments
    )

    # overwrite the '0' values with old width values
    for uvk, data in G.edges.items():
        if data['width'] == '0':
            data['width'] = data['_width_old']
        data['width'] = utils.safe_float(data['width'])
        data['_width_old'] = utils.safe_float(data['_width_old'])


def _match_parking_generic(
        G, parking_spots,
        source_column,
        agg,
        parking_space_length=7, parking_space_length_for_one_lane=24,
        remove_previous_parking=True,
        parking_type='car'
):

    parking_count_key = '_n_parking_spots'

    column_configs = [
        {'source_column': source_column, 'target_column': parking_count_key, 'agg': agg}
    ]

    # create a copy of the street graph that only contains ground-level edges
    H = copy.deepcopy(G)
    for uvk, data in G.edges.items():
        if data.get('layer') != 0:
            H.remove_edge(*uvk)

    if parking_type == 'car':
        modes = [MODE_PRIVATE_CARS]
    elif parking_type == 'bicycle':
        modes = [MODE_CYCLING]

    # match on the ground-level street graph
    match_linestrings(
        H, parking_spots, column_configs, remove_short_overlaps=False,
        modes=modes,
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
            if parking_type == 'car':
                remove_modes = {MODE_CAR_PARKING}
            elif parking_type == 'bicycle':
                remove_modes = {MODE_BICYCLE_PARKING}
            data[KEY_LANES_DESCRIPTION] = space_allocation.filter_lanes_by_modes(
                data[KEY_LANES_DESCRIPTION],
                MODES.difference(remove_modes)
            )

        for i in range(n_parking_lanes):
            if parking_type == 'car':
                lanetype = LANETYPE_PARKING_PARALLEL
            elif parking_type == 'bicycle':
                lanetype = LANETYPE_BICYCLE_PARKING
            data[KEY_LANES_DESCRIPTION].append(space_allocation.Lane(lanetype, DIRECTION_FORWARD))


def match_parking_spots_lines(
        G, parking_spots,
        capacity_column='capacity',
        **kwargs
):
    """
    Not reliable, due to double counting of parking spots
    """

    _match_parking_generic(
        G, parking_spots,
        source_column=capacity_column,
        agg='sum',
        **kwargs
    )


def match_parking_spots_polygons(
        G, parking_spots,
        capacity_column='capacity',
        **kwargs
):
    """
    Not reliable, due to double counting of parking spots
    """

    parking_spots = copy.deepcopy(parking_spots)
    parking_spots.geometry = parking_spots.geometry.apply(
        lambda poly: geometry_tools.get_polygon_axis(poly)
    )

    _match_parking_generic(
        G, parking_spots,
        source_column=capacity_column,
        agg='sum',
        **kwargs
    )


def match_parking_spots_polygons_points_sampling(
        G, parking_spots,
        capacity_column='capacity',
        **kwargs
):

    parking_spots = copy.deepcopy(parking_spots)
    parking_spots.geometry = parking_spots.apply(
        lambda row: geometry_tools.random_points_in_polygon(row['geometry'], row[capacity_column]),
        axis=1
    )

    parking_spots = parking_spots.explode().reset_index(drop=True)
    parking_spots['id'] = parking_spots.index

    match_parking_spots(
        G, parking_spots,
        **kwargs
    )


def match_parking_spots(
        G, parking_spots,
        id_column='id',
        **kwargs
):
    """
    Match on-street parking spots provided as points onto the street graph as added parking lanes.
    The resulting density of parking spots will be automatically converted to a number of parking lanes.
    No parking lanes will be added if the density os very low.

    Parameters
    ----------
    G: street_graph.StreetGraph
    parking_spots: gpd.GeoDataFrame
        a points dataset with parking spots
    parking_space_length: float
        typical parking space length
    parking_space_length_for_one_lane: float
        typical parking space length for a single parking lane
    remove_previous_parking: bool
        if true, all previously existing parking lanes will be overwritten
    """

    # copy the parking spaces dataset and convert the points into zero-length linestrings
    # that can be matched onto the edges
    parking_spots = copy.deepcopy(parking_spots.reset_index())
    parking_spots['id'] = parking_spots.index
    parking_spots.geometry = parking_spots.geometry.apply(lambda x: shp.LineString([x, x]))

    _match_parking_generic(
        G, parking_spots,
        source_column=id_column,
        agg='count',
        **kwargs
    )

def match_public_transit_zvv(
        G, pt_routes,
        route_number_column='LINIENNUMM',
        direction_column='RICHTUNG',
        **matcher_kwargs
):
    """
    Match public transit routes onto the street graph using mapmatching.
    The geodata must follow the same structure as the ZVV dataset of Zurich:
    https://www.stadt-zuerich.ch/geodaten/download/108
     * Every route direction is represented by a separate feature with proper directions
     * The features can be linestrings or polylinestrings. In case of the latter they are automatically dissolved into
       single linestrings

    Parameters
    ----------
    G : nx.MultiDiGraph
        street graph
    pt_routes : gpd.GeoDataFrame
        transit routes
    max_dist : int
        see settings of the leuven mapmatcher
    max_lattice_width : int
        see settings of the leuven mapmatcher
    route_number_column : str
        which column in the routes dataset holds the route number
    direction_column : str
        which column in the routes dataset holds the route direction id
    """

    # convert multilinestrings into linestrings
    pt_routes = pt_routes.explode(ignore_index=True)

    column_configs = [
        {'source_column': route_number_column, 'target_column': '_pt_routes', 'agg': 'list'},
        {'source_column': direction_column, 'target_column': '_pt_directions', 'agg': 'list'},
    ]

    match_linestrings(
        G, pt_routes, column_configs, remove_short_overlaps=True,
        modes=(MODE_TRANSIT, MODE_PRIVATE_CARS),
        **matcher_kwargs
    )

    for uvk, data in G.edges.items():
        data['pt_forward'] = data['_pt_routes_forward'] != '[]'
        data['pt_backward'] = data['_pt_routes_backward'] != '[]'
        # for backward compatibility, to be removed later
        data['pt_bus'] = data['pt_forward'] or data['pt_backward']


def match_public_transit_simple(
        G, pt_routes,
        max_dist=400,
        max_dist_init=500,
        max_lattice_width=5,
        route_number_column='route',
        is_bidirectional_column='both_directions'
):
    """
    Matches public transit routes onto the network. The format must be as follows:
     * Each route section is represented as a separate linestring
     * An attribute defines whether the section is bidirectional or not

    This function is suited for manually created public transit routes datasets.

    Parameters
    ----------
    G : nx.MultiDiGraph
        street graph
    pt_routes : gpd.GeoDataFrame
        transit routes
    max_dist : int
        see settings of the leuven mapmatcher
    max_lattice_width : int
        see settings of the leuven mapmatcher
    route_number_column : str
        which column in the routes dataset holds the route number
    is_bidirectional_column : str
        which column says whether the route section is bidirectional (True) or not (False)
    """

    # create opposite geometry for each bidirectional segment
    selected_rows = pt_routes[pt_routes['both_directions'] == True]
    selected_rows['geometry'] = selected_rows['geometry'].apply(lambda geom: geom.reverse())
    selected_rows['direction'] = 2

    # append the opposite geometries to the original ones
    pt_routes = pd.concat([pt_routes, selected_rows], ignore_index=True)
    pt_routes = gpd.GeoDataFrame(pt_routes)
    pt_routes['direction'] = pt_routes['direction'].fillna(1)

    # use the zvv function for matching
    match_public_transit_zvv(
            G, pt_routes,
            max_dist=max_dist,
            max_dist_init=max_dist_init,
            max_lattice_width=max_lattice_width,
            route_number_column='route',
            direction_column='direction'
    )


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


def match_traffic_counts_npvm(G, traffic_counts_gpd):
    """
    Map match traffic counts from the Swiss NPVM to the street graph.

    Parameters
    ----------
    G: nx.MultiDiGraph
        street graph
    traffic_counts_gpd: gpd.GeoDataFrame
        dataset with the traffic flows, must follow the data structure of the Swiss NPVM 2017

    Returns
    -------

    """

    # remove links with zero traffic (otherwise they will distort the averages on the matched links)
    traffic_counts_gpd = traffic_counts_gpd[traffic_counts_gpd['DTV_ALLE'] > 0]

    # match
    enrichment.match_linestrings(G, traffic_counts_gpd, [
        {'source_column': 'DTV_ALLE',   'target_column': 'adt_avg',         'agg': 'avg' },
        {'source_column': 'DTV_ALLE',   'target_column': 'adt_max',         'agg': 'max' },
        {'source_column': 'FROMNODENO', 'target_column': 'npvm_fromnodeno', 'agg': 'list'},
        {'source_column': 'TONODENO',   'target_column': 'npvm_tonodeno',   'agg': 'list'}
    ])
