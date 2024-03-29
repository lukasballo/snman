from . import osmnx_customized as oxc
from . import geometry_tools, space_allocation, utils, street_graph, street_graph_edge, street_graph_node, graph
from .constants import *
import geopandas as gpd
import pandas as pd
import pyproj
import shapely.ops
import shapely.geometry
import shapely
import xml.etree.ElementTree as ET
import itertools
import networkx as nx
import copy
import numpy as np
import json


# - LOADING BASIC DATASETS ---------------------------------------------------------------------------------------------

def load_street_graph(edges_path, nodes_path, crs=DEFAULT_CRS, recreate_iterables=True):
    """
    Load a pre-generated street graph that has been saved as geofile (shp, gpkg, etc.)

    Parameters
    ----------
    edges_path : string
        path to the file containing the edges
    nodes_path : string
        path to the file containing the nodes
    crs : int
        target coordinate reference system of the imported street graph

    Returns
    -------
    G : nx.MultiGraph
        street graph
    """

    edges_gdf = import_geofile_to_gdf(edges_path, index=['u', 'v', 'key'], crs=crs).replace(np.nan, None)
    nodes_gdf = import_geofile_to_gdf(nodes_path, index='osmid', crs=crs).replace(np.nan, None)

    if recreate_iterables:
        _iterable_columns_from_strings(edges_gdf, {'ln_desc', 'ln_desc_after', 'given_lanes'}, separator=' | ')
        _iterable_columns_from_strings(edges_gdf, {'sensors_forward', 'sensors_backward', '_intermediate_nodes'},
                                       method='str')
        _iterable_columns_from_strings(nodes_gdf, {'layers'}, separator=',')

    G = nx.MultiDiGraph(crs=crs)
    nodes_gdf.apply(lambda n: G.add_node(n.name, **n), axis=1)
    edges_gdf.apply(lambda e: G.add_edge(*e.name, **e), axis=1)

    return G


def import_geofile_to_gdf(file_path, crs=DEFAULT_CRS, index=None, filter_index=None, perimeter=None):
    """
    Import a geofile (shp, gpkg, etc.) as a GeoDataFrame

    Parameters
    ----------
    file_path : string
    crs : int
        target coordinate reference system of the imported geodataframe
    index : str or list
        which column(s) should be used as index
    filter_index : list
        which rows should be included, by index values
    perimeter : polygon
        will be used to crop the imported geometries

    Returns
    -------
    gdf : gpd.GeoDataFrame
        resulting GeoDataFrame
    """

    gdf = gpd.read_file(file_path).to_crs(crs)
    if index is not None:
        gdf = gdf.set_index(index)

    if filter_index is not None:
        gdf = gdf.filter(items=filter_index, axis=0)

    if perimeter is not None:
        gdf = gdf.overlay(perimeter, how='intersection')

    return gdf


def load_perimeters(path, filter=None, crs=DEFAULT_CRS):
    """
    Load a geofile (shp, gpkg, etc.) with network perimeters. These will be used to download the data from OSM
    and prepare the simplified street graph

    The geofile must contain following columns:
        * geometry: polygon
        * id: str

    Parameters
    ----------
    path : str
    filter : list
        which perimeters should be loaded, e.g. ['zollikerberg']

    Returns
    -------
    perimeters : gpd.GeoDataFrame
    """

    perimeters = import_geofile_to_gdf(path, index='id', filter_index=filter, crs=crs)
    return perimeters


def load_regions(path, default_tolerance=None, street_graph=None, crs=DEFAULT_CRS):
    """
    Load a geofile (shp, gpkg, etc.) with regions. These will be used to override the network simplification settings
    in specified areas

    The geofile must contain following columns:
        * geometry: polygon
        * tolerance: integer (what radius should be applied to the intersections)

    Parameters
    ----------
    path : str
    default_tolerance : int
    street_graph : nx.MultiGraph or nx.MultiDiGraph

    Returns
    -------
    regions : gpd.GeoDataFrame
    """

    regions = import_geofile_to_gdf(path, crs=crs)

    # create a new region containing all points that don't belong to a region yet
    if default_tolerance and street_graph is not None:
        # create a convex hull around all node geometries + some buffer to be safe
        nodes_gdf = oxc.utils_graph.graph_to_gdfs(street_graph, edges=False)
        polygon = nodes_gdf.geometry.unary_union.convex_hull.buffer(1000)
        # merge all other polygons
        other_regions_polygons = geometry_tools.ensure_multipolygon(regions['geometry'].unary_union)
        # cut them out of this default polygon
        polygon = polygon.difference(other_regions_polygons)
        new_row = gpd.GeoDataFrame({'geometry': [polygon], 'tolerance': [default_tolerance]}, geometry='geometry')
        regions = pd.concat([regions, new_row], ignore_index=True)

    return regions


def load_intersections(path, crs=DEFAULT_CRS):
    """
    Load an geofile (shp, gpkg, etc.) with intersection polygons. These will be used to override the automatically
    detected intersections. It is useful in cases where the auto-detection delivers unsatisfactory results

    The geofile must contain following columns:
        * geometry: polygon
        * id: int

    Parameters
    ----------
    path : str

    Returns
    -------
    intersections : gpd.GeoDataFrame
    """

    polygons = import_geofile_to_gdf(path, crs=crs)
    # Duplicate the point geometries so that they get taken along after the spatial join
    intersections = polygons
    return intersections


def load_rebuilding_regions(path, crs=DEFAULT_CRS, only_active=True, filter_ids=None):
    """
    Load a geofile (shp, gpkg, etc.) with rebuilding regions. These will be used to define areas where the streets
    should be rebuilt. In each rebuilding region, you can also specify which street hierarchies should be considered
    and which should be kept the way they are.

    The geofile must contain following columns:
        * geometry: polygon
        * fid: int
        * hierarchies_to_include: str, separated by comma (which street hierarchies should be considered
          in the rebuilding process, example: '1_main_road,2_local_road')
        * hierarchies_to_fix: str, separated by comma (which street hierarchies should not be changed,
          example: '1_main_road,2_local_road')
        * keep_all_streets: bool
            (True: enforce that all streets are passable for cars, False: only all nodes will be accessible for cars)
        * order: int (in which order should the rebuilding steps be applied)
        * active: bool (whether each step should be executed or ignored)

    Parameters
    ----------
    path : str

    Returns
    -------
    gpd.GeoDataFrame
    """

    rebuilding_regions = import_geofile_to_gdf(path, crs=crs, index='id')

    if only_active:
        rebuilding_regions = rebuilding_regions.query(f"active == True")

    if filter_ids:
        rebuilding_regions = rebuilding_regions.loc[filter_ids]

    # convert strings into lists
    rebuilding_regions['hierarchies_to_include'] = \
        rebuilding_regions['hierarchies_to_include'].apply(
            lambda x: set(x.split(',')) if isinstance(x, str) and len(x) > 0 else set()
        )
    rebuilding_regions['hierarchies_to_fix'] = \
        rebuilding_regions['hierarchies_to_fix'].apply(
            lambda x: set(x.split(',')) if isinstance(x, str) and len(x) > 0 else set()
        )
    rebuilding_regions = rebuilding_regions.sort_values(['order'])
    return rebuilding_regions


def load_measurement_regions(path, crs=DEFAULT_CRS, only_active=True, filter_names=None, set_filter=None):
    measurement_regions = import_geofile_to_gdf(path, crs=crs, index='name')
    if filter_names:
        measurement_regions = measurement_regions.loc[filter_names]
    if only_active:
        measurement_regions = measurement_regions[measurement_regions['active'] == True]
    if set_filter:
        measurement_regions = measurement_regions[measurement_regions['set'] == set_filter]
    measurement_regions['area'] = measurement_regions.geometry.area
    return measurement_regions


def load_poi(path, perimeter=None, crs=DEFAULT_CRS):
    poi = import_geofile_to_gdf(path, perimeter=perimeter, crs=crs)
    return poi


def load_sensors(path):
    sensors = pd.read_csv(path).set_index('id')
    return sensors


# - LOADING ENRICHMENT DATASETS ----------------------------------------------------------------------------------------


def infer_parking_orientation(parking_spots):
    res = utils.join_nearest_points(parking_spots, parking_spots)
    res['orientation'] = res.apply(
        lambda row:
            LANETYPE_PARKING_PERPENDICULAR if row['distance_to_next_space'] < 3.25 else
            LANETYPE_PARKING_PARALLEL if row['distance_to_next_space'] < 12.00 else None
        , axis=1)
    return res


def load_parking_spots(path, crs=DEFAULT_CRS):
    """
    Load a geofile with parking spaces. Each parking spot needs to be represented as a point.
    For Zurich, see this dataset:
    https://data.stadt-zuerich.ch/dataset/geo_oeffentlich_zugaengliche_strassenparkplaetze_ogd

    Parameters
    ----------
    path : str
    crs : int

    Returns
    -------
    gpd.GeoDataFrame
    """

    parking_spots = import_geofile_to_gdf(path, crs=crs)
    parking_spots = infer_parking_orientation(parking_spots)
    return parking_spots


def load_lane_edits(path, crs=DEFAULT_CRS):
    lane_edits = import_geofile_to_gdf(path, crs=crs)
    # convert the lists encoded as strings into lists
    lane_edits['add_lanes'] = (lane_edits['add_lanes']
                               .apply(lambda x: x.split(',') if isinstance(x, str) else []))
    lane_edits['replace_lanes'] = (lane_edits['replace_lanes']
                                   .apply(lambda x: x.split(',') if isinstance(x, str) else []))
    lane_edits['fid'] = lane_edits.index
    return lane_edits


def load_public_transit_routes(path, perimeter=None, crs=DEFAULT_CRS):
    """
    Load a geofile with public transit routes. Each route is represented by a (Multi)LineString.
    For Zurich, see this dataset: 'Linien des öffentlichen Verkehrs OGD (kantonaler Datensatz)'
    https://data.stadt-zuerich.ch/dataset/ktzh_linien_des_oeffentlichen_verkehrs__ogd_

    Parameters
    ----------
    path : str
    crs : int

    Returns
    -------
    gpd.GeoDataFrame
    """

    pt_routes = import_geofile_to_gdf(path, crs=crs)

    if perimeter:
        pt_routes.geometry = pt_routes.geometry.apply(lambda route: shapely.intersection(route, perimeter))
        pt_routes = pt_routes[~shapely.is_empty(pt_routes.geometry)]

    return pt_routes


# - OTHER --------------------------------------------------------------------------------------------------------------


def _get_nodes_within_polygon(G, polygon):
    """
    Return nodes of a graph that fit into a polygon

    TODO: Merge with a corresponding function of osmnx

    Parameters
    ----------
    G : nx.MultiDiGraph or nx.MultiGraph
    polygon : shapely.geometry.Polygon

    Returns
    -------
    nodes : list
    """

    nodes_gdf = oxc.graph_to_gdfs(G, edges=False)
    nodes_gdf = nodes_gdf[nodes_gdf.within(polygon)]
    return set(nodes_gdf.index.values)


def export_street_graph(
        G, path_edges, path_nodes, edge_columns=None, node_columns=None, cast_attr=True, crs=None,
        lane_keys=(KEY_LANES_DESCRIPTION, KEY_LANES_DESCRIPTION_AFTER, KEY_GIVEN_LANES_DESCRIPTION)
):
    """
    Export street graph as a geofile (shp, gpkg, etc.)

    Parameters
    ----------
    G : nx.MultiGraph or nx.MultiDiGraph
        street graph
    path_edges : str
        where should the edges file be saved
    path_nodes : str
        where should the nodes file be saved
    edge_columns : list
        which columns should be included in the edges file (None -> include all columns)
    node_columns : list
        dito for nodes file
    crs : int

    Returns
    -------
    None
    """

    H = copy.deepcopy(G)

    if cast_attr:
        graph.apply_function_to_each_edge(
            H, lambda H, uvk: street_graph_edge.cast_attributes_for_export(H, uvk, lane_keys)
        )
        graph.apply_function_to_each_node(
            H, street_graph_node.cast_attributes_for_export
        )

    nodes, edges = oxc.graph_to_gdfs(H)

    if edge_columns:
        edges = edges[list(set(set(edges.columns) & set(edge_columns)).union({'geometry'}))]

    if node_columns:
        nodes = edges[list(set(set(nodes.columns) & set(node_columns)).union({'geometry'}))]

    # Limit every attribute to 10 characters to match the SHP format restrictions
    if path_edges.split()[-1] == 'shp':
        edges.columns = [column[0:10] for column in edges.columns]

    # write files
    export_gdf(edges, path_edges, crs=crs)
    export_gdf(nodes, path_nodes, crs=crs)


def export_street_graph_with_lanes(G, lanes_attributes, path, scaling=1, crs=None):
    """
    Export a geofile with individual lane geometries. This is helpful for visualization purposes.

    Parameters
    ----------
    G : nx.MultiGraph or nx.MultiDiGraph
        street graph
    lanes_attributes : str, list
        which attribute should be used as a source of the lane configuration of each edge
    path : str
        where should the file be saved
    scaling : float
        optional scaling of the lane width, for an optimized visualization
    crs : int

    Returns
    -------
    None
    """
    # TODO: Build on top of a lane graph

    if type(lanes_attributes) == str:
        lanes_attributes = [lanes_attributes]

    # Create empty list for the lanes
    lanes_list = []

    for lanes_attribute in lanes_attributes:

        for uvk, data in G.edges.items():
            given_total_width = 0

            # Reconstruct total width of given lanes
            for lane in data.get(lanes_attribute, []):
                lane_properties = space_allocation._lane_properties(lane)
                given_total_width += lane_properties.width

            offset = -given_total_width / 2
            for lane in data.get(lanes_attribute, []):
                lane_properties = space_allocation._lane_properties(lane)

                centerline_offset = offset + lane_properties.width / 2
                offset += lane_properties.width
                geom = data.get('geometry')

                if geom and round(centerline_offset, 1) != 0:
                    geom = geom.parallel_offset(centerline_offset * scaling, 'right')
                    # the above function reverses direction when the offset is positive, this steps reverses it back
                    if centerline_offset > 0 and not geom.is_empty:
                        # geom.coords = list(geom.coords)[::-1]
                        shapely.ops.substring(geom, 1, 0, normalized=True)
                        pass

                lanes_list.append({
                    'lanes_key': lanes_attribute,
                    'type': lane_properties.lanetype,
                    'direction': lane_properties.direction,
                    'descr': lane,
                    'width_m': lane_properties.width * scaling,
                    'layer': data.get('layer'),
                    'length': data.get('length'),
                    'geometry': geom
                })

    lanes_gdf = gpd.GeoDataFrame(
        lanes_list,
        columns=['lanes_key', 'type', 'direction', 'descr', 'width_m', 'layer', 'length', 'geometry'],
        geometry='geometry',
        crs=G.graph['crs']
    )

    export_gdf(lanes_gdf, path, crs=crs)


def export_lane_geometries(L, path, scaling=1, crs=None):
    nodes, edges = oxc.graph_to_gdfs(L)
    print(nodes.head())
    print(edges.head())


def export_gdf(gdf, path, columns=[], crs=None):
    """
    Export a GeoDataFrame into a geofile.
    A wrapper around GeoDataFrame.to_file() but with an optional column filter.

    Parameters
    ----------
    gdf : gpd.GeoDataFrame
    path : str
    columns : list
        which columns should be included ([] -> include all columns)
    crs : int

    Returns
    -------
    None
    """

    # convert crs
    if crs:
        gdf = gdf.to_crs(crs)

    if columns == []:
        gdf.to_file(path)
    else:
        gdf[columns].to_file(path)


def export_osm_xml(
        G, path, tags,
        uv_tags=False,
        lanes_tag=False,
        tag_all_nodes=False,
        key_lanes_description=KEY_LANES_DESCRIPTION,
        as_oneway_links=False,
        modes=(MODE_TRANSIT, MODE_CYCLING, MODE_PRIVATE_CARS)
):
    """
    Generates an OSM file from the street graph

    Parameters
    ----------
    G : nx.MultiDiGraph
        street graph
    path : str
    tags : list
        which OSM tags should be included
    uv_tags : bool
        include special tags for the start and end node of each edge (for debugging)
    tag_all_nodes : bool
        add a special tag to each node so that all nodes appear as points when imported in QGIS (for debugging)
    key_lanes_description : str
        which attribute should be used for the lanes
    as_oneway_links : bool
        if true, every link with bidirectional traffic will be exported as two one-way links

    Returns
    -------
    None
    """

    # keep only the relevant modes
    H = street_graph.filter_lanes_by_modes(
        G,
        modes,
        lane_description_key=key_lanes_description
    )

    if as_oneway_links:
        H = street_graph.separate_edges_for_lane_directions(H, lanes_key=key_lanes_description)
    else:
        H = copy.deepcopy(H)

    # initial ID value for new OSM objects, avoid duplicity with graph node ids
    max_node_id = max(list(H.nodes))
    osm_id = itertools.count(max_node_id * 100)

    # ensure right edge directions, tags and crs
    street_graph.organize_edge_directions(H, method='by_osm_convention', key_lanes_description=key_lanes_description)
    space_allocation.update_osm_tags(H, lanes_description_key=key_lanes_description)
    street_graph.convert_crs(H, 'epsg:4326')

    # create the overall structure of the XML
    tree = ET.ElementTree('tree')

    # add OSM metadata
    osm = ET.Element('osm', attrib={
        'version': '0.6',
        'generator': 'osmium/1.14.0'
    })

    # add bounds
    lats = [data.get('y') for id, data in H.nodes.items()]
    lons = [data.get('x') for id, data in H.nodes.items()]
    ET.SubElement(osm, 'bounds', attrib={
        'minlat': str(min(lats)),
        'minlon': str(min(lons)),
        'maxlat': str(max(lats)),
        'maxlon': str(max(lons))
    })

    # add ways
    ways = []
    node_points = {}
    for uvk, data in H.edges.items():

        way = ET.Element('way', attrib={
            'id': str(next(osm_id)),  # assign a new unique osm id to each way
            'version': '1',
            'timestamp': '2000-01-01T00:00:00Z'
        })

        ways.append(way)

        # way nodes
        for i in [0, 'intermediary_nodes', 1]:

            # first and last node
            if i in {0, 1}:
                node_id = uvk[i]
                # avoid id=0
                node_points[node_id + 1] = shapely.Point(H.nodes[node_id]['x'], H.nodes[node_id]['y'])
                # add node along the way, use the existing node id
                ET.SubElement(way, 'nd', attrib={'ref': str(node_id + 1)})

            # intermediary nodes
            elif i == 'intermediary_nodes':
                linestring = data.get('geometry')
                # iterate over all points along the linestring geometry but exclude the first and last one
                for point in linestring.coords[1:-1]:
                    node_id = next(osm_id)
                    node_points[node_id] = shapely.geometry.Point(point)
                    # add node along the way, assign a new unique id
                    ET.SubElement(way, 'nd', attrib={'ref': str(node_id)})

        # way tags
        for tag in tags:
            # skip if the tag is not defined for this way
            if data.get(tag, None) is None:
                continue
            ET.SubElement(way, 'tag', attrib={
                'k': tag,
                'v': str(data.get(tag, ''))
            })

        if uv_tags:
            ET.SubElement(way, 'tag', attrib={'k': '_u', 'v': str(uvk[0])})
            ET.SubElement(way, 'tag', attrib={'k': '_v', 'v': str(uvk[1])})
            ET.SubElement(way, 'tag', attrib={'k': '_key', 'v': str(uvk[2])})

        if lanes_tag:
            ET.SubElement(way, 'tag', attrib={'k': 'lanes', 'v': str(data.get(key_lanes_description))})

        for mode in modes:
            for direction in [DIRECTION_BACKWARD, DIRECTION_FORWARD]:
                tag = 'cost_' + mode + '_' + direction
                value = str(data.get('cost_' + key_lanes_description + '_' + mode + '_' + direction, ''))
                ET.SubElement(way, 'tag', attrib={'k': tag, 'v': value})

    # add nodes
    for node_id, point in node_points.items():
        node = ET.SubElement(osm, 'node', attrib={
            'id': str(node_id),
            'version': '1',
            'timestamp': '2000-01-01T00:00:00Z',
            'lat': str(point.y),
            'lon': str(point.x)
        })

        if tag_all_nodes:
            ET.SubElement(node, 'tag', attrib={'k': '_node', 'v': 'true'})

    # append ways
    for way in ways:
        osm.append(way)

    # Save into xml file
    tree._setroot(osm)
    ET.indent(tree)
    tree.write(path, encoding='UTF-8', xml_declaration=True)


def _iterable_columns_from_strings(df, columns, method='separator', separator=','):
    """
    Converts selected columns of a DataFrame from strings into iterables

    Parameters
    ----------
    df : pd.DataFrame
    columns : list
        column names
    method : str
        how should the string be interpreted
            * 'separator' -> list from separated values
            * 'str' -> auto-interpret stringified iterable (not implemented yet)
    separator : str
        only used if method='separator'

    Returns
    -------
    None
    """
    for column in columns:
        if column in df:
            if method == 'separator':
                df[column] = df[column].apply(lambda x: x.split(separator))
            elif method == 'str':
                df[column] = df[column].apply(lambda x: json.loads(x) if x != 'nan' else [])
