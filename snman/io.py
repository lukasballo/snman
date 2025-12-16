from . import osmnx_customized as oxc
from . import geometry_tools, space_allocation, utils, street_graph, lane_graph
from .constants import *
import geopandas as gpd
import pandas as pd
import shapely.ops
import shapely.geometry
import shapely
import xml.etree.ElementTree as ET
import itertools
import networkx as nx
import copy
import numpy as np
import json
import ast
import os
import importlib.resources as pkg_resources
import subprocess


# - CREATING BASIC DATASETS --------------------------------------------------------------------------------------------


def create_street_graph_from_OSM(perimeter_geometry, CRS_internal, return_raw=False):
    """
    Download OSM data and create a street graph from it. A typical SNMan project will start with this function.

    Parameters
    ----------
    perimeter_geometry: shapely.Polygon
    CRS_internal: int
        Which projected crs should be used for the street graph
    return_raw: bool
        If true this function will return a tuple of (G_raw, G) where G_raw is the raw street graph directly
        after the OSM import and G is the processed street graph


    Returns
    -------
    G: street_graph.StreetGraph
    """

    G = oxc.graph_from_polygon(
        perimeter_geometry,
        custom_filter=OSM_FILTER,
        simplify=True, simplify_strict=False, retain_all=True, one_edge_per_direction=False
    )

    G = street_graph.StreetGraph(G)

    if return_raw:
        G_raw = copy.deepcopy(G)
    else:
        G_raw = None

    # operations to prepare the imported street graph
    for uvk, data in G.edges.items():
        # ensure consistent data types: maxspeed
        maxspeed = data.get('maxspeed', '')
        data['maxspeed'] = int(maxspeed) if maxspeed.isdigit() else -1

        # ensure consistent data types: layer
        layer = data.get('layer', '')
        # isdigit only supports positive numbers, so we need to remove any '-' first
        data['layer'] = int(layer) if layer.lstrip('-').isdigit() else 0

    for i, data in G.nodes.items():
        # prepare traffic_signals attribute
        data['traffic_signals'] = 1 * (data.get('highway') == 'traffic_signals')

    hierarchy.add_hierarchy(G)
    street_graph.convert_crs(G, CRS_internal)
    street_graph.fill_wrong_edge_geometries(G)
    space_allocation.generate_lanes(G)
    street_graph.filter_lanes_by_modes(G, {MODE_TRANSIT, MODE_PRIVATE_CARS, MODE_CYCLING, MODE_CAR_PARKING})

    if return_raw:
        return G_raw, G
    else:
        return G


# - LOADING BASIC DATASETS ---------------------------------------------------------------------------------------------


def load_street_graph(
        edges_path, nodes_path, crs=DEFAULT_CRS,
        unstringify_attributes={
            KEY_LANES_DESCRIPTION: space_allocation.space_allocation_from_string,
            KEY_LANES_DESCRIPTION_AFTER: space_allocation.space_allocation_from_string,
            KEY_GIVEN_LANES_DESCRIPTION: space_allocation.space_allocation_from_string,
            KEY_FORCED_GIVEN_LANES_DESCRIPTION: space_allocation.space_allocation_from_string,
            'layers': ast.literal_eval
        }
):
    """
    Load a pre-generated street graph that has been saved as geofile (shp, gpkg, etc.)

    Parameters
    ----------
    edges_path : string
        path to the file containing the edges
    nodes_path : string
        path to the file containing the nodes
    unstringify_attributes : dict
        defines which function should be applied to each attribute to unstringify it
    crs : int
        target coordinate reference system of the imported street graph

    Returns
    -------
    G : street_graph.StreetGraph
    """

    edges_gdf = import_geofile_to_gdf(edges_path, crs=crs).replace(np.nan, None)
    edges_gdf['u'] = edges_gdf['u'].astype(int)
    edges_gdf['v'] = edges_gdf['v'].astype(int)
    edges_gdf['key'] = edges_gdf['key'].astype(int)
    edges_gdf.set_index(['u', 'v', 'key'], inplace=True)

    nodes_gdf = import_geofile_to_gdf(nodes_path, crs=crs).replace(np.nan, None)
    nodes_gdf['osmid'] = nodes_gdf['osmid'].astype(int)
    nodes_gdf.set_index('osmid', inplace=True)

    for key, fn in unstringify_attributes.items():
        if key in edges_gdf:
            edges_gdf[key] = edges_gdf[key].apply(fn)
        if key in nodes_gdf:
            nodes_gdf[key] = nodes_gdf[key].apply(fn)

    G = street_graph.street_graph_from_gdf(nodes_gdf, edges_gdf, crs=crs)

    return G


def load_lane_graph(
        edges_path, nodes_path, crs=DEFAULT_CRS,
        unstringify_attributes={
            'lane': space_allocation.lane_from_string,
            'layers': ast.literal_eval
        }
):
    """
    Load a pre-generated lane graph that has been saved as geofile (shp, gpkg, etc.)

    Parameters
    ----------
    edges_path : string
        path to the file containing the edges
    nodes_path : string
        path to the file containing the nodes
    unstringify_attributes : dict
        defines which function should be applied to each attribute to unstringify it
    crs : int
        target coordinate reference system of the imported street graph

    Returns
    -------
    L : lane_graph.LaneGraph
    """

    edges_gdf = import_geofile_to_gdf(edges_path, crs=crs).replace(np.nan, None)
    edges_gdf['u'] = edges_gdf['u'].astype(int)
    edges_gdf['v'] = edges_gdf['v'].astype(int)
    edges_gdf['key'] = edges_gdf['key'].astype(str)
    edges_gdf.set_index(['u', 'v', 'key'], inplace=True)

    nodes_gdf = import_geofile_to_gdf(nodes_path, crs=crs).replace(np.nan, None)
    nodes_gdf['osmid'] = nodes_gdf['osmid'].astype(int)
    nodes_gdf.set_index('osmid', inplace=True)

    for key, fn in unstringify_attributes.items():
        if key in edges_gdf:
            edges_gdf[key] = edges_gdf[key].apply(fn)
        if key in nodes_gdf:
            nodes_gdf[key] = nodes_gdf[key].apply(fn)

    L = lane_graph.LaneGraph(crs=crs)
    nodes_gdf.apply(lambda n: L.add_node(n.name, **n), axis=1)
    edges_gdf.apply(lambda e: L.add_edge(*e.name, **e), axis=1)

    return L


def import_geofile_to_gdf(file_path, crs=DEFAULT_CRS, index=None, filter_index=None, perimeter=None):
    """
    Import a geofile (shp, gpkg, etc.) or a parquet file (gzip) as a GeoDataFrame

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

    if file_path[-5:] == '.gzip':
        gdf = gpd.read_parquet(file_path)
    else:
        try:
            gdf = gpd.read_file(file_path).to_crs(crs)
        except Exception as e:
            # Handle invalid geometries that may occur with different geopandas versions
            # The default engine (pyogrio) is stricter about geometry validation
            # Try using fiona engine which is more lenient with invalid geometries
            try:
                gdf = gpd.read_file(file_path, engine='fiona').to_crs(crs)
            except Exception:
                # If fiona also fails, try to fix invalid geometries
                try:
                    gdf = gpd.read_file(file_path, engine='fiona')
                    # Fix invalid geometries using buffer(0) trick (works for most cases)
                    invalid_mask = ~gdf.geometry.is_valid
                    if invalid_mask.any():
                        gdf.loc[invalid_mask, 'geometry'] = gdf.loc[invalid_mask, 'geometry'].buffer(0)
                    gdf = gdf.to_crs(crs)
                except Exception:
                    # Last resort: re-raise the original exception
                    raise e

    if index is not None:
        gdf = gdf.set_index(index)

    if filter_index is not None:
        gdf = gdf.filter(items=filter_index, axis=0)

    if perimeter is not None:
        gdf.geometry = gdf.geometry.apply(lambda route: shapely.intersection(route, perimeter))
        gdf = gdf[~shapely.is_empty(gdf.geometry)]

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
    intersections['simplify'] = intersections['simplify'].fillna(1)
    return intersections


def load_rebuilding_regions(path, crs=DEFAULT_CRS, projects=None, only_active=False, filter_ids=None):
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

    if projects:
        rebuilding_regions = rebuilding_regions.query(f"project in {projects}")

    if only_active:
        rebuilding_regions = rebuilding_regions.query(f"active == True")

    if filter_ids:
        rebuilding_regions = rebuilding_regions.loc[filter_ids]

    # convert strings into lists
    if 'hierarchies_to_include' in rebuilding_regions.columns:
        rebuilding_regions['hierarchies_to_include'] = \
            rebuilding_regions['hierarchies_to_include'].apply(
                lambda x: set(x.split(',')) if isinstance(x, str) and len(x) > 0 else set()
            )
    if 'hierarchies_to_fix' in rebuilding_regions.columns:
        rebuilding_regions['hierarchies_to_fix'] = \
            rebuilding_regions['hierarchies_to_fix'].apply(
                lambda x: set(x.split(',')) if isinstance(x, str) and len(x) > 0 else set()
            )
    if 'hierarchies_with_cycling_lanes' in rebuilding_regions.columns:
        rebuilding_regions['hierarchies_with_cycling_lanes'] = \
            rebuilding_regions['hierarchies_with_cycling_lanes'].apply(
                lambda x: set(x.split(',')) if isinstance(x, str) and len(x) > 0 else set()
            )
    rebuilding_regions = rebuilding_regions.sort_values(['order'])
    return rebuilding_regions


def load_measurement_regions(
        path, crs=DEFAULT_CRS, only_active=True, filter_names=None, set_filter=None, name_column='name'
):
    measurement_regions = import_geofile_to_gdf(path, crs=crs, index=name_column)
    if filter_names:
        measurement_regions = measurement_regions.loc[filter_names]
    if only_active:
        measurement_regions = measurement_regions[measurement_regions['active'] == True]
    if set_filter:
        measurement_regions = measurement_regions[measurement_regions['set'] == set_filter]
    measurement_regions['area'] = measurement_regions.geometry.area
    return measurement_regions


def load_poi(path, perimeter=None, crs=DEFAULT_CRS):
    """
    Load points of interest from a geofile.

    Parameters
    ----------
    path : str
        Path to geofile
    perimeter : shapely.Polygon, optional
        Polygon to crop geometries to
    crs : int, optional
        Coordinate reference system (default: DEFAULT_CRS)

    Returns
    -------
    gpd.GeoDataFrame
        GeoDataFrame with points of interest
    """
    poi = import_geofile_to_gdf(path, perimeter=perimeter, crs=crs)
    return poi


def load_sensors(path):
    """
    Load traffic sensor data from CSV file.

    Parameters
    ----------
    path : str
        Path to CSV file containing sensor data

    Returns
    -------
    pd.DataFrame
        DataFrame with sensors indexed by 'id'
    """
    sensors = pd.read_csv(path).set_index('id')
    return sensors


# - LOADING ENRICHMENT DATASETS ----------------------------------------------------------------------------------------


def infer_parking_orientation(parking_spots):
    """
    Infer parking orientation based on distance to nearest parking spot.

    Parameters
    ----------
    parking_spots : gpd.GeoDataFrame
        GeoDataFrame with parking spot points

    Returns
    -------
    gpd.GeoDataFrame
        GeoDataFrame with added 'orientation' column indicating parking type
        (perpendicular, parallel, or None)
    """
    res = utils.join_nearest_points(parking_spots, parking_spots)
    res['orientation'] = res.apply(
        lambda row:
            LANETYPE_PARKING_PERPENDICULAR if row['distance_to_next_space'] < 3.25 else
            LANETYPE_PARKING_PARALLEL if row['distance_to_next_space'] < 12.00 else None
        , axis=1)
    return res


def load_parking_spots(path, crs=DEFAULT_CRS, perimeter=None):
    """
    Load a geofile with parking spaces. Each parking spot needs to be represented as a point.
    For Zurich, see this dataset:
    https://data.stadt-zuerich.ch/dataset/geo_oeffentlich_zugaengliche_strassenparkplaetze_ogd

    Parameters
    ----------
    path : str
    crs : int
    perimeter : shp.Polygon
        crop geometries to a perimeter polygon

    Returns
    -------
    gpd.GeoDataFrame
    """

    parking_spots = import_geofile_to_gdf(path, crs=crs, perimeter=perimeter)
    parking_spots = infer_parking_orientation(parking_spots)
    return parking_spots


def load_lane_edits(path, lane_columns=('lanes',), crs=DEFAULT_CRS):
    """
    Loads a geofile with manual lane edits

    Parameters
    ----------
    path: str
    lane_columns: list
        which columns contain lane lists, these will be converted into SpaceAllocation
    crs: int

    Returns
    -------
    gpd.GeoDataFrame
    """

    lane_edits = import_geofile_to_gdf(path, crs=crs)

    # convert strings to SpaceAllocation
    for column in lane_columns:
        lane_edits[column] = lane_edits[column].apply(
            lambda x: space_allocation.space_allocation_from_string(x)
        )

    return lane_edits


def load_hierarchy_edits(path, crs=DEFAULT_CRS):
    """
    Loads a geofile with manual hierarchy edits

    Parameters
    ----------
    path: str
    crs: int

    Returns
    -------
    gpd.GeoDataFrame
    """

    lane_edits = import_geofile_to_gdf(path, crs=crs)
    return lane_edits


def load_public_transit_routes_zvv(path, perimeter=None, crs=DEFAULT_CRS):
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


def load_traffic_counts_npvm(path, crs=DEFAULT_CRS):
    """
    Load traffic counts from NPVM dataset.

    Parameters
    ----------
    path : str
        Path to geofile with traffic counts
    crs : int, optional
        Coordinate reference system (default: DEFAULT_CRS)

    Returns
    -------
    gpd.GeoDataFrame
        GeoDataFrame with traffic count data
    """
    return gpd.read_file(path).to_crs(crs)


def load_street_polygons_zurich(path, crs=DEFAULT_CRS, only_roads=True):
    """
    Load the street polygon geometries from the official land surveying data of Canton Zurich:
    https://www.stadt-zuerich.ch/geodaten/download/10016,
    layer 'Bodenbedeckung_BoFlaeche_Area'

    Parameters
    ----------
    path: str
    crs: int
    only_roads: bool
        if False, all polygon geometries will be loaded (not recommended)

    Returns
    -------
    gpd.GeoDataFrame
    """

    street_polygons = gpd.read_file(path).to_crs(crs).set_index('OBJECTID')
    if only_roads:
        street_polygons = street_polygons.query('Art == 8')
    return street_polygons


def load_access_needs(path, crs=DEFAULT_CRS):
    """
    Load the dataset with access_needs.

    Parameters
    ----------
    path: str
    crs: int

    Returns
    -------
    gpd.GeoDataFrame
    """

    access_needs = gpd.read_file(path).to_crs(crs)
    return access_needs


def load_statpop(path, crs=DEFAULT_CRS):
    """
    Load the statpop dataset

    Parameters
    ----------
    path: str
    crs: int

    Returns
    -------
    gpd.GeoDataFrame
    """

    statpop = gpd.read_file(path).to_crs(crs)
    return statpop



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


def export_graph(
        G, path_edges, path_nodes, edge_columns=None, node_columns=None,
        stringify_attributes=[], reset_index=False, crs=4326
):
    """
    Export a graph as geofiles of nodes and edges

    Parameters
    ----------
    G: nx.Graph, nx.MultiGraph, nx.DiGraph, nx.MultiDiGraph
    path_edges: str
    path_nodes: str
    edge_columns: list
        which columns should be included, if None, all will be included
    node_columns: list
        dito
    stringify_attributes: list
        The attributes in this list will be converted to strings both in edges and nodes.
        This is useful if they have types that are incompatible with geofiles, such as lists.
    crs: int
    """

    if len(G.edges) == 0:
        return

    H = copy.deepcopy(G)

    # stringify node ids
    nx.relabel_nodes(H, {i: str(i) for i in list(H.nodes)}, copy=False)

    # convert the graph to geodataframes
    nodes, edges = oxc.graph_to_gdfs(H)

    # remove columns
    if edge_columns:
        edges = edges[list(set(set(edges.columns) & set(edge_columns)).union({'geometry'}))]
    if node_columns:
        nodes = edges[list(set(set(nodes.columns) & set(node_columns)).union({'geometry'}))]

    if reset_index:
        edges.reset_index(inplace=True)
        nodes.reset_index(inplace=True)

    # stringify attributes
    for attr in stringify_attributes:
        if attr in edges:
            edges[attr] = edges[attr].apply(lambda x: str(x))
        if attr in nodes:
            nodes[attr] = nodes[attr].apply(lambda x: str(x))

    # write files
    export_gdf(edges, path_edges, crs=crs)
    export_gdf(nodes, path_nodes, crs=crs)


def export_street_graph(
    G, path_edges, path_nodes, edge_columns=None, node_columns=None,
    stringify_attributes=(
            KEY_LANES_DESCRIPTION, KEY_LANES_DESCRIPTION_AFTER,
            KEY_GIVEN_LANES_DESCRIPTION, KEY_FORCED_GIVEN_LANES_DESCRIPTION,
            KEY_SENSORS_FORWARD, KEY_SENSORS_BACKWARD,
            '_extension', '_intermediary_nodes', 'highway', 'layers'
    ),
    stringify_additional_attributes=(),
    crs=4326
):
    """
    Exports a street graph to geofile. See `export_graph()` for details.
    Parameters
    ----------
    G: street_graph.StreetGraph
    path_edges: str
    path_nodes: str
    edge_columns: list
    node_columns: list
    stringify_attributes: list
    crs: int
    """
    export_graph(
        G, path_edges, path_nodes, edge_columns=edge_columns, node_columns=node_columns,
        stringify_attributes=list(stringify_attributes) + list(stringify_additional_attributes),
        crs=crs
    )


def export_lane_graph(
    L, path_edges, path_nodes, edge_columns=None, node_columns=None,
    stringify_attributes=('lane', 'highway', 'layers'),
    crs=4326
):
    """
    Exports a lane graph to geofile. See `export_graph()` for details.
    Parameters
    ----------
    L: lane_graph.LaneGraph
    path_edges: str
    path_nodes: str
    edge_columns: list
    node_columns: list
    stringify_attributes: list
    crs: int
    """
    export_graph(
        L, path_edges, path_nodes, edge_columns=edge_columns, node_columns=node_columns,
        stringify_attributes=stringify_attributes, crs=crs
    )


def export_access_graph(A, path_edges, path_nodes, stringify_attributes=[]):
    """
    Export access graph to files.

    Parameters
    ----------
    A : AccessGraph
        Access graph to export
    path_edges : str
        Path for edges file
    path_nodes : str
        Path for nodes file
    stringify_attributes : list, optional
        List of attributes to convert to strings (default: [])

    Returns
    -------
    None
    """
    pass


def export_lane_geometries(L, path_edges, path_nodes, scaling=1, include_opposite_direction=True, crs=4326):
    """
    Export a geofile with individual lane geometries. This is helpful for visualization purposes.

    Parameters
    ----------
    L: lane_graph.LaneGraph
    path_edges: str
    path_nodes: str
    scaling: float
    include_opposite_direction: bool
        if True bidirectional lanes will be included with both directions
    crs: int
    """

    M = copy.deepcopy(L)

    # remove opposite direction edges of bidirectional lanes
    for uvk, data in L.edges.items():
        if data['instance'] != 1 and not include_opposite_direction:
            M.remove_edge(*uvk)

    for uvk, data in M.edges.items():

        lane = data['lane']

        # create attributes for easy visualization
        data['width'] = lane.width
        data['lanetype'] = lane.lanetype
        data['direction'] = lane.direction
        data['status'] = lane.status
        data['horizontal_position'] = lane_graph.get_horizontal_position_of_lane(M, *uvk)

        # scale widths
        data['width_scaled'] = data['width'] * scaling
        data['horizontal_position_scaled'] = data['horizontal_position'] * scaling

        if data['horizontal_position_scaled'] != 0:
            if data['backward'] == 1:
                offset = data['horizontal_position_scaled']
            else:
                offset = -data['horizontal_position_scaled']
            if -np.inf < offset and offset < np.inf:
                try:
                    data['geometry'] = data['geometry'].offset_curve(offset)

                    # In older geos versions, the linestring gets reversed when the offset is negative,
                    # see here: https://shapely.readthedocs.io/en/stable/reference/shapely.offset_curve.html
                    # In that case, we have to reverse it back:
                    if offset < 0 and shapely.geos_version < (3,11,0):
                        data['geometry'] = shapely.reverse(data['geometry'])

                except shapely.errors.GEOSException:
                    print(uvk, 'geometry offset failed')
                except AttributeError:
                    print(uvk, 'geometry offset failed')


    export_lane_graph(M, path_edges, path_nodes, crs=crs)


def export_HLA(path, step, H=None, L=None, A=None, B=None, C=None, scaling_factor=1, export_crs=4326):
    if H is not None:
        # Export street graph
        export_street_graph(
            H,
            os.path.join(path, f'{step}_H_edges.gpkg'),
            os.path.join(path, f'{step}_H_nodes.gpkg'),
            crs=export_crs
        )

    if L is not None:
        # Export lane geometries
        export_lane_geometries(
            L,
            os.path.join(path, f'{step}_L_edges.gpkg'),
            os.path.join(path, f'{step}_L_nodes.gpkg'),
            scaling=scaling_factor,
            crs=export_crs
        )

    if A is not None:
        # Save graph
        export_graph(
            A,
            os.path.join(path, f'{step}_A_edges.gpkg'),
            os.path.join(path, f'{step}_A_nodes.gpkg'),
            reset_index=True,
            stringify_attributes=['osmid', 'u', 'v', 'has_parking_spots'],
            crs=export_crs
        )

    if B is not None:
        # Save graph
        export_graph(
            B,
            os.path.join(path, f'{step}_B_edges.gpkg'),
            os.path.join(path, f'{step}_B_nodes.gpkg'),
            reset_index=True,
            stringify_attributes=['osmid', 'u', 'v', 'has_parking_spots'],
            crs=export_crs
        )

    if C is not None:
        # Save graph
        export_graph(
            C,
            os.path.join(path, f'{step}_C_edges.gpkg'),
            os.path.join(path, f'{step}_C_nodes.gpkg'),
            reset_index=True,
            stringify_attributes=['osmid', 'u', 'v', 'has_parking_spots'],
            crs=export_crs
        )


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
        modes=(MODE_TRANSIT, MODE_CYCLING, MODE_PRIVATE_CARS),
        overwrite_highway=False,
        dont_overwrite_highway=(),
        set_maxspeed_by_cost=False,
        floor_maxspeed=0.5,
        ceil_maxspeed=120,
        neutral_speed_kmh=None,
        raw_graph=False
):
    """
    Generates an OSM file from the street graph

    Parameters
    ----------
    G : street_graph.StreetGraph
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
    modes : list
        which modes should be included
    overwrite_highway : bool or str
        set all highway tag values to this value, do nothing if False,
        this setting is needed to encode a cycling network as travel lanes for routing with link costs in R5
    dont_overwrite_highway : list
        which highway tag values should not be overwritten
    set_maxspeed_by_cost : bool or str
        set maxspeed of links according to their length and a cost attribute defined here,
        do nothing if False
    floor_maxspeed : float
        the maxspeed will be no lower than this number
    ceil_maxspeed : float
        the maxspeed will be no higher than this number
    neutral_speed_kmh : int
        only valid if set_maxspeed_by_cost is set
    raw_graph : bool
        If True, ignore any lanes and use the graph as is, i.e. every edge will be converted into one osm way.
        Please note that filtering by mode and creating one-way links will not work.

    Returns
    -------
    None
    """

    H = copy.deepcopy(G)

    if not raw_graph:

        # keep only the relevant modes
        street_graph.filter_lanes_by_modes(
            H,
            modes,
            operator='or',
            lane_description_key=key_lanes_description
        )

        if as_oneway_links:
            H = street_graph.separate_edges_for_lane_directions(H, lanes_key=key_lanes_description)
        else:
            H = copy.deepcopy(H)

        # ensure right edge directions, tags and crs
        street_graph.organize_edge_directions(H, method='by_osm_convention', key_lanes_description=key_lanes_description)
        street_graph.add_edge_costs(H, lanes_description=key_lanes_description)
        space_allocation.update_osm_tags(H, lanes_description_key=key_lanes_description)

    street_graph.convert_crs(H, 'epsg:4326')

    # initial ID value for new OSM objects, avoid duplicity with graph node ids
    max_node_id = max(list(H.nodes))
    osm_id = itertools.count(max_node_id * 100)


    if overwrite_highway:
        for uvk, data in H.edges.items():
            if data['highway'] not in dont_overwrite_highway:
                data['highway'] = overwrite_highway

    # remove all maxspeed tags to avoid duplicities
    if set_maxspeed_by_cost:
        for uvk, data in H.edges.items():
            if 'maxspeed' in data.keys():
                del data['maxspeed']

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
                node_points[node_id+1] = shapely.Point(H.nodes[node_id]['x'], H.nodes[node_id]['y'])
                # add node along the way, use the existing node id
                ET.SubElement(way, 'nd', attrib={'ref': str(node_id+1)})

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
            if tag == 'maxspeed' and set_maxspeed_by_cost:
                if neutral_speed_kmh:
                    neutral_speed_ms = neutral_speed_kmh / 3.6
                else:
                    neutral_speed_ms = 1
                tag_value = 3.6 * neutral_speed_ms * data['length'] / data[set_maxspeed_by_cost]
                # avoid 0 speeds
                tag_value = max([tag_value, floor_maxspeed])
                tag_value = min([tag_value, ceil_maxspeed])
                tag_value = str(tag_value)
            elif data.get(tag, None) is not None:
                tag_value = str(data.get(tag, ''))
            else:
                continue

            ET.SubElement(way, 'tag', attrib={
                'k': tag,
                'v': tag_value
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
                value = str(data.get('cost_' + key_lanes_description + '_' + mode + '_' + direction))
                ET.SubElement(way, 'tag', attrib={'k': tag, 'v': value})

    # add nodes
    node_points = dict(sorted(node_points.items()))
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


def osm_to_pbf(osm_file, dry_run=False):
    """
    Converts an osm file to pbf, using the osmconvert.exe utility.

    Using the wizzard of the utility:
        1. copy the osm file into this directory
        2. run osmconvert64-0.8.8p.exe¨
        3. press 'a' and enter
        4. enter the osm file name and press enter
        5. press 1
        6. press 3

    Alternatively, you can use the following command:

        osmconvert before_oneway_links.osm --out-pbf -o=before_oneway_links_01.pbf
        Parameters
    ----------
    osm_file: str
        Path to the osm file. The resulting pbf file will be written to the same directory
    dry_run: bool
        if true, output the command string instead of running osmconvert
    """

    current_file_path = os.path.abspath(__file__)
    current_directory = os.path.dirname(current_file_path)

    osmconvert = os.path.join(current_directory, 'osmconvert.exe')
    command = f'"{osmconvert}" "{osm_file}" --out-pbf -o="{osm_file}.pbf"'
    if dry_run:
        return command
    else:
        subprocess.run(command, capture_output=True, text=True)

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
