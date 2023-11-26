import shapely

from . import osmnx_customized as oxc
from .constants import *
import pyproj
import functools
import pandas as pd
import geopandas as gpd
import networkx as nx

from . import constants
from . import io
from . import street_graph
from . import lane_graph
from . import hierarchy
from . import enrichment
from . import simplification
from . import merge_edges
from . import graph
from . import distribution
from . import rebuilding
from . import space_allocation
from . import stats
from . import street_graph_node
from . import street_graph_edge


def get_street_graph(
        perimeter,
        crs,
        given_intersections_gdf,
        sensors_df=None,
        simplification_radius_for_edge_geometries=35,
        intersection_tolerance=10,
        simplification_iterations=3,
        osm_filter=OSM_FILTER,
        perimeter_crs=None
):
    """
    Create a snman street graph from OpenStreetMap

    Parameters
    ----------
    perimeter : shapely.Polygon
        which area should be used
    crs : int
        select a geometric projection that provides acceptable accuracy in the chosen perimeter
    given_intersections_gdf : gpd.GeoDataFrame
        manually defined intersections for cases where the automatic detection fails
    sensors_df : pd.DataFrame
        a list of traffic sensors to be matched onto the raw street graph before simplification
    simplification_radius_for_edge_geometries : int
    intersection_tolerance : int
    simplification_iterations : int
    osm_filter : list
    perimeter_crs : int
        crs of the provided perimeter, if None, than the main crs will be used

    Returns
    -------
    G : nx.MultiDiGraph

    """

    # transform perimeter to 4326 if necessary
    perimeter_crs = crs if perimeter_crs is None else perimeter_crs
    if perimeter_crs != 4326:
        transformer = pyproj.Transformer.from_crs(perimeter_crs, 4326, always_xy=True)
        perimeter = shapely.ops.transform(transformer.transform, perimeter)

    print('Get data from OSM server')
    # At this step, simplification means only removing degree=2 edges
    G = oxc.graph_from_polygon(
        # set the perimeter here
        perimeter,
        custom_filter=osm_filter,
        simplify=True, simplify_strict=False, retain_all=True, one_edge_per_direction=False
    )

    print('Prepare graph')
    street_graph.prepare_graph(G)

    print('Convert CRS of street graph to 2056')
    street_graph.convert_crs(G, crs)

    print('Identify hierarchy')
    # split the edges into hierarchy categories, such as main roads, local roads, etc.
    hierarchy.add_hierarchy(G)

    print('Generate lanes')
    # interpreting the OSM tags into a collection of lanes on each edge
    space_allocation.generate_lanes(G)

    if sensors_df is not None:
        print('Load sensors and assign them to edges in the raw street graph')
        enrichment.match_sensors(G, sensors_df)

    print('Simplify edge geometries (1/2)')
    simplification.simplify_edge_geometries(G, 5)

    for i in range(simplification_iterations):
        print('ITERATION', i)

        print('Detect intersections')
        intersections_gdf = simplification.merge_nodes_geometric(
            G, intersection_tolerance,
            given_intersections_gdf=given_intersections_gdf
        )

        print('Split through edges in intersections')
        simplification.split_through_edges_in_intersections(G, intersections_gdf)

        print('Detect intersections (repeat to ensure that no points are outside of intersections)')
        intersections_gdf = simplification.merge_nodes_geometric(
            G, intersection_tolerance,
            given_intersections_gdf=given_intersections_gdf
        )

        print('Add layers to nodes')
        simplification.add_layers_to_nodes(G)

        print('Add connections between components in intersections')
        simplification.connect_components_in_intersections(G, intersections_gdf, separate_layers=True)

        print('Consolidate intersections')
        G = simplification.consolidate_intersections(
            G, intersections_gdf,
            reconnect_edges=True
        )

        print('Merge consecutive edges')
        merge_edges.merge_consecutive_edges(G)

        print('Merge parallel edges')
        merge_edges.merge_parallel_edges(G)

        print('Update precalculated attributes')
        street_graph.update_precalculated_attributes(G)

    print('Simplify edge geometries (2/2)')
    simplification.simplify_edge_geometries(G, 35)

    print('Keep only the largest weakly connected component')
    G = graph.keep_only_the_largest_connected_component(G, weak=True)

    print('Add lane stats to edges')
    space_allocation.generate_lane_stats(G)

    print('Update OSM tags')
    # to match the simplified and merged edges
    space_allocation.update_osm_tags(G)

    print('Update street counts per node')
    spn = oxc.stats.count_streets_per_node(G, nodes=G.nodes)
    nx.set_node_attributes(G, values=spn, name="street_count")

    return G


def add_elevation(G, raster, cpus=1):
    print('Add elevation')
    G = oxc.elevation.add_node_elevations_raster(G, raster, cpus=1)
    G = oxc.elevation.add_edge_grades(G, add_absolute=False)
    return G


def add_public_transport(G, pt_network):
    enrichment.match_pt(G, pt_network)


def save_street_graph(G, path):
    io.export_street_graph(G, path + '_edges.gpkg', path + '_nodes.gpkg')


def save_street_graph_as_osm(G, path, key_lanes='before', osm_tags=EXPORT_OSM_TAGS):

    if key_lanes == 'before':
        key_lanes = KEY_LANES_DESCRIPTION
    elif key_lanes == 'after':
        key_lanes = KEY_LANES_DESCRIPTION_AFTER

    io.export_osm_xml(G, path + '.osm', osm_tags, uv_tags=True, tag_all_nodes=False, key_lanes_description=key_lanes)


def save_lane_geometries(G, path, scaling=1, key_lanes='before'):

    if key_lanes == 'before':
        key_lanes = KEY_LANES_DESCRIPTION
    elif key_lanes == 'after':
        key_lanes = KEY_LANES_DESCRIPTION_AFTER

    io.export_street_graph_with_lanes(
        G,
        key_lanes,
        path + '.shp',
        scaling=scaling
    )


