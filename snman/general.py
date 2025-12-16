from snman import street_graph, io, hierarchy, space_allocation, simplification, enrichment, merge_edges, graph
from snman.constants import *
from snman import osmnx_customized as oxc

import shapely
import pyproj
import networkx as nx


def get_street_graph(
        perimeter,
        crs,
        given_intersections_gdf,
        sensors_df=None,
        simplification_radius_for_edge_geometries=35,
        intersection_tolerance=10,
        simplification_iterations=3,
        osm_filter=OSM_FILTER,
        perimeter_crs=None,
    elevation_file=None,
        export_raw_streetgraph=None,
    ):
    """
    Create a snman street graph from OpenStreetMap.

    Parameters
    ----------
    perimeter : shapely.Polygon
        Area to extract street network from
    crs : int
        Geometric projection CRS code that provides acceptable accuracy in the chosen perimeter
    given_intersections_gdf : gpd.GeoDataFrame
        Manually defined intersections for cases where the automatic detection fails
    sensors_df : pd.DataFrame, optional
        List of traffic sensors to be matched onto the raw street graph before simplification
    simplification_radius_for_edge_geometries : int, optional
        Radius for edge geometry simplification (default: 35)
    intersection_tolerance : int, optional
        Tolerance for intersection detection in meters (default: 10)
    simplification_iterations : int, optional
        Number of simplification iterations to perform (default: 3)
    osm_filter : list, optional
        OSM filter to use for data extraction (default: OSM_FILTER)
    perimeter_crs : int, optional
        CRS of the provided perimeter. If None, the main crs will be used
    elevation_file : str, optional
        Path to elevation raster file
    export_raw_streetgraph : str, optional
        Path prefix for exporting raw street graph (edges and nodes)

    Returns
    -------
    G : nx.MultiDiGraph
        Processed street graph
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

    if export_raw_streetgraph:
        print('Export raw street graph')
        # each street is one edge, the lanes are saved as an attribute
        io.export_street_graph(
            G,
            export_raw_streetgraph + '_edges.gpkg',
            export_raw_streetgraph + '_nodes.gpkg',
            crs=CRS_for_export
        )

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

    street_graph.surrogate_missing_edge_geometries(G)

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

    if elevation_file:
        print('Add elevation')
        G = add_elevation(G, elevation_file)

    return G


def add_elevation(G, raster, cpus=1):
    """
    Add elevation data to nodes and calculate edge grades.

    Parameters
    ----------
    G : nx.MultiDiGraph
        Street graph
    raster : str
        Path to elevation raster file
    cpus : int, optional
        Number of CPUs to use for processing (default: 1)

    Returns
    -------
    nx.MultiDiGraph
        Street graph with elevation attributes added
    """
    print('Add elevation')
    G = oxc.elevation.add_node_elevations_raster(G, raster, cpus=1)
    G = oxc.elevation.add_edge_grades(G, add_absolute=False)
    return G


def add_public_transport(G, pt_network):
    """
    Add public transport network data to the street graph.

    Parameters
    ----------
    G : nx.MultiDiGraph
        Street graph
    pt_network : Any
        Public transport network data
    """
    enrichment.match_pt(G, pt_network)


def save_street_graph(G, path):
    """
    Export street graph to GeoPackage files.

    Parameters
    ----------
    G : nx.MultiDiGraph
        Street graph to export
    path : str
        Base path for output files (will append '_edges.gpkg' and '_nodes.gpkg')
    """
    io.export_street_graph(G, path + '_edges.gpkg', path + '_nodes.gpkg')


def save_street_graph_as_osm(G, path, key_lanes='before', osm_tags=EXPORT_OSM_TAGS):
    """
    Export street graph to OSM XML format.

    Parameters
    ----------
    G : nx.MultiDiGraph
        Street graph to export
    path : str
        Base path for output file (will append '.osm')
    key_lanes : str, optional
        Which lane configuration to export: 'before' or 'after' (default: 'before')
    osm_tags : set, optional
        Set of OSM tags to include in export (default: EXPORT_OSM_TAGS)
    """
    if key_lanes == 'before':
        key_lanes = KEY_LANES_DESCRIPTION
    elif key_lanes == 'after':
        key_lanes = KEY_LANES_DESCRIPTION_AFTER

    io.export_osm_xml(G, path + '.osm', osm_tags, uv_tags=True, tag_all_nodes=False, key_lanes_description=key_lanes)


def save_lane_geometries(G, path, scaling=1, key_lanes='before'):
    """
    Export lane geometries to shapefile.

    Parameters
    ----------
    G : nx.MultiDiGraph
        Street graph to export
    path : str
        Base path for output file (will append '.shp')
    scaling : float, optional
        Scaling factor for lane geometries (default: 1)
    key_lanes : str, optional
        Which lane configuration to export: 'before' or 'after' (default: 'before')
    """
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

