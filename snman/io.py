from . import osmnx_customized as oxc
from . import geometry_tools, lanes
import geopandas as gpd
import pandas as pd
import pyproj
import shapely.ops
import shapely.geometry
import shapely
import xml.etree.ElementTree as ET
import itertools
import networkx as nx


def load_street_graph(edges_path, nodes_path, crs=2056):
    """
    Load a pre-generated street graph that has been saved as geofile (shp, gpkg, etc.)

    Parameters
    ----------
    edges_path : string
        path to the file containing the edges
    nodes_path : string
        path to the file containing the nodes
    crs : int
        number of the coordinate reference system

    Returns
    -------
    G : nx.MultiGraph
        street graph
    """

    edges_gdf = import_geofile_to_gdf(edges_path, index=['u', 'v', 'key'], crs=crs)
    nodes_gdf = import_geofile_to_gdf(nodes_path, index='osmid', crs=crs)
    _iterable_columns_from_strings(edges_gdf, {'ln_desc', 'ln_desc_after', 'given_lanes'}, separator=' | ')

    G = nx.MultiGraph(crs=crs)
    nodes_gdf.apply(lambda n: G.add_node(n.name, **n), axis=1)
    edges_gdf.apply(lambda e: G.add_edge(*e.name, **e), axis=1)

    return G


def import_geofile_to_gdf(file_path, crs=2056, index=None):
    """
    Import a geofile (shp, gpkg, etc.) as a GeoDataFrame

    Parameters
    ----------
    file_path : string
    crs : int
    index : str or list
        which column(s) should be used as index

    Returns
    -------
    gdf : gpd.GeoDataFrame
        resulting GeoDataFrame
    """

    gdf = gpd.read_file(file_path).to_crs(crs)
    if index:
        gdf = gdf.set_index(index)
    return gdf


def load_perimeters(path):
    """
    Load a geofile (shp, gpkg, etc.) with network perimeters. These will be used to download the data from OSM
    and prepare the simplified street graph

    The geofile must contain following columns:
        * geometry: polygon
        * id: str

    Parameters
    ----------
    path : str

    Returns
    -------
    perimeters : gpd.GeoDataFrame
    """

    perimeters = import_geofile_to_gdf(path, crs=4326, index='id')
    return perimeters


def load_regions(path, default_tolerance=None, street_graph=None):
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

    regions = import_geofile_to_gdf(path)

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


def load_intersections(path):
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

    polygons = import_geofile_to_gdf(path)
    # Duplicate the point geometries so that they get taken along after the spatial join
    intersections = polygons
    return intersections


def load_rebuilding_regions(path):
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
        * order: int (in which order should the rebuilding steps be applied)
        * active: bool (whether each step should be executed or ignored)

    Parameters
    ----------
    path : str

    Returns
    -------
    rebuilding_regions : gpd.GeoDataFrame
    """

    rebuilding_regions = import_geofile_to_gdf(path)
    # convert strings into lists
    rebuilding_regions['hierarchies_to_include'] = \
        rebuilding_regions['hierarchies_to_include'].apply(
            lambda x: set(x.split(',')) if isinstance(x, str) and len(x) > 0 else {}
        )
    rebuilding_regions['hierarchies_to_fix'] = \
        rebuilding_regions['hierarchies_to_fix'].apply(
            lambda x: set(x.split(',')) if isinstance(x, str) and len(x) > 0 else {}
        )
    rebuilding_regions = rebuilding_regions.sort_values(['order'])
    return rebuilding_regions


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


def export_streetgraph(G, path_edges, path_nodes, edge_columns=None, node_columns=None):
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

    Returns
    -------
    None
    """

    nodes, edges = oxc.graph_to_gdfs(G)

    if edge_columns:
        edges = edges[list(set(set(edges.columns) & set(edge_columns)).union({'geometry'}))]

    if node_columns:
        nodes = edges[list(set(set(nodes.columns) & set(node_columns)).union({'geometry'}))]

    # Limit every attribute to 10 characters to match the SHP format restrictions
    if path_edges.split()[-1] == 'shp':
        edges.columns = [column[0:10] for column in edges.columns]

    # stringify iterable columns
    _stringify_iterable_columns(edges, {'ln_desc', 'ln_desc_after', 'given_lanes'}, separator=' | ')
    _stringify_iterable_columns(nodes, {'layers'}, method='str')

    # write files
    export_gdf(edges, path_edges)
    export_gdf(nodes, path_nodes)


def export_streetgraph_with_lanes(G, lanes_attribute, path, scaling=1):
    """
    Export a geofile with individual lane geometries. This is helpful for visualization purposes.

    Parameters
    ----------
    G : nx.MultiGraph or nx.MultiDiGraph
        street graph
    lanes_attribute : str
        which attribute should be used as a source of the lane configuration of each edge
    path : str
        where should the file be saved
    scaling : float
        optional scaling of the lane width, for an optimized visualization

    Returns
    -------
    None
    """

    # Create empty list for the lanes
    lanes_list = []

    for id, data in G.edges.items():
        given_total_width = 0

        # Reconstruct total width of given lanes
        for lane in data.get(lanes_attribute, []):
            lane_properties = lanes._lane_properties(lane)
            given_total_width += lane_properties.width

        offset = -given_total_width / 2
        for lane in data.get(lanes_attribute, []):
            lane_properties = lanes._lane_properties(lane)

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
                'type': lane_properties.lanetype,
                'direction': lane_properties.direction,
                'descr': lane,
                'width_m': lane_properties.width * scaling,
                'layer': data.get('layer'),
                'geometry': geom
            })

    lanes_gdf = gpd.GeoDataFrame(
        lanes_list,
        columns=['type', 'direction', 'descr', 'width_m', 'layer', 'geometry'],
        geometry='geometry'
    )

    export_gdf(lanes_gdf, path)


def export_gdf(gdf, path, columns=[]):
    """
    Export a GeoDataFrame into a geofile.
    A wrapper around GeoDataFrame.to_file() but with an optional column filter.

    Parameters
    ----------
    gdf : gpd.GeoDataFrame
    path : str
    columns : list
        which columns should be included ([] -> include all columns)

    Returns
    -------
    None
    """

    if columns == []:
        gdf.to_file(path)
    else:
        gdf[columns].to_file(path)


def convert_crs_of_street_graph(G, to_crs):
    """
    Convert the coordinate reference system of the geometries in a graph.
    The graph must follow the convention of networkx:
        * edges have a 'geometry' attribute with a shapely geometry
        * nodes have 'x' and 'y' attributes holding the coordinates as floats
        * the graph has a 'crs' attribute

    Parameters
    ----------
    G : nx.MultiGraph or nx.MultiDiGraph
        street graph
    to_crs : int
        target crs

    Returns
    -------
    None
    """

    # Initialize the CRS transformer
    from_crs = pyproj.CRS(G.graph['crs'])
    to_crs = pyproj.CRS(to_crs)
    project = pyproj.Transformer.from_crs(from_crs, to_crs, always_xy=True).transform

    # Update the street_graph's metadata
    G.graph["crs"] = to_crs

    # Transform the geometry of all edges
    for edge in G.edges(data=True, keys=True):
        if "geometry" in edge[3]:
            edge[3]["geometry"] = shapely.ops.transform(project, edge[3]["geometry"])

    # Transform the geometry of all nodes
    for id, data in G.nodes.items():
        geom = shapely.geometry.Point(data.get('x'), data.get('y'))
        geom = shapely.ops.transform(project, geom)
        data['x'] = geom.x
        data['y'] = geom.y


def export_osm_xml(G, path, tags, uv_tags=False, tag_all_nodes=False):
    """
    Generates an OSM file from the street graph

    Parameters
    ----------
    G : nx.MultiGraph or nx.MultiDiGraph
        street graph
    path : str
    tags : list
        which OSM tags should be included
    uv_tags : bool
        include special tags for the start and end node of each edge (for debugging)
    tag_all_nodes : bool
        add a special tag to each node so that all nodes appear as points when imported in QGIS (for debugging)

    Returns
    -------
    None
    """

    # initial ID value for new OSM objects, avoid duplicity with graph node ids
    max_node_id = max(list(G.nodes))
    osm_id = itertools.count(max_node_id * 100)

    # make a copy of the original graph and convert it to pseudo mercator which is the official OSM crs
    G = G.copy(G)
    convert_crs_of_street_graph(G, 'epsg:4326')

    # here, we build the XML tree...
    # create the overall structure of the XML
    tree = ET.ElementTree('tree')

    # add OSM metadata
    osm = ET.Element('osm', attrib={
        'version': '0.6',
        'generator': 'osmium/1.14.0'
    })

    # bounds
    lats = [data.get('y') for id, data in G.nodes.items()]
    lons = [data.get('x') for id, data in G.nodes.items()]
    ET.SubElement(osm, 'bounds', attrib={
        'minlat': str(min(lats)),
        'minlon': str(min(lons)),
        'maxlat': str(max(lats)),
        'maxlon': str(max(lons))
    })

    # nodes dictionary, will be filled while creating ways
    nodes = {}

    # ways
    for id, data in G.edges.items():

        way = ET.SubElement(osm, 'way', attrib={
            'id': str(next(osm_id)),  # assign a new unique osm id to each way
            'version': '1',
            'timestamp': '2000-01-01T00:00:00Z'
        })

        # way nodes
        for i in [0, 0.5, 1]:

            # first and last node
            if i in [0, 1]:
                this_id = id[i]
                this_graph_node = G.nodes[this_id]
                nodes[this_id] = shapely.geometry.Point(this_graph_node['x'], this_graph_node['y'])
                # add node along the way, use the existing node id
                ET.SubElement(way, 'nd', attrib={'ref': str(this_id)})

            # intermediary nodes
            elif i == 0.5:
                linestring = data.get('geometry')
                # iterate over all points along the linestring geometry but exclude the first ans last one
                for point in linestring.coords[1:-1]:
                    this_new_osm_id = next(osm_id)
                    nodes[this_new_osm_id] = shapely.geometry.Point(point)
                    # add node along the way, assign a new unique id
                    ET.SubElement(way, 'nd', attrib={'ref': str(this_new_osm_id)})

        # way tags
        for tag in tags:
            # skip if the tag is note defined for this way
            if data.get(tag, None) is None:
                continue
            ET.SubElement(way, 'tag', attrib={
                'k': tag,
                'v': str(data.get(tag, ''))
            })

        if uv_tags:
            ET.SubElement(way, 'tag', attrib={'k': '_u', 'v': str(id[0])})
            ET.SubElement(way, 'tag', attrib={'k': '_v', 'v': str(id[1])})
            ET.SubElement(way, 'tag', attrib={'k': '_key', 'v': str(id[2])})

    # nodes
    for id, point in nodes.items():
        node = ET.SubElement(osm, 'node', attrib={
            'id': str(id),
            'version': '1',
            'timestamp': '2000-01-01T00:00:00Z',
            'lat': str(point.y),
            'lon': str(point.x)
        })

        if tag_all_nodes:
            ET.SubElement(node, 'tag', attrib={'k': '_node', 'v': 'true'})

        # TODO: add other tags

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
                # TODO: to be implemented
                pass


def _stringify_iterable_columns(df, columns, method='separator', separator=','):
    """
    Convert iterables in selected columns into strings

    Parameters
    ----------
    df : pd.DataFrame
    columns : list
        column names
    method : 'str'
        how should the iterable be stringified
            * 'separator' -> join the values with a separator
            * 'str' -> stringify the iterable using str(x)
    separator : 'str'
        only used if method='separator'

    Returns
    -------
    None
    """

    for column in columns:
        if column in df:
            if method == 'separator':
                df[column] = df[column].apply(lambda x: separator.join(x))
            elif method == 'str':
                df[column] = df[column].apply(lambda x: str(x))
