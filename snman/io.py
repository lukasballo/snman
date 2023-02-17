from . import osmnx_customized as oxc
from . import geometry_tools, lanes
import geopandas as gpd
import pandas as pd
import pyproj
import shapely.ops
import shapely.geometry
import shapely
import xml.etree.ElementTree as ET
import numpy as np
import copy
import itertools
import networkx as nx

def export_streetgraph(street_graph, file_name_edges, file_name_nodes, edge_columns=None, node_columns=None):
    """
    Exports street graph as a shape file

    Params
    ------
    street_graph : nx.MultiGraph
    """
    nodes, edges = oxc.graph_to_gdfs(street_graph)

    if edge_columns:
        edges = edges[list(set(set(edges.columns) & set(edge_columns)).union({'geometry'}))]

    if node_columns:
        nodes = edges[list(set(set(nodes.columns) & set(node_columns)).union({'geometry'}))]

    # Limit every attribute to 10 characters to match the SHP format restrictions
    if file_name_edges.split()[-1] == 'shp':
        edges.columns = [column[0:10] for column in edges.columns]

    #TODO: Detect lists and convert them to strings automatically

    # Convert list attributes to strings
    if 'ln_desc' in edges:
        edges['ln_desc'] = edges['ln_desc'].apply(lambda ln_desc: ' | '.join(ln_desc))
    if 'given_lanes' in edges:
        edges['given_lanes'] = edges['given_lanes'].apply(lambda ln_desc: ' | '.join(ln_desc))
    if 'ln_desc_after' in edges:
            edges['ln_desc_after'] = edges['ln_desc_after'].apply(lambda ln_desc: ' | '.join(ln_desc))
    if 'layers' in nodes:
        nodes['layers'] = nodes['layers'].apply(lambda layers: str(layers))

    export_gdf(edges, file_name_edges)
    export_gdf(nodes, file_name_nodes)

def export_streetgraph_with_lanes(street_graph, lanes_attribute, file_name, scaling=1):
    """
    Exports street graph as a shape file with a polyline for each lane

    Params
    ------
    street_graph : nx.MultiGraph
    scaling : a scaling factor for the width, for visualization purposes
    """

    # Create empty list for the lanes
    lanes_list = []

    for id, data in street_graph.edges.items():
        given_total_width = 0

        # Reconstruct total width of given lanes
        for lane in data.get(lanes_attribute, []):
            lane_properties = lanes._lane_properties(lane)
            given_total_width += lane_properties.width

        offset = -given_total_width / 2
        for lane in data.get(lanes_attribute, []):
            lane_properties = lanes._lane_properties(lane)

            centerline_offset = offset + lane_properties.width/2
            offset += lane_properties.width
            geom = data.get('geometry')

            if geom and round(centerline_offset,1) != 0:
                geom = geom.parallel_offset(centerline_offset * scaling, 'right')
                # the above function reverses direction when the offset is positive, this steps reverses it back
                if centerline_offset > 0 and not geom.is_empty:
                    #geom.coords = list(geom.coords)[::-1]
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

    export_gdf(lanes_gdf, file_name)


def export_gdf(gdf, file_path, columns=[]):
    if columns == []:
        gdf.to_file(file_path)
    else:
        gdf[columns].to_file(file_path)



def import_shp_to_gdf(file_path, crs=2056, index=None):
    gdf = gpd.read_file(file_path).to_crs(crs)
    if index:
        gdf=gdf.set_index(index)
    return gdf


def convert_crs_of_street_graph(street_graph, to_crs):

    # Initialize the CRS transformer
    from_crs = pyproj.CRS(street_graph.graph['crs'])
    to_crs = pyproj.CRS(to_crs)
    project = pyproj.Transformer.from_crs(from_crs, to_crs, always_xy=True).transform

    # Update the street_graph's metadata
    street_graph.graph["crs"] = to_crs

    # Transform the geometry of all edges
    for edge in street_graph.edges(data=True, keys=True):
        if "geometry" in edge[3]:
            edge[3]["geometry"] = shapely.ops.transform(project, edge[3]["geometry"])

    # Transform the geometry of all nodes
    for id, data in street_graph.nodes.items():
        geom = shapely.geometry.Point(data.get('x'), data.get('y'))
        geom = shapely.ops.transform(project, geom)
        data['x'] = geom.x
        data['y'] = geom.y


def export_osm_xml(G:nx.MultiGraph, file_name:str, tags:set, uv_tags:bool=False, tag_all_nodes:bool=False) -> None:
    """
    Generates an OSM file from the street graph
        :param G: street graph
        :param file_name: where to save
        :param tags: which tags should be included
        :param uv_tags: include tags indicating the origin (_u) and destination (_v) nodes, as well as the key (_key),
            for debugging purposes
        :tag_all_nodes: give every node a tag '_node' so that it is visible in qgis import, for debugging purposes
    """

    # initial ID value for new OSM objects, avoid duplicity with graph node ids
    max_node_id = max(list(G.nodes))
    osm_id = itertools.count(max_node_id*100)

    # make a copy of the original graph and concert it to pseudo mercator which is the official OSM crs
    G = copy.copy(G)
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
            'id': str(next(osm_id)),    # assign a new unique osm id to each way
            'version': '1',
            'timestamp': '2000-01-01T00:00:00Z'
        })

        # way nodes
        for i in [0, 0.5, 1]:

            # first and last node
            if i in [0,1]:
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
            ET.SubElement(way, 'tag', attrib={'k': '_u',   'v': str(id[0])})
            ET.SubElement(way, 'tag', attrib={'k': '_v',   'v': str(id[1])})
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

        #TODO: add other tags

    # Save into xml file
    tree._setroot(osm)
    ET.indent(tree)
    tree.write(file_name, encoding='UTF-8', xml_declaration=True)




def export_matsim_xml(street_graph, file_name):
    # TODO: Currently just a preview, not functional

    # Create the overall structure
    tree = ET.ElementTree('tree')
    network = ET.Element('network')
    attributes = ET.SubElement(network, 'attributes')
    ET.SubElement(attributes, 'attribute', attrib={
        'name': 'coordinateReferenceSystem',
        'class': 'java.lang.string'
    }).text = 'Atlantis'
    nodes = ET.SubElement(network, 'nodes')
    links = ET.SubElement(network, 'links', {
        'capperiod': '01:00:00',
        'effectivecellsize': '7.5',
        'effectivelanewidth': '3.75'
    })

    # Add Nodes
    ET.SubElement(nodes, 'node', {'id': '1', 'x': '1000000', 'y': '2000000'})
    ET.SubElement(nodes, 'node', {'id': '2', 'x': '1000500', 'y': '2000500'})

    # Add Links
    link = ET. SubElement(links, 'link', {
        'id': '10',
        'from': '1',
        'to': '2',
        'length': '123',
        'freespeed': '6.9',
        'capacity': '600',
        'permlanes': '1.0',
        'oneway': '1',
        'modes': 'car_passenger,truck,car'
    })

    link_attr = ET.SubElement(link, 'attributes')
    ET.SubElement(link_attr, 'attribute').text = 'unclassified'

    # Save into xml file
    tree._setroot(network)
    ET.indent(tree)
    tree.write(file_name, encoding='UTF-8', xml_declaration=True)
    #TODO: Add DOCTYPE to the file: <!DOCTYPE network SYSTEM "http://www.matsim.org/files/dtd/network_v2.dtd">


def load_perimeters(path):
    perimeters = import_shp_to_gdf(path, crs=4326, index='id')
    return perimeters


def load_regions(path, default_tolerance=None, street_graph=None):
    regions = import_shp_to_gdf(path)

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


def load_intersections(path_polygons):
    polygons = import_shp_to_gdf(path_polygons)
    # Duplicate the point geometries so that they get taken along after the spatial join
    intersections = polygons
    return intersections


def _get_nodes_within_polygon(street_graph, polygon):
    nodes_gdf = oxc.graph_to_gdfs(street_graph, edges=False)
    nodes_gdf = nodes_gdf[nodes_gdf.within(polygon)]
    return set(nodes_gdf.index.values)