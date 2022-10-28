import networkx as nx
from . import osmnx as ox
from . import config, graph_tools, lanes
import geopandas as gpd
import pyproj
import shapely.ops
import shapely
import momepy
import xml.etree.ElementTree as ET

def export_streetgraph(street_graph, file_name):
    """
    Exports street graph as a shape file

    Params
    ------
    street_graph : nx.MultiGraph
    """
    nodes, edges = ox.graph_to_gdfs(street_graph)

    # Limit every attribute to 10 characters to match the SHP format restrictions
    if file_name.split()[-1] == 'shp':
        edges.columns = [column[0:10] for column in edges.columns]

    # Convert list attributes to strings
    if 'ln_desc' in edges:
        edges['ln_desc'] = edges['ln_desc'].apply(lambda ln_desc: ' | '.join(ln_desc))
    if 'given_lanes' in edges:
        edges['given_lanes'] = edges['given_lanes'].apply(lambda ln_desc: ' | '.join(ln_desc))

    export_gdf(edges, file_name)

def export_streetgraph_with_lanes(street_graph, lanes_attribute, file_name):
    """
    Exports street graph as a shape file with a polyline for each lane

    Params
    ------
    street_graph : nx.MultiGraph
    """

    # Create empty list for the lanes
    lanes_list = []

    for id, data in street_graph.edges.items():
        given_total_width = 0

        # Reconstruct total width of given lanes
        for lane in data.get(lanes_attribute, []):
            lane_properties = lanes._get_lane_properties(lane)
            given_total_width += lane_properties['width']

        offset = -given_total_width / 2
        for lane in data.get(lanes_attribute, []):
            lane_properties = lanes._get_lane_properties(lane)

            centerline_offset = offset + lane_properties['width']/2
            offset += lane_properties['width']
            geom = data.get('geometry')

            if geom and round(centerline_offset,1) != 0:
                geom = geom.parallel_offset(centerline_offset, 'right')
                # the above function reverses direction when the offset is positive, this steps reverses it back
                if centerline_offset > 0:
                    #geom.coords = list(geom.coords)[::-1]
                    shapely.ops.substring(geom, 1, 0, normalized=True)
                    pass

            lanes_list.append({
                'type': lane_properties['type'],
                'direction': lane_properties['direction'],
                'descr': lane,
                'width_m': lane_properties['width'],
                'geometry': geom
            })

    lanes_gdf = gpd.GeoDataFrame(
        lanes_list,
        columns=['type', 'direction', 'descr', 'width_m', 'geometry'],
        geometry='geometry'
    )

    export_gdf(lanes_gdf, file_name)


def export_gdf(gdf, file_path):
    gdf.to_file(file_path)


def import_shp_to_gdf(file_path):
    edges = gpd.read_file(file_path)
    edges = edges.to_crs(2056)
    return edges


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

def export_matsim_xml(street_graph, file_name):

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
