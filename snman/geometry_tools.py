from . import lanes
import networkx as nx
import shapely as shp


def remove_multipart_geometries(street_graph):

    for id, data in street_graph.edges.items():
        geom = data.get('geometry', None)
        if not isinstance(geom, shp.geometry.linestring.LineString):
            print(type(geom))
            start_node = street_graph.nodes.items()[id[0]]
            end_node = street_graph.nodes.items()[id[1]]
            simple_line = shp.geometry.LineString(
                shp.geometry.Point(
                    start_node.get('x',0),
                    start_node.get('y',0)
                ),
                shp.geometry.Point(
                    end_node.get('x', 0),
                    end_node.get('y', 0)
                )
            )
            data['geometry'] = simple_line
