from . import lanes
import networkx as nx
import shapely as shp
import numpy as np
import math


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


def _offset_distance(linestrings):

    # Start and end point of axis based on the first geometry
    u = shp.ops.Point(linestrings[0].coords[0])
    v = shp.ops.Point(linestrings[0].coords[-1])
    axis = shp.ops.LineString([u,v])

    # check if the axis vector is in quadrant 2 or 3 (in these cases, we need to add 180deg to the value of arctan)
    quadrants_23 = v.x - u.x < 0

    # Points at which we will measure the distance between axis and each geometry
    measurement_points = [axis.interpolate(step, normalized=True) for step in np.arange(0,1,0.1)]
    axis_angle = math.degrees(np.arctan((v.y - u.y) / (v.x - u.x)))

    offsets = []

    for linestring in linestrings:
        nearest_points = [shp.ops.nearest_points(measurement_point, linestring)[1] for measurement_point in measurement_points]
        vectors = [
            (points_pair[1].x - points_pair[0].x, points_pair[1].y - points_pair[0].y)
            for points_pair
            in list(zip(measurement_points, nearest_points))
        ]
        diff_vector = np.mean(vectors, axis=0)
        # rotate the difference vector based on the axis
        rotated_diff_vector = shp.affinity.rotate(
            shp.ops.LineString([(0,0), diff_vector]),
            -(axis_angle + (quadrants_23 * 180)),
            origin=[0,0]
        )
        # take the y coordinate as offset
        offset = rotated_diff_vector.coords[1][1]
        offsets.append(offset)

    return offsets

def ensure_multipolygon(geometry):
    if isinstance(geometry, shp.geometry.MultiPolygon):
        return geometry
    else:
        return shp.geometry.MultiPolygon([geometry])
