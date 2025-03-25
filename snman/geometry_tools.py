import networkx as nx
import shapely as shp
import numpy as np
import math

from . import utils


def remove_multipart_geometries(G):
    """
    Replace all multiline geometries in a street graph by simple start-end lines based in the node coordinates.
    Normally, this should not be necessary but there are rare cases when some edges have multipart geometries,
    leading to errors in geometry operations.

    Parameters
    ----------
    G : nx.MultiGraph
        street graph

    Returns
    -------
    None
    """

    for id, data in G.edges.items():
        geom = data.get('geometry', None)
        if not isinstance(geom, shp.geometry.linestring.LineString):
            start_node = G.nodes.items()[id[0]]
            end_node = G.nodes.items()[id[1]]
            simple_line = shp.geometry.LineString(
                shp.geometry.Point(
                    start_node.get('x', 0),
                    start_node.get('y', 0)
                ),
                shp.geometry.Point(
                    end_node.get('x', 0),
                    end_node.get('y', 0)
                )
            )
            data['geometry'] = simple_line


def _offset_distance(linestrings):
    """
    Takes a set of linestrings and calculates the average offset for each of them.
    This is useful for detecting the order of approximately parallel edges geometries that should be merged into
    a single street.

    Parameters
    ----------
    linestrings : list
        a list of approximately parallel linestrings

    Returns
    -------
    offsets : list
        a list of the calculated offsets
    """

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
    """
    Converts polygons into multipolygons if necessary

    Parameters
    ----------
    geometry : shp.Polygon, shp.MultiPolygon

    Returns
    -------
    geometry : shp.MultiPolygon
    """

    if isinstance(geometry, shp.geometry.MultiPolygon):
        return geometry
    else:
        return shp.geometry.MultiPolygon([geometry])


def ensure_multilinestring(geometry):
    """
    Converts linestrings into multilinestrings if necessary

    Parameters
    ----------
    geometry : shp.LineString, shp.MultiLineString

    Returns
    -------
    geometry : shp.MultiLineString
    """

    if isinstance(geometry, shp.geometry.MultiLineString):
        return geometry
    else:
        return shp.geometry.MultiLineString([geometry])


def reverse_linestring(geometry):
    return shp.ops.substring(geometry, 1, 0, normalized=True)


def get_polygon_axis(polygon):

    # Get the exterior coordinates of the rectangle
    rectangle = polygon.minimum_rotated_rectangle
    coords = list(rectangle.exterior.coords)

    # Identify the edges (as LineStrings)
    edges = [shp.LineString([coords[i], coords[i + 1]]) for i in range(len(coords) - 1)]

    # Sort the edges based on their length
    sorted_edges = sorted(edges, key=lambda x: x.length)

    # The first two in the sorted list will be the shorter edges
    shorter_edges = sorted_edges[:4]

    axis = shp.LineString([
        shp.line_interpolate_point(shorter_edges[0], 0.5, normalized=True),
        shp.line_interpolate_point(shorter_edges[1], 0.5, normalized=True)
    ])

    return axis


def random_points_in_polygon(polygon, number):
    """
    Adapted from https://www.matecdev.com/posts/random-points-in-polygon.html

    Parameters
    ----------
    polygon: shp.Polygon
    number: int
    """

    points = []
    minx, miny, maxx, maxy = polygon.bounds
    while len(points) < number:
        pnt = shp.Point(np.random.uniform(minx, maxx), np.random.uniform(miny, maxy))
        if polygon.contains(pnt):
            points.append(pnt)
    return shp.MultiPoint(points)

def linestring_angle(geometry):
    dx = geometry.coords[-1][0] - geometry.coords[0][0]
    dy = geometry.coords[-1][1] - geometry.coords[0][1]
    angle = np.degrees(np.arctan(
        utils.safe_division(dy,dx)
    ))
    if dx < 0:
        angle += 180
    if angle < 0:
        angle += 360
    return angle
