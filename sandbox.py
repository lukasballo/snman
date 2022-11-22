import snman
import shapely as shp
import numpy as np
import geopandas as gpd
import math


def _offset_distance(linestrings):

    # Start and end point of axis based on the first geometry
    u = shp.ops.Point(linestrings[0].coords[0])
    v = shp.ops.Point(linestrings[0].coords[-1])
    axis = shp.ops.LineString([u,v])

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
        rotated_diff_vector = shp.affinity.rotate(shp.ops.LineString([(0,0), diff_vector]), -axis_angle, origin=[0,0])
        # take the y coordinate as offset
        offset = rotated_diff_vector.coords[1][1]
        offsets.append(offset)

    return offsets

offsets = _offset_distance([
    shp.ops.LineString([(0,0), (5,6), (10,10)]),
    shp.ops.LineString([(0,0), (5,5), (10,10)]),
    shp.ops.LineString([(0,0), (5,4), (10,10)])
])