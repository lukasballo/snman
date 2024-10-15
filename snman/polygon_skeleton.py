# Adapted from https://github.com/plouvart/marrow/

__version__ = "0.0.1"

from shapely.geometry import Polygon, Point, LineString, MultiLineString, MultiPoint, LinearRing
from shapely.affinity import translate, rotate
from shapely.ops import nearest_points, unary_union, linemerge, substring
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import time


def disp(poly, *args, ax=None, **kwargs):
    res = gpd.GeoDataFrame(geometry=[poly]).plot(ax=ax, *args, **kwargs)
    if ax is None:
        return res
    else:
        return ax


def cut(line, distance):
    # Cuts a line in two at a distance from its starting point
    if distance <= 0.0 or distance >= line.length:
        return [LineString(line)]
    coords = list(line.coords)
    for i, p in enumerate(coords):
        pd = line.project(Point(p))
        if pd == distance:
            return [
                LineString(coords[:i + 1]),
                LineString(coords[i:])]
        if pd > distance:
            cp = line.interpolate(distance)
            return [
                LineString(coords[:i] + [(cp.x, cp.y)]),
                LineString([(cp.x, cp.y)] + coords[i:])]


def create_ballpoint(p0, p1, p2, v, contour, min_dist, max_dist):
    p = Point(p1)
    corner = LineString(
        [p0, p1, p2]
    )
    max_length = p.distance(
        (back_point := nearest_points(
            p,
            LineString(
                [
                    p1 + v * min_dist,
                    p1 + v * max_dist,
                ]
            ).intersection(contour)
        )[1])
    )
    right_l = max_length
    left_l = min_dist
    good_l = None
    while 1:
        mid_l = (right_l + left_l) / 2

        # ax = disp(contour)
        # ax = disp(LineString([p, back_point]), color="yellow", ax=ax)
        # ax = disp(Point(p1), color="b", ax=ax)
        # ax = disp(Point(p1 + v * right_l), color="r", ax=ax)
        # ax = disp(Point(p1 + v * left_l), color="g", ax=ax)
        # ax = disp(Point(p1 + v * mid_l), color="purple", ax=ax)
        # ax = disp(corner, color="yellow", ax=ax)
        # ax = disp(
        # nearest_points(
        # Point(p1 + v * mid_l),
        # contour,
        # )[1], ax=ax, color="orange"
        # )
        # plt.show()

        if corner.distance(
                nearest_points(
                    Point(p1 + v * mid_l),
                    contour,
                )[1]
        ) < 0.00001:
            good_l = mid_l
            right_l = right_l
            left_l = mid_l
        else:
            right_l = mid_l
            left_l = left_l

        if ((right_l - left_l) < min_dist) and good_l is not None:
            break

    return Point(p1 + v * good_l)


def create_ballpoint_interior(p0, p1, p2, v, interior, contour, min_dist, max_dist):
    p = Point(p1)
    corner = LineString(
        [p0, p1, p2]
    )
    max_length = p.distance(
        (back_point := nearest_points(
            p,
            LineString(
                [
                    p1 + v * min_dist,
                    p1 + v * max_dist,
                ]
            ).intersection(contour)
        )[1])
    )
    right_l = max_length
    left_l = min_dist
    good_l = None
    while 1:
        mid_l = (right_l + left_l) / 2

        # ax = disp(interior)
        # ax = disp(LineString([p, back_point]), color="yellow", ax=ax)
        # ax = disp(Point(p1), color="b", ax=ax)
        # ax = disp(Point(p1 + v * right_l), color="r", ax=ax)
        # ax = disp(Point(p1 + v * left_l), color="g", ax=ax)
        # ax = disp(Point(p1 + v * mid_l), color="purple", ax=ax)
        # ax = disp(corner, color="yellow", ax=ax)
        # ax = disp(
        # nearest_points(
        # Point(p1 + v * mid_l),
        # interior,
        # )[1], ax=ax, color="orange"
        # )
        # plt.show()

        if corner.distance(
                nearest_points(
                    Point(p1 + v * mid_l),
                    contour,
                )[1]
        ) < 0.00001:
            good_l = mid_l
            right_l = right_l
            left_l = mid_l
        else:
            right_l = mid_l
            left_l = left_l

        if ((right_l - left_l) < min_dist) and good_l is not None:
            break

    return Point(p1 + v * good_l)


def compute_corners(contour):
    c1 = contour.coords[:-1]
    c0 = c1[-1:] + c1[:-1]
    c2 = c1[1:] + c1[:1]
    v = np.array(list(zip(c0, c1, c2)))
    v1 = v[:, 1, :] - v[:, 0, :]
    v2 = v[:, 2, :] - v[:, 1, :]
    v1 = v1 / (np.sum(v1 ** 2, axis=1) ** .5).reshape(-1, 1)
    v2 = v2 / (np.sum(v2 ** 2, axis=1) ** .5).reshape(-1, 1)

    normals = np.array([
        s * (v - u) if (s := (np.sign(np.arctan2(np.cross(u, v), sum(u ** 2) ** .5 * sum(v ** 2) ** .5)))) != 0
        else np.array([[0, -1], [1, 0]]).dot(v)
        for u, v in zip(v1, v2)
    ])
    normals = normals / (np.sum(normals ** 2, axis=1) ** .5).reshape(-1, 1)

    return c0, c1, c2, normals


def simple_polygon_skeleton(poly):
    t1 = time.time()
    ax = None
    contour = poly.exterior
    if not contour.is_ccw:
        contour = LinearRing(contour.coords[::-1])

    min_dist = 0.125 / 20  # Should be < min distance between all points / 2
    max_dist = 10  # Should be > max distance between all points

    c0, c1, c2, vv = compute_corners(contour)
    # ax = disp(poly, alpha=0.4)
    # for p,q in zip(c1,c2):
    # ax = disp(LineString([p,q]), ax=ax, color="r")
    # plt.show()

    balls = []
    ax = None
    for i, (p0, p1, p2, v) in enumerate(zip(c0, c1, c2, vv)):
        b = create_ballpoint(p0, p1, p2, v, contour, min_dist, max_dist)
        # ax = disp(b.buffer(b.distance(contour)), ax=ax, alpha=0.5)
        # ax = disp(LineString([p1,b]), ax=ax)
        balls.append(b)
    # ax = disp(contour, color="r", ax=ax)
    # plt.show()

    t = time.time()
    dists = np.ones((len(balls), len(balls))) * np.inf
    inds = list(range(len(balls)))
    for i in inds:
        for j in inds[i + 1:]:
            if not LineString([balls[i], balls[j]]).intersects(contour):
                dists[i, j] = balls[i].distance(balls[j])
    skeleton = MultiLineString([])
    used_points = [inds.pop(0)]
    unused_points = inds

    min_used = None
    min_unused = None
    while len(unused_points):
        min_d = np.inf
        for i1 in used_points:
            for i2 in unused_points:
                if i1 > i2:
                    if (d := dists[i2, i1]) < min_d:
                        min_used = i1
                        min_unused = i2
                        min_d = d
                else:
                    if (d := dists[i1, i2]) < min_d:
                        min_used = i1
                        min_unused = i2
                        min_d = d
        used_points.append(min_unused)
        unused_points.remove(min_unused)
        skeleton = skeleton.union(LineString([balls[min_unused], balls[min_used]]))
    # print("Time ellapsed: ", time.time() - t1)

    # ax = disp(contour, color="r")
    # disp(skeleton, color="g", ax=ax)
    # plt.show()
    return skeleton


def compute_interior_corners(interior):
    c1 = interior.coords[1:][::-1]
    c0 = c1[-1:] + c1[:-1]
    c2 = c1[1:] + c1[:1]
    v = np.array(list(zip(c0, c1, c2)))
    v1 = v[:, 1, :] - v[:, 0, :]
    v2 = v[:, 2, :] - v[:, 1, :]
    v1 = v1 / (np.sum(v1 ** 2, axis=1) ** .5).reshape(-1, 1)
    v2 = v2 / (np.sum(v2 ** 2, axis=1) ** .5).reshape(-1, 1)

    normals = np.array([
        s * (v - u) if (s := (np.sign(np.arctan2(np.cross(u, v), sum(u ** 2) ** .5 * sum(v ** 2) ** .5)))) != 0
        else np.array([[0, -1], [1, 0]]).dot(v)
        for u, v in zip(v1, v2)
    ])
    normals = normals / (np.sum(normals ** 2, axis=1) ** .5).reshape(-1, 1)

    print(c0, c1, c2)
    return c0, c1, c2, normals


def complex_polygon_skeleton(poly):
    t1 = time.time()
    ax = None
    contour = MultiLineString(
        [poly.exterior] + list(poly.interiors)
    )
    # if not contour.is_ccw:
    # contour = LinearRing(contour.coords[::-1])

    min_dist = 0.125 / 20  # Should be < min distance between all points / 2
    max_dist = 10  # Should be > max distance between all points

    balls = []
    ax = None

    for interior in contour.geoms:
        c0, c1, c2, vv = compute_interior_corners(interior)
        # ax = disp(poly, alpha=0.4)
        # for p,q in zip(c1,c2):
        # ax = disp(LineString([p,q]), ax=ax, color="r")
        # plt.show()

        ax = None
        for i, (p0, p1, p2, v) in enumerate(zip(c0, c1, c2, vv)):
            b = create_ballpoint_interior(p0, p1, p2, v, interior, contour, min_dist, max_dist)
            # ax = disp(b.buffer(b.distance(interior)), ax=ax, alpha=0.5)
            # ax = disp(LineString([p1,b]), ax=ax)
            balls.append(b)
        # ax = disp(contour, color="r", ax=ax)
        # plt.show()

    t = time.time()
    dists = np.ones((len(balls), len(balls))) * np.inf
    inds = list(range(len(balls)))
    for i in inds:
        for j in inds[i + 1:]:
            if LineString([balls[i], balls[j]]).within(poly):
                dists[i, j] = balls[i].distance(balls[j])
    skeleton = MultiLineString([])
    used_points = [inds.pop(0)]
    unused_points = inds

    min_used = None
    min_unused = None
    while len(unused_points):
        min_d = np.inf
        for i1 in used_points:
            for i2 in unused_points:
                if i1 > i2:
                    if (d := dists[i2, i1]) < min_d:
                        min_used = i1
                        min_unused = i2
                        min_d = d
                else:
                    if (d := dists[i1, i2]) < min_d:
                        min_used = i1
                        min_unused = i2
                        min_d = d
        used_points.append(min_unused)
        unused_points.remove(min_unused)
        skeleton = skeleton.union(LineString([balls[min_unused], balls[min_used]]))
    print("Time ellapsed: ", time.time() - t1)

    ax = disp(poly, color="r", alpha=0.3)
    ax = disp(skeleton, color="g", ax=ax)

    plt.show()
    return skeleton
