from typing import Iterable
import json

import shapely
import shapely as shp
import numpy as np
from scipy.spatial import cKDTree
import pandas as pd
import geopandas as gpd


def flatten_list(items):
    """
    A simple function for flattening nested lists, from this post:
    https://stackoverflow.com/questions/952914/how-do-i-make-a-flat-list-out-of-a-list-of-lists/40857703?r=Saves_UserSavesList#40857703

    Parameters
    ----------
    items : list
        nested list

    Returns
    -------
    list
        flattened list
    """
    for x in items:
        if isinstance(x, Iterable) and not isinstance(x, (str, bytes)):
            for sub_x in flatten_list(x):
                yield sub_x
        else:
            yield x


def convert_list_items_to_strings(items):
    return [str(item) for item in items]


def set_last_or_append(given_list, value):
    if not given_list:
        given_list.append(value)
    else:
        given_list[-1] = value


def is_convertible_to_int(s):
    try:
        int(s)
        return True
    except ValueError:
        return False


def safe_int(s, fallback_value=None):
    if is_convertible_to_int(s):
        return int(s)
    else:
        return fallback_value


def _is_json_serializable(obj):
    try:
        json.dumps(obj)
        return True
    except (TypeError, OverflowError):
        return False


def safe_dumps(data):
    if _is_json_serializable(data):
        return json.dumps(data)
    else:
        return 'NULL'


def multilinestring_to_linestring(geom, how='only_first'):
    """
    Converts a shapely.MultiLineString to shapely.LineString,
    using the selected method.

    Parameters
    ----------
    geom : shp.MultiLineString
        the multilinestring geometry
    how : str
        - only_first: only first element will be kept
        - merge: all elements will be merged into one

    Returns
    -------

    """
    if isinstance(geom, shp.geometry.MultiLineString):
        if how == 'only_first':
            return geom.geoms[0]
        elif how == 'merge':
            points = []
            for single_geom in list(geom.geoms):
                points += list(single_geom.coords)
            return shapely.LineString(points)
        else:
            return None
    else:
        return geom



def safe_division(dividend, divisor, if_division_by_zero=float('inf')):
    return dividend / divisor if divisor != 0 else if_division_by_zero


def join_nearest_points(gdA, gdB):
    """
    Joins two GeoDataFrames with points such that each point gets joined with its nearest neighbor

    Copied from here:
    https://gis.stackexchange.com/questions/222315/finding-nearest-point-in-other-geodataframe-using-geopandas

    Parameters
    ----------
    gdA : gpd.GeoDataFrame
    gdB : gpd.GeoDataFrame

    Returns
    -------
    gpd.GeoDataFrame

    """
    nA = np.array(list(gdA.geometry.apply(lambda x: (x.x, x.y))))
    nB = np.array(list(gdB.geometry.apply(lambda x: (x.x, x.y))))
    btree = cKDTree(nB)
    dist, idx = btree.query(nA, k=2)
    dist = np.array(list(map(lambda x: x[1], dist)))
    idx = np.array(list(map(lambda x: x[1], idx)))
    gdB_nearest = gdB.iloc[idx].drop(columns="geometry").reset_index(drop=True)
    gdf = pd.concat(
        [
            gdA.reset_index(drop=True),
            pd.Series(dist, name='distance_to_next_space')
        ],
        axis=1)

    return gdf


def object_from_string(string):
    return json.loads(string) if string != 'nan' else []


def get_nth_element_of_list(my_list, n):
    """
    A safe way to get nth element of a list. Returns None if that element does not exist.

    Parameters
    ----------
    my_list: list
    n: int

    Returns
    -------
    Any

    """
    if n < len(my_list):
        return my_list[n]
    else:
        return None
