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
    """
    Convert all items in a list to strings.

    Parameters
    ----------
    items : list
        List of items to convert

    Returns
    -------
    list
        List of string representations of items
    """
    return [str(item) for item in items]


def set_last_or_append(given_list, value):
    """
    Set the last element of a list to a value, or append if list is empty.

    Parameters
    ----------
    given_list : list
        List to modify
    value : Any
        Value to set as last element or append
    """
    if not given_list:
        given_list.append(value)
    else:
        given_list[-1] = value


def is_convertible_to_int(s):
    """
    Check if a value can be converted to an integer.

    Parameters
    ----------
    s : Any
        Value to check

    Returns
    -------
    bool
        True if value can be converted to int, False otherwise
    """
    try:
        int(s)
        return True
    except ValueError:
        return False


def safe_int(s, fallback_value=None):
    """
    Safely convert a value to integer, returning fallback if conversion fails.

    Parameters
    ----------
    s : Any
        Value to convert to int
    fallback_value : Any, optional
        Value to return if conversion fails (default: None)

    Returns
    -------
    int or Any
        Converted integer value or fallback_value if conversion fails
    """
    if is_convertible_to_int(s):
        return int(s)
    else:
        return fallback_value


def safe_float(s, fallback_value=None):
    """
    Safely convert a value to float, returning fallback if conversion fails.

    Parameters
    ----------
    s : Any
        Value to convert to float
    fallback_value : Any, optional
        Value to return if conversion fails (default: None)

    Returns
    -------
    float or Any
        Converted float value or fallback_value if conversion fails
    """
    try:
        return float(s)
    except (ValueError, TypeError):
        return fallback_value


def _is_json_serializable(obj):
    """
    Check if an object can be serialized to JSON.

    Parameters
    ----------
    obj : Any
        Object to check

    Returns
    -------
    bool
        True if object is JSON serializable, False otherwise
    """
    try:
        json.dumps(obj)
        return True
    except (TypeError, OverflowError):
        return False


def safe_dumps(data):
    """
    Safely serialize data to JSON string, returning 'NULL' if serialization fails.

    Parameters
    ----------
    data : Any
        Data to serialize

    Returns
    -------
    str
        JSON string representation or 'NULL' if serialization fails
    """
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
    shp.LineString or None
        Converted LineString geometry or None if conversion fails
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
    """
    Safely divide two numbers, returning a default value if divisor is zero.

    Parameters
    ----------
    dividend : float
        Number to be divided
    divisor : float
        Number to divide by
    if_division_by_zero : float, optional
        Value to return if divisor is zero (default: inf)

    Returns
    -------
    float
        Result of division or if_division_by_zero if divisor is zero
    """
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
    """
    Parse a JSON string to a Python object, returning empty list for 'nan'.

    Parameters
    ----------
    string : str
        JSON string to parse

    Returns
    -------
    Any
        Parsed object or empty list if string is 'nan'
    """
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


def isnan(val):
    """
    Check if a value is NaN (Not a Number).

    Parameters
    ----------
    val : Any
        Value to check

    Returns
    -------
    bool
        True if value is a float and NaN, False otherwise
    """
    return type(val) == float and np.isnan(val)


def merge_dicts(dicts, ignored_values=(None, '', [], '[]', np.nan, 'nan')):
    """
    Merges a list of dictionaries using the *update()* method, while ignoring a list of values

    Parameters
    ----------
    dicts: list
        the list of dictionaries
    ignored_values: list
        which values should be ignored

    Returns
    -------
    dict
    """

    merged_dict = {}
    for d in dicts:
        d = {
            k: v for k, v in d.items()
            if not any(v is val or isnan(v) for val in ignored_values)
        }
        merged_dict.update(d)
    return merged_dict


def is_in_list(obj, lst):
    """
    Check if an object is in a list using identity (is) rather than equality (==).

    Parameters
    ----------
    obj : Any
        Object to search for
    lst : list
        List to search in

    Returns
    -------
    bool
        True if object is found in list (by identity), False otherwise
    """
    return any(obj is item for item in lst)


class IncrementingVariable:
    """
    An iterator that returns incrementing integer values.

    Each call to next() returns the current value and then increments it.
    """
    def __init__(self, start=0):
        """
        Initialize the incrementing variable.

        Parameters
        ----------
        start : int, optional
            Starting value (default: 0)
        """
        self.value = start

    def __iter__(self):
        return self

    def __next__(self):
        """
        Return the current value and increment.

        Returns
        -------
        int
            Current value before incrementing
        """
        self.value += 1
        return self.value - 1  # Return the previous value before incrementing
