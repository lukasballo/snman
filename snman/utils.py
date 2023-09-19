from typing import Iterable
import json
import shapely as shp


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


def multilinestring_to_linestring(geom):
    if isinstance(geom, shp.geometry.MultiLineString):
        geom = geom.geoms[0]
    return geom


def safe_division(dividend, divisor, if_division_by_zero=float('inf')):
    return dividend / divisor if divisor != 0 else if_division_by_zero
