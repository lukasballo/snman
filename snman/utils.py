from typing import Iterable


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
