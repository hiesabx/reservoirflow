#%%
import warnings
from functools import lru_cache
from numpy import ndarray
import numpy as np


def _lru_cache(maxsize=None):
    def wrapper_cache(func):
        _func = lru_cache(maxsize=maxsize)(func)
        return _func
        # @wraps(func)
        # def wrapped_func(*args, **kwargs):
        #     return func(*args, **kwargs)
        # return wrapped_func

    return wrapper_cache


#%%
def get_boundary_str(boundary):
    return "with boundary" if boundary else "without boundary"


def get_points_str(points):
    return " (as points)" if points else ""


def get_fshape_str(fshape):
    return "with fshape" if fshape else "without fshape"


def get_verbose_str(boundary, fshape):
    return get_boundary_str(boundary), get_fshape_str(fshape)


def isin(x: tuple, in_data):
    """Check if x in or not.

    This function checks if x is in data. If x itself is a python
    iterable or array (e.g. coords), x will be converted to a tuple to
    allow for the check.

    Parameters
    ----------
    x : int, tuple, list, tuple, ndarray
        int (e.g. id) or list-like (e.g. coords) of len 3.
    data : list, set, tuple, ndarray
        list-like of int or tuples of len 3.

    Returns
    -------
    Boolean
        True is x in data, otherwise False

    ToDo
    ----
    test if data contains points data and x is i,j,k tuple.
    - must use sets to check instead of fore loops.
    """
    if not isinstance(x, tuple) and not isinstance(x, (np.integer, int)):
        x = tuple(x)
    if isinstance(in_data, np.ndarray):
        if len(in_data.shape) > 1:
            for a in in_data:
                if tuple(a) == x:
                    return True
            return False
    return x in in_data


def intersection(array_x: ndarray, array_y: ndarray, fmt="array"):
    """Find common tuples between two arrays.

    arrays must be flatten

    Parameters
    ----------
    array_x : ndarray
        array of tuples or arrays of len 3.
    array_y : ndarray
        array of tuples or arrays of len 3.
    fmt : str, optional
        output format as str in ['array', 'list', 'tuple'].

    Returns
    -------
    ndarray, list
        common tuples between two arrays.
    """
    argmin = np.argmin([np.max(array_x.shape), np.max(array_y.shape)])
    if argmin == 0:
        xy = [tuple(a) for a in array_x if isin(a, array_y)]
    else:
        xy = [tuple(a) for a in array_y if isin(a, array_x)]

    if fmt in ("tuple", "list"):
        return xy
    elif fmt == "array":
        return np.array(xy)
    else:
        raise ValueError("fmt argument is unknown.")


def fshape_warn(class_unify, func_unify):
    if class_unify != func_unify:
        warnings.warn("Inconsistent argument was used.")
        print(
            "[WARNING]: "
            f"Class was initiated with unify option set to {class_unify}. "
            f"Setting unify argument to {func_unify} may cause some errors."
        )


def issametype(in_data, fmt):
    if isinstance(in_data, np.ndarray) and fmt == "array":
        return in_data
    elif isinstance(in_data, [tuple, list]) and fmt in ["tuple", "list"]:
        return in_data
    elif isinstance(in_data, set) and fmt == "set":
        return in_data
    else:
        return False


def ispoints(in_data):
    """Check if data contains points or scaler.

    Parameters
    ----------
    in_data : _type_
        data must be flatten.

    Returns
    -------
    Boolean
        True if data is tuple of 3 (i.e. points). Otherwise, False.
    """
    in_data = np.array(in_data)
    shape_check = len(in_data.shape) > 1 and in_data.shape[-1] == 3
    if isinstance(in_data, np.ndarray) and shape_check:
        return True
    return False


def reformat(in_data, fmt="tuple"):
    """Reformat input data.

    Parameters
    ----------
    in_data : _type_
        input data to be reformated. If data is ndarray, it must be
        flatten before reformated.
    points : bool, optional, by default False
        _description_
    fmt : str, optional, by default "tuple"
        output format as str from ['array', 'list', 'tuple', 'set']. For
        a better performance, use 'set' to check if an item is in a list
        or not. Use tuples to iterate through items. When option 'array'
        is used, utils.isin() must be used to check if a tuple of 3 is
        in the array.

    Returns
    -------
    list, tuple, set, array, dict
        iterable format of data based on fmt argument.

    Raises
    ------
    ValueError
        fmt is unknown.
    """
    if isinstance(in_data, dict):
        if fmt == "dict":
            return in_data
        in_data = sum(in_data.values(), [])

    points = ispoints(in_data)

    if fmt == "list":
        if points:
            return [tuple(a) for a in in_data]
        return list(in_data)
    elif fmt == "tuple":
        if points:
            return tuple(tuple(a) for a in in_data)
        return tuple(in_data)
    elif fmt == "set":
        if points:
            return set(tuple(a) for a in in_data)
        return set(in_data)
    elif fmt == "array":
        return np.array(in_data)
    else:
        raise ValueError("fmt is unknown or can't reformat to dict.")


def shape_error(in_shape, fshape):
    msg = (
        "boundaries are not included or "
        + "points argument must be correctly assigned.\n"
        + "     - shape info >> "
        + f"in_data shape: {in_shape} - "
        + f"required shape: {fshape}"
    )
    raise ValueError(msg)


def remove_diag(arr):
    # reference: https://stackoverflow.com/a/63531391/11549398
    return arr[~np.eye(len(arr), dtype=bool)].reshape(len(arr), -1)
