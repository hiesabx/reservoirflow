import warnings
from functools import lru_cache
from xmlrpc.client import Boolean
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


def get_boundary_str(boundary):
    return "with boundary" if boundary else "without boundary"


def get_fshape_str(fshape):
    return "with fshape" if fshape else "without fshape"


def get_verbose_str(boundary, fshape):
    return get_boundary_str(boundary), get_fshape_str(fshape)


def isin(x: tuple, array: ndarray):
    """Check if x in array.

    x must be tuple or array of len 3.
    
    Parameters
    ----------
    x : tuple
        tuple or array of len 3.
    array : ndarray
        array of tuples or arrays of len 3.

    Returns
    -------
    Boolean
        True is x in array, otherwise False
    """
    if not isinstance(x, tuple):
        x = tuple(x)
    for a in array:
        if tuple(a) == x:
            return True
    return False


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
