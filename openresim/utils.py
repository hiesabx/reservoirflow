import warnings
from functools import lru_cache
from numpy import ndarray


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
    for a in array:
        if tuple(a) == x:
            return True
    return False


def fshape_warn(class_unify, func_unify):
    if class_unify != func_unify:
        warnings.warn("Inconsistent argument was used.")
        print(
            "[WARNING]: "
            f"Class was initiated with unify option set to {class_unify}. "
            f"Setting unify argument to {func_unify} may cause some errors."
        )
