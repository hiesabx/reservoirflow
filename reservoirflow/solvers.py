import warnings
from reservoirflow.utils import _lru_cache
import scipy.linalg as sl
import scipy.sparse as ss
import scipy.sparse.linalg as ssl


def get_dsolver(name):
    pass


@_lru_cache(maxsize=1)
def get_isolver(name):
    """Returns an iterative solver (isolver).

    Iterative solvers for linear systems in sparse matrices.

    Parameters
    ----------
    name : str, optional, by default "cgs"
        name of th eiterative solver. Available solvers are
        ["bicg", "bicgstab", "cg", "cgs", "gmres", "lgmres",
        "minres", "qmr", "gcrotmk", "tfqmr"].
        If None, direct solver is used. Only relevant when argument
        sparse=True. Option "cgs" is recommended to increase
        performance while option "minres" is not recommended due to
        high MB error. For more information check [1][2].

    Returns
    -------
    solver
        iterative solver for sparse matrices

    Raises
    ------
    ValueError
        isolver name is unknown.

    References
    ----------
    [1] https://docs.scipy.org/doc/scipy/reference/sparse.linalg.html#solving-linear-problems
    [2] https://scipy-lectures.org/advanced/scipy_sparse/solvers.html#iterative-solvers
    """

    if name == "bicg":
        solver = ssl.bicg
    elif name == "bicgstab":
        solver = ssl.bicgstab
    elif name == "cg":
        solver = ssl.cg
    elif name == "cgs":
        solver = ssl.cgs
    elif name == "gmres":
        solver = ssl.gmres
    elif name == "lgmres":
        solver = ssl.lgmres
    elif name == "minres":
        solver = ssl.minres
        warnings.warn("option isolver='minres' is not recommended.")
    elif name == "qmr":
        solver = ssl.qmr
    elif name == "gcrotmk":
        solver = ssl.gcrotmk
    elif name == "tfqmr":
        solver = ssl.tfqmr
    else:
        raise ValueError("isolver name is unknown.")

    return solver
