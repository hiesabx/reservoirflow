"""
solvers
=======

This module contains direct and iterative solvers to solve the system
of linear equations provided by numerical methods.
"""

import warnings

import scipy.sparse.linalg as ssl

from reservoirflow.utils.helpers import _lru_cache


def get_dsolver(name):
    pass


@_lru_cache(maxsize=1)
def get_isolver(name):
    """Returns an iterative solver (isolver).

    Iterative solvers for linear systems in sparse matrices using SciPy.
    Available iterative solvers are:

    .. highlight:: python
    .. code-block:: python

        isolvers = [
                "bicg",
                "bicgstab",
                "cg",
                "cgs",
                "gmres",
                "lgmres",
                "minres",
                "qmr",
                "gcrotmk",
                "tfqmr",
            ]

    Parameters
    ----------
    name : str, optional
        name of the iterative solver. If None, direct solver is used.
        Only relevant when argument sparse=True.
        Option "cgs" is recommended to increase performance while option
        "minres" is not recommended due to high MB error.
        For more information, check References section below.

    Returns
    -------
    isolver
        iterative solver for sparse matrices

    Raises
    ------
    ValueError
        isolver name is unknown.

    References
    ----------
    - SciPy: `Solving linear problems <https://docs.scipy.org/doc/scipy/reference/sparse.linalg.html#solving-linear-problems>`_.
    - SciPy: `Iterative Solvers <https://scipy-lectures.org/advanced/scipy_sparse/solvers.html#iterative-solvers>`_.
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
