"""
backends
========

This module is used to offer multiple backends to allow using and 
benchmarking different computing tools and frameworks such as ``NumPy``,
``JAX``, ``PyTorch``, ``TensorFlow``, and others.

.. attention::
    This module is not ready.

Information:
    - design pattern: inheritance, abstraction
    - base class: `Backend </api/reservoirflow.backends.Backend.html>`_
    - base class type: ABS (abstract)
"""

__all__ = [
    "Backend",
    "NumPy",
    "PyTorch",
    "TensorFlow",
    "JAX",
]

from .backend import Backend
from .jax import JAX
from .numpy import NumPy
from .pytorch import PyTorch
from .tensorflow import TensorFlow


def set_backend(name: str):
    """Set backend.

    .. attention::
        This function is not ready.

    Parameters
    ----------
    name : str
        backend name.
    """
    pass
