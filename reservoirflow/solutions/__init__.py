"""
solutions
=========

This module provides two important classes which are ``Compiler`` and
``Solution``. The ``Compiler`` is used by a model from ``models``
module (e.g. using a ``model.compile()`` method) to build a solution 
object for the corresponding model. The solution object has a class 
inherited from ``Solution`` which is an abstract class used to unify 
the interface of all solution classes available under this module.

There are three submodules available under this module, each have its
own classes based on the solution type. These submodules are:
``analytical``, ``numerical``, and ``neurical``.

.. warning::
    Only classes ``Compiler`` and ``Solution`` are loaded by default.
    The submodules under ``solutions`` are not loaded by default in 
    ``reservoirflow``. As a result, using ``rf.solutions.numerical.FDM``
    after ``import reservoirflow as rf`` will fail.
    
.. tip:: 
    Specific submodules such as ``numerical`` will only be available 
    once at least a single solution from that module (e.g. ``FDM``)
    was used to compile a model.

For more information, check the corresponding documentation of each
submodule and class.
"""

__all__ = [
    "Compiler",
    "Solution",
    "analytical",
    "neurical",
    "numerical",
]

from .compiler import Compiler
from .solution import Solution
