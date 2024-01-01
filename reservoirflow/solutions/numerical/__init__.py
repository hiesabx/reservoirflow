"""
numerical
=========

This module provides solutions using the scientific computing approach 
based on discretized :term:`PDE` where a system of linear equations is
built. Additionally, this module includes ``solvers`` which is used to 
solve the system  of linear equations.

List of numerical solutions: 
    * ``FDM``: Finite-Difference-Method
    * ``FVM``: Finite-Volume-Method (not available)
    * ``FEM``: Finite-Element-Method (not available)

Information:
    - design pattern: inheritance, abstraction
    - base class: `Solution </api/reservoirflow.solutions.Solution.html>`_
    - base class type: ABS (abstract)
"""

__all__ = [
    "FDM",
    "FVM",
    "FEM",
    "solvers",
]


from . import solvers
from .fdm import FDM
from .fem import FEM
from .fvm import FVM
