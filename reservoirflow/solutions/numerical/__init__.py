"""
numerical
=========

This module provides solutions using the scientific computing approach
based on discretized :term:`PDE` where a system of linear equations is
built. Additionally, this module includes ``solvers`` which is used to
solve the system  of linear equations.

List of numerical solutions:
    * ``FdmSolution``: Finite-Difference-Method
    * ``FvmSolution``: Finite-Volume-Method (not available)
    * ``FemSolution``: Finite-Element-Method (not available)

Information:
    - design pattern: inheritance, abstraction
    - base class: `Solution </api/reservoirflow.solutions.Solution.html>`_
    - base class type: ABS (abstract)
"""

__all__ = [
    "FdmSolution",
    "FemSolution",
    "FvmSolution",
    "solvers",
]


from . import solvers
from .fdm_solution import FdmSolution
from .fem_solution import FemSolution
from .fvm_solution import FvmSolution
