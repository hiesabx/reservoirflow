"""
analytical
==========

This module provides specific analytical solutions for basic problems
where analytical solutions based on the continues form of :term:`PDE`
can be established. Usually, these problems are simpler and the
analytical solution can be considered as an exact solution which can be
used as a reference solution to evaluate and confirm other solutions (
e.g. ``numerical`` and ``neurical``) before these solutions are used in
more complex problems.

Analytical solutions are a problem specific and the list below shows the
plan to work on these solutions in the future.

List of analytical solutions:
    * ``D1P1Solution``: 1-Dimension-1-Phase
    * ``D1P2Solution``: 1-Dimension-2-Phase (not available)
    * ``D1P3Solution``: 1-Dimension-3-Phase (not available)
    * ``D2P1Solution``: 2-Dimension-1-Phase (not available)
    * ``D2P2Solution``: 2-Dimension-2-Phase (not available)
    * ``D2P3Solution``: 2-Dimension-3-Phase (not available)
    * ``D3P1Solution``: 3-Dimension-1-Phase (not available)
    * ``D3P2Solution``: 3-Dimension-2-Phase (not available)
    * ``D3P3Solution``: 3-Dimension-3-Phase (not available)

.. note::
    This list is rather a very optimistic one. For example, finding the
    analytical for a ``D3P3Solution`` might be impossible. However, constructing
    a solution based on simpler cases (e.g. ``D1P3Solution``) is a topic we
    would like to investigate.

Information:
    - design pattern: inheritance, abstraction
    - base class: `Solution </api/reservoirflow.solutions.Solution.html>`_
    - base class type: ABS (abstract)
"""

__all__ = [
    "D1P1Solution",
    "D1P2Solution",
    "D1P3Solution",
    "D2P1Solution",
    "D2P2Solution",
    "D2P3Solution",
    "D3P1Solution",
    "D3P2Solution",
    "D3P3Solution",
]

from .d1p1_solution import D1P1Solution
from .d1p2_solution import D1P2Solution
from .d1p3_solution import D1P3Solution
from .d2p1_solution import D2P1Solution
from .d2p2_solution import D2P2Solution
from .d2p3_solution import D2P3Solution
from .d3p1_solution import D3P1Solution
from .d3p2_solution import D3P2Solution
from .d3p3_solution import D3P3Solution
