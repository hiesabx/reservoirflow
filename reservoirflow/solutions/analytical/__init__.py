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
    * ``D1P1``: 1-Dimension-1-Phase
    * ``D1P2``: 1-Dimension-2-Phase (not available)
    * ``D1P3``: 1-Dimension-3-Phase (not available)
    * ``D2P1``: 2-Dimension-1-Phase (not available)
    * ``D2P2``: 2-Dimension-2-Phase (not available)
    * ``D2P3``: 2-Dimension-3-Phase (not available)
    * ``D3P1``: 3-Dimension-1-Phase (not available)
    * ``D3P2``: 3-Dimension-2-Phase (not available)
    * ``D3P3``: 3-Dimension-3-Phase (not available)

.. note::
    This list is rather a very optimistic one. For example, finding the 
    analytical for a ``D3P3`` might be impossible. However, constructing 
    a solution based on simpler cases (e.g. ``D1P3``) is a topic we
    would like to investigate.

Information:
    - design pattern: inheritance, abstraction
    - base class: `Solution </api/reservoirflow.solutions.Solution.html>`_
    - base class type: ABS (abstract)
"""

__all__ = [
    "D1P1",
    "D1P2",
    "D1P3",
    "D2P1",
    "D2P2",
    "D2P3",
    "D3P1",
    "D3P2",
    "D3P3",
]

from .d1p1 import D1P1
from .d1p2 import D1P2
from .d1p3 import D1P3
from .d2p1 import D2P1
from .d2p2 import D2P2
from .d2p3 import D2P3
from .d3p1 import D3P1
from .d3p2 import D3P2
from .d3p3 import D3P3
