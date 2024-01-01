"""
wells
=====

This module can be used to create wells.

Information:
    - design pattern: inheritance, abstraction
    - base class: `Well </api/reservoirflow.wells.Well.html>`_
    - base class type: ABS (abstract)
"""

__all__ = [
    "Well",
    "SingleCell",
    "MultiCell",
    "Directional",
]

from .directional import Directional
from .multi_cell import MultiCell
from .single_cell import SingleCell
from .well import Well
