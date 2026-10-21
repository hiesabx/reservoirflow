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
    "DirectionalWell",
    "MultiCellWell",
    "SingleCellWell",
    "Well",
]

from .directional_well import DirectionalWell
from .multi_cell_well import MultiCellWell
from .single_cell_well import SingleCellWell
from .well import Well
