"""
grids
=====


.. note::
    Flow dimension is defined based on grid shape.
    
Information:
    - design pattern: inheritance, abstraction
    - base class: `Grid </api/reservoirflow.grids.Grid.html>`_
    - base class type: ABS (abstract)
"""

__all__ = [
    "Grid",
    "RegularCartesian",
    "Radial",
    "IrregularCartesian",
]

from .grid import Grid
from .irregular_cartesian import IrregularCartesian
from .radial import Radial
from .regular_cartesian import RegularCartesian

# __all_exports = [RegularCartesian]

# for e in __all_exports:
#     e.__module__ = __name__

# __all__ = [e.__name__ for e in __all_exports]
