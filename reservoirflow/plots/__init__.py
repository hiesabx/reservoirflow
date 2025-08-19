"""
plots
=====
This module contains plotting classes for visualizing reservoir simulation data.

Information:
    - design pattern: inheritance, abstraction
    - base class: `Plot </api/reservoirflow.plots.Plot.html>`_
    - base class type: ABS (abstract)
"""

__all__ = [
    "Plot",
    "Plot1D",
    "Plot2D",
    "Contour1D",
]

from .plot import Plot
from .plot_1d import Plot1D
from .plot_2d import Plot2D
from .contour_1d import Contour1D
