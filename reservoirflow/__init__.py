"""
ReservoirFlow `<reservoirflow>`
===============================

a Petroleum Reservoir Simulation Library in Python developed by 
`Zakariya Abugrin`.

Available subpackages
---------------------
fluids
    Fluid classes (required to build a Model).
grids
    Grid classes (required to build a Model).
models
    Model classes to combine Grid and Fluid to build a simulation model.
plots
    Plot classes to plot data i
scalers
    Scaler classes to scale model data.
utils
    Utilities functions for profiling, solvers, etc.
visuals
    Visual classes to visualize a Model in 3D.
wells
    Wells classes to add a well to a Model.
    
How to Use
----------

First import:
The docstring examples assume that `reservoirflow` has been imported as `rf`::

  >>> import reservoirflow as rf
    
Print the current version::

  >>> print(rf.__version__)
    
Create a `Cartesian` grid from `grids` module::

  >>> grid = rf.grids.Cartesian()
    
Create a `SinglePhase` fluid from `fluids` module::

  >>> fluid = rf.fluids.SinglePhase()
    
Construct a `BlackOil` model from `models` module::

  >>> model = rf.models.BlackOil(
        grid=grid,
        fluid=fluid,
      )
"""

__version__ = "0.0.1"

from ._base import UNITS, FACTORS
from . import fluids, grids, models, scalers, utils
