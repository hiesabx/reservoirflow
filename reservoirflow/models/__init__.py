"""
models
======

Information:
    - design pattern: inheritance, abstraction
    - base class: `Model </api/reservoirflow.models.Model.html>`_
    - base class type: ABS (abstract)
"""

__all__ = [
    "BlackOilModel",
    "CompositionalModel",
    "Model",
    "ThermalModel",
]


from .black_oil_model import BlackOilModel
from .compositional_model import CompositionalModel
from .model import Model
from .thermal_model import ThermalModel

# __all_exports = [BlackOilModel]

# for e in __all_exports:
#     e.__module__ = __name__

# __all__ = [e.__name__ for e in __all_exports]
