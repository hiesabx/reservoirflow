"""
models
======

Information:
    - design pattern: inheritance, abstraction
    - base class: `Model </api/reservoirflow.models.Model.html>`_
    - base class type: ABS (abstract)
"""

__all__ = [
    "Model",
    "BlackOil",
    "Compositional",
    "Thermal",
]


from .black_oil import BlackOil
from .compositional import Compositional
from .model import Model
from .thermal import Thermal

# __all_exports = [BlackOil]

# for e in __all_exports:
#     e.__module__ = __name__

# __all__ = [e.__name__ for e in __all_exports]
