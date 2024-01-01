"""
fluids
======


Information:
    - design pattern: inheritance, abstraction
    - base class: `Fluid </api/reservoirflow.fluids.Fluid.html>`_
    - base class type: ABS (abstract)
"""

__all__ = [
    "Fluid",
    "SinglePhase",
    "TwoPhase",
    "ThreePhase",
    "MultiPhase",
]


from .fluid import Fluid
from .multi_phase import MultiPhase
from .single_phase import SinglePhase
from .three_phase import ThreePhase
from .two_phase import TwoPhase

# __all_exports = [SinglePhase]

# for e in __all_exports:
#     e.__module__ = __name__

# __all__ = [e.__name__ for e in __all_exports]
