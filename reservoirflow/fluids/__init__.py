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
    "MultiPhaseFluid",
    "SinglePhaseFluid",
    "ThreePhaseFluid",
    "TwoPhaseFluid",
]


from .fluid import Fluid
from .multi_phase_fluid import MultiPhaseFluid
from .single_phase_fluid import SinglePhaseFluid
from .three_phase_fluid import ThreePhaseFluid
from .two_phase_fluid import TwoPhaseFluid

# __all_exports = [SinglePhaseFluid]

# for e in __all_exports:
#     e.__module__ = __name__

# __all__ = [e.__name__ for e in __all_exports]
