"""
fluids
======
"""

__all__ = [
    "Fluid",
    "SinglePhase",
    "MultiPhase",
]


from .fluid import Fluid
from .multi_phase import MultiPhase
from .single_phase import SinglePhase

# __all_exports = [SinglePhase]

# for e in __all_exports:
#     e.__module__ = __name__

# __all__ = [e.__name__ for e in __all_exports]
