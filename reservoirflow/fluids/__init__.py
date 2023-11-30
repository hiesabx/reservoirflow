"""
fluids
======
"""
# .. autosummary::
#    :toctree:
#    :template: class_custom.rst
#    :recursive:

#    SinglePhase
#    MultiPhase

__all__ = [
    "SinglePhase",
    "MultiPhase",
]


from ._multiphase import MultiPhase

# from ._fluid import _Fluid
from ._single_phase import SinglePhase

# __all_exports = [SinglePhase]

# for e in __all_exports:
#     e.__module__ = __name__

# __all__ = [e.__name__ for e in __all_exports]
