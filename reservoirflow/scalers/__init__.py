"""
scalers
=======

This module provides scaler classes which are used to transform (i.e.
scale) and/or inverse transform (i.g. descale) input data in column-wise
fashion. These scalers are used to scale simulation data (e.g. between 1
and -1) based on the selected activation functions used in neural-
networks to achieve successful and efficient training process.

Information:
    - design pattern: inheritance, abstraction
    - base class: `Scaler </api/reservoirflow.scalers.Scaler.html>`_
    - base class type: ABS (abstract)
"""

__all__ = [
    "Scaler",
    # "ScalerDC",
    "Dummy",
    "MinMax",
]

from .dummy import Dummy
from .minmax import MinMax
from .scaler import Scaler

# from .scalerdc import ScalerDC
