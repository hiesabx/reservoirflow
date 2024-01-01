"""
neurical
========

This module provides solutions using the scientific machine learning 
approach (i.e. mostly based on neural networks of deep learning) where
a specific loss function based on :term:`PDE` with specific boundary 
conditions are used to force the laws of physics within the domain of 
interest. In contrast to numerical and analytical solutions, solution 
classes based on neurical solution usually requires a training 
(i.e. ``model.fit()``) after compiling (i.e. ``model.compile()``).


.. hint:: 
    This module might be very interesting for you if you are learning 
    about scientific machine learning.

List of neurical solutions: 
    * ``PINN``: Physics-Informed-Neural-Network
    * ``DeepONet``: Deep-Operator-Network (not available)
    
.. attention::
    More efficient and state-of-the-art solutions will be added in the 
    future based on the research progress in the field of scientific 
    machine learning. 
    
Information:
    - design pattern: inheritance, abstraction
    - base class: `Solution </api/reservoirflow.solutions.Solution.html>`_
    - base class type: ABS (abstract)
"""

__all__ = [
    "PINN",
    "DeepONet",
]

from .deeponet import DeepONet
from .pinn import PINN
