"""
Grid classes for reservoir simulation models.

This module contains all grid classes that are required to build the 
Model class. Grid class represents both the rock geometry and the rock 
properties which are required for the fluid-flow in porous-media 
calculations.
"""
from abc import ABC, abstractmethod

# from ..base import Base
from reservoirflow.base import Base


class Grid(ABC, Base):
    """Abstract grid class.

    Grid class represents both the rock geometry and the rock properties
    using numpy arrays including pyvista object for visualization.

    .. attention::

        This is an abstract class and can't be instantiated. This class
        is only used as a parent for other classes of ``grids`` module.

    Returns
    -------
    Grid
        Grid object.
    """

    name = "Grid"

    def __init__(self, unit, dtype, verbose, unify):
        """Construct grid object.

        Parameters
        ----------
        dtype : str or `np.dtype`, optional
            data type used in all arrays. Numpy dtype such as
            `np.single` or `np.double` can be used.
        unit : str ('field', 'metric', 'lab'), optional
            unit used in input and output. Both `units` and `factors`
            attributes will be updated based on the selected `unit` and
            can be accessed directly from this class.
        unify : bool, optional
            unify shape to be always tuple of 3 when set to True. When
            set to False, shape includes only the number of girds in
            flow direction as tuple. This option is only relevant in
            case of 1D or 2D flow. This option may be required to make
            1D and 2D shapes shapes of this class more consistent with
            each other or with 3D shape.
        verbose : bool, optional
            print information for debugging.
        """
        super().__init__(unit, dtype, verbose)
        self.unify = unify
        props_keys = ["kx", "ky", "kz", "phi", "z", "comp"]
        self.k = {}
        self.d = {}
        self.A = {}
        self.__props__ = dict.fromkeys(props_keys)

    def set_comp(self, comp: float):
        """Set grid compressibility.

        Parameters
        ----------
        comp : float
            grid compressibility.

        Raises
        ------
        ValueError
            Compressibility smaller than zero is not allowed.
        """
        self.comp = comp
        if comp == 0:
            self.comp_type = "incompressible"
        elif comp > 0:
            self.comp_type = "compressible"
        else:
            self.comp_type = None
            raise ValueError("Compressibility smaller than zero is not allowed.")

    # -------------------------------------------------------------------------
    # Synonyms:
    # -------------------------------------------------------------------------

    def allow_synonyms(self):
        """Allow full descriptions.

        This function maps functions as following:

        .. highlight:: python
        .. code-block:: python

            self.set_compressibility = self.set_comp
            self.compressibility = self.comp
            self.compressibility_type = self.comp_type

        """
        self.set_compressibility = self.set_comp
        self.compressibility = self.comp
        self.compressibility_type = self.comp_type

    # -------------------------------------------------------------------------
    # End
    # -------------------------------------------------------------------------


if __name__ == "__main__":
    dtype = "double"
    unit = "field"
    verbose = False
    unify = True
    grid = Grid(unit, dtype, verbose, unify)
    print(grid)
