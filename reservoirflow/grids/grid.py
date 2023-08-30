"""
Grid classes for reservoir simulation models.

This module contains all grid classes that are required to build the 
Model class. Grid class represents both the rock geometry and the rock 
properties which are required for the fluid-flow in porous-media 
calculations.
"""
# from ..base import Base
from reservoirflow.base import Base


class Grid(Base):
    """Base Grid class.

    Grid class represents both the rock geometry and the rock properties
    using numpy arrays including pyvista object for visualization.

    Parameters
    ----------
    Base : class
        Base class with universal settings.
    """

    name = "Grid"

    def __init__(self, unit, dtype, verbose, unify):
        """Returns parent Grid class.

        Parameters
        ----------
        dtype : str or `np.dtype`, optional, by default 'double'
            data type used in all arrays. Numpy dtype such as
            `np.single` or `np.double` can be used.
        unit : str ('field', 'metric'), optional, by default 'field'
            units used in input and output. Parameters can be defined as
            `unit='field'` (default) or `unit='metric'`. `units`
            attribute can be accessed from this class using
            (`Cartesian.units`) or from the base class (`Grid.units`).
        unify : bool, optional, by default False
            unify shape to be always tuple of 3 when set to True. When
            set to False, shape includes only the number of girds in
            flow direction as tuple. This option is only relevant in
            case of 1D or 2D flow. This option may be required to make
            1D and 2D shapes shapes of this class more consistent with
            each other or with 3D shape.
        verbose : bool, optional, by default False
            print information for debugging.
        """
        super().__init__(unit, dtype, verbose)
        self.unify = unify
        props_keys = ["kx", "ky", "kz", "phi", "z", "comp"]
        self.k = {}
        self.d = {}
        self.A = {}
        self.__props__ = dict.fromkeys(props_keys)


if __name__ == "__main__":
    dtype = "double"
    unit = "field"
    verbose = False
    unify = True
    grid = Grid(unit, dtype, verbose, unify)
    print(grid)
