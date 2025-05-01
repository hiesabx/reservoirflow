from abc import ABC, abstractmethod

from reservoirflow.base import Base


class Fluid(ABC, Base):
    """Abstract fluid class.

    .. attention::

        This is an abstract class and can't be instantiated. This class
        is only used as a parent for other classes of ``fluids`` module.

    Returns
    -------
    Fluid
        Fluid object.
    """

    name = "Fluid"

    def __init__(
        self, 
        unit,
        dtype,
        verbose,
    ):
        """Construct fluid object.

        Parameters
        ----------
        unit : str ('field', 'metric', 'lab'), optional
            unit used in input and output. Both `units` and `factors`
            attributes will be updated based on the selected `unit` and
            can be accessed directly from this class.
        dtype : str or `np.dtype`, optional
            data type used in all arrays. Numpy dtype such as
            `np.single` or `np.double` can be used.
        verbose : bool, optional
            print information for debugging.
        """
        super().__init__(unit, dtype, verbose)

    def set_comp(self, comp: float):
        """Set fluid compressibility.

        Parameters
        ----------
        comp : float
            fluid compressibility.

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
    fluid = Fluid(unit, dtype, verbose)
    print(fluid)
