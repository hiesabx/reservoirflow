from abc import ABC, abstractmethod

from reservoirflow.base import Base
from reservoirflow.solutions.compiler import Compiler


class Model(ABC, Base):
    """Abstract model class.

    Model class used to create a reservoir simulation model. Model class
    represents the fluid flow process in a reservoir due to pressure
    change cause by boundary conditions and/or by (production or
    injection) wells.

    .. attention::

        This is an abstract class and can't be instantiated. This class
        is only used as a parent for other classes of ``models`` module.

    Returns
    -------
    Model
        Model object.
    """

    name = "Model"

    def __init__(
        self, 
        unit,
        dtype,
        verbose,
    ):
        """Construct model object.

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
        self.solution = None
        self.solutions = {}

    def set_comp(self, comp: float):
        """Set model compressibility

        Parameters
        ----------
        comp : float
            model compressibility which is the total compressibility
            of Grid and Fluid classes.

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

    def compile(
        self,
        stype: str,
        method: str,
        sparse: bool=True,
    ):
        """Build a solution (equation system) for the model.

        This function will add ``model.solution`` and ``model.compiler``
        which are defined based on the parameters used in this function.
        In addition, methods ``model.solve()`` and ``model.run()`` are
        actually mapped from ``model.solution``.

        Note that methods ``solve()`` and ``run()`` in addition
        to many other methods can be accessed using the solution object
        (e.g. ``model.solution.run()``). For more information about the
        assigned solution object, check the
        `documentation </api/reservoirflow.solutions.html#>`_.

        Parameters
        ----------
        stype : str
            solution type in ``['numerical', 'analytical', 'neurical']``.
        method : str
            solution method as following:
            
            - 'numerical' methods: ``['FDM', 'FVM', 'FEM']``.
            - 'analytical' methods: ``['1D1P', '1D2P', etc.]``.
            - 'neurical' methods: ``['PINN', 'DeepONet', etc.]``.
            
        sparse : bool, optional, default: True
            using sparse computing for a better performance.
        """
        self.compiler = Compiler(self, stype, method, sparse)
        self.solutions[self.compiler.method] = self.solution
        self.solve = self.solution.solve
        self.run = self.solution.run
        
    
    def set_solution(self, method):
        if method in self.solutions:
            self.solution = self.solutions[method]
            self.pressures = self.solution.pressures
            self.rates = self.solution.rates
            self.solve = self.solution.solve
            self.run = self.solution.run
        else:
            raise ValueError("Solution method was not compiled.")
        
    def get_solutions(self):
        return list(self.solutions.keys())

    def solve(self, **kwargs):
        """Solve a single timestep.

        .. attention::
            This method is not available until the model is compiled
            using ``model.compile()``.

        Once the model is compiled, the documentation of the assigned
        solution can be accessed using one of the following methods:

        >>> help(model.solve) # or help(model.solution.solve)
        >>> print(model.solve.__doc__) # or print(model.solution.solve.__doc__)
        """
        print(
            "The model is not compiled.",
            "Use model.compile() to add solve() and run() methods.",
        )
        return None

    def run(self, **kwargs):
        """Solve multiple timesteps.

        .. attention::
            This method is not available until the model is compiled
            using ``model.compile()``.

        Once the model is compiled, the documentation of the assigned
        solution can be accessed using one of the following methods:

        >>> help(model.run) # or help(model.solution.run)
        >>> print(model.run.__doc__) # or print(model.solution.run.__doc__)
        """
        self.solve()

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
    model = Model(unit, dtype, verbose)
    print(model)
