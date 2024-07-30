"""
SinglePhase
===========
"""

from reservoirflow.fluids.fluid import Fluid


class SinglePhase(Fluid):
    """SinglePhase fluid class.

    Returns
    -------
    Fluid
        SinglePhase fluid object.
    """

    name = "SinglePhase"

    def __init__(
        self,
        mu: float = None,
        B: float = None,
        rho: float = None,
        comp: float = None,
        unit: str = "field",
        dtype="double",
        verbose: bool = False,
    ):
        """Create SinglePhase Fluid.

        Parameters
        ----------
        mu : float, optional
            fluid viscosity.
        B : float, optional
            fluid formation volume factor.
        rho : float, optional
            fluid density.
        comp : float, optional
            fluid compressibility.
        unit : str ('field', 'metric', 'lab'), optional
            unit used in input and output. Both `units` and `factors`
            attributes will be updated based on the selected `unit` and
            can be accessed directly from this class.
        dtype : str or `np.dtype`, optional
            data type used in all arrays. Numpy dtype such as
            `np.single` or `np.double` can be used.
        verbose : bool, optional
            print information for debugging.

        Notes
        -----
        .. note::
            Both attributes units and factors are defined based on `unit`
            argument, for more details, check
            `Units & Factors </user_guide/units_factors/units_factors.html>`_.
            For definitions, check
            `Glossary </user_guide/glossary/glossary.html>`_.

        Examples
        --------

        .. testcode::

            >>> import reservoirflow as rf
            >>> fluid = rf.fluids.SinglePhase(
            ...         mu=0.5,
            ...         B=1,
            ...         rho=50,
            ...         comp=1e-5,
            ...         unit="field",
            ...         )
            >>> print(fluid)
            
        .. code-block:: python
            :linenos:
            
            import reservoirflow as rf
            fluid = rf.fluids.SinglePhase(
                    mu=0.5,
                    B=1,
                    rho=50,
                    comp=1e-5,
                    unit="field",
                    )
            print(fluid)
        """
        super().__init__(unit, dtype, verbose)
        self.set_props(mu, B, rho, comp)

    def set_mu(self, mu):
        """Set fluid viscosity.

        Parameters
        ----------
        mu : float
            fluid viscosity.
        """
        self.mu = mu  #: Fluid viscosity.

    def set_B(self, B: float):
        """Set fluid formation volume factor (FVF).

        Parameters
        ----------
        B : float
            fluid formation volume factor.
        """
        self.B = B  # Fluid formation volume factor (FVF).

    def set_rho(self, rho: float):
        """Set fluid density.

        Parameters
        ----------
        rho : float
            fluid density.
        """
        self.rho = rho  #: Fluid density.
        self.g = (
            self.factors["gravity conversion"]
            * self.rho
            * self.factors["gravitational acceleration"]
        )

    def set_props(
        self,
        mu: float = None,
        B: float = None,
        rho: float = None,
        comp: float = None,
    ):
        """Set fluid properties.

        This function allows to set/update fluid properties at once.

        Parameters
        ----------
        mu : float, optional
            fluid viscosity.
        B : float, optional
            fluid formation volume factor.
        rho : float, optional
            fluid density.
        comp : float, optional
            fluid compressibility.
        """
        if mu is not None:
            self.set_mu(mu)
        if B is not None:
            self.set_B(B)
        if rho is not None:
            self.set_rho(rho)
        if comp is not None:
            self.set_comp(comp)
        if not hasattr(self, "rho"):
            self.set_rho(0)
        if not hasattr(self, "comp"):
            self.set_comp(0)

    # -------------------------------------------------------------------------
    # Synonyms:
    # -------------------------------------------------------------------------

    def allow_synonyms(self):
        """Allow full descriptions.

        This function maps functions as following:

        .. code-block:: python

            self.set_viscosity = self.set_mu
            self.viscosity = self.mu
            self.set_density = self.set_rho
            self.density = self.rho
            self.gravity = self.g
            self.set_formation_volume_factor = self.set_FVF = self.set_B
            self.formation_volume_factor = self.FVF = self.B
            self.set_properties = self.set_props

        """
        self.set_viscosity = self.set_mu
        self.viscosity = self.mu
        self.set_density = self.set_rho
        self.density = self.rho
        self.gravity = self.g
        self.set_formation_volume_factor = self.set_FVF = self.set_B
        self.formation_volume_factor = self.FVF = self.B
        self.set_properties = self.set_props

    # -------------------------------------------------------------------------
    # End
    # -------------------------------------------------------------------------


if __name__ == "__main__":
    fluid = SinglePhase(mu=0.5, B=1, rho=1, unit="metric")
    fluid.set_units("metric")
    fluid.set_rho(10)
    print(fluid)
