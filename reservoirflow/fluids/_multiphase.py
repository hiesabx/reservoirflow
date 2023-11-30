"""
MultiPhase
-----------

This module is used to create a single phase fluid objects.
"""

from reservoirflow.fluids._fluid import _Fluid


class MultiPhase(_Fluid):
    """Single phase fluid class.

    Returns
    -------
    Fluid
        SinglePhase fluid object.
    """

    name = "Single Phase Fluid"

    def __init__(
        self,
        mu: float = None,
        B: float = None,
        rho: float = None,
        comp: float = None,
        dtype="double",
        unit: str = "field",
        verbose: bool = False,
    ):
        """_summary_

        Parameters
        ----------
        mu : float, optional
            fluid viscosity.
        B : _type_, optional
            fluid formation volume factor.
        rho : _type_, optional
            fluid density.
        comp : _type_, optional
            fluid compressibility.
        dtype : str or `np.dtype`, optional
            data type used in all arrays. Numpy dtype such as
            `np.single` or `np.double` can be used.
        unit : str ('field', 'metric', 'lab'), optional
            units used in input and output. Parameters can be defined as
            `unit='field'` (default), `unit='metric'`, or `unit='lab'`.
            `units`attribute can be accessed from this class using
            (`Model.units`).
        verbose : bool, optional
            print information for debugging.

        Notes
        -----
        .. note::
            Units are defined based on `unit` argument, for more
            details, check
            `Units & Factors </user_guide/units_factors/units_factors.html>`_.
            For definitions, check
            `Glossary </user_guide/glossary/glossary.html>`_.
        """
        super().__init__(unit, dtype, verbose)
        self.set_props(mu, B, rho, comp)

    def set_mu(self, mu):
        self.mu = mu

    def set_B(self, B):
        self.B = B

    def set_rho(self, rho):
        self.rho = rho
        self.g = (
            self.factors["gravity conversion"]
            * self.rho
            * self.factors["gravitational acceleration"]
        )

    def set_props(self, mu=None, B=None, rho=None, comp=None):
        if mu != None:
            self.set_mu(mu)
        if B != None:
            self.set_B(B)
        if rho != None:
            self.set_rho(rho)
        if comp != None:
            self.set_comp(comp)
        if not hasattr(self, "rho"):
            self.set_rho(0)
        if not hasattr(self, "comp"):
            self.set_comp(0)

    # -------------------------------------------------------------------------
    # Synonyms:
    # -------------------------------------------------------------------------

    def allow_synonyms(self):
        self.set_viscosity = self.set_mu
        self.viscosity = self.mu
        self.set_rho = self.set_rho
        self.density = self.rho
        self.gravity = self.g
        self.set_formation_volume_factor = self.set_FVF = self.set_B
        self.formation_volume_factor = self.FVF = self.B

        self.set_props = self.set_props

    # -------------------------------------------------------------------------
    # End
    # -------------------------------------------------------------------------


if __name__ == "__main__":
    fluid = SinglePhase(mu=0.5, B=1, rho=1, unit="metric")
    fluid.set_units("metric")
    fluid.set_rho(10)
    print(fluid)
