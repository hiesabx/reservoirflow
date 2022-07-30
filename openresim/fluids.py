#%% 1. Import Statements:
from openresim.base import Base
import numpy as np


#%% 2. Fluid Class:
class Fluid(Base):
    def __init__(self, unit, dtype, verbose):
        super().__init__(unit, dtype, verbose)


class SinglePhase(Fluid):
    """
    Fluid class to create a fluid.
    """

    name = "Single Phase Fluid"

    def __init__(
        self,
        mu=None,
        B=None,
        rho=None,
        comp=None,
        comp_type=None,
        dtype="double",
        unit="field",
        verbose=False,
    ):
        super().__init__(unit, dtype, verbose)
        self.set_props(mu, B, rho, comp, comp_type)

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

    def set_props(self, mu=None, B=None, rho=None, comp=None, comp_type=None):
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
