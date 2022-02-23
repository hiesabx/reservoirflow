#%% 1. Import Statements: 
from openresim.base import Base
import numpy as np


#%% 2. Fluid Class: 
class Fluid(Base):
    '''
    Fluid class to create a fluid.
    '''
    name = 'Single Phase Fluid'
    def __init__(self, mu=None, B=None, rho=None, comp=None, comp_type=None, dtype='double', unit='field'):
        # super().__init__(unit)
        self.set_properties(mu, B, rho, comp, comp_type)


    def set_viscosity(self, mu):
        self.viscosity = self.mu = mu
    set_mu = set_viscosity


    def set_formation_volume_factor(self, B):
        self.formation_volume_factor = self.FVF = self.B = B
    set_B = set_FVF = set_formation_volume_factor


    def set_density(self, rho):
        self.density = self.rho = rho
        self.gravity = self.g = self.factors['gravity conversion'] * self.density * self.factors['gravitational acceleration']
    set_rho = set_density


    def set_properties(self, mu=None, B=None, rho=None, comp=None, comp_type=None):
        if mu != None:
            self.set_viscosity(mu)
        if B != None:
            self.set_formation_volume_factor(B)
        if rho != None:
            self.set_density(rho)
        if comp != None:
            self.set_compressibility(comp)
        if not hasattr(self, 'density'):
            self.set_density(0)
        if not hasattr(self, 'compressibility'):
            self.set_compressibility(0)
    set_props = set_properties


if __name__ == '__main__':
    fluid = Fluid(mu=0.5 , B=1, rho=1, unit='metric')
    fluid.set_units('metric')
    fluid.set_density(10)
    print(fluid)

# %%
