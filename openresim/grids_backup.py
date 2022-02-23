#%% 1. Import Statements:
from pyrsistent import ny
from openresim.base import Base
import numpy as np


#%% 2. Grid Class: 
class Grid(Base):
    '''
    Grid class to create a grid using numpy arrays. Note that the following conventions are used: 
        1. rows for dx. 1D grids are stored as 1 column with multiple rows.  
        2. columns for dy. 2D grids are stored as a table where rows refer to x-direction while columns refere to y-direction. 
        3. layers for dz. 3D grids are stored as cubes where rows refer to x-direction, columns to y-direction, and layers for z-direction.
    '''
    def __init__(self, shape, dx, dy, dz):
        self.dim = list(dim for dim in shape if dim > 1)
        self.dims = (dx, dy, dz) # lengths: dx, dy, dz
        self.dtype = np.single # np.single, np.double
        if len(self.dim) == 1:
            self.grid_type = '1D'
        elif len(self.dim) == 2:
            self.grid_type = '2D'
        elif len(self.dim) == 3:
            self.grid_type = '3D'
        else:
            assert len(self.dim) <= 3, 'Geometry higher than 3 dimensions is not allowed!'


    def set_tops(self, values):
        self.z = self.blocks * values


    def set_rock_props(self, phi, k):
        self.phi = self.blocks * phi # np.zeros(self.shape) + phi
        self.k = self.blocks * k # np.zeros(self.shape) + k


    def set_system(self, sys_type, units):
        self.sys_units = units
        if sys_type in ['incompressible']:
            self.RHS = 0
            self.sys_type = 'incompressible'
        else:
            self.RHS = 0
            print("System type is unknown!. RHS is set to 0 'incompressible'")


    def get_neighbors(self, index):
        if len(self.dim) == 1:
            assert index > 0 and index <= self.shape-2, 'grid index is out of range'
            if index == 1:
                return [index+1]
            if index == self.shape-2:
                return [index-1]
            else:
                return [index-1, index+1]


    def get_boundaries(self, index):
        if len(self.dim) == 1:
            assert index > 0 and index <= self.shape-2, 'grid index is out of range'
            if index == 1:
                return [index-1]
            if index == self.shape-2:
                return [index+1]
            else:
                return []