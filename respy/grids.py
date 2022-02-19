#%% 1. Import Statements:
from respy.base import Base
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


#%% 2. 1D Grid Class: 
class Grid1D(Base):
    '''
    Grid class to create a grid using numpy arrays. 1D grids are stored as 1 column with multiple rows which represent x-direction.  
    '''
    name = '1D Grid'
    def __init__(self, shape, dx, dy, dz, z=None, phi=None, k=None, comp=None, comp_type=None, dtype='double', unit='us'):
        """
        
        """
        # super().__init__(unit)
        self.dim = list(dim for dim in shape if dim > 1)
        assert len(self.dim) == 1, 'Geometry higher than 1D is not allowed using Grid1D class.'
        self.dims = (dx, dy, dz) # lengths: dx, dy, dz
        self.dtype = dtype # np.single, np.double
        self.type = '1D'
        self.shape = np.max(shape) + 2 # + 2 boundaries in 1D
        self.order = np.arange(self.shape) # natural order: 0, -1 are boundaries.
        self.blocks = np.ones(self.shape, dtype='int')
        self.iBlocks = self.blocks.copy()
        self.iBlocks[[0, -1]] = 0
        self.bBlocks = np.zeros(self.shape, dtype='int') #ss.lil_matrix(self.shape, dtype='int')
        self.bBlocks[[0, -1]] = 1            
        self.boundaries = np.argwhere(self.bBlocks != 0)
        self.dx = self.blocks * dx
        self.dy = self.blocks * dy
        self.dz = self.blocks * dz
        self.area = self.dy * self.dz
        self.set_properties(phi, k, z, comp, comp_type)


    def set_porosity(self, phi, i=None):
        """
        
        """
        if not i:
            self.porosity = self.phi = self.blocks * phi # np.zeros(self.shape) + phi
        else:
            self.porosity[i] = phi
    set_phi = set_porosity


    def set_permeability(self, k, i=None):
        """
        
        """
        if not i:
            self.permeability = self.k = self.blocks * k # np.zeros(self.shape) + k
        else:
            self.permeability[i] = k  
    set_k = set_permeability

    
    def set_tops(self, z=None, i=None):
        """
        
        """
        if not i:
            self.tops = self.z = self.blocks * z
        else:
            self.tops[i] = z
    set_z = set_tops


    def set_properties(self, phi=None, k=None, z=None, comp=None, comp_type=None, i=None):
        """
        
        """
        if phi != None:
            self.set_porosity(phi, i)
        if k != None:
            self.set_permeability(k, i)
        if z != None:
            self.set_tops(z, i)
        if comp != None:
            self.set_compressibility(comp)
        if not hasattr(self, 'tops'):
            self.set_tops(0)
        if not hasattr(self, 'compressibility'):
            self.set_compressibility(0)
    set_props = set_properties


    # def set_boundaries_dict(self, loc, d):
    #     self.boundaries_dict[loc] = 
        
            
    # set_b_dict = set_boundaries_dict


    def get_boundaries(self, i=None):
        """
        Arguments:
            - i: index in x-direction.
        """
        if i == None:
            return [0, self.shape-1]
        assert i > 0 and i <= self.shape-2, 'grid index is out of range(1, {}).'.format(self.shape-1)
        if i == 1:
            if i == self.shape-2:
                return [i-1, i+1]
            return [i-1]
        if i == self.shape-2:
            return [i+1]
        else:
            return []
    get_b = get_boundaries


    def get_neighbors(self, i):
        """
        Arguments:
            - i: index in x-direction.
        """
        assert i > 0 and i <= self.shape-2, 'grid index is out of range.'
        if i == 1:
            return [i+1]
        if i == self.shape-2:
            return [i-1]
        else:
            return [i-1, i+1]
    get_n = get_neighbors


#%%
if __name__ == '__main__':
    # Define:
    grid = Grid1D(shape=(4, 1, 1), dx=300, dy=350, dz=40, phi=0.27, k=270)
    # Doc:
    print(grid.__doc__)
    # Setters:
    grid.set_comp(1e-6)
    # grid.set_props(0.3, 11)
    # grid.set_phi(0.2)
    # grid.set_k(10)
    # grid.set_tops(10)
    # Getters:
    # print(grid.get_boundaries())
    print(*grid.get_boundaries(1), '1', *grid.get_neighbors(1))
    print(grid) # or grid.report()
    # print(grid.set_boundaries(loc=[0, -1], conditions=['pressure', 'rate'], values=[4000, 0]))