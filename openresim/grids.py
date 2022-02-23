#%% 1. Import Statements:
from openresim.base import Base
import numpy as np
import pyvista as pv


#%% 2. Grid Class: 
class Grid(Base):
    '''
    Grid class to create a grid using numpy arrays. Note that the following conventions are used: 
        1. rows for dx. 1D grids are stored as 1 column with multiple rows.  
        2. columns for dy. 2D grids are stored as a table where rows refer to x-direction while columns refere to y-direction. 
        3. layers for dz. 3D grids are stored as cubes where rows refer to x-direction, columns to y-direction, and layers for z-direction.
    '''
    def __init__(self):
        pass


#%% 2. 1D Grid Class: 
class Grid1D(Base):
    '''
    Grid class to create a grid using numpy arrays. 1D grids are stored as 1 column with multiple rows which represent x-direction.  
    '''
    name = '1D Grid'
    def __init__(self, nx, ny, nz, dx, dy, dz, z=None, phi=None, k=None, comp=None, dtype='double', unit='field'):
        """
        
        """
        # super().__init__(unit)
        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.get_dimension() # self.dimension = self.D
        self.get_shape() # self.shape
        self.get_order(type='natural')
        assert self.D <= 1, 'Geometry higher than 1D is not allowed using Grid1D class.'
        self.dtype = dtype # np.single, np.double
        self.type = 'cartesian'
        self.blocks = np.ones(self.nx + 2, dtype='int')
        self.i_blocks = self.blocks.copy()
        self.i_blocks[[0, -1]] = 0
        self.b_blocks = np.zeros(self.nx + 2, dtype='int') #ss.lil_matrix(self.shape, dtype='int')
        self.b_blocks[[0, -1]] = 1            
        self.boundaries = self.get_boundaries()
        self.dx = self.blocks * dx
        self.dy = self.blocks * dy
        self.dz = self.blocks * dz
        self.meshgrid = np.meshgrid(self.dx[:self.nx], self.dy[:self.ny], self.dz[:self.nz])
        # self.get_pv_grid()
        self.get_area() # self.area
        self.set_properties(phi, k, z, comp)
        self.get_G() # self.G


    def get_corners(self, show_boundary, verbose=False):
        xcorn = np.insert(self.dx.cumsum(), 0, 0) # or np.append([0],arr) or np.concatenate([0],arr)
        if show_boundary:
            i = 2
        else:
            i = 0
            xcorn = xcorn[1:-1]
        xcorn = np.repeat(xcorn, 2)
        xcorn = xcorn[1:-1]
        xcorn = np.tile(xcorn, 4*self.ny*self.nz)
        
        ycorn = np.arange(0, (self.ny+1)*self.dy[1], self.dy[1])
        ycorn = np.repeat(ycorn, 2)
        ycorn = ycorn[1:-1]
        ycorn = np.tile(ycorn, (2*(self.nx+i), 2*self.nz))
        ycorn = np.transpose(ycorn)
        ycorn = ycorn.flatten()

        zcorn = np.arange(0, (self.nz+1)*self.dz[1], self.dz[1])
        zcorn = np.repeat(zcorn, 2)
        zcorn = zcorn[1:-1]
        zcorn = np.repeat(zcorn, (4*(self.nx+i)*self.ny))
        if verbose: print(xcorn.shape, ycorn.shape, zcorn.shape)
        corners = np.stack((xcorn, ycorn, zcorn))
        self.corners = corners.transpose()
        return self.corners

    
    def get_pv_grid(self, show_boundary, method=2, verbose=False):
        if show_boundary:
            dims = np.array((self.nx+2,self.ny,self.nz))+1
        else:
            dims = np.array((self.nx,self.ny,self.nz))+1
        if method == 1:
            self.pv_grid = pv.UniformGrid(dims=dims)
            self.pv_grid.spacing = (self.dx[0], self.dy[0], self.dz[0])
        else:
            self.pv_grid = pv.ExplicitStructuredGrid(dims, self.get_corners(show_boundary))
            #self.pv_grid = self.pv_grid.compute_connectivity()
        if verbose: print(self.pv_grid)
        return self.pv_grid


    def get_is_homogeneous(self):
        if  all(self.dx[1:-1] == self.dx[1]) & \
            all(self.dy[1:-1] == self.dy[1]) & \
            all(self.dz[1:-1] == self.dz[1]) & \
            all(self.k[1:-1] == self.k[1]) & \
            all(self.porosity[1:-1] == self.porosity[1]):
            self.is_homogeneous = True # homogeneous
        else:
            self.is_homogeneous = False # heterogeneous


    def mean(self, property, type='geometric'):
        if self.is_homogeneous:
            return property[1:]
        else:
            if type == 'geometric':
                return (property[:-1] + property[1:])/2
            else:
                return (property[:-1] + property[1:])/2


    def get_G(self):
        if self.is_homogeneous:
            self.G = self.factors['transmissibility conversion'] * \
                self.mean(self.k) * self.mean(self.area) / self.mean(self.dx)
        else:
            self.G = 2 * self.factors['transmissibility conversion'] / \
                (
                    (self.dx[:-1] / (self.area[:-1] * self.k[:-1])) + \
                    (self.dx[1:] / (self.area[1:] * self.k[1:]))
                )


    def set_porosity(self, phi, i=None):
        """
        
        """
        if not i:
            self.porosity = self.phi = self.blocks * phi # np.zeros(self.shape) + phi
        else:
            self.porosity[i] = phi
            self.get_is_homogeneous()
    set_phi = set_porosity


    def set_permeability(self, k, i=None):
        """
        
        """
        if not i:
            self.permeability = self.k = self.blocks * k # np.zeros(self.shape) + k
        else:
            self.permeability[i] = k
            self.get_is_homogeneous()
    set_k = set_permeability

    
    def set_tops(self, z=None, i=None):
        """
        
        """
        if not i:
            self.tops = self.z = self.blocks * z
        else:
            self.tops[i] = z
    set_z = set_tops


    def set_properties(self, phi=None, k=None, z=None, comp=None, i=None):
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
        # Set default values if not defined:
        if not hasattr(self, 'tops'):
            self.set_tops(0)
        if not hasattr(self, 'compressibility'):
            self.set_compressibility(0)
        if not hasattr(self, 'is_homogeneous'):
            self.get_is_homogeneous()
    set_props = set_properties


    def get_area(self):
        self.area = self.A = self.dy * self.dz
    get_A = get_area


    def get_dimension(self):
        self.dimension = self.D = sum([1 if n>1 else 0 for n in (self.nx, self.ny, self.nz)])
    get_D = get_dimension

    
    def get_shape(self):
        self.shape = np.array((self.nx, self.ny, self.nz), dtype='int')
    get_s = get_shape


    def get_order(self, type='natural'):
        if type == 'natural':
            self.order = np.arange(self.nx + 2) # natural order: 0, -1 are boundaries.
        else:
            raise ValueError("""Order type is not supported.\n
                            Supported order types: ['natural']""")
    
    def get_boundaries(self, i=None):
        """
        Arguments:
            - i: index in x-direction.
        """
        if i == None:
            return [0, self.nx + 1] # np.argwhere(self.b_blocks == 1)
        assert i > 0 and i <= self.nx, 'grid index is out of range(1, {}).'.format(self.nx + 1)
        if self.nx == 1:
            return [0, 2]
        else:
            if i == 1:
                return [0]
            if i == self.nx:
                return [self.nx + 1]
            else:
                return []
    get_b = get_boundaries


    def get_neighbors(self, i):
        """
        Arguments:
            - i: index in x-direction.
        """
        assert i > 0 and i <= self.nx, 'grid index is out of range.'
        if self.nx == 1:
            return []
        else:
            if i == 1:
                return [i+1]
            if i == self.nx:
                return [i-1]
            else:
                return [i-1, i+1]
    get_n = get_neighbors


#%%
if __name__ == '__main__':
    # Define:
    grid = Grid1D(nx=10, ny=1, nz=1, dx=300, dy=350, dz=40, phi=0.27, k=270)
    # Doc:
    print(grid.__doc__)
    # Setters:
    grid.set_comp(1e-6)
    # grid.set_permeability(200, 1)
    # grid.set_props(0.3, 11)
    # grid.set_phi(0.2)
    # grid.set_k(10)
    # grid.set_tops(10)
    # Getters:
    # print(grid.get_boundaries())
    print(*grid.get_boundaries(1), '1', *grid.get_neighbors(1))
    print(grid) # or grid.report()
    # print(grid.set_boundaries(loc=[0, -1], conditions=['pressure', 'rate'], values=[4000, 0]))