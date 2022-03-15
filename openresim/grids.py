#%% Import Statements:
from openresim.base import Base
import numpy as np
import pyvista as pv
import math


#%% Grid Class: 
class Grid(Base):
    '''
    Grid class to create a grid using numpy arrays. Note that the following conventions are used: 
        1. rows for dx. 1D grids are stored as 1 column with multiple rows.  
        2. columns for dy. 2D grids are stored as a table where rows refer to x-direction while columns refere to y-direction. 
        3. layers for dz. 3D grids are stored as cubes where rows refer to x-direction, columns to y-direction, and layers for z-direction.
    '''
    def __init__(self, nx, ny, nz, dtype, unit):
        super().__init__(unit)
        self.type = 'cartesian'
        self.dtype = dtype # np.single, np.double
        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.get_dimension() # > self.dimension, self.D
        
        
    def get_dimension(self):
        self.dimension = self.D = sum([1 if n>1 else 0 for n in (self.nx, self.ny, self.nz)])
        return self.dimension
    get_D = get_dimension
    
    
    def get_shape(self):
        self.shape = np.array((self.nx, self.ny, self.nz), dtype='int')
        return self.shape
    get_s = get_shape
    
    
    def get_area(self):
        if self.D == 1:
            self.area = self.A = self.dy * self.dz
            return self.area
        else:
            return None
    get_A = get_area


    def get_volume(self):
        if self.D == 1:
            self.volume = self.V = self.dx * self.dy * self.dz
            return self.volume
        else:
            return None
    get_V = get_volume
    
    


        # if cell_id==0:
        #     return 0,0,0
        # k=math.ceil(cell_id/(self.nx*self.ny))-1
        # j=math.ceil((cell_id-(self.nx*self.ny)*k)/self.nx)-1
        # i=math.ceil(cell_id-(self.nx*self.ny*k)-self.nx*j)
        # return i,j,k


#%% CartGrid Class: 
class CartGrid(Grid):
    '''CartGrid class to create a explicit structured grid.
     
    1D grids are stored as 1 column with multiple rows which represent x-direction. Flow only allowed in x-direction.
    Parameters
    nx: int
        number of grids in x-direction.
    ny: int
        number of grids in y-direction.
    nz: int
        number of grids in z-direction.
    dx: int, list, or array
        grid width in x-direction.
    dy: int, list, or array
        grid width in y-direction.
    dz: int, list, or array
        grid width in z-direction.
    z: int, list, or array (Default: 0)
    phi: 
    k: 
    comp:
    dtype: data type str or np.dtype  (Default: double)
        array data type.
    unit: str of ['field', 'metric'] (Default: field)
        grid properties' unit.
    '''
    name = 'CartGrid'
    
    def __init__(self, nx, ny, nz, dx, dy, dz, z=None, phi=None, k=None, comp=None, dtype='double', unit='field'):
        super().__init__(nx, ny, nz, dtype, unit)
        self.get_shape()
        self.get_flow_direction()
        
        if self.D == 1:
            if self.flowdir == 'x':
                self.nx_b = self.nx + 2
                self.ny_b = self.ny
                self.nz_b = self.nz
                flowshape = self.nx_b
            elif self.flowdir == 'y':
                self.nx_b = self.nx
                self.ny_b = self.ny + 2
                self.nz_b = self.nz
                flowshape = self.ny_b
            elif self.flowdir == 'z':
                self.nx_b = self.nx
                self.ny_b = self.ny
                self.nz_b = self.nz
                flowshape = self.nz_b
        elif self.D == 2:
            if self.flowdir == 'xy':
                self.nx_b = self.nx + 2
                self.ny_b = self.ny + 2
                self.nz_b = self.nz
                flowshape = self.nx_b
            elif self.flowdir == 'xz':
                self.nx_b = self.nx + 2
                self.ny_b = self.ny + 2 # wrong: we do not need boundary
                self.nz_b = self.nz
                flowshape = self.nx_b
            elif self.flowdir == 'yz':
                self.nx_b = self.nx
                self.ny_b = self.ny
                self.nz_b = self.nz
                flowshape = self.nx_b
        elif self.D == 3:
            self.nx_b = self.nx + 2
            self.ny_b = self.ny + 2
            self.nz_b = self.nz # no boundaries in z-direction
            flowshape = self.nx_b
        else:
            raise ValueError('Geometry higher than 3D is not allowed using CartGrid class.')
        
        self.blocks = np.ones(flowshape, dtype='int')
        self.i_blocks = self.blocks.copy()
        self.i_blocks[[0, -1]] = 0
        self.b_blocks = np.zeros(flowshape, dtype='int') #ss.lil_matrix(self.shape, dtype='int')
        self.b_blocks[[0, -1]] = 1
        
        
        
        print(self.flowdir)
        if self.D == 1:
            self.dx = np.ones(flowshape, dtype='int') * dx
            self.dy = np.ones(flowshape, dtype='int') * dy
            self.dz = np.ones(flowshape, dtype='int') * dz
        elif self.D == 2:
            self.dx = np.ones(self.nx_b, dtype='int') * dx
            self.dy = np.ones(self.ny_b, dtype='int') * dy
            if self.flowdir == 'xy':
                self.dz = np.ones(self.nx_b, dtype='int') * dz
            elif self.flowdir == 'xz':
                self.dz = np.ones(self.nz_b, dtype='int') * dz
        elif self.D == 3:  
            self.dx = np.ones(self.nx_b, dtype='int') * dx
            self.dy = np.ones(self.ny_b, dtype='int') * dy
            self.dz = np.ones(self.nz_b, dtype='int') * dz
        else:
            raise ValueError('Geometry higher than 3D is not allowed using CartGrid class.')
            
        self.set_properties(phi, k, z, comp)
        self.__get_properties__()


    def set_properties(self, phi=None, k=None, z=None, comp=None, i=None):
        """
        Set properties
        """
        # Set user defined properties:
        if phi is not None:
            self.set_porosity(phi, i)
        if k is not None:
            self.set_permeability(k, i)
        if z is not None:
            self.set_tops(z, i)
        if comp is not None:
            self.set_compressibility(comp)
        # Set default values if not defined:
        if not hasattr(self, 'tops'):
            self.set_tops(0)
        if not hasattr(self, 'compressibility'):
            self.set_compressibility(0)
        if not hasattr(self, 'is_homogeneous'):
            self.get_is_homogeneous()
    set_props = set_properties


    def __get_properties__(self):
        # self.get_shape() # > self.shape
        self.get_order(type='natural') # > self.order
        self.get_boundaries() # > self.boundaries, self.b
        self.get_pv_grid(show_boundary=True) # > pv_grid
        self.get_area() # self.area
        self.get_volume() # self.volume
        self.get_G() # > self.G
        self.get_is_homogeneous() # > self.is_homogeneous
        self.get_centers()
        self.get_flow_direction()


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
        '''
        
        '''
        if not i:
            self.tops = self.z = self.blocks * z
        else:
            self.tops[i] = z
    set_z = set_tops


    def get_is_homogeneous(self):
        if  all(self.dx[1:-1] == self.dx[0]) & \
            all(self.dy[1:-1] == self.dy[0]) & \
            all(self.dz[1:-1] == self.dz[0]) & \
            all(self.k[1:-1] == self.k[0]) & \
            all(self.porosity[1:-1] == self.porosity[0]):
            self.is_homogeneous = True # homogeneous
        else:
            self.is_homogeneous = False # heterogeneous
        return self.is_homogeneous


    def get_order(self, type='natural'):
        if type == 'natural':
            self.order = np.arange(self.nx_b) # natural order: 0, -1 are boundaries.
        else:
            raise ValueError("""Order type is not supported.\n
                            Supported order types: ['natural']""")
        return self.order
    

    def get_boundaries(self, i=None, from_pv_grid=False):
        '''
        Arguments:
            - i: index in x-direction.
        '''
        if from_pv_grid:
            pv_grid = self.get_pv_grid(show_boundary=True)
            shape = self.get_flowdir_shape()
            self.cells_ids = np.arange(pv_grid.n_cells).reshape(shape)
            if self.D == 1:
                self.boundaries = self.cells_ids[[0,-1]].flatten()
            elif self.D == 2:
                self.boundaries = np.sort(np.concatenate([
                    self.cells_ids[[0,-1],:].flatten(),
                    self.cells_ids[1:-1, [0,-1]].flatten()
                ]))
            elif self.D == 3:
                print(self.cells_ids)
                self.boundaries = np.sort(np.concatenate([
                    self.cells_ids[:,[0,-1],:].flatten(),
                    self.cells_ids[:, 1:-1, [0,-1]].flatten(),
                ]))
            if i == None:
                return self.boundaries
            else:
                print('not finished')
                
        if self.D == 1:
            if i == None:
                self.boundaries = self.b = [0, self.nx + 1]
                return self.boundaries
            else:
                if i < 0:
                    i = self.nx_b + i
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
        elif self.D == 2:
            print('not ready')
        elif self.D == 3:
            print('not ready')
    get_b = get_boundaries


    def get_neighbors(self, i, from_pv_grid=True):
        '''
        Arguments:
            - i: index in x-direction.
        '''
        if from_pv_grid:
                lst = self.pv_grid.neighbors(i, rel='connectivity') # topological, connectivity, geometric
                if self.D == 1:
                    if self.flowdir == 'x' and i == self.nx:
                        lst.remove(self.nx+1)
                    elif self.flowdir == 'y' and i == self.ny:
                        lst.remove(self.ny+1)
                    elif self.flowdir == 'z' and i == self.nz:
                        lst.remove(self.nz+1)
                return lst
        else:
            if i < 0:
                i = self.nx_b + i
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


    def get_G(self):
        '''
        Grid geometry factor. 
        '''
        if self.D == 1:
            if self.is_homogeneous:
                if 'x' in self.flowdir:
                    d = self.dx
                elif 'y' in self.flowdir:
                    d = self.dy
                elif 'z' in self.flowdir:
                    d = self.dz
                self.G = self.factors['transmissibility conversion'] * \
                    self.mean(self.k) * self.mean(self.area) / self.mean(d)
            else:
                self.G = 2 * self.factors['transmissibility conversion'] / \
                    (
                        (self.dx[:-1] / (self.area[:-1] * self.k[:-1])) + \
                        (self.dx[1:] / (self.area[1:] * self.k[1:]))
                    )
            return self.G
        else:
            return None

    
    def mean(self, property, type='geometric'):
        '''
        '''
        if self.is_homogeneous:
            return property[1:]
        else:
            if type == 'geometric':
                return (property[:-1] + property[1:])/2
            else:
                raise ValueError('Unknown mean type')

    
    def get_pv_grid(self, show_boundary=True, verbose=False):
        '''
        https://docs.pyvista.org/api/core/_autosummary/pyvista.ExplicitStructuredGrid.html
        '''
        if show_boundary:
            dims = np.array((self.nx_b,self.ny_b,self.nz_b))+1
        else:
            dims = np.array((self.nx,self.ny,self.nz))+1

        corners = self.get_corners(show_boundary)
        print('dims:',dims)
        print('corners.shape:',corners.shape)
        self.pv_grid = pv.ExplicitStructuredGrid(dims, corners)

        if verbose: print(self.pv_grid)
        return self.pv_grid


    def get_corners(self, show_boundary=True, verbose=True):
        '''
        
        '''
        # Cumulative sum: 
        # also: np.append([0],arr), np.concatenate([0],arr)
        if 'x' in self.flowdir:
            xcorn = np.insert(self.dx.cumsum(), 0, 0)
        else:
            xcorn = np.arange(0, (self.nx+1)*self.dx[0], self.dx[0])
        
        if self.D == 1:
            if self.flowdir == 'x':
                ycorn = np.arange(0, (self.ny+1)*self.dy[1], self.dy[1])
                zcorn = np.arange(0, (self.nz+1)*self.dz[1], self.dz[1])
            else:
                ycorn = np.arange(0, (self.ny+1)*self.dy[0], self.dy[0])
                zcorn = np.arange(0, (self.nz+1)*self.dz[0], self.dz[0])
        elif self.D == 2:
            ycorn = np.insert(self.dy.cumsum(), 0, 0)
            zcorn = np.arange(0, (self.nz+1)*self.dz[1], self.dz[1])
        elif self.D == 3:
            ycorn = np.insert(self.dy.cumsum(), 0, 0)
            zcorn = np.insert(self.dz.cumsum(), 0, 0)
        
        # Show boundary:
        if show_boundary:
            if 'x' in self.flowdir:
                ix = 2
            else:
                ix = 0
            if self.D == 1:
                if self.flowdir == 'y':
                    iy = 2
                else:
                    iy = 0
            else:
                iy = ix
        else:
            ix = 0
            iy = 0
            xcorn = xcorn[1:-1]
            if self.D >= 2:
                ycorn = ycorn[1:-1]
            
        # X corners:
        xcorn = np.repeat(xcorn, 2)
        xcorn = xcorn[1:-1]
        xcorn = np.tile(xcorn, 4*(self.ny+iy)*self.nz)
        
        # Y corners:
        ycorn = np.repeat(ycorn, 2)
        ycorn = ycorn[1:-1]
        ycorn = np.tile(ycorn, (2*(self.nx+ix), 2*(self.nz)))
        ycorn = np.transpose(ycorn)
        ycorn = ycorn.flatten()

        # Z corners:
        zcorn = np.repeat(zcorn, 2)
        zcorn = zcorn[1:-1]
        zcorn = np.repeat(zcorn, 4*(self.nx+ix)*(self.ny+iy))

        if verbose: print(xcorn.shape, ycorn.shape, zcorn.shape)

        # Combine corners:
        corners = np.stack((xcorn, ycorn, zcorn))
        corners = corners.transpose()
        
        return corners


    def get_flow_direction(self):
        if self.D == 1:
            flowdir_i = np.argmax(self.shape)
            if flowdir_i == 0:
                self.flowdir = 'x'
            elif flowdir_i == 1:
                self.flowdir = 'y'
            elif flowdir_i == 2:
                self.flowdir = 'z'
        elif self.D == 2:
            flowdir_i = np.argmin(self.shape)
            if flowdir_i == 2:
                self.flowdir = 'xy'
            elif flowdir_i == 1:
                self.flowdir = 'xz'
            elif flowdir_i == 0:
                self.flowdir = 'yz'
        else: # elif self.D == 3:
            self.flowdir = 'xyz'
        return self.flowdir
    

    def get_flowdir_shape(self):
        flowdir = self.get_flow_direction()
        if flowdir == 'x':
            self.flowdir_shape = (self.nx_b,)
        elif flowdir == 'y':
            self.flowdir_shape = (self.ny_b,)
        elif flowdir == 'z':
            self.flowdir_shape = (self.nz_b,)
        elif flowdir == 'xy':
            self.flowdir_shape = (self.ny_b, self.nx_b)
        elif flowdir == 'xz':
            self.flowdir_shape = (self.nz_b, self.nx_b, 3)
        elif flowdir == 'yz':
            self.flowdir_shape = (self.ny_b, self.nz_b, 3)
        elif flowdir == 'xyz':
            self.flowdir_shape = (self.nz_b, self.ny_b, self.nx_b)
        return self.flowdir_shape
    
    
    def get_centers(self, with_boundary=True):
        pv_grid = self.get_pv_grid(show_boundary=True)
        shape = self.get_flowdir_shape() + (3,)
        self.centers = pv_grid.cell_centers().points.reshape(shape)
        if with_boundary:
            return self.centers.flatten()
        else:
            if self.D == 1:
                return self.centers[1:-1,:].flatten()
            elif self.D == 2:
                return self.centers[1:-1,1:-1,:].flatten()
            elif self.D == 3:
                return self.centers[:,1:-1,1:-1,:].flatten()
    
    
    def get_cell_id(self, i,j,k):
        return self.pv_grid.cell_id((i,j,k))


    def get_cell_ijk(self, cell_id):
        return self.pv_grid.cell_coords(cell_id)

    
    def get_cells_ids(self, with_boundary=True):
        pv_grid = self.get_pv_grid(show_boundary=True)
        shape = self.get_flowdir_shape()
        self.cells_ids = np.arange(pv_grid.n_cells).reshape(shape)
        if with_boundary:
            return self.cells_ids.flatten()
        else:
            if self.D == 1:
                return self.cells_ids[1:-1].flatten()
            elif self.D == 2:
                return self.cells_ids[1:-1,1:-1].flatten()
            elif self.D == 3:
                return self.cells_ids[:, 1:-1,1:-1].flatten()
    
    '''Todo:
    def set_d(self, d, n):
        if isinstance(d, (list, tuple, np.ndarray)):
            if len(d) == n + 2:
                return self.blocks * d
            elif len(d) == n:
                d_array = np.insert(d, 0, d[0])
                d_array = np.append(d, d[-1])
                return d_array
        elif isinstance(d, (int, float)):
            return self.blocks * d
        else:
            pass
        
    def get_geometry(self, dx, dy, dz):        
        
        for c in range(self.pv_grid.n_cells):
            bounds = np.array(self.pv_grid.cell_bounds(c))
            dx, dy, dz = bounds[1::2] - bounds[::2]
            self.dx[c] = dx
            self.dy[c] = dy
            self.dz[c] = dz
    # todo: copy functionality
    # def copy(self):
    #     return self
    '''
    
    def show(self, show_lables=False, show_centers=True, show_boundary=True):
        pv_grid = self.get_pv_grid(show_boundary)
        pl = pv.Plotter()
        pl.add_mesh(
            pv_grid, 
            show_edges=True, 
            color="tan",
            opacity=0.7,
        )
        
        if show_lables:
            points = pv_grid.points
            mask = points[:,1] == 0
            pl.add_point_labels(
                points[mask], 
                points[mask].tolist(), 
                point_size=10, 
                font_size=10,
            )
        
        if show_centers:
            labels = self.get_cells_ids(show_boundary)
            centers = self.get_centers(show_boundary)
            pl.add_point_labels(
                points=centers,
                labels=labels,
                point_size=10,
                font_size=10,
            )
        
        pl.add_camera_orientation_widget()
        pl.enable_fly_to_right_click()
        pl.show_axes()
        pl.camera_position = 'xy'
        pl.set_background('black', top='gray')
        pl.show(full_screen=True)
        
    
#%%
if __name__ == '__main__':
    # Canvas:
    dx = 100 # np.array([20,40,60,10,70,80,10])
    dy = 100 # np.array([40,60,10,70])
    dz = 100 # np.array([20,40])
    grid = CartGrid(nx=3, ny=1, nz=2, dx=dx, dy=dy, dz=dz, phi=0.27, k=270)
    
    # indecies:
    # print(grid.cell_id(2,1,0))
    
    # cells:
    print('cell_id(1,1,1):', grid.get_cell_id(1,1,1))
    print('cell_ijk(36):', grid.get_cell_ijk(36))
    # print(grid.get_cells_ids(with_boundary=True))
    # print(grid.get_cells_ids())
    print(grid.get_boundaries(from_pv_grid=True))
    
    # Neighbors:
    # print('Neighbors:')
    # print('1:', grid.get_neighbors(1))
    # print(grid.get_neighbors(1, from_pv_grid=True)) # connectivity, geometric, topological
    # print('2:', grid.get_neighbors(2))
    # print(grid.get_neighbors(2, from_pv_grid=True))
    # print('5:', grid.get_neighbors(5))
    # print(grid.get_neighbors(5, from_pv_grid=True))
    
    
    # Boundaries: pyvista is not suitable
    # print('Boundaries:')
    # print('1:', grid.get_boundaries(1))
    # print(grid.pv_grid.neighbors(1)) # connectivity, geometric, topological
    # print('2:', grid.get_boundaries(2))
    # print(grid.pv_grid.neighbors(2))
    # print('5:', grid.get_boundaries(5))
    # print(grid.pv_grid.neighbors(5))
    
    # Doc and reporting:
    # print(grid.__doc__)
    # print(grid) # or grid.report()
    
    # Setters:
    # grid.set_comp(1e-6)
    # grid.set_permeability(200, 1)
    # grid.set_props(0.3, 11)
    # grid.set_phi(0.2)
    # grid.set_k(10)
    # grid.set_tops(10)

    # Getters:
    # grid.get_pv_grid(show_boundary=True)
    # print(grid.get_boundaries())
    # print(grid.pv_grid.neighbors(1, rel='connectivity')) # connectivity, geometric
    
    # Show:
    grid.show(show_boundary=True)
    # grid.show(show_boundary=False)