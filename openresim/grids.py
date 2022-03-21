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
        self.get_shape() # > self.shape
        

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
        self.get_flow_direction() # > self.flowdir
        
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
            self.nz_b = self.nz # + 2 # new
            flowshape = self.nz_b
        elif self.flowdir == 'xy':
            self.nx_b = self.nx + 2
            self.ny_b = self.ny + 2
            self.nz_b = self.nz
            flowshape = self.nx_b
        elif self.flowdir == 'xz':
            self.nx_b = self.nx + 2
            self.ny_b = self.ny
            self.nz_b = self.nz # + 2 # new
            flowshape = self.nx_b
        elif self.flowdir == 'yz':
            self.nx_b = self.nx
            self.ny_b = self.ny + 2
            self.nz_b = self.nz # + 2 # new
            flowshape = self.nz_b
        elif self.flowdir == 'xz+':
            self.nx_b = self.nx + 2
            self.ny_b = self.ny + 2
            self.nz_b = self.nz # + 2 # new
            flowshape = self.nx_b
        elif self.flowdir == 'yz+':
            self.nx_b = self.nx + 2
            self.ny_b = self.ny + 2
            self.nz_b = self.nz # + 2 # new
            flowshape = self.nz_b
        elif self.flowdir == 'xyz':
            self.nx_b = self.nx + 2
            self.ny_b = self.ny + 2
            self.nz_b = self.nz # + 2 # new
            flowshape = self.nx_b
        else:
            raise ValueError('Geometry higher than 3D is not allowed using CartGrid class.')
        self.get_flow_shape() # > self.flow_shape
        
        self.blocks = np.ones(flowshape, dtype='int')
        self.i_blocks = self.blocks.copy()
        self.i_blocks[[0, -1]] = 0
        self.b_blocks = np.zeros(flowshape, dtype='int') #ss.lil_matrix(self.shape, dtype='int')
        self.b_blocks[[0, -1]] = 1
        
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
            elif self.flowdir == 'yz':
                self.dz = np.ones(self.nz_b, dtype='int') * dz               
        elif self.D == 3:  
            self.dx = np.ones(self.nx_b, dtype='int') * dx
            self.dy = np.ones(self.ny_b, dtype='int') * dy
            self.dz = np.ones(self.nz_b, dtype='int') * dz
        else:
            raise ValueError('Geometry higher than 3D is not allowed using CartGrid class.')
            
        self.pv_grid_b = self.get_pv_grid(boundary=True)
        self.pv_grid = self.get_pv_grid(boundary=False)
        
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
        # self.get_dimension() # > self.dimension, self.D
        # self.get_shape() # > self.shape
        self.get_order(type='natural') # > self.order
        self.get_boundaries() # > self.boundaries, self.b
        # self.get_pv_grid(boundary=True) # > pv_grid
        self.get_area() # self.area
        self.get_volume() # self.volume
        self.get_G() # > self.G
        self.get_is_homogeneous() # > self.is_homogeneous
        self.get_centers()
        # self.get_flow_direction()


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
            self.cells_ids = np.arange(self.pv_grid_b.n_cells).reshape(self.flow_shape)
            if self.D == 1:
                self.boundaries = self.cells_ids[[0,-1]].flatten()
            elif self.D == 2:
                self.boundaries = np.sort(np.concatenate([
                    self.cells_ids[[0,-1],:].flatten(),
                    self.cells_ids[1:-1, [0,-1]].flatten()
                ]))
            elif self.D == 3:
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

    
    def get_pv_grid(self, boundary=True, verbose=False):
        '''
        https://docs.pyvista.org/api/core/_autosummary/pyvista.ExplicitStructuredGrid.html
        '''
        if boundary:
            dims = np.array((self.nx_b,self.ny_b,self.nz_b))+1
        else:
            dims = np.array((self.nx,self.ny,self.nz))+1

        corners = self.get_corners(boundary, verbose)
        pv_grid = pv.ExplicitStructuredGrid(dims, corners)

        s = 'with boundary' if boundary else 'without boundary'
        print(f'- pv_grid {s} was created.')
        if verbose: print(pv_grid)
        return pv_grid


    def get_corners(self, boundary=True, verbose=True):
        '''
        
        '''

        if self.D == 1:
            if self.flowdir == 'x':
                xcorn = np.insert(self.dx.cumsum(), 0, 0)
                ycorn = np.arange(0, (self.ny+1)*self.dy[1], self.dy[1])
                zcorn = np.arange(0, (self.nz+1)*self.dz[1], self.dz[1])
            if self.flowdir == 'y':
                xcorn = np.arange(0, (self.nx+1)*self.dx[0], self.dx[0])
                ycorn = np.insert(self.dy.cumsum(), 0, 0)
                zcorn = np.arange(0, (self.nz+1)*self.dz[0], self.dz[0])
            if self.flowdir == 'z':
                xcorn = np.arange(0, (self.nx+1)*self.dx[0], self.dx[0])
                ycorn = np.arange(0, (self.ny+1)*self.dy[1], self.dy[1])
                zcorn = np.insert(self.dz.cumsum(), 0, 0)
        elif self.D == 2:
            xcorn = np.insert(self.dx.cumsum(), 0, 0)
            ycorn = np.insert(self.dy.cumsum(), 0, 0)
            zcorn = np.arange(0, (self.nz+1)*self.dz[1], self.dz[1])
        elif self.D == 3:
            xcorn = np.insert(self.dx.cumsum(), 0, 0)
            ycorn = np.insert(self.dy.cumsum(), 0, 0)
            zcorn = np.insert(self.dz.cumsum(), 0, 0)
            
        # Show boundary:
        if boundary:
            if 'x' in self.flowdir: 
                ix = 2
            else:
                ix = 0
                
            if 'y' in self.flowdir:
                iy = 2
                # if self.D > 1:
                #     ix = 2
            else:
                iy = 0
                
            # if self.D > 1 and '+' in self.flowdir:
            #     iy = ix 
        else:
            ix = 0
            iy = 0
           
            if 'x' in self.flowdir:
                xcorn = xcorn[1:-1]
                  
            if 'y' in self.flowdir:
                ycorn = ycorn[1:-1]
            # if self.D > 1 and '+' in self.flowdir:
            #     xcorn = xcorn[1:-1]
            #     ycorn = ycorn[1:-1]
            
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
            flowdir_id = np.argmax(self.shape)
            if flowdir_id == 0:
                self.flowdir = 'x'
            elif flowdir_id == 1:
                self.flowdir = 'y'
            elif flowdir_id == 2:
                self.flowdir = 'z'
        elif self.D == 2:
            flowdir_id = np.argmin(self.shape)
            if flowdir_id == 2:
                self.flowdir = 'xy'
            elif flowdir_id == 1:
                self.flowdir = 'xz'
            elif flowdir_id == 0:
                self.flowdir = 'yz'
        else: # elif self.D == 3:
            self.flowdir = 'xyz'
        print(f'- flowdir is set to {self.flowdir}.')
        return self.flowdir
    

    def get_flow_shape(self):
        if self.flowdir == 'x':
            self.flow_shape = (self.nx_b,)
        elif self.flowdir == 'y':
            self.flow_shape = (self.ny_b,)
        elif self.flowdir == 'z':
            self.flow_shape = (self.nz_b,)
        elif self.flowdir == 'xy':
            self.flow_shape = (self.ny_b, self.nx_b)
        elif self.flowdir == 'xz':
            self.flow_shape = (self.nz_b, self.nx_b)
        elif self.flowdir == 'yz':
            self.flow_shape = (self.nz_b, self.ny_b)
        elif self.flowdir == 'xz+':
            self.flow_shape = (self.nz_b, self.nx_b, 3)
        elif self.flowdir == 'yz+':
            self.flow_shape = (self.nz_b, self.ny_b, 3)
        elif self.flowdir == 'xyz':
            self.flow_shape = (self.nz_b, self.ny_b, self.nx_b)
        print(f'- flow_shape is set to {self.flow_shape}')
        return self.flow_shape
    
    
    def get_centers(self, boundary=True):
        flow_shape = self.flow_shape + (3,)
        self.centers = self.pv_grid_b.cell_centers().points.reshape(flow_shape)
        if boundary:
            return self.centers.flatten()
        else:
            if self.D == 1 and self.flowdir != 'z':
                return self.centers[1:-1,:].flatten()
            elif self.D == 1 and self.flowdir == 'z':
                return self.centers.flatten()           
            elif self.D == 2 and 'z' not in self.flowdir:
                return self.centers[1:-1,1:-1,:].flatten()
            elif self.D == 2 and self.flowdir == 'xz':
                return self.centers[:,1:-1].flatten()
            elif self.D == 2 and self.flowdir == 'yz':
                return self.centers[:,1:-1].flatten()
            elif self.D == 2 and self.flowdir == 'xz+':
                return self.centers[:,2:-2,:].flatten()
            elif self.D == 2 and self.flowdir == 'yz+':
                return self.centers[:,1:-1,1:-1].flatten()
            elif self.D == 3:
                return self.centers[:,1:-1,1:-1,:].flatten()

    
    def get_cells_ids(self, boundary=True):     
        self.cells_ids = np.arange(self.pv_grid_b.n_cells).reshape(self.flow_shape)
        if boundary:
            return self.cells_ids.flatten()
        else:
            if self.D == 1 and self.flowdir != 'z':
                return self.cells_ids[1:-1].flatten()
            elif self.D == 1 and self.flowdir == 'z':
                return self.cells_ids.flatten()
            elif self.D == 2 and 'z' not in self.flowdir:
                return self.cells_ids[1:-1,1:-1].flatten()
            elif self.D == 2 and self.flowdir == 'xz':
                return self.cells_ids[:,1:-1].flatten()
            elif self.D == 2 and self.flowdir == 'yz':
                return self.cells_ids[:,1:-1].flatten()
            elif self.D == 2 and self.flowdir == 'xz+':
                return self.cells_ids[:,2:-2,:].flatten()
            elif self.D == 2 and self.flowdir == 'yz+':
                return self.cells_ids[:,1:-1,1:-1].flatten()
            elif self.D == 3:
                return self.cells_ids[:, 1:-1,1:-1].flatten()
    
    
    def get_cell_id(self, i,j,k):
        return self.pv_grid_b.cell_id((i,j,k))


    def get_cell_ijk(self, cell_id):
        return self.pv_grid_b.cell_coords(cell_id)
    
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
    
    def show(self, show_lables=False, show_centers=True, boundary=True):
        if boundary:
            pv_grid = self.pv_grid_b
        else:
            pv_grid = self.pv_grid
            
        pl = pv.Plotter()
        pl.add_mesh(
            pv_grid, 
            show_edges=True, 
            color="white",
            opacity=0.8,
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
            labels = self.get_cells_ids(boundary)
            centers = self.get_centers(boundary)
            pl.add_point_labels(
                points=centers,
                labels=labels,
                point_size=10,
                font_size=10,
            )
        
        s = 'with boundary' if boundary else 'without boundary'
        title = f'OpenReSim: {self.D}D model (flow at {self.flowdir}-direction {s})'
        pl.add_title(title, 
                     font='courier', 
                     color='white',
                     font_size=8)
        pl.add_camera_orientation_widget()
        pl.enable_fly_to_right_click()
        pl.show_axes()
        pl.camera_position = 'xz'
        pl.set_background('black', top='gray')
        pl.show(title=title, window_size=[1000,1000], full_screen=True)
        
    
#%%
if __name__ == '__main__':
    # Canvas:
    dx = 10 # np.array([60,10,70,80,10])
    dy = 10 # np.array([40,40,60,10,70])
    dz = 10 # np.array([20,40, 60])
    
    # grid = CartGrid(nx=2, ny=1, nz=1, dx=dx, dy=dy, dz=dz, phi=0.27, k=270)
    # grid.show(boundary=True)
    # grid.show(boundary=False)
    # grid = CartGrid(nx=1, ny=2, nz=1, dx=dx, dy=dy, dz=dz, phi=0.27, k=270)
    # grid.show(boundary=True)
    # grid.show(boundary=False)
    # grid = CartGrid(nx=1, ny=1, nz=2, dx=dx, dy=dy, dz=dz, phi=0.27, k=270)
    # grid.show(boundary=True)
    # grid.show(boundary=False)
    # grid = CartGrid(nx=3, ny=2, nz=1, dx=dx, dy=dy, dz=dz, phi=0.27, k=270)
    # grid.show(boundary=True)
    # grid.show(boundary=False)
    # grid = CartGrid(nx=1, ny=2, nz=2, dx=dx, dy=dy, dz=dz, phi=0.27, k=270)
    # grid.show(boundary=True)
    # grid.show(boundary=False)
    # grid = CartGrid(nx=2, ny=1, nz=2, dx=dx, dy=dy, dz=dz, phi=0.27, k=270)
    # grid.show(boundary=True)
    # grid.show(boundary=False)
    # grid = CartGrid(nx=2, ny=2, nz=2, dx=dx, dy=dy, dz=dz, phi=0.27, k=270)
    # grid.show(boundary=True)
    # grid.show(boundary=False)
    
    grid = CartGrid(nx=3, ny=1, nz=1, dx=dx, dy=dy, dz=dz, phi=0.27, k=270)
    grid.show(boundary=True)
    grid.show(boundary=False)
    grid = CartGrid(nx=1, ny=3, nz=1, dx=dx, dy=dy, dz=dz, phi=0.27, k=270)
    grid.show(boundary=True)
    grid.show(boundary=False)
    grid = CartGrid(nx=1, ny=1, nz=3, dx=dx, dy=dy, dz=dz, phi=0.27, k=270)
    grid.show(boundary=True)
    grid.show(boundary=False)
    grid = CartGrid(nx=3, ny=3, nz=1, dx=dx, dy=dy, dz=dz, phi=0.27, k=270)
    grid.show(boundary=True)
    grid.show(boundary=False)
    grid = CartGrid(nx=1, ny=3, nz=3, dx=dx, dy=dy, dz=dz, phi=0.27, k=270)
    grid.show(boundary=True)
    grid.show(boundary=False)
    grid = CartGrid(nx=3, ny=1, nz=3, dx=dx, dy=dy, dz=dz, phi=0.27, k=270)
    grid.show(boundary=True)
    grid.show(boundary=False)
    grid = CartGrid(nx=3, ny=2, nz=3, dx=dx, dy=dy, dz=dz, phi=0.27, k=270)
    grid.show(boundary=True)
    grid.show(boundary=False)
    
    # indecies:
    # print(grid.cell_id(2,1,0))
    
    # cells:
    # print('cell_id(1,1,1):', grid.get_cell_id(1,1,1))
    # print('cell_ijk(36):', grid.get_cell_ijk(36))
    # print(grid.get_cells_ids(with_boundary=True))
    # print(grid.get_cells_ids())
    # print('- boundaries:', grid.get_boundaries(from_pv_grid=True))
    
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