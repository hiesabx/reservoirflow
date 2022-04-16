#%% Import Statements:
from openresim.base import Base
import numpy as np
import pyvista as pv
from functools import lru_cache


#%% Grid Class: 
class Grid(Base):
    '''Grid class to create a grid using numpy arrays and pyvista grids.
    '''
    def __init__(self, dtype, unit):
        super().__init__(unit)
        self.dtype = dtype # np.single, np.double


#%% CartGrid Class: 
class CartGrid(Grid):
    '''CartGrid class to create a explicit structured grid.
    
    Note that the following conventions are used: 
        1. rows for dx. 1D grids are stored as 1 column with multiple rows.  
        2. columns for dy. 2D grids are stored as a table where rows refer to x-direction while columns refere to y-direction. 
        3. layers for dz. 3D grids are stored as cubes where rows refer to x-direction, columns to y-direction, and layers for z-direction.
    
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
        super().__init__(dtype, unit)
        self.type = 'cartesian'
        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.get_dimension() # > self.D
        self.get_shape() # > self.shape
        self.get_flow_dir() # > self.flow_dir
        self.get_n_b() # > self.nx_b, self.ny_b, self.n_zb
        self.get_flow_shape() # > self.flow_shape
        self.get_n_cells() # > self.n_cells_b, self.n_cells
        self.get_order(type='natural') # > self.order_b, self.order
        self.max_nb = max(self.nx_b, self.ny_b, self.nz_b)
        self.get_cells_dims(dx, dy, dz)
        self.get_cells_area()
        self.get_cells_volume_()
        
        self.pyvista_grid_b = self.get_pyvista_grid(boundary=True)
        self.pyvista_grid = self.get_pyvista_grid(boundary=False)
        self.get_cells_coords() # > cells_coords_b, cells_coords
        
        self.set_properties(phi, k, z, comp)
        self.get_boundaries() # > self.boundaries_id, self.boundaries_coords
        self.get_volume() # self.volume
        self.get_G() # > self.G
        self.get_is_homogeneous() # > self.is_homogeneous
        self.get_cells_center()


    def get_dimension(self):
        if not hasattr(self, 'D'):
            self.D = sum([1 if n>1 else 0 for n in (self.nx, self.ny, self.nz)])
            print(f'- D (dimension) is set to {self.D}.')
        return self.D
    
    
    def get_shape(self):
        if not hasattr(self, 'shape'):
            self.shape = np.array((self.nx, self.ny, self.nz), dtype='int')
            print(f'- shape is set to {self.shape}.')
        return self.shape
    

    def get_flow_dir(self):
        if not hasattr(self, 'flow_dir'):
            if self.D == 0:
                self.flow_dir = '-'
            elif self.D == 1:
                flowdir_id = np.argmax(self.shape)
                if flowdir_id == 0:
                    self.flow_dir = 'x'
                elif flowdir_id == 1:
                    self.flow_dir = 'y'
                elif flowdir_id == 2:
                    self.flow_dir = 'z'
            elif self.D == 2:
                flowdir_id = np.argmin(self.shape)
                if flowdir_id == 2:
                    self.flow_dir = 'xy'
                elif flowdir_id == 1:
                    self.flow_dir = 'xz'
                elif flowdir_id == 0:
                    self.flow_dir = 'yz'
            elif self.D == 3:
                self.flow_dir = 'xyz'
            print(f'- flow_dir is set to {self.flow_dir}.')
        return self.flow_dir
    
    
    def get_n_b(self):
        if 'x' in self.flow_dir:
            self.nx_b = self.nx + 2
        else:
            self.nx_b = self.nx
        if 'y' in self.flow_dir:
            self.ny_b = self.ny + 2
        else:
            self.ny_b = self.ny
        if 'z' in self.flow_dir:
            self.nz_b = self.nz + 2
        else:
            self.nz_b = self.nz
        return (self.nx_b, self.ny_b, self.nz_b)


    def get_flow_shape(self):
        if self.flow_dir == '-':
            self.flow_shape = (1,)
        elif self.flow_dir == 'x':
            self.flow_shape = (self.nx_b,)
        elif self.flow_dir == 'y':
            self.flow_shape = (self.ny_b,)
        elif self.flow_dir == 'z':
            self.flow_shape = (self.nz_b,)
        elif self.flow_dir == 'xy':
            self.flow_shape = (self.ny_b, self.nx_b)
        elif self.flow_dir == 'xz':
            self.flow_shape = (self.nz_b, self.nx_b)
        elif self.flow_dir == 'yz':
            self.flow_shape = (self.nz_b, self.ny_b)
        # elif self.flowdir == 'xz+':
        #     self.flow_shape = (self.nz_b, self.nx_b, 3)
        # elif self.flowdir == 'yz+':
        #     self.flow_shape = (self.nz_b, self.ny_b, 3)
        elif self.flow_dir == 'xyz':
            self.flow_shape = (self.nz_b, self.ny_b, self.nx_b)
        print(f'- flow_shape is set to {self.flow_shape}')
        return self.flow_shape


    def get_n_cells(self, boundary=True):
        if not hasattr(self, 'n_cells_b'):
            self.n_cells_b = self.nx_b * self.ny_b * self.nz_b
        if not hasattr(self, 'n_cells'):
            self.n_cells = self.nx * self.ny * self.nz
        if boundary:
            return self.n_cells_b
        else:
            return self.n_cells
    
    
    def get_order(self, type='natural', boundary=True):
        if type == 'natural':
            if not hasattr(self, 'order_b'):
                self.order_b = np.arange(self.n_cells_b) # .reshape(self.flow_shape)
            if not hasattr(self, 'order'):
                self.order = np.arange(self.n_cells) # .reshape(self.flow_shape)
        else:
            raise ValueError("""Order type is not supported.\n
                            Supported order types: ['natural']""")
        if boundary:
            return self.order_b
        else:
            return self.order
    
    
    def get_cells_dims(self, dx, dy, dz):
        self.blocks = np.ones(self.max_nb, dtype='int')
        self.i_blocks = self.blocks.copy()
        self.i_blocks[[0, -1]] = 0
        self.b_blocks = np.zeros(self.max_nb, dtype='int') #ss.lil_matrix(self.shape, dtype='int')
        self.b_blocks[[0, -1]] = 1
        
        m_list = []
        if 'x' in self.flow_dir:
            self.dxx = np.ones(self.nx_b, dtype='int') * dx
            m_list.append(self.dxx)
        else:
            self.dxx = np.ones(self.max_nb, dtype='int') * dx
            m_list.append(dx)
        if 'y' in self.flow_dir:
            self.dyy = np.ones(self.ny_b, dtype='int') * dy
            m_list.append(self.dyy)
        else:
            self.dyy = np.ones(self.max_nb, dtype='int') * dy
            m_list.append(dy)
        if 'z' in self.flow_dir:
            self.dzz = np.ones(self.nz_b, dtype='int') * dz
            m_list.append(self.dzz)
        else:
            self.dzz = np.ones(self.max_nb, dtype='int') * dz
            m_list.append(dz)
        
        self.dx, self.dy, self.dz = np.meshgrid(*m_list, copy=False)
        self.dx = np.transpose(self.dx, axes=(0,2,1))
        self.dy = np.transpose(self.dy, axes=(2,0,1))
        self.dz = np.transpose(self.dz, axes=(2,1,0))
        
        
    def get_cells_area(self, dir=None, boundary=True):
        if not hasattr(self, 'area_x'):
            self.area_x_b = (self.dy.flatten() * self.dz.flatten()).reshape(self.flow_shape)
            self.area_x = self.remove_boundaries(self.area_x_b)
        if not hasattr(self, 'area_y'):
            self.area_y_b = (self.dx.flatten() * self.dz.flatten()).reshape(self.flow_shape)
            self.area_y = self.remove_boundaries(self.area_y_b)
        if not hasattr(self, 'area_z'):
            self.area_z_b = (self.dx.flatten() * self.dy.flatten()).reshape(self.flow_shape)
            self.area_z = self.remove_boundaries(self.area_z_b)
        if not hasattr(self, 'area'):
            self.area_b = {'x': self.area_x_b, 'y': self.area_y_b, 'z': self.area_z_b}
            self.area = {'x': self.area_x, 'y': self.area_y, 'z': self.area_z}
        if dir == 'x':
            if boundary:
                return self.area_x_b
            else:
                return self.area_x
        elif dir == 'y':
            if boundary:
                return self.area_y_b
            else:
                return self.area_y
        elif dir == 'z':
            if boundary:
                return self.area_z_b
            else:
                return self.area_z
        else:
            if boundary:
                return self.area_b
            else:
                return self.area
    
    
    def get_cells_area_x(self, boundary=True):
        return self.get_cells_area(dir='x', boundary=boundary)
    def get_cells_area_y(self, boundary=True):
        return self.get_cells_area(dir='y', boundary=boundary)
    def get_cells_area_z(self, boundary=True):
        return self.get_cells_area(dir='z', boundary=boundary)
    
    
    def get_cells_volume_(self, boundary=True):
        if not hasattr(self, 'cells_volume_b'):
            self.cells_volume_b = (self.dx.flatten() * 
                                   self.dy.flatten() * 
                                   self.dz.flatten()).reshape(self.flow_shape)
        if not hasattr(self, 'cells_volume'):
            self.cells_volume = self.remove_boundaries(self.cells_volume_b)
        if boundary:
            return self.cells_volume_b
        else:
            return self.cells_volume
        
    
    def fix_boundary(self, d):
        if isinstance(d, [list, tuple, np.array]):
            if len(d) == self.max_nb - 2:
                d = np.concatenate((d[0], d, d[-1]))
        return d
    
    
    def set_properties(self, phi=None, k=None, z=None, comp=None, id=None):
        """
        Set properties
        """
        # Set user defined properties:
        if phi is not None:
            self.set_porosity(phi, id)
        if k is not None:
            self.set_permeability(k, id)
        if z is not None:
            self.set_tops(z, id)
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
    

    def set_porosity(self, phi, id=None):
        """
        
        """
        if not id:
            self.porosity = self.phi = self.blocks * phi # np.zeros(self.shape) + phi
        else:
            self.porosity[id] = phi
            self.get_is_homogeneous()
    set_phi = set_porosity


    def set_permeability(self, k, id=None):
        """
        
        """
        if not id:
            self.permeability = self.k = self.blocks * k # np.zeros(self.shape) + k
        else:
            self.permeability[id] = k
            self.get_is_homogeneous()
    set_k = set_permeability

    
    def set_tops(self, z=None, id=None):
        '''
        
        '''
        if not id:
            self.tops = self.z = self.blocks * z
        else:
            self.tops[id] = z
    set_z = set_tops


    def get_is_homogeneous(self):
        if  all(self.dxx[1:-1] == self.dxx[0]) & \
            all(self.dyy[1:-1] == self.dyy[0]) & \
            all(self.dzz[1:-1] == self.dzz[0]) & \
            all(self.k[1:-1] == self.k[0]) & \
            all(self.porosity[1:-1] == self.porosity[0]):
            self.is_homogeneous = True # homogeneous
        else:
            self.is_homogeneous = False # heterogeneous
        return self.is_homogeneous
    

    def get_G(self):
        '''
        Grid geometry factor. 
        '''
        if self.D == 1:
            if self.is_homogeneous:
                if 'x' in self.flow_dir:
                    d = self.dxx
                elif 'y' in self.flow_dir:
                    d = self.dyy
                elif 'z' in self.flow_dir:
                    d = self.dzz
                self.G = self.factors['transmissibility conversion'] * \
                    self.get_mean(self.k) * self.get_mean(self.area) / self.get_mean(d)
            else:
                self.G = 2 * self.factors['transmissibility conversion'] / \
                    (
                        (self.dxx[:-1] / (self.area[:-1] * self.k[:-1])) + \
                        (self.dxx[1:] / (self.area[1:] * self.k[1:]))
                    )
            return self.G
        else:
            return None

    
    def get_mean(self, property, type='geometric'):
        '''
        '''
        if self.is_homogeneous:
            return property[1:]
        else:
            if type == 'geometric':
                return (property[:-1] + property[1:])/2
            else:
                raise ValueError('Unknown mean type')

    
    def get_pyvista_grid(self, boundary=True, verbose=False):
        '''
        https://docs.pyvista.org/api/core/_autosummary/pyvista.ExplicitStructuredGrid.html
        '''
        if boundary:
            dims = np.array((self.nx_b,self.ny_b,self.nz_b)) + 1
        else:
            dims = np.array((self.nx,self.ny,self.nz)) + 1

        corners = self.get_corners(boundary, verbose)
        pyvista_grid = pv.ExplicitStructuredGrid(dims, corners)

        s = 'with boundary' if boundary else 'without boundary'
        print(f'- pv_grid {s} was created.')
        if verbose: print(pyvista_grid)
        return pyvista_grid


    def get_corners(self, boundary=True, verbose=True):
        '''
        
        (Reference: https://docs.pyvista.org/examples/00-load/create-explicit-structured-grid.html)
        '''

        if 'x' in self.flow_dir:
            xcorn = np.insert(self.dxx.cumsum(), 0, 0)
        else:
            xcorn = np.arange(0, (self.nx+1)*self.dxx[0], self.dxx[0])
            
        if 'y' in self.flow_dir:
            ycorn = np.insert(self.dyy.cumsum(), 0, 0)
        else:
            ycorn = np.arange(0, (self.ny+1)*self.dyy[0], self.dyy[0])
            
        if 'z' in self.flow_dir:
            zcorn = np.insert(self.dzz.cumsum(), 0, 0)
        else:
            zcorn = np.arange(0, (self.nz+1)*self.dzz[0], self.dzz[0])
            
        # Show boundary:
        if boundary:
            ix = 2 if 'x' in self.flow_dir else 0
            iy = 2 if 'y' in self.flow_dir else 0
            iz = 2 if 'z' in self.flow_dir else 0
        else:
            ix = 0
            iy = 0
            iz = 0
            if 'x' in self.flow_dir:
                xcorn = xcorn[1:-1]
            if 'y' in self.flow_dir:
                ycorn = ycorn[1:-1]
            if 'z' in self.flow_dir:
                zcorn = zcorn[1:-1]

        # X corners:
        xcorn = np.repeat(xcorn, 2)
        xcorn = xcorn[1:-1]
        xcorn = np.tile(xcorn, 4*(self.ny+iy)*(self.nz+iz))
        
        # Y corners:
        ycorn = np.repeat(ycorn, 2)
        ycorn = ycorn[1:-1]
        ycorn = np.tile(ycorn, (2*(self.nx+ix), 2*(self.nz+iz)))
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


    def get_cells_center(self, boundary=True):
        if not hasattr(self, 'cells_center_b'):
            flow_shape = self.flow_shape + (3,)
            self.cells_center_b = self.pyvista_grid_b.cell_centers().points.reshape(flow_shape)
        if not hasattr(self, 'cells_center'):
            self.cells_center = self.remove_boundaries(self.cells_center_b)
        if boundary:
            return self.cells_center_b
        else:
            return self.cells_center


    def get_volume(self, boundary=True):
        if not hasattr(self, 'volume_b'):
            self.volume_b = self.pyvista_grid_b.volume
        if not hasattr(self, 'volume'):
            self.volume = self.pyvista_grid.volume
        if boundary:
            return self.pyvista_grid_b.volume
        else:
            return self.pyvista_grid.volume
        
        
    def get_cells_d(self, dir, boundary=True):
        if dir == 'x':
            cells_d = self.dx
        elif dir == 'y':
            cells_d = self.dy
        elif dir == 'z':
            cells_d = self.dz
        if boundary:
            return cells_d
        else:
            return self.remove_boundaries(cells_d)
    
    
    def get_cells_dx(self, boundary=True):
        return self.get_cells_d(dir='x', boundary=boundary)
    def get_cells_dy(self, boundary=True):
        return self.get_cells_d(dir='y', boundary=boundary)
    def get_cells_dz(self, boundary=True):
        return self.get_cells_d(dir='z', boundary=boundary)
    
    
    def get_cell_d(self, dir, id=None, coords=None):
        cells_d = self.get_cells_d(dir=dir, boundary=True)
        if id is not None:
            return cells_d.flatten()[id]
        elif coords is not None:
            return cells_d[coords[2],coords[1],coords[0]]
        else:
            raise ValueError('at least id or coords argument must be defined.')
        
        
    def get_cell_dx(self, id=None, coords=None):
        return self.get_cell_d(dir='x', id=id, coords=coords)
    def get_cell_dy(self, id=None, coords=None):
        return self.get_cell_d(dir='y', id=id, coords=coords)
    def get_cell_dz(self, id=None, coords=None):
        return self.get_cell_d(dir='z', id=id, coords=coords)
        
        
    def get_cell_area(self, dir, id=None, coords=None):
        cells_area = self.get_cells_area(dir=dir, boundary=True)  
        if id is not None:
            return cells_area.flatten()[id]
        elif coords is not None:
            return cells_area[coords[2],coords[1],coords[0]]
        else:
            raise ValueError('at least id or coords argument must be defined.')
        

    def get_cell_area_x(self, id=None, coords=None):
        return self.get_cell_area(dir='x', id=id, coords=coords)
    def get_cell_area_y(self, id=None, coords=None):
        return self.get_cell_area(dir='y', id=id, coords=coords)
    def get_cell_area_z(self, id=None, coords=None):
        return self.get_cell_area(dir='z', id=id, coords=coords)
    
        
    def get_cell_volume(self, id=None, coords=None):   
        cells_volume = self.get_cells_volume_(boundary=True)  
        if id is not None:
            return cells_volume.flatten()[id]
        elif coords is not None:
            return cells_volume[coords[2],coords[1],coords[0]]
        else:
            raise ValueError('at least id or coords argument must be defined.')
            
            
    def get_cells_volume(self, boundary=True, round_int=2):
        if not hasattr(self, 'cells_volume_b'):
            self.cells_volume_b = self.pyvista_grid_b.compute_cell_sizes()['Volume']
            self.cells_volume_b = self.cells_volume_b.round(round_int).reshape(self.flow_shape)
        if not hasattr(self, 'cells_volume'):
            self.cells_volume = self.remove_boundaries(self.cells_volume_b)
        if boundary:
            return self.cells_volume_b
        else:
            return self.cells_volume
    
    
    # @lru_cache(maxsize=None)
    def get_cell_id(self, coords, boundary=True):
        if boundary:
            return self.pyvista_grid_b.cell_id(coords)
        else:
            return self.pyvista_grid.cell_id(coords)


    def isin(self, arr, coords):
        for x in arr:
            if tuple(x) == coords:
                return True
        return False
    
    
    def remove_boundaries(self, in_array):
        if self.D == 0:
            return in_array
        elif self.D == 1:
            return in_array[1:-1]
        elif self.D == 2:
            return in_array[1:-1,1:-1]
        elif self.D == 3:
            return in_array[1:-1,1:-1,1:-1]
        # if self.D == 0:
        #     self.cells_center = self.cells_center_b
        # elif self.D == 1: # and self.flowdir != 'z':
        #     self.cells_center = self.cells_center_b[1:-1,:]    
        # elif self.D == 2 and 'z' not in self.flow_dir:
        #     self.cells_center = self.cells_center_b[1:-1,1:-1,:]
        # elif self.D == 2 and self.flow_dir == 'xz':
        #     self.cells_center = self.cells_center_b[1:-1,1:-1]
        # elif self.D == 2 and self.flow_dir == 'yz':
        #     self.cells_center = self.cells_center_b[1:-1,1:-1]
        # elif self.D == 2 and self.flow_dir == 'xz+':
        #     self.cells_center = self.cells_center_b[1:-1,2:-2,:]
        # elif self.D == 2 and self.flow_dir == 'yz+':
        #     self.cells_center = self.cells_center_b[1:-1,1:-1,1:-1]
        # elif self.D == 3:
        #     self.cells_center = self.cells_center_b[1:-1,1:-1,1:-1,:]
        
        
    def get_cells_id(self, boundary=True):
        if not hasattr(self, 'cells_id_b'):
            # self.cells_id_b = np.arange(self.pyvista_grid_b.n_cells).reshape(self.flow_shape)
            self.cells_id_b = self.order_b.reshape(self.flow_shape)
        if not hasattr(self, 'cells_id'):
            self.cells_id = self.remove_boundaries(self.cells_id_b)
        if boundary:
            return self.cells_id_b
        else:
            return self.cells_id


    def get_cell_coords(self, id, boundary=True):
        if isinstance(id, (list, tuple, np.ndarray)):
            if boundary:
                return [tuple(x) for x in self.pyvista_grid_b.cell_coords(id)]
            else:
                return [tuple(x) for x in self.pyvista_grid.cell_coords(id)]
        else:
            if boundary:
                return tuple(self.pyvista_grid_b.cell_coords(id))
            else:
                return tuple(self.pyvista_grid.cell_coords(id))


    def get_cells_coords(self, boundary=True, asarray=True):
        if not hasattr(self, 'cells_coords_b'):
            cells_id_b = self.get_cells_id(boundary=True).flatten()
            self.cells_coords_b = [tuple(x) for x in self.pyvista_grid_b.cell_coords(cells_id_b)]
            # needs reshape!
            # arr = np.array(self.cells_coords_b).reshape(self.flow_shape+(3,))
        if not hasattr(self, 'cells_coords'):
            cells_id = self.get_cells_id(boundary=False).flatten()
            self.cells_coords = [tuple(x) for x in self.pyvista_grid_b.cell_coords(cells_id)]
        if boundary:
            return self.cells_coords_b
        else:
            return self.cells_coords


    def get_cell_neighbors(self, id=None, coords=None, boundary=False, fmt='dict'):
        cell_neighbors = {'x': [], 'y': [], 'z': []}
        if id is not None:
            # coords = self.get_cell_coords(id, boundary=True)
            cells_id = self.get_cells_id(boundary).flatten()
            assert id in cells_id, f'cell id is out of range {cells_id}.'
            if self.D >= 1:
                neighbors = [n for n in [id-1, id+1] if n in cells_id]
                cell_neighbors[self.flow_dir[0]] =  neighbors
            if self.D >= 2:
                neighbors = [n for n in [id-self.max_nb, id+self.max_nb] if n in cells_id]
                cell_neighbors[self.flow_dir[1]] = neighbors
            if self.D >= 3:
                neighbors =  [n for n in [id-(self.nx_b*self.ny_b), id+(self.nx_b*self.ny_b)] if n in cells_id]
                cell_neighbors[self.flow_dir[2]] = neighbors
        elif coords is not None:
            cells_coords = self.get_cells_coords(boundary)
            assert coords in cells_coords, f'cell coords are out of range {cells_coords}.'
            i,j,k = coords
            if 'x' in self.flow_dir:
                neighbors = [n for n in [(i-1,j,k), (i+1,j,k)] if n in cells_coords]
                cell_neighbors['x'] = neighbors
            if 'y' in self.flow_dir:
                neighbors = [n for n in [(i,j-1,k), (i,j+1,k)] if n in cells_coords]
                cell_neighbors['y'] = neighbors
            if 'z' in self.flow_dir:
                neighbors = [n for n in [(i,j,k-1), (i,j,k+1)] if n in cells_coords]
                cell_neighbors['z'] = neighbors
        else:
            raise ValueError('at least id or coords argument must be defined.')
        
        if fmt == 'dict':
            return cell_neighbors
        elif fmt == 'list':
            return sum(cell_neighbors.values(), [])
        else: 
            raise ValueError("format argument must be either 'dict' or 'list'.")
        
        
    # @lru_cache(maxsize=None)
    # def get_cell_neighbors(self, id=None, coords=None, boundary=False):
    #     cell_neighbors = []
    #     if id is not None:
    #         cells_id = self.get_cells_id(boundary).flatten()
    #         assert id in cells_id, f'cell id is out of range {cells_id}.'
    #         if self.D >= 1:
    #             cell_neighbors = cell_neighbors + [id-1, id+1]
    #         if self.D >= 2:
    #             cell_neighbors = cell_neighbors + [id-self.max_nb, id+self.max_nb]
    #         if self.D >= 3:
    #             cell_neighbors = cell_neighbors + [id-(self.nx_b*self.ny_b), id+(self.nx_b*self.ny_b)]
    #         return [n for n in cell_neighbors if n in cells_id]
    #     elif coords is not None:
    #         cells_coords = self.get_cells_coords(boundary)
    #         assert coords in cells_coords, f'cell coords are out of range {cells_coords}.'
    #         i,j,k = coords
    #         if 'x' in self.flow_dir:
    #             cell_neighbors = cell_neighbors + [(i-1,j,k), (i+1,j,k)]
    #         if 'y' in self.flow_dir:
    #             cell_neighbors = cell_neighbors + [(i,j-1,k), (i,j+1,k)]
    #         if 'z' in self.flow_dir:
    #             cell_neighbors = cell_neighbors + [(i,j,k-1), (i,j,k+1)]
    #         return [n for n in cell_neighbors if n in cells_coords]
    #     else:
    #         raise ValueError('at least id or coords argument must be defined.')


    # @lru_cache(maxsize=None)
    def get_cell_boundaries(self, id=None, coords=None, fmt='dict'):
        if id is not None:
            boundaries = self.boundaries_id
            cells_id_b = self.get_cells_id(boundary=True).flatten()
            assert id in cells_id_b, f'cell id is out of range {cells_id_b}.'
            cell_neighbors = self.get_cell_neighbors(id=id, boundary=True, fmt=fmt)
        elif coords is not None:
            boundaries = self.boundaries_coords
            cells_coords = self.get_cells_coords(boundary=True)
            assert coords in cells_coords, f'cell coords are out of range {cells_coords}.'
            cell_neighbors = self.get_cell_neighbors(coords=coords, boundary=True, fmt=fmt)
        else:
            raise ValueError('at least id or coords argument must be defined.')
        
        if fmt == 'dict':
            cell_boundaries = {}
            cell_boundaries['x'] = list(set(cell_neighbors['x']).intersection(set(boundaries)))
            cell_boundaries['y'] = list(set(cell_neighbors['y']).intersection(set(boundaries)))
            cell_boundaries['z'] = list(set(cell_neighbors['z']).intersection(set(boundaries)))
            return cell_boundaries
        elif fmt == 'list':
            return list(set(cell_neighbors).intersection(set(boundaries)))
        else:
            raise ValueError("format argument must be either 'dict' or 'list'.")
         
    
    def get_boundaries(self, ids=True):
        '''
        Arguments:
            - i: index in x-direction.
        '''
        if not hasattr(self, 'boundaries_id'):
            cells_id_b = self.get_cells_id(boundary=True) # or self.cells_id_b
            if self.D == 0:
                self.boundaries_id = cells_id_b.flatten()
            elif self.D == 1:
                self.boundaries_id = cells_id_b[[0,-1]].flatten()
            elif self.D == 2:
                self.boundaries_id = np.sort(
                    np.concatenate([
                        cells_id_b[[0,-1],:].flatten(),
                        cells_id_b[1:-1, [0,-1]].flatten()
                ]))
            elif self.D == 3:
                self.boundaries_id = np.sort(np.concatenate([
                    cells_id_b[:,[0,-1],:].flatten(),
                    cells_id_b[:, 1:-1, [0,-1]].flatten(),
                    cells_id_b[[0,-1], 1:-1, 1:-1].flatten(),
                ]))
        if not hasattr(self, 'boundaries_coords'):
            self.boundaries_coords = self.get_cell_coords(self.boundaries_id, boundary=True)
        if ids:
            return self.boundaries_id
        else:
            return self.boundaries_coords
    
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
    
    def show(self, 
            boundary=False,
            points=False, 
            centers=None, # 'coords' or 'id',
            ):
        '''
        - centers_label: str ('coords', 'id')
        '''
        if boundary:
            pv_grid = self.pyvista_grid_b
        else:
            pv_grid = self.pyvista_grid
        
        if self.max_nb > 12:
            opacity = 1
        else:
            opacity = 0.8
            
        pl = pv.Plotter()
        pl.add_mesh(
            pv_grid, 
            show_edges=True, 
            color="white",
            opacity=opacity,
        )
        
        if points:
            points_lst = pv_grid.points
            mask = points_lst[:,1] == 0
            pl.add_point_labels(
                points_lst[mask], 
                points_lst[mask].tolist(), 
                point_size=10, 
                font_size=10,
            )
        
        if centers is not None:
            if centers == 'coords':
                labels = self.get_cells_coords(boundary)
            elif centers == 'id':
                labels = self.get_cells_id(boundary).flatten()
            elif centers == 'volume':
                labels = self.get_cells_volume(boundary).flatten()
            elif centers == 'vol':
                labels = self.get_cells_volume_(boundary).flatten()
            elif centers == 'dx':
                labels = self.get_cells_dx(boundary).flatten()
            elif centers == 'dy':
                labels = self.get_cells_dy(boundary).flatten()
            elif centers == 'dz':
                labels = self.get_cells_dz(boundary).flatten()
            elif centers == 'area_x':
                labels = self.get_cells_area_x(boundary).flatten()
            elif centers == 'area_y':
                labels = self.get_cells_area_y(boundary).flatten()
            elif centers == 'area_z':
                labels = self.get_cells_area_z(boundary).flatten()
            points = self.get_cells_center(boundary).flatten()
            pl.add_point_labels(
                points=points,
                labels=labels,
                point_size=10,
                font_size=10,
            )
        
        s = 'with boundary' if boundary else 'without boundary'
        title = f'{centers} for {self.D}D model (flow at {self.flow_dir}-direction {s})'
        pl.add_title(title, 
                     font='courier', 
                     color='white',
                     font_size=8)
        pl.add_camera_orientation_widget()
        pl.enable_fly_to_right_click()
        pl.show_axes()
        pl.camera_position = 'xz'
        pl.set_background('black', top='gray')
        pl.show(title='openresim 3D show', full_screen=True)


#%%
if __name__ == '__main__':
    # Canvas:
    dx = [11,21,31,41]
    dy = [12,22,32,42]
    dz = [13,23,33,43]
    grid = CartGrid(nx=2, ny=2, nz=2, dx=dx, dy=dy, dz=dz, phi=0.27, k=270)
    print('(1,1,1):',grid.get_cell_d(dir='y', coords=(1,1,1)))
    print('(2,1,1):',grid.get_cell_d(dir='y', coords=(2,1,1)))
    print('(1,2,1):',grid.get_cell_d(dir='y', coords=(1,2,1)))
    print('(1,1,2):',grid.get_cell_d(dir='y', coords=(1,1,2)))
    print('(2,2,2):',grid.get_cell_d(dir='y', coords=(2,2,2)))
    # grid.show(boundary=True, centers='dx')
    # grid.show(boundary=True, centers='dy')
    # grid.show(boundary=True, centers='dz')
    # grid.show(boundary=True, centers='area_x')
    # grid.show(boundary=True, centers='area_y')
    # grid.show(boundary=True, centers='area_z')
    # grid.show(boundary=True, centers='coords')
    grid.show(boundary=True, centers='area_x')
    # grid.show(boundary=True, centers='volume')
    
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