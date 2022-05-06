#%% Import Statements:
from multiprocessing.sharedctypes import Value
from turtle import shape

from sympy import true
from openresim.base import Base
import numpy as np
import pyvista as pv
from functools import lru_cache


#%% Grid Class:
class Grid(Base):
    """Grid class to create a grid using numpy arrays and pyvista grids."""

    def __init__(self, dtype, unit):
        super().__init__(unit)
        self.dtype = dtype  # np.single, np.double


#%% CartGrid Class:
class CartGrid(Grid):
    """Create a cartesian explicit structured grid.

    Note that the following conventions are used:
        1. rows for dx. 1D grids are stored as 1 column with multiple rows.
        2. columns for dy. 2D grids are stored as a table where rows refer to x-direction while columns refere to y-direction.
        3. layers for dz. 3D grids are stored as cubes where rows refer to x-direction, columns to y-direction, and layers for z-direction.

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
    """

    name = "CartGrid"
    type = "cartesian"

    def __init__(
        self,
        nx, ny, nz,
        dx, dy, dz,
        kx=None, ky=None, kz=None,
        phi=None, z=None, comp=None,
        dtype="double", unit="field",
    ):
        super().__init__(dtype, unit)
        self.nx, self.ny, self.nz = nx, ny, nz
        self.get_D()  # > self.D
        self.get_shape()  # > self.shape
        self.get_flow_dir()  # > self.flowdir
        self.get_n_b()  # > self.nx_b, self.ny_b, self.n_zb
        self.get_fshape()  # > self.fshape
        self.get_n_cells()  # > self.n_cells
        self.get_order(type="natural")  # > self.order
        self.get_blocks()  # > self.blocks
        self.max_nb = max(self.nx_b, self.ny_b, self.nz_b)
        self.get_cells_dims(dx, dy, dz)
        self.get_cells_area_x()
        self.get_cells_area_y()
        self.get_cells_area_z()
        self.get_cells_volume()
        self.get_cells_coords()  # > cells_coords
        self.set_props(phi, kx, ky, kz, z, comp)
        self.get_boundaries()  # > self.boundaries_id, self.boundaries_coords
        self.get_volume()  # self.volume
        self.get_Gx()  # > self.Gx
        self.get_Gy()  # > self.Gy
        self.get_Gz()  # > self.Gz
        self.get_cells_center()

    @lru_cache(maxsize=1)
    def get_D(self):
        self.D = sum([1 if n > 1 else 0 for n in (self.nx, self.ny, self.nz)])
        print(f"- Dimension (D) is set to {self.D}.")
        # self.dimension = self.D
        return self.D
    # get_dimension = get_D

    @lru_cache(maxsize=2)
    def get_shape(self, boundary=False):
        if boundary:
            self.get_n_b()
            self.shape = np.array((self.nx_b, self.ny_b, self.nz_b),
                                  dtype="int")
        else:
            self.shape = np.array((self.nx, self.ny, self.nz), dtype="int")
        print(f"- Shape is set to {self.shape}.")
        return self.shape

    @lru_cache(maxsize=1)
    def get_flow_dir(self):
        self.get_D()
        if self.D == 0:
            self.flow_dir = "-"
        elif self.D == 1:
            flow_dir_id = np.argmax(self.shape)
            if flow_dir_id == 0:
                self.flow_dir = "x"
            elif flow_dir_id == 1:
                self.flow_dir = "y"
            elif flow_dir_id == 2:
                self.flow_dir = "z"
        elif self.D == 2:
            flow_dir_id = np.argmin(self.shape)
            if flow_dir_id == 2:
                self.flow_dir = "xy"
            elif flow_dir_id == 1:
                self.flow_dir = "xz"
            elif flow_dir_id == 0:
                self.flow_dir = "yz"
        elif self.D == 3:
            self.flow_dir = "xyz"
        print(f"- Flow direction (flowdir) is set to {self.flow_dir}.")
        return self.flow_dir

    @lru_cache(maxsize=1)
    def get_n_b(self):
        self.get_flow_dir()
        if "x" in self.flow_dir:
            self.nx_b = self.nx + 2
        else:
            self.nx_b = self.nx
        if "y" in self.flow_dir:
            self.ny_b = self.ny + 2
        else:
            self.ny_b = self.ny
        if "z" in self.flow_dir:
            self.nz_b = self.nz + 2
        else:
            self.nz_b = self.nz
        return (self.nx_b, self.ny_b, self.nz_b)

    @lru_cache(maxsize=1)
    def get_fshape(self):
        self.get_flow_dir()
        if self.flow_dir == "-":
            self.fshape = (1,)
        elif self.flow_dir == "x":
            self.fshape = (self.nx_b,)
        elif self.flow_dir == "y":
            self.fshape = (self.ny_b,)
        elif self.flow_dir == "z":
            self.fshape = (self.nz_b,)
        elif self.flow_dir == "xy":
            self.fshape = (self.ny_b, self.nx_b)
        elif self.flow_dir == "xz":
            self.fshape = (self.nz_b, self.nx_b)
        elif self.flow_dir == "yz":
            self.fshape = (self.nz_b, self.ny_b)
        elif self.flow_dir == "xyz":
            self.fshape = (self.nz_b, self.ny_b, self.nx_b)
        print(f"- Flow shape (fshape) is set to {self.fshape}")
        """
        - in case of unified 3d shape:
        # if self.flow_dir == "-":
        #     self.fshape = (1, 1, 1)
        # elif self.flow_dir == "x":
        #     self.fshape = (1, 1, self.nx_b)
        # elif self.flow_dir == "y":
        #     self.fshape = (1, self.ny_b, 1)
        # elif self.flow_dir == "z":
        #     self.fshape = (self.nz_b, 1, 1)
        # elif self.flow_dir == "xy":
        #     self.fshape = (1, self.ny_b, self.nx_b)
        # elif self.flow_dir == "xz":
        #     self.fshape = (self.nz_b, 1, self.nx_b)
        # elif self.flow_dir == "yz":
        #     self.fshape = (self.nz_b, self.ny_b, 1)
        
        - in case of 2D but also with boundary conditions:
            elif self.flowdir == 'xz+':
                self.fshape = (self.nz_b, self.nx_b, 3)
            elif self.flowdir == 'yz+':
                self.fshape = (self.nz_b, self.ny_b, 3)
        """
        # self.flow_shape = self.fshape
        return self.fshape
    # get_flow_shape = get_fshape

    @lru_cache(maxsize=2)
    def get_n_cells(self, boundary=True):
        if boundary:
            self.get_n_b()
            self.n_cells = self.nx_b * self.ny_b * self.nz_b
        else:
            self.n_cells = self.nx * self.ny * self.nz
        return self.n_cells

    @lru_cache(maxsize=4)
    def get_order(
        self, 
        type="natural", 
        boundary=True, 
        fshape=False, 
        verbose=True
        ):

        if type == "natural":
            self.order = np.arange(self.get_n_cells(True))
        else:
            raise ValueError(
                "Order type is not supported or unknown. "
                "Supported order types: ['natural']"
            )

        if fshape:
            self.order = self.order.reshape(self.fshape)

        if not boundary:
            self.order = self.remove_boundaries(self.order)

        if verbose:
            s1, s2 = self.get_verbose_str(boundary, fshape)
            print(f"- Cells order (order) was computed ({s1} - {s2}).")
        return self.order

    def get_blocks(self):
        self.blocks = np.ones(self.fshape, dtype="int")
        # self.i_blocks = self.blocks.copy()
        # self.i_blocks[[0, -1]] = 0
        # self.b_blocks = np.zeros(self.max_nb, dtype='int') #ss.lil_matrix(self.shape, dtype='int')
        # self.b_blocks[[0, -1]] = 1
        
            
    """
    Cells id and coordinates:
    """

    @lru_cache(maxsize=None)
    def get_cell_id(self, coords, boundary=True):
        pyvista_grid = self.get_pyvista_grid(boundary)
        return pyvista_grid.cell_id(coords)

    @lru_cache(maxsize=4)
    def get_cells_id(self, boundary=True, fshape=True, verbose=True):
        self.cells_id = self.get_order("natural", boundary, fshape, False)
        if verbose:
            s1, s2 = self.get_verbose_str(boundary, fshape)
            print(f"- Cells id (cells_id) was computed ({s1} - {s2}).")
        return self.cells_id

    # @lru_cache(maxsize=None): does not work
    def get_cell_coords(self, id, boundary=True):
        pyvista_grid = self.get_pyvista_grid(boundary)
        if isinstance(id, (list, tuple, np.ndarray)):
            cell_coords = [tuple(x) for x in 
                            pyvista_grid.cell_coords(id)]
        else:
            cell_coords = tuple(pyvista_grid.cell_coords(id))
        return cell_coords

    @lru_cache(maxsize=4)
    def get_cells_coords(self, boundary=True, fshape=False, verbose=True):
        cells_id = self.get_cells_id(boundary, False, False)
        pyvista_grid = self.get_pyvista_grid(True)
        self.cells_coords = [tuple(x) for x in 
                             pyvista_grid.cell_coords(cells_id)]

        if fshape:
            self.cells_coords = self.cells_coords.reshape(self.fshape)

        if verbose:
            s1, s2 = self.get_verbose_str(boundary, fshape)
            print(f"- Cells coords (cells_coords) was computed ({s1} - {s2}).")

        return self.cells_coords

    def get_cell_icoords(self, coords):
        if self.D <= 2:
            icoords = tuple(c for c in coords[::-1] if c > 0)
        else:
            icoords = tuple(c for c in coords[::-1])
        return icoords

    def add_boundary(self, d):
        if isinstance(d, [list, tuple, np.array]):
            if len(d) == self.max_nb - 2:
                d = np.concatenate((d[0], d, d[-1]))
        return d

    """
    Properties:
    """

    def set_props(
        self,
        phi=None,
        kx=None,
        ky=None,
        kz=None,
        z=None,
        comp=None,
        id=None,
        coords=None,
    ):
        """
        Set properties
        """
        # Set user defined properties:
        if phi is not None:
            self.set_phi(phi, id, coords)
        if kx is not None:
            self.set_kx(kx, id, coords)
        if ky is not None:
            self.set_ky(ky, id, coords)
        if kz is not None:
            self.set_kz(kz, id, coords)
        if z is not None:
            self.set_z(z, id, coords)
        if comp is not None:
            self.set_compressibility(comp)
        # Set default values if not defined:
        if not hasattr(self, "tops"):
            self.set_z(0)
        if not hasattr(self, "compressibility"):
            self.set_compressibility(0)
        if not hasattr(self, "is_homogeneous"):
            self.get_is_homogeneous()

    # set_properties = set_props

    def set_phi(self, phi, id=None, coords=None):
        """ """
        if id is None and coords is None:
            self.phi = self.blocks * phi
            # self.porosity = self.phi
            # np.zeros(self.shape) + phi
            s = "all cells"
        else:
            if id is not None:
                coords = self.get_cell_coords(id, True)
            self.phi[coords[2], coords[1], coords[0]] = phi
            s = "cell " + str(coords)
            # self.get_is_homogeneous()
        print(f"- Porosity (phi) is set to {phi} for {s}.")
    # set_porosity = set_phi

    def set_kx(self, kx, id=None, coords=None):
        if id is None and coords is None:
            self.kx = self.blocks * kx
            # self.permeability_x = self.kx
            s = "all cells"
        else:
            if id is not None:
                kx_f = self.kx.flatten()
                kx_f[id] = kx
                self.kx = kx_f.reshape(self.fshape)
            if coords is not None:
                icoords = self.get_cell_icoords(coords)
                self.kx[icoords] = kx
            s = "cell " + str(coords)
            self.get_is_homogeneous()
        print(f"- Permeability at x-direction (kx) is set to {kx} for {s}.")
    # set_permeability_x = set_kx

    def set_ky(self, ky, id=None, coords=None):
        if id is None and coords is None:
            self.ky = self.blocks * ky
            # self.permeability_y = self.ky
            s = "all cells"
        else:
            if id is not None:
                ky_f = self.ky.flatten()
                ky_f[id] = ky
                self.ky = ky_f.reshape(self.fshape)
            if coords is not None:
                icoords = self.get_cell_icoords(coords)
                self.ky[icoords] = ky
            s = "cell " + str(coords)
            self.get_is_homogeneous()
        print(f"- Permeability at y direction (ky) is set to {ky} for {s}.")
    # set_permeability_y = set_ky

    def set_kz(self, kz, id=None, coords=None):
        if id is None and coords is None:
            self.kz = self.blocks * kz
            # self.permeability_z = self.kz
            s = "all cells"
        else:
            if id is not None:
                kz_f = self.kz.flatten()
                kz_f[id] = kz
                self.kz = kz_f.reshape(self.fshape)
            if coords is not None:
                icoords = self.get_cell_icoords(coords)
                self.kz[icoords] = kz
            s = "cell " + str(coords)
            self.get_is_homogeneous()
        print(f"- Permeability at z direction (kz) is set to {kz} for {s}.")
    # set_permeability_z = set_kz

    @lru_cache(maxsize=6)
    def get_k(self, dir, boundary=True):
        if dir == "x":
            k = self.kx
        elif dir == "y":
            k = self.ky
        elif dir == "z":
            k = self.kz
        elif dir in [None, "all"]:
            k = {"x": self.kx, "y": self.ky, "z": self.kz}
        else:
            raise ValueError("dir is unknown!")

        if boundary:
            return k
        else:
            return self.remove_boundaries(k)
    # get_permeability = get_k

    def set_z(self, z=None, id=None, coords=None):
        """ """
        if id is None and coords is None:
            self.z = self.blocks * z
            # self.tops = self.z
            s = "all cells"
        else:
            if id is not None:
                z_f = self.z.flatten()
                z_f[id] = z
                self.z = z_f.reshape(self.fshape)
            if coords is not None:
                icoords = self.get_cell_icoords(coords)
                self.z[icoords] = z
            s = "cell " + str(coords)
        print(f"- Tops (z) is set to {z} for {s}.")
    set_tops = set_z

    def get_is_homogeneous(self):
        if (
            np.all(self.kx[1:-1] == self.kx[0])
            & np.all(self.ky[1:-1] == self.kx[0])
            & np.all(self.kz[1:-1] == self.kx[0])
            & np.all(self.phi[1:-1] == self.phi[0])
        ):
            self.is_homogeneous = True
        else:
            self.is_homogeneous = False
        print(f"- Flag is_homogeneous is set to {self.is_homogeneous}.")
        return self.is_homogeneous

    def get_G(self, dir):
        k = self.get_k(dir=dir, boundary=True)
        area = self.get_cells_area(dir=dir, boundary=True)
        d = self.get_cells_d(dir=dir, boundary=True)
        if self.is_homogeneous:
            G = (
                self.factors["transmissibility conversion"]
                * self.get_G_homo_mean(k)
                * self.get_G_homo_mean(area)
                / self.get_G_homo_mean(d)
            )
        else:
            G = (
                2
                * self.factors["transmissibility conversion"]
                / self.get_G_hetro_denom(d, area, k)
            )
        return G

    def get_Gx(self):
        """
        Grid geometry factor at x-direction.
        """
        self.Gx = self.get_G(dir="x")
        return self.Gx

    def get_Gy(self):
        """
        Grid geometry factor at y-direction.
        """
        self.Gy = self.get_G(dir="y")
        return self.Gy

    def get_Gz(self):
        """
        Grid geometry factor at z-direction.
        """
        self.Gz = self.get_G(dir="z")
        return self.Gz

    def get_G_hetro_denom(self, dx, area, k):
        if self.D == 0:
            return (dx / (area * k))
        elif self.D == 1:
            return  (dx[:-1] / (area[:-1] * k[:-1])) + \
                    (dx[1:]  / (area[1:]  * k[1:]))
        elif self.D == 2:
            return  (dx[:-1, :-1] / (area[:-1, :-1] * k[:-1, :-1])) + \
                    (dx[1:, 1:]   / (area[1:, 1:]   * k[1:, 1:]))
        elif self.D == 3:
            return (dx[:-1,:-1,:-1]/(area[:-1,:-1,:-1]*k[:-1,:-1,:-1])) + \
                   (dx[1:, 1:, 1:] /(area[1:, 1:, 1:] *k[1:, 1:, 1:]))

    def get_G_homo_mean(self, property, type="geometric"):
        """ """
        if self.is_homogeneous:
            if self.D == 0:
                return property
            elif self.D == 1:
                return property[1:]
            elif self.D == 2:
                return property[1:, 1:]
            elif self.D == 3:
                return property[1:, 1:, 1:]
        else:
            if type == "geometric":
                if self.D == 1:
                    return (property[:-1] + property[1:]) / 2
                elif self.D == 2:
                    return (property[:-1, :-1] + property[1:, 1:]) / 2
                elif self.D == 2:
                    return (property[:-1, :-1, :-1] + property[1:, 1:, 1:]) / 2
            else:
                raise ValueError("Unknown mean type")

    
    """
    Geometry and dimensions
    """
    
    def get_cells_dims(self, dx, dy, dz, verbose=True):
        self.get_n_b()
        m_list = []
        if "x" in self.flow_dir:
            self.dxx = np.ones(self.nx_b, dtype="int") * dx
            m_list.append(self.dxx)
        else:
            self.dxx = np.ones(self.max_nb, dtype="int") * dx
            m_list.append(dx)
        if "y" in self.flow_dir:
            self.dyy = np.ones(self.ny_b, dtype="int") * dy
            m_list.append(self.dyy)
        else:
            self.dyy = np.ones(self.max_nb, dtype="int") * dy
            m_list.append(dy)
        if "z" in self.flow_dir:
            self.dzz = np.ones(self.nz_b, dtype="int") * dz
            m_list.append(self.dzz)
        else:
            self.dzz = np.ones(self.max_nb, dtype="int") * dz
            m_list.append(dz)

        self.dx, self.dy, self.dz = np.meshgrid(*m_list, copy=False)
        self.dx = np.transpose(self.dx, axes=(0, 2, 1)).reshape(self.fshape)
        self.dy = np.transpose(self.dy, axes=(2, 0, 1)).reshape(self.fshape)
        self.dz = np.transpose(self.dz, axes=(2, 1, 0)).reshape(self.fshape)

        if verbose:
            print(f"- Cells dims (dx, dy, dz) were computed.")
            
    @lru_cache(maxsize=2)
    def get_pyvista_grid(self, boundary=True, verbose=True):
        """
        https://docs.pyvista.org/api/core/_autosummary/pyvista.ExplicitStructuredGrid.html
        """
        if boundary:
            dims = np.array((self.nx_b, self.ny_b, self.nz_b)) + 1
        else:
            dims = np.array((self.nx, self.ny, self.nz)) + 1

        corners = self.get_corners(boundary, verbose)
        pyvista_grid = pv.ExplicitStructuredGrid(dims, corners)

        if verbose:
            s = "with boundary" if boundary else "without boundary"
            print(f"- Pyvista grid (pyvista_grid) {s} was created.")
        
        return pyvista_grid

    @lru_cache(maxsize=2)
    def get_corners(self, boundary=True, verbose=True):
        """

        (Reference: https://docs.pyvista.org/examples/00-load/create-explicit-structured-grid.html)
        """

        if "x" in self.flow_dir:
            xcorn = np.insert(self.dxx.cumsum(), 0, 0)
        else:
            xcorn = np.arange(0, (self.nx + 1) * self.dxx[0], self.dxx[0])

        if "y" in self.flow_dir:
            ycorn = np.insert(self.dyy.cumsum(), 0, 0)
        else:
            ycorn = np.arange(0, (self.ny + 1) * self.dyy[0], self.dyy[0])

        if "z" in self.flow_dir:
            zcorn = np.insert(self.dzz.cumsum(), 0, 0)
        else:
            zcorn = np.arange(0, (self.nz + 1) * self.dzz[0], self.dzz[0])

        # Boundary:
        if boundary:
            ix = 2 if "x" in self.flow_dir else 0
            iy = 2 if "y" in self.flow_dir else 0
            iz = 2 if "z" in self.flow_dir else 0
        else:
            ix = 0
            iy = 0
            iz = 0
            if "x" in self.flow_dir:
                xcorn = xcorn[1:-1]
            if "y" in self.flow_dir:
                ycorn = ycorn[1:-1]
            if "z" in self.flow_dir:
                zcorn = zcorn[1:-1]

        # X corners:
        xcorn = np.repeat(xcorn, 2)
        xcorn = xcorn[1:-1]
        xcorn = np.tile(xcorn, 4 * (self.ny + iy) * (self.nz + iz))

        # Y corners:
        ycorn = np.repeat(ycorn, 2)
        ycorn = ycorn[1:-1]
        ycorn = np.tile(ycorn, (2 * (self.nx + ix), 2 * (self.nz + iz)))
        ycorn = np.transpose(ycorn)
        ycorn = ycorn.flatten()

        # Z corners:
        zcorn = np.repeat(zcorn, 2)
        zcorn = zcorn[1:-1]
        zcorn = np.repeat(zcorn, 4 * (self.nx + ix) * (self.ny + iy))

        if verbose:
            s = "with boundary" if boundary else "without boundary"
            print(f"- Grid corners {s} were calculated.")
            print("    - xcorn shape:", xcorn.shape,
                "- ycorn shape:", ycorn.shape,
                "- zcorn shape:", zcorn.shape
            )

        # Combine corners:
        corners = np.stack((xcorn, ycorn, zcorn))
        corners = corners.transpose()

        return corners

    def get_verbose_str(self, boundary, fshape):
        s1 = "with fshape" if fshape else "without fshape"
        s2 = "with boundary" if boundary else "without boundary"
        return s1, s2

    @lru_cache(maxsize=4)
    def get_cells_center(self, boundary=True, fshape=False, verbose=True):
        pyvista_grid = self.get_pyvista_grid(True)
        self.cells_center = pyvista_grid.cell_centers().points

        if fshape:
            flow_shape = self.fshape + (3,)
            self.cells_center = self.cells_center.reshape(flow_shape)

        if not boundary:
            self.cells_center = self.remove_boundaries(self.cells_center, True)

        if verbose:
            s1, s2 = self.get_verbose_str(boundary, fshape)
            print(f"- Cells center (cells_center) was computed ({s1} - {s2}).")

        return self.cells_center

    """
    Volume Calculations:
    """

    @lru_cache(maxsize=2)
    def get_volume(self, boundary=True):
        pyvista_grid = self.get_pyvista_grid(boundary)
        self.V = pyvista_grid.volume
        return self.V

    @lru_cache(maxsize=6)
    def get_cells_volume(self, boundary=True, fshape=False, pyvista=False):

        if pyvista:
            pyvista_grid = self.get_pyvista_grid(True)
            self.cells_V = pyvista_grid.compute_cell_sizes()["Volume"]
            self.cells_V = self.cells_V.round(2)
        else:
            self.cells_V = self.dx.flatten() * self.dy.flatten() * self.dz.flatten()
        if fshape:
            self.cells_V = self.cells_V.reshape(self.fshape)
        if not boundary:
            self.cells_V = self.remove_boundaries(self.cells_V)
        print("- Cells volumes (cells_V) was computed.")
        return self.cells_V

    @lru_cache(maxsize=None)
    def get_cell_volume(self, id=None, coords=None):
        if id is not None:
            cells_volume = self.get_cells_volume(True, False)
            return cells_volume.flatten()[id]
        elif coords is not None:
            cells_volume = self.get_cells_volume(True, True)
            return cells_volume[coords[2], coords[1], coords[0]]
        else:
            raise ValueError("at least id or coords argument must be defined.")

    """
    Cells Dimensions (dx, dy, dz) Calculations:
    """

    @lru_cache(maxsize=5)
    def get_cells_d(self, dir, boundary=True, fshape=True):
        if dir == "x":
            cells_d = self.dx
        elif dir == "y":
            cells_d = self.dy
        elif dir == "z":
            cells_d = self.dz

        if fshape:
            cells_d = cells_d.reshape(self.fshape)
        if not boundary:
            cells_d = self.remove_boundaries(cells_d)
        return cells_d

    @lru_cache(maxsize=2)
    def get_cells_dx(self, boundary=True, fshape=True):
        return self.get_cells_d("x", boundary, fshape)

    @lru_cache(maxsize=2)
    def get_cells_dy(self, boundary=True, fshape=True):
        return self.get_cells_d("y", boundary, fshape)

    @lru_cache(maxsize=2)
    def get_cells_dz(self, boundary=True, fshape=True):
        return self.get_cells_d("z", boundary, fshape)

    @lru_cache(maxsize=None)
    def get_cell_d(self, dir, id=None, coords=None):
        cells_d = self.get_cells_d(dir=dir, boundary=True, fshape=True)
        if id is not None:
            return cells_d.flatten()[id]
        elif coords is not None:
            icoords = self.get_cell_icoords(coords)
            return cells_d[icoords]
        else:
            raise ValueError("at least id or coords argument must be defined.")

    @lru_cache(maxsize=None)
    def get_cell_dx(self, id=None, coords=None):
        return self.get_cell_d("x", id, coords)

    @lru_cache(maxsize=None)
    def get_cell_dy(self, id=None, coords=None):
        return self.get_cell_d("y", id, coords)

    @lru_cache(maxsize=None)
    def get_cell_dz(self, id=None, coords=None):
        return self.get_cell_d("z", id, coords)

    """
    Area calculations:
    """

    @lru_cache(maxsize=4)
    def get_cells_area_x(self, boundary=True, fshape=True):
        self.area_x = self.dy.flatten() * self.dz.flatten()
        if fshape:
            self.area_x = self.area_x.reshape(self.fshape)
        if not boundary:
            self.area_x = self.remove_boundaries(self.area_x)
        return self.area_x

    @lru_cache(maxsize=4)
    def get_cells_area_y(self, boundary=True, fshape=True):
        self.area_y = self.dx.flatten() * self.dz.flatten()
        if fshape:
            self.area_y = self.area_y.reshape(self.fshape)
        if not boundary:
            self.area_y = self.remove_boundaries(self.area_y)
        return self.area_y

    @lru_cache(maxsize=4)
    def get_cells_area_z(self, boundary=True, fshape=True):
        self.area_z = self.dx.flatten() * self.dy.flatten()
        if fshape:
            self.area_z = self.area_z.reshape(self.fshape)
        if not boundary:
            self.area_z = self.remove_boundaries(self.area_z)
        return self.area_z

    @lru_cache(maxsize=8)
    def get_cells_area(self, dir=None, boundary=True, fshape=True):
        if dir == "x":
            return self.get_cells_area_x(boundary, fshape)
        elif dir == "y":
            return self.get_cells_area_y(boundary, fshape)
        elif dir == "z":
            return self.get_cells_area_z(boundary, fshape)
        elif dir in [None, "all"]:
            self.area = {
                "x": self.get_cells_area_x(boundary, fshape),
                "y": self.get_cells_area_y(boundary, fshape),
                "z": self.get_cells_area_z(boundary, fshape),
            }
            return self.area
        else:
            raise ValueError(
                "Direction (dir) is unknown! "
                "Argument dir can be one of the following: "
                "['x','y','z','all', None]"
            )

    @lru_cache(maxsize=None)
    def get_cell_area(self, dir, id=None, coords=None):
        if id is not None:
            cells_area = self.get_cells_area(dir, True, False)
            return cells_area[id]
        elif coords is not None:
            cells_area = self.get_cells_area(dir, True, True)
            return cells_area[coords[2], coords[1], coords[0]]
        else:
            raise ValueError("At least id or coords argument must be defined.")

    @lru_cache(maxsize=None)
    def get_cell_area_x(self, id=None, coords=None):
        return self.get_cell_area("x", id, coords)

    @lru_cache(maxsize=None)
    def get_cell_area_y(self, id=None, coords=None):
        return self.get_cell_area("y", id, coords)

    @lru_cache(maxsize=None)
    def get_cell_area_z(self, id=None, coords=None):
        return self.get_cell_area("z", id, coords)

    def isin(self, arr, coords):
        for x in arr:
            if tuple(x) == coords:
                return True
        return False

    def remove_boundaries(self, in_data, points=False):
        if not isinstance(in_data, dict):
            if points:
                fshape = self.fshape + (3,)
            else:
                fshape = self.fshape

            if in_data.shape != fshape:
                try:
                    in_data = in_data.reshape(fshape)
                    flatten = True
                except:
                    print(
                        "in_data shape:",
                        in_data.shape,
                        "- required shape:",
                        self.fshape,
                    )
                    raise ValueError("array must have all cells to remove boundaries!")
            else:
                flatten = False
            if self.D == 0:
                out_data = in_data
            elif self.D == 1:
                out_data = in_data[1:-1]
            elif self.D == 2:
                out_data = in_data[1:-1, 1:-1]
            elif self.D == 3:
                out_data = in_data[1:-1, 1:-1, 1:-1]
            else:
                out_data = in_data[1:-1, 1:-1, 1:-1]
            if flatten and not points:
                out_data = out_data.flatten()
            if flatten and points:
                out_data = out_data.reshape((-1, 3))
            return out_data
        elif isinstance(in_data, dict):
            for k, v in in_data.items():
                in_data[k] = self.remove_boundaries(v)
            return in_data
        else:
            raise ValueError("dtype is unknown.")
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

    @lru_cache(maxsize=None)
    def get_cell_neighbors(self, id=None, coords=None, boundary=False, fmt="dict"):
        cell_neighbors = {"x": [], "y": [], "z": []}
        if id is not None:
            # coords = self.get_cell_coords(id, boundary=True)
            cells_id = self.get_cells_id(boundary).flatten()
            assert id in cells_id, f"cell id is out of range {cells_id}."
            if self.D >= 1:
                neighbors = [n for n in [id - 1, id + 1] if n in cells_id]
                cell_neighbors[self.flow_dir[0]] = neighbors
            if self.D >= 2:
                neighbors = [
                    n for n in [id - self.max_nb, id + self.max_nb] if n in cells_id
                ]
                cell_neighbors[self.flow_dir[1]] = neighbors
            if self.D >= 3:
                neighbors = [
                    n for n in [
                        id - (self.nx_b * self.ny_b),
                        id + (self.nx_b * self.ny_b),
                    ] if n in cells_id
                ]
                cell_neighbors[self.flow_dir[2]] = neighbors
        elif coords is not None:
            cells_coords = self.get_cells_coords(boundary=boundary, fshape=False)
            assert (
                coords in cells_coords
            ), f"cell coords are out of range {cells_coords}."
            i, j, k = coords
            if "x" in self.flow_dir:
                neighbors = [
                    n for n in [
                        (i - 1, j, k), 
                        (i + 1, j, k)
                    ] if n in cells_coords
                ]
                cell_neighbors["x"] = neighbors
            if "y" in self.flow_dir:
                neighbors = [
                    n for n in [
                        (i, j - 1, k), 
                        (i, j + 1, k)
                    ] if n in cells_coords
                ]
                cell_neighbors["y"] = neighbors
            if "z" in self.flow_dir:
                neighbors = [
                    n for n in [
                        (i, j, k - 1), 
                        (i, j, k + 1)
                    ] if n in cells_coords
                ]
                cell_neighbors["z"] = neighbors
        else:
            raise ValueError("at least id or coords argument must be defined.")

        if fmt == "dict":
            return cell_neighbors
        elif fmt == "list":
            return sum(cell_neighbors.values(), [])
        else:
            raise ValueError("format argument must be either 'dict' or 'list'.")

    # @lru_cache(maxsize=None)
    def get_cell_boundaries(self, id=None, coords=None, fmt="dict"):
        if id is not None:
            boundaries = self.boundaries_id
            cells_id_b = self.get_cells_id(boundary=True).flatten()
            assert id in cells_id_b, f"cell id is out of range {cells_id_b}."
            cell_neighbors = self.get_cell_neighbors(id=id, boundary=True, fmt=fmt)
        elif coords is not None:
            boundaries = self.boundaries_coords
            cells_coords = self.get_cells_coords(boundary=True, fshape=False)
            assert (
                coords in cells_coords
            ), f"cell coords are out of range {cells_coords}."
            cell_neighbors = self.get_cell_neighbors(
                coords=coords, boundary=True, fmt=fmt
            )
        else:
            raise ValueError("at least id or coords argument must be defined.")

        if fmt == "dict":
            cell_boundaries = {}
            cell_boundaries["x"] = list(
                set(cell_neighbors["x"]).intersection(set(boundaries))
            )
            cell_boundaries["y"] = list(
                set(cell_neighbors["y"]).intersection(set(boundaries))
            )
            cell_boundaries["z"] = list(
                set(cell_neighbors["z"]).intersection(set(boundaries))
            )
            return cell_boundaries
        elif fmt == "list":
            return list(set(cell_neighbors).intersection(set(boundaries)))
        else:
            raise ValueError("format argument must be either 'dict' or 'list'.")

    @lru_cache(maxsize=2)
    def get_boundaries(self, ids=True):
        """
        Arguments:
            - i: index in x-direction.
        """
        cells_id = self.get_cells_id(True, True)  # or self.cells_id
        if self.D == 0:
            self.boundaries_id = cells_id.flatten()
        elif self.D == 1:
            self.boundaries_id = cells_id[[0, -1]].flatten()
        elif self.D == 2:
            self.boundaries_id = np.sort(
                np.concatenate(
                    [
                        cells_id[[0, -1], :].flatten(),
                        cells_id[1:-1, [0, -1]].flatten(),
                    ]
                )
            )
        elif self.D == 3:
            self.boundaries_id = np.sort(
                np.concatenate(
                    [
                        cells_id[:, [0, -1], :].flatten(),
                        cells_id[:, 1:-1, [0, -1]].flatten(),
                        cells_id[[0, -1], 1:-1, 1:-1].flatten(),
                    ]
                )
            )
        self.boundaries_coords = self.get_cell_coords(
            self.boundaries_id, boundary=True
        )
        if ids:
            return self.boundaries_id
        else:
            return self.boundaries_coords

    def show(
        self,
        boundary=False,
        corners=False,
        label=None,  # 'coords' or 'id',
    ):
        """
        - centers_label: str ('coords', 'id')
        """
        pyvista_grid = self.get_pyvista_grid(boundary)

        if self.max_nb > 12:
            opacity = 1
        else:
            opacity = 0.8

        pl = pv.Plotter()
        pl.add_mesh(
            pyvista_grid,
            show_edges=True,
            color="white",
            opacity=opacity,
        )

        if corners:
            points_lst = self.get_pyvista_grid(True).points
            if boundary:
                mask = points_lst[:, 1] == 0
            else:
                id = self.get_order(boundary=False)[0]
                mask = points_lst[:, 1] == self.get_cell_dy(id=id)
            pl.add_point_labels(
                points_lst[mask],
                points_lst[mask].tolist(),
                point_size=10,
                font_size=10,
            )

        if label is not None:
            if label == "coords":
                labels = self.get_cells_coords(boundary, False)
            elif label == "id":
                labels = self.get_cells_id(boundary, False)
            elif label == "volume":
                labels = self.get_cells_volume(boundary, False)
            elif label == "center":
                labels = self.get_cells_center(boundary, False)
            elif label == "dx":
                labels = self.get_cells_dx(boundary, False)
            elif label == "dy":
                labels = self.get_cells_dy(boundary, False)
            elif label == "dz":
                labels = self.get_cells_dz(boundary, False)
            elif label == "area_x":
                labels = self.get_cells_area_x(boundary, False)
            elif label == "area_y":
                labels = self.get_cells_area_y(boundary, False)
            elif label == "area_z":
                labels = self.get_cells_area_z(boundary, False)
            else:
                raise ValueError(f"{label} can't be used!")
            points = self.get_cells_center(boundary, False)
            pl.add_point_labels(
                points=points,
                labels=labels,
                point_size=10,
                font_size=10,
            )

        s = "with boundary" if boundary else "without boundary"
        title = f"{self.D}D model by {label}" + \
                f"(flow at {self.flow_dir}-direction {s})"
        pl.add_title(title, font="courier", color="white", font_size=8)
        pl.add_camera_orientation_widget()
        pl.enable_fly_to_right_click()
        pl.show_axes()
        pl.camera_position = "xz"
        pl.set_background("black", top="gray")
        pl.show(title="openresim 3D show", full_screen=True)


#%%
if __name__ == "__main__":
    # Canvas:
    dx = 11  # [11,21,31,41]
    dy = 12  # [12,22,32,42]
    dz = 13  # [13,23,33,43]
    grid = CartGrid(
        nx=2, ny=2, nz=1,
        dx=dx, dy=dy, dz=dz,
        kx=270, ky=10, kz=20,
        phi=0.27,
    )

    id = grid.get_cell_coords(3, True)
    # print(id)
    
    # grid.set_kx(kx=10, coords=(1,1,1))
    # print(grid.kx)
    # print(grid.get_cells_id())
    # print(grid.kx.shape)
    # print(grid.cells_volume_b.shape)
    grid.show(boundary=True, label="id", corners=False)
    grid.show(boundary=False, label="id", corners=False)
    # grid.show(boundary=True, label="dy")
    # grid.show(boundary=True, label="dz")
    # grid.show(boundary=True, label="area_x")
    # grid.show(boundary=True, label="area_y")
    # grid.show(boundary=True, label="area_z")

    # grid.show(boundary=True, centers='coords')
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
