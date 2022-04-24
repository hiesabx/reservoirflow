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
    """CartGrid class to create a explicit structured grid.

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
    """

    name = "CartGrid"

    def __init__(
        self,
        nx,
        ny,
        nz,
        dx,
        dy,
        dz,
        kx=None,
        ky=None,
        kz=None,
        phi=None,
        z=None,
        comp=None,
        dtype="double",
        unit="field",
    ):
        super().__init__(dtype, unit)
        self.type = "cartesian"
        self.nx, self.ny, self.nz = nx, ny, nz
        self.get_D()  # > self.D
        self.get_shape()  # > self.shape
        self.get_flowdir()  # > self.flowdir
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

        self.pyvista_grid_b = self.get_pyvista_grid(boundary=True)
        self.pyvista_grid = self.get_pyvista_grid(boundary=False)
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
            self.shape = np.array((self.nx_b, self.ny_b, self.nz_b), dtype="int")
        else:
            self.shape = np.array((self.nx, self.ny, self.nz), dtype="int")
        print(f"- Shape is set to {self.shape}.")
        return self.shape

    @lru_cache(maxsize=1)
    def get_flowdir(self):
        if self.D == 0:
            self.flowdir = "-"
        elif self.D == 1:
            flowdir_id = np.argmax(self.shape)
            if flowdir_id == 0:
                self.flowdir = "x"
            elif flowdir_id == 1:
                self.flowdir = "y"
            elif flowdir_id == 2:
                self.flowdir = "z"
        elif self.D == 2:
            flowdir_id = np.argmin(self.shape)
            if flowdir_id == 2:
                self.flowdir = "xy"
            elif flowdir_id == 1:
                self.flowdir = "xz"
            elif flowdir_id == 0:
                self.flowdir = "yz"
        elif self.D == 3:
            self.flowdir = "xyz"
        print(f"- flow direction (flowdir) is set to {self.flowdir}.")
        return self.flowdir

    @lru_cache(maxsize=1)
    def get_n_b(self):
        if "x" in self.flowdir:
            self.nx_b = self.nx + 2
        else:
            self.nx_b = self.nx
        if "y" in self.flowdir:
            self.ny_b = self.ny + 2
        else:
            self.ny_b = self.ny
        if "z" in self.flowdir:
            self.nz_b = self.nz + 2
        else:
            self.nz_b = self.nz
        return (self.nx_b, self.ny_b, self.nz_b)

    @lru_cache(maxsize=1)
    def get_fshape(self):
        if self.flowdir == "-":
            self.fshape = (1,)
        elif self.flowdir == "x":
            self.fshape = (self.nx_b,)
        elif self.flowdir == "y":
            self.fshape = (self.ny_b,)
        elif self.flowdir == "z":
            self.fshape = (self.nz_b,)
        elif self.flowdir == "xy":
            self.fshape = (self.ny_b, self.nx_b)
        elif self.flowdir == "xz":
            self.fshape = (self.nz_b, self.nx_b)
        elif self.flowdir == "yz":
            self.fshape = (self.nz_b, self.ny_b)
        elif self.flowdir == "xyz":
            self.fshape = (self.nz_b, self.ny_b, self.nx_b)
        print(f"- flow shape (fshape) is set to {self.fshape}")
        # self.flow_shape = self.fshape
        """
        - in case of 2D but also with boundary conditions:
            elif self.flowdir == 'xz+':
                self.fshape = (self.nz_b, self.nx_b, 3)
            elif self.flowdir == 'yz+':
                self.fshape = (self.nz_b, self.ny_b, 3)
        """
        return self.fshape

    # get_flow_shape = get_fshape

    @lru_cache(maxsize=2)
    def get_n_cells(self, boundary=True):
        if boundary:
            self.n_cells = self.nx_b * self.ny_b * self.nz_b
        else:
            self.n_cells = self.nx * self.ny * self.nz
        return self.n_cells

    @lru_cache(maxsize=4)
    def get_order(self, type="natural", boundary=True, fshape=False, verbose=True):

        if type == "natural":
            self.order = np.arange(self.get_n_cells(boundary=True))
        else:
            raise ValueError(
                """Order type is not supported.\n
                            Supported order types: ['natural']"""
            )

        if fshape:
            self.order = self.order.reshape(self.fshape)

        if not boundary:
            self.order = self.remove_boundaries(self.order)

        if verbose:
            s1 = "with fshape" if fshape else "without fshape"
            s2 = "with boundary" if boundary else "without boundary"
            print(f"- Cells order (order) was computed ({s1} - {s2}).")
        return self.order

    def get_blocks(self):
        self.blocks = np.ones(self.fshape, dtype="int")
        # self.i_blocks = self.blocks.copy()
        # self.i_blocks[[0, -1]] = 0
        # self.b_blocks = np.zeros(self.max_nb, dtype='int') #ss.lil_matrix(self.shape, dtype='int')
        # self.b_blocks[[0, -1]] = 1

    """
    
    """

    def get_cells_dims(self, dx, dy, dz):
        m_list = []
        if "x" in self.flowdir:
            self.dxx = np.ones(self.nx_b, dtype="int") * dx
            m_list.append(self.dxx)
        else:
            self.dxx = np.ones(self.max_nb, dtype="int") * dx
            m_list.append(dx)
        if "y" in self.flowdir:
            self.dyy = np.ones(self.ny_b, dtype="int") * dy
            m_list.append(self.dyy)
        else:
            self.dyy = np.ones(self.max_nb, dtype="int") * dy
            m_list.append(dy)
        if "z" in self.flowdir:
            self.dzz = np.ones(self.nz_b, dtype="int") * dz
            m_list.append(self.dzz)
        else:
            self.dzz = np.ones(self.max_nb, dtype="int") * dz
            m_list.append(dz)

        self.dx, self.dy, self.dz = np.meshgrid(*m_list, copy=False)
        self.dx = np.transpose(self.dx, axes=(0, 2, 1)).reshape(self.fshape)
        self.dy = np.transpose(self.dy, axes=(2, 0, 1)).reshape(self.fshape)
        self.dz = np.transpose(self.dz, axes=(2, 1, 0)).reshape(self.fshape)

    @lru_cache(maxsize=None)
    def get_cell_id(self, coords, boundary=True):
        if boundary:
            return self.pyvista_grid_b.cell_id(coords)
        else:
            return self.pyvista_grid.cell_id(coords)

    @lru_cache(maxsize=4)
    def get_cells_id(self, boundary=True, fshape=True, verbose=True):
        self.cells_id = self.get_order(
            type="natural", boundary=boundary, fshape=fshape, verbose=False
        )
        if verbose:
            s1 = "with fshape" if fshape else "without fshape"
            s2 = "with boundary" if boundary else "without boundary"
            print(f"- Cells id (cells_id) was computed ({s1} - {s2}).")
        return self.cells_id

    # @lru_cache(maxsize=None): does not work
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

    @lru_cache(maxsize=2)
    def get_cells_coords(self, boundary=True, fshape=False, verbose=True):
        cells_id = self.get_cells_id(boundary=boundary, fshape=False)
        self.cells_coords = [
            tuple(x) for x in self.pyvista_grid_b.cell_coords(cells_id)
        ]

        if fshape:
            self.cells_coords = self.cells_coords.reshape(self.fshape)

        if verbose:
            s1 = "with fshape" if fshape else "without fshape"
            s2 = "with boundary" if boundary else "without boundary"
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
                coords = self.get_cell_coords(id)
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
        print(f"- is_homogeneous is set to {self.is_homogeneous}.")
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
        if self.D == 1:
            return (dx[:-1] / (area[:-1] * k[:-1])) + (dx[1:] / (area[1:] * k[1:]))
        elif self.D == 2:
            return (dx[:-1, :-1] / (area[:-1, :-1] * k[:-1, :-1])) + (
                dx[1:, 1:] / (area[1:, 1:] * k[1:, 1:])
            )
        elif self.D == 3:
            return (dx[:-1, :-1, :-1] / (area[:-1, :-1, :-1] * k[:-1, :-1, :-1])) + (
                dx[1:, 1:, 1:] / (area[1:, 1:, 1:] * k[1:, 1:, 1:])
            )

    def get_G_homo_mean(self, property, type="geometric"):
        """ """
        if self.is_homogeneous:
            if self.D == 1:
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

    def get_pyvista_grid(self, boundary=True, verbose=False):
        """
        https://docs.pyvista.org/api/core/_autosummary/pyvista.ExplicitStructuredGrid.html
        """
        if boundary:
            dims = np.array((self.nx_b, self.ny_b, self.nz_b)) + 1
        else:
            dims = np.array((self.nx, self.ny, self.nz)) + 1

        corners = self.get_corners(boundary, verbose)
        pyvista_grid = pv.ExplicitStructuredGrid(dims, corners)

        s = "with boundary" if boundary else "without boundary"
        print(f"- pv_grid {s} was created.")
        if verbose:
            print(pyvista_grid)
        return pyvista_grid

    def get_corners(self, boundary=True, verbose=True):
        """

        (Reference: https://docs.pyvista.org/examples/00-load/create-explicit-structured-grid.html)
        """

        if "x" in self.flowdir:
            xcorn = np.insert(self.dxx.cumsum(), 0, 0)
        else:
            xcorn = np.arange(0, (self.nx + 1) * self.dxx[0], self.dxx[0])

        if "y" in self.flowdir:
            ycorn = np.insert(self.dyy.cumsum(), 0, 0)
        else:
            ycorn = np.arange(0, (self.ny + 1) * self.dyy[0], self.dyy[0])

        if "z" in self.flowdir:
            zcorn = np.insert(self.dzz.cumsum(), 0, 0)
        else:
            zcorn = np.arange(0, (self.nz + 1) * self.dzz[0], self.dzz[0])

        # Boundary:
        if boundary:
            ix = 2 if "x" in self.flowdir else 0
            iy = 2 if "y" in self.flowdir else 0
            iz = 2 if "z" in self.flowdir else 0
        else:
            ix = 0
            iy = 0
            iz = 0
            if "x" in self.flowdir:
                xcorn = xcorn[1:-1]
            if "y" in self.flowdir:
                ycorn = ycorn[1:-1]
            if "z" in self.flowdir:
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
            print(xcorn.shape, ycorn.shape, zcorn.shape)

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
        self.cells_center = self.pyvista_grid_b.cell_centers().points

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
        if boundary:
            self.V = self.pyvista_grid_b.volume
        else:
            self.V = self.pyvista_grid.volume
        return self.V

    @lru_cache(maxsize=6)
    def get_cells_volume(self, boundary=True, fshape=False, pyvista=False):

        if pyvista:
            self.cells_V = self.pyvista_grid_b.compute_cell_sizes()["Volume"]
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
        cells_volume = self.get_cells_volume(boundary=True)
        if id is not None:
            return cells_volume.flatten()[id]
        elif coords is not None:
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
        cells_d = self.get_cells_d(dir=dir, boundary=True)
        if id is not None:
            return cells_d.flatten()[id]
        elif coords is not None:
            return cells_d[coords[2], coords[1], coords[0]]
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
                cell_neighbors[self.flowdir[0]] = neighbors
            if self.D >= 2:
                neighbors = [
                    n for n in [id - self.max_nb, id + self.max_nb] if n in cells_id
                ]
                cell_neighbors[self.flowdir[1]] = neighbors
            if self.D >= 3:
                neighbors = [
                    n
                    for n in [
                        id - (self.nx_b * self.ny_b),
                        id + (self.nx_b * self.ny_b),
                    ]
                    if n in cells_id
                ]
                cell_neighbors[self.flowdir[2]] = neighbors
        elif coords is not None:
            cells_coords = self.get_cells_coords(boundary=boundary, fshape=False)
            assert (
                coords in cells_coords
            ), f"cell coords are out of range {cells_coords}."
            i, j, k = coords
            if "x" in self.flowdir:
                neighbors = [
                    n for n in [(i - 1, j, k), (i + 1, j, k)] if n in cells_coords
                ]
                cell_neighbors["x"] = neighbors
            if "y" in self.flowdir:
                neighbors = [
                    n for n in [(i, j - 1, k), (i, j + 1, k)] if n in cells_coords
                ]
                cell_neighbors["y"] = neighbors
            if "z" in self.flowdir:
                neighbors = [
                    n for n in [(i, j, k - 1), (i, j, k + 1)] if n in cells_coords
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

    def get_boundaries(self, ids=True):
        """
        Arguments:
            - i: index in x-direction.
        """
        if not hasattr(self, "boundaries_id"):
            cells_id_b = self.get_cells_id(boundary=True)  # or self.cells_id
            if self.D == 0:
                self.boundaries_id = cells_id_b.flatten()
            elif self.D == 1:
                self.boundaries_id = cells_id_b[[0, -1]].flatten()
            elif self.D == 2:
                self.boundaries_id = np.sort(
                    np.concatenate(
                        [
                            cells_id_b[[0, -1], :].flatten(),
                            cells_id_b[1:-1, [0, -1]].flatten(),
                        ]
                    )
                )
            elif self.D == 3:
                self.boundaries_id = np.sort(
                    np.concatenate(
                        [
                            cells_id_b[:, [0, -1], :].flatten(),
                            cells_id_b[:, 1:-1, [0, -1]].flatten(),
                            cells_id_b[[0, -1], 1:-1, 1:-1].flatten(),
                        ]
                    )
                )
        if not hasattr(self, "boundaries_coords"):
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

        if corners:
            points_lst = pv_grid.points
            mask = points_lst[:, 1] == 0
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
        title = f"{label} for {self.D}D model (flow at {self.flowdir}-direction {s})"
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
        nx=2,
        ny=2,
        nz=1,
        dx=dx,
        dy=dy,
        dz=dz,
        kx=270,
        ky=10,
        kz=20,
        phi=0.27,
    )

    print(grid.cells_V.shape)
    V = grid.cells_id.reshape((1, 4, 4))
    print(V)
    print(grid.cells_V.shape)
    # grid.set_kx(kx=10, coords=(1,1,1))
    # print(grid.kx)
    # print(grid.get_cells_id())
    # print(grid.kx.shape)
    # print(grid.cells_volume_b.shape)
    # grid.show(boundary=True, label="center")

    # grid.show(boundary=False, label="id", corners=True)
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
