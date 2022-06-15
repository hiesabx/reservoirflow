"""
Grid classes for reservoir simulation models.

This module contains all grid classes that are required to build the 
Model class. Grid class represents both the rock geometry and the rock 
properties that are required for fluid flow calculations.

"""
import warnings
from openresim.base import Base
import numpy as np
import scipy.sparse as ss
import pyvista as pv
from openresim.utils import _lru_cache
from openresim import utils


class Grid(Base):
    """Create a Grid class.

    Grid class holds grid and rock properties using numpy arrays
    including pyvista object for visualization.

    Parameters
    ----------
    Base : class
        Base class with universal settings.
    """

    def __init__(self, dtype, unit, unify, verbose):
        super().__init__(unit)
        self.dtype = dtype  # np.single, np.double
        self.unify = unify
        self.verbose = verbose


class CartGrid(Grid):
    """Cartesian grid class with explicit structure. Parameters can be
    defined as `unit='field'` (default) or `unit='metric'`. `units`
    can be accessed from this class or base class `CartGrid.units` or
    `Grid.units`.

    Parameters
    ----------
    Grid : class
        parent grid class with universal grid settings.

    Returns
    -------
    CartGrid
        cartesian grid object.

    Raises
    ------
    ValueError
        _description_
    """

    name = "CartGrid"
    type = "cartesian"

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
        z=0,
        comp=None,
        dtype="double",
        unit="field",
        unify=True,
        verbose=False,
    ):
        """Return cartesian grid object with explicit structure.
        Parameters can be defined as `unit='field'` (by default) or
        `unit='metric'`. `units` can be accessed from this class with
        `CartGrid.units` or from or the base class with `Grid.units`.

        Parameters
        ----------
        nx : int
            number of grids in x-direction and must be >= 1.
        ny : int
            number of grids in y-direction and must be >= 1.
        nz : int
            number of grids in z-direction and must be >= 1.
        dx : int, float, array-like
            grid dimension in x-direction. In case of a list or array,
            the length should be equal to nx+2 for all cells including
            boundary cells. Vales should be in natural order (i.e. from
            left to right).
        dy : int, float, array-like
            grid dimension in y-direction. In case of a list or array,
            the length should be equal to ny+2 for all cells including
            boundary cells. Vales should be in natural order (i.g. from
            front to back).
        dz : int, float, array-like
            grid dimension in z-direction. In case of a list or array,
            the length should be equal to nz+2 for all cells including
            boundary cells. Vales should be in natural order (i.g. from
            down to up).
        kx : int, float, array-like, optional, by default None
            permeability in x-direction (relevant only if 'x' was in
            fluid flow direction). In case of a list or array,
            the length should be equal to ny+2 for all cells including
            boundary cells. Vales should be in natural order (i.g. from
            left to right).
        ky : int, float, array-like, optional, by default None
            permeability in x-direction (relevant only if 'y' was in
            fluid flow direction). In case of a list or array,
            the length should be equal to ny+2 for all cells including
            boundary cells. Vales should be in natural order (i.g.from
            front to back).
        kz : int, float, array-like, optional, by default None
            permeability in x-direction (relevant only if 'z' was in
            fluid flow direction). In case of a list or array,
            the length should be equal to nz+2 for all cells including
            boundary cells. Vales should be in natural order (i.g. from
            down to up).
        phi : float, array-like, optional, by default None
            porosity. In case of an array, the shape should be equal to
            grid.shape with boundaries. Vales should be in natural
            order.
        z : int, float, array-like, optional, by default 0.
            depth of grid tops (NOT FULLY IMPLEMENTED).
        comp : float, optional, by default None
            compressibility.
        dtype : str or `np.dtype`, optional, by default 'double'
            data type used in all arrays.
        unit : str ('field', 'metric'), optional, by default 'field'
            units used in input and output.
        unify : bool, optional, by default False
            unify shape to be always tuple of 3 when set to True. When
            set to False, shape includes only the number of girds in
            flow direction as tuple. This option is only relevant in
            case of 1D or 2D flow. This option may be required to make
            1D and 2D shapes shapes of this class more consistent with
            each other or with 3D shape. Warning: True option is not
            yet fully compatible.
        verbose : bool, optional, by default False
            print information for debugging.

        ToDo
        ----
        Complete unify feature to be compatible with all class
        components.
        """
        super().__init__(dtype, unit, unify, verbose)
        assert nx >= 1, "nx must be 1 or larger."
        assert ny >= 1, "ny must be 1 or larger."
        assert nz >= 1, "nz must be 1 or larger."
        self.nx, self.ny, self.nz = nx, ny, nz
        # self.get_D()  # > self.D
        # self.get_shape()  # > self.shape
        # self.get_fdir()  # > self.flowdir
        # self.get_n(True)  # > self.nx_b, self.ny_b, self.n_zb
        # self.get_fshape(True, False, False)  # > self.fshape
        # self.get_n_cells()  # > self.n_cells
        # self.get_order(type="natural")  # > self.order
        # self.get_ones(True, False)  # > self.blocks
        # self.get_n_max(True)  #
        self.get_cells_dims(dx, dy, dz)
        self.get_cells_area_x()
        self.get_cells_area_y()
        self.get_cells_area_z()
        self.get_cells_volume()
        self.get_cells_coords()  # > cells_coords
        self.set_props(phi, kx, ky, kz, z, comp)
        self.get_boundaries()  # > self.boundaries_id, self.boundaries_coords
        self.get_volume()  # self.volume
        # self.get_Gx()  # > self.Gx
        # self.get_Gy()  # > self.Gy
        # self.get_Gz()  # > self.Gz
        self.get_cells_center()

    # -------------------------------------------------------------------------
    # Basic:
    # -------------------------------------------------------------------------

    @_lru_cache(maxsize=1)
    def get_D(self):
        """Return the grid dimension as int.

        Returns
        -------
        int
            number of dimensions higher than 1.
        """
        self.D = sum([1 if n > 1 else 0 for n in (self.nx, self.ny, self.nz)])

        if self.verbose:
            print(f"- Dimension (D) is set to {self.D}.")

        return self.D

    @_lru_cache(maxsize=2)
    def get_shape(self, boundary=False):
        """Return the grid shape in x, y, z as ndarray.

        Parameters
        ----------
        boundary : bool, optional, by default False
            values with boundary (True) or without boundary (False).

        Returns
        -------
        ndarray
            number of cells as np.array([nx, ny, nz]).
        """
        n = self.get_n(boundary)
        self.shape = np.array(n, dtype="int")

        if self.verbose:
            s = utils.get_boundary_str(boundary)
            print(f"- Shape {s} is set to {self.shape}.")

        return self.shape

    @_lru_cache(maxsize=1)
    def get_fdir(self):
        """Return the flow direction as str.

        Returns
        -------
        str
            contains one or more of ('-','x','y','z') based on the grid
            dimensions that are higher than 1.
        """
        self.get_D()
        self.get_shape(False)

        if self.D == 0:
            self.fdir = "-"
        elif self.D == 1:
            flow_dir_id = np.argmax(self.shape)
            if flow_dir_id == 0:
                self.fdir = "x"
            elif flow_dir_id == 1:
                self.fdir = "y"
            elif flow_dir_id == 2:
                self.fdir = "z"
        elif self.D == 2:
            flow_dir_id = np.argmin(self.shape)
            if flow_dir_id == 2:
                self.fdir = "xy"
            elif flow_dir_id == 1:
                self.fdir = "xz"
            elif flow_dir_id == 0:
                self.fdir = "yz"
        elif self.D == 3:
            self.fdir = "xyz"

        if self.verbose:
            print(f"- Flow direction (fdir) is set to {self.fdir}.")

        return self.fdir

    @_lru_cache(maxsize=2)
    def get_n(self, boundary=False):
        """Return the number of grids in x, y, and z as tuple.

        Parameters
        ----------
        boundary : bool, optional, by default False
            values with boundary (True) or without boundary (False).

        Returns
        -------
        tuple
            the number of grids as (nx, ny, nz).
        """
        if boundary:
            self.get_fdir()
            if "x" in self.fdir:
                self.nx_b = self.nx + 2
            else:
                self.nx_b = self.nx
            if "y" in self.fdir:
                self.ny_b = self.ny + 2
            else:
                self.ny_b = self.ny
            if "z" in self.fdir:
                self.nz_b = self.nz + 2
            else:
                self.nz_b = self.nz
            self.n = (self.nx_b, self.ny_b, self.nz_b)
        else:
            self.n = (self.nx, self.ny, self.nz)

        if self.verbose:
            s = utils.get_boundary_str(boundary)
            print(f"- N-cells {s} is {self.n}.")

        return self.n

    @_lru_cache(maxsize=2)
    def get_n_max(self, boundary=True):
        """Return the maximum number of grids as int.

        Parameters
        ----------
        boundary : bool, optional, by default True
            values with boundary (True) or without boundary (False).

        Returns
        -------
        int
            maximum number of grids as max(nx, ny, nz).
        """
        n = self.get_n(boundary)
        self.n_max = max(n)

        if self.verbose:
            s = utils.get_boundary_str(boundary)
            print(f"- Maximum number of cells {s} is {self.n_max}.")

        return self.n_max

    @_lru_cache(maxsize=2)
    def get_fshape(self, boundary=True, points=False, unify=False):
        """Return flow shape as tuple.

        Parameters
        ----------
        boundary : bool, optional, by default True
            values with boundary (True) or without boundary (False).
            Warning: True option is required.
        points : bool, optional, by default False
            True for points (i.e. tuples of len 3 like coords, icoords)
            and False for scaler values (e.g. id).
        unify : bool, optional, by default False
            unify shape to be always tuple of 3 as (nz,ny,nx) when set
            to True. When set to False, shape includes only the number
            of girds in flow direction as tuple. This option is only
            relevant in case of 1D or 2D flow. Warning: True option is
            not yet fully compatible.

        Returns
        -------
        tuple
            number of grids at flow directions in z, y, x order.

        ToDo
        ----
        add usage with boundary=False.
        add optional boundaries for 2D models:
            if self.flowdir == 'xz+':
                self.fshape = (self.nz_b, self.nx_b, 3)
            if self.flowdir == 'yz+':
                self.fshape = (self.nz_b, self.ny_b, 3)
        """
        assert boundary == True, "False boundary is not permitted."
        utils.fshape_warn(self.unify, unify)

        nx, ny, nz = self.get_n(boundary)  # Includes self.get_fdir()

        if self.fdir == "xyz":
            self.fshape = (nz, ny, nx)
        else:
            if not self.unify:
                if self.fdir == "-":
                    self.fshape = (1,)
                elif self.fdir == "x":
                    self.fshape = (nx,)
                elif self.fdir == "y":
                    self.fshape = (ny,)
                elif self.fdir == "z":
                    self.fshape = (nz,)
                elif self.fdir == "xy":
                    self.fshape = (ny, nx)
                elif self.fdir == "xz":
                    self.fshape = (nz, nx)
                elif self.fdir == "yz":
                    self.fshape = (nz, ny)
                else:
                    raise ValueError("unknown fdir value.")
            else:
                if self.fdir == "-":
                    self.fshape = (1, 1, 1)
                elif self.fdir == "x":
                    self.fshape = (1, 1, nx)
                elif self.fdir == "y":
                    self.fshape = (1, ny, 1)
                elif self.fdir == "z":
                    self.fshape = (nz, 1, 1)
                elif self.fdir == "xy":
                    self.fshape = (1, ny, nx)
                elif self.fdir == "xz":
                    self.fshape = (nz, 1, nx)
                elif self.fdir == "yz":
                    self.fshape = (nz, ny, 1)
                else:
                    raise ValueError("unknown fdir value.")

        if points:
            self.fshape = self.fshape + (3,)

        if self.verbose:
            print(f"- Flow shape (fshape) is set to {self.fshape}.")

        return self.fshape

    @_lru_cache(maxsize=2)
    def get_n_cells(self, boundary=True):
        """Return total number of grid cells as int.

        Parameters
        ----------
        boundary : bool, optional, by default True
            values with boundary (True) or without boundary (False).

        Returns
        -------
        int
            total number of cells.
        """
        n = self.get_n(boundary)
        self.n_cells = np.prod(n)

        if self.verbose:
            s = utils.get_boundary_str(boundary)
            print(f"- Number of cells {s} is {self.n_cells}.")

        return self.n_cells

    @_lru_cache(maxsize=4)
    def get_order(self, type="natural", boundary=True, fshape=False):
        """Return grid order as ndarray.

        Parameters
        ----------
        type : str, optional, by default "natural"
            order type in which grids are numbered.
        boundary : bool, optional, by default True
            values with boundary (True) or without boundary (False).
        fshape : bool, optional, by default False
            values in flow shape (True) or flatten (False).

        Returns
        -------
        ndarray
            gird order as array.
        """
        if type == "natural":
            self.order = np.arange(self.get_n_cells(True))
        else:
            raise ValueError(
                "Order type is not supported or unknown. "
                "Supported order types: ['natural']"
            )

        if fshape:
            shape = self.get_fshape(boundary, False, False)
            self.order = self.order.reshape(shape)

        if not boundary:
            self.order = self.remove_boundaries(self.order)

        if self.verbose:
            s1, s2 = self.get_verbose_str(boundary, fshape)
            print(f"- Cells order (order) was computed ({s1} - {s2}).")

        return self.order

    @_lru_cache(maxsize=1)
    def get_ones(self, boundary=True, sparse=False):
        """Return array of ones in flow shape.

        Parameters
        ----------
        boundary : bool, optional, by default True
            values with boundary (True) or without boundary (False).
        sparse : bool, optional, by default False
            values as sparse matrix (True) or as ndarray (False).

        Returns
        -------
        array
            array in flow shape filled with ones.
        """

        assert boundary == True, "False boundary option is not allowed."
        assert sparse == False, "True sparse option is not allowed."

        fshape = self.get_fshape(boundary, False, False)

        if not sparse:
            self.ones = np.ones(fshape, dtype="int")
        else:
            self.ones = ss.lil_matrix(fshape, dtype="int")

        if self.verbose:
            print(f"- Ones array (ones) was computed.")

        return self.ones

    # -------------------------------------------------------------------------
    # Cells id and coordinates:
    # -------------------------------------------------------------------------

    @_lru_cache(maxsize=None)
    def get_cell_id(self, coords=[], boundary=True):
        """Return cell/cells id based on natural as int/list.

        Parameters
        ----------
        coords : tuple of int, tuple of tuples of int
            cell coordinates (i,j,k) as a tuple of int. For multiple
            cells, tuple of tuples of int as ((i,j,k),(i,j,k),..).
        boundary : bool, optional, by default True
            values with boundary (True) or without boundary (False).

        Returns
        -------
        int/list
            cell id based on natural order as int for a single cell
            coords or as list of int for multiple cells coords.
        """
        pyvista_grid = self.get_pyvista_grid(boundary)

        if all(isinstance(c, tuple) for c in coords):
            return [pyvista_grid.cell_id(c) for c in coords]
        else:
            return pyvista_grid.cell_id(coords)

    @_lru_cache(maxsize=4)
    def get_cells_id(self, boundary=True, fshape=False, fmt="tuple"):
        """Return all cells id based on natural order as ndarray.

        Parameters
        ----------
        boundary : bool, optional, by default False
            values with boundary (True) or without boundary (False).
        fshape : bool, optional, by default True
            values in flow shape (True) or flatten (False). If set to
            True, fmt argument will be ignored.
        fmt : str, optional, by default "tuple"
            output format as str from ['array', 'list', 'tuple', 'set'].
            This argument is ignored if fshape argument is set to True.
            For a better performance, use 'set' to check if an item is
            in a list or not. Use tuples to iterate through items. When
            option 'array' is used, utils.isin() must be used to check
            if a tuple of 3 is in the array.

        Returns
        -------
        ndarray
            cells id in natural order as array.

        ToDo
        ----
        can be a generator instead.
        """
        self.cells_id = self.get_order("natural", boundary, fshape)

        if not fshape:
            self.cells_id = utils.reformat(self.cells_id, fmt)

        if self.verbose:
            s1, s2 = self.get_verbose_str(boundary, fshape)
            print(f"- Cells id (cells_id) was computed ({s1} - {s2}).")

        return self.cells_id

    @_lru_cache(maxsize=None)
    def get_cell_coords(self, id, boundary=True):
        """Return cell/cells coordinates as tuple/list of tuples.

        Parameters
        ----------
        id : int, tuple of int
            cell id based on natural order as int. For multiple cells,
            list of int [id,id,..] or tuple of int (id,id,...). Warning:
            passing list or arrays instead of tuples does not work with
            cache decorator used in this method since lists and ndarray
            are both unhashable.
        boundary : bool, optional, by default True
            values in flow shape (True) or flatten (False).

        Returns
        -------
        tuple/list of tuples
            cell/cells coordinates as tuple/list of tuples.
        """
        pyvista_grid = self.get_pyvista_grid(boundary)

        if isinstance(id, (list, tuple, np.ndarray)):
            return [tuple(x) for x in pyvista_grid.cell_coords(id)]
        else:
            return tuple(pyvista_grid.cell_coords(id))

    @_lru_cache(maxsize=4)
    def get_cells_coords(self, boundary=True, fshape=False, fmt="tuple"):
        """Return all cells coords based on (i,j,k) as ndarray.

        Parameters
        ----------
        boundary : bool, optional, by default True
            values with boundary (True) or without boundary (False).
        fshape : bool, optional, by default False
            values in flow shape (True) or flatten (False). If set to
            True, fmt argument will be ignored.
        fmt : str, optional, by default 'tuple'
            output format as str from ['array', 'list', 'tuple', 'set'].
            This argument is ignored if fshape argument is set to True.
            For a better performance, use 'set' to check if an item is
            in a list or not. Use tuples to iterate through items. When
            option 'array' is used, utils.isin() must be used to check
            if a tuple of 3 is in the array.

        Returns
        -------
        ndarray
            cells coords in (i,j,k) as array.

        ToDo
        ----
        can be a generator instead.
        """
        cells_id = self.get_cells_id(boundary, False, "array")
        pyvista_grid = self.get_pyvista_grid(True)
        self.cells_coords = pyvista_grid.cell_coords(cells_id)

        if fshape:
            coords_fshape = self.get_fshape(boundary, True, False)
            self.cells_coords = self.cells_coords.reshape(coords_fshape)
        else:
            self.cells_coords = utils.reformat(self.cells_coords, fmt)

        if self.verbose:
            s1, s2 = self.get_verbose_str(boundary, fshape)
            print(f"- Cells coords (cells_coords) was computed ({s1} - {s2}).")

        return self.cells_coords

    @_lru_cache(maxsize=None)
    def get_cell_icoords(self, coords, unify=False):
        """Convert `coords` from `(i,j,k)` into `(k,j,i)`.

        This method is required to create `icoords` based on `(k,j,i)`
        which can be used to access ndarrays in this class. icoords is
        not compatible with pyvista grid which use `coords` based on
        `(i,j,k)`.

        Parameters
        ----------
        coords : tuple of int, tuple of tuples of int
            cell coordinates (i,j,k) as a tuple of int. For
            multiple cells, tuple of tuples of int as
            ((i,j,k),(i,j,k),..).
        unify : bool, optional, by default False
            unify shape to be always tuple of 3 as (k,j,i) when set to
            True. When set to False, shape includes only the number
            of girds in flow direction as tuple. This option is only
            relevant in case of 1D or 2D flow. Warning: True option is
            not yet fully compatible.

        Returns
        -------
        tuple/list of tuples
            internal coords (icoords)

        ToDo
        ----
        add tuple of tuples for coords as input.
        set unify as None.
        """
        cells_coords = self.get_cells_coords(True, False)
        shape_bool = cells_coords.shape == (self.get_n_cells(True), 3)
        assert shape_bool, "coords should include boundary and be flatten."
        assert utils.isin(coords, cells_coords), "coords are out of range."
        utils.fshape_warn(self.unify, unify)

        if not self.unify and self.D <= 2:
            icoords = tuple(c for c in coords[::-1] if c > 0)
            assert len(icoords) == self.get_D(), "icoords is not compatible"
        else:
            icoords = tuple(c for c in coords[::-1])

        return icoords

    # -------------------------------------------------------------------------
    # Neighbors and Boundaries:
    # -------------------------------------------------------------------------

    @_lru_cache(maxsize=None)
    def get_cell_neighbors(
        self,
        id=None,
        coords=None,
        boundary=False,
        fmt="dict",
    ):
        """Return cell neighbors.

        This method returns cell neighbors by id or coords. If
        neighbors are desired by id, then id argument should be used.
        The same applies coords argument. This method will raise
        ValueError if none of id or coords arguments were defined or if
        undefined fmt argument was used. Boundary cells are not allowed.

        Warning: passing ndarray of len(shape) > 1 (e.g. coords ndarray)
        causes a TypeError due to the cache decorator used in this
        method since multi-dim ndarray is unhashable.

        Parameters
        ----------
        id : int, iterable of int, by default None
            cell id based on natural order as int. For multiple cells,
            list of int [id,id,..] or tuple of int (id,id,...).
        coords : iterable of int, iterable of tuples of int, by default
            None cell coordinates (i,j,k) as a tuple of int. For
            multiple cells, tuple of tuples of int as
            ((i,j,k),(i,j,k),..).
        boundary : bool, optional, by default True
            values in flow shape (True) or flatten (False).
        fmt : str, optional, by default 'dict'
            output format as str from ['list', 'tuple', 'set', 'dict'].
            Use 'dict' to output neighbors in x,y,z directions as keys.
            Use 'tuple' or 'list' for list of neighbors when directions
            are not needed.

        Returns
        -------
        iterable
            cell neighbors.

        Raises
        ------
        ValueError
            None of id or coords arguments are not defined.
        ValueError
            Unknown str value was used for fmt argument.

        ToDo
        ----
        nothing to do.
        """
        cell_neighbors = {"x": [], "y": [], "z": []}

        if id is not None:
            assert not isinstance(id, np.ndarray), "block"
            boundaries = self.get_boundaries("id", "set")
            isin_boundary = utils.isin(id, boundaries)
            assert not isin_boundary, "boundary cells are not allowed."
            cells_id = self.get_cells_id(boundary, False, "set")
            isin_cells_id = utils.isin(id, cells_id)
            assert isin_cells_id, f"id is out of range {cells_id}."
            if self.D >= 1:
                n_lst = [id - 1, id + 1]
                neighbors = [i for i in n_lst if i in cells_id]
                cell_neighbors[self.fdir[0]] = neighbors
            if self.D >= 2:
                nx, ny, _ = self.get_n(True)
                if "x" in self.fdir:
                    n_lst = [id - nx, id + nx]
                elif "y" in self.fdir:
                    n_lst = [id - ny, id + ny]
                neighbors = [i for i in n_lst if i in cells_id]
                cell_neighbors[self.fdir[1]] = neighbors
            if self.D >= 3:
                nx_ny_b = self.nx_b * self.ny_b
                n_lst = [id - nx_ny_b, id + nx_ny_b]
                neighbors = [i for i in n_lst if i in cells_id]
                cell_neighbors[self.fdir[2]] = neighbors
        elif coords is not None:
            boundaries = self.get_boundaries("coords", "set")
            isin_boundary = utils.isin(coords, boundaries)
            assert not isin_boundary, "boundary cells are not allowed."
            cells_coords = self.get_cells_coords(boundary, False, "set")
            isin_cells_coords = utils.isin(coords, cells_coords)
            assert isin_cells_coords, f"coords are out of range {cells_coords}."
            i, j, k = coords
            if "x" in self.fdir:
                n_lst = [(i - 1, j, k), (i + 1, j, k)]
                neighbors = [c for c in n_lst if utils.isin(c, cells_coords)]
                cell_neighbors["x"] = neighbors
            if "y" in self.fdir:
                n_lst = [(i, j - 1, k), (i, j + 1, k)]
                neighbors = [c for c in n_lst if utils.isin(c, cells_coords)]
                cell_neighbors["y"] = neighbors
            if "z" in self.fdir:
                n_lst = [(i, j, k - 1), (i, j, k + 1)]
                neighbors = [c for c in n_lst if utils.isin(c, cells_coords)]
                cell_neighbors["z"] = neighbors
        else:
            raise ValueError("at least id or coords argument must be defined.")

        return utils.reformat(cell_neighbors, fmt=fmt)

    def remove_boundaries(self, in_data, points: bool = None):
        """Remove boundary cells from ndarray.

        Parameters
        ----------
        in_data : ndarray, dict of ndarray
            input data where boundaries need to be removed. Input data
            must be an ndarray with boundaries. Input data as dict with
            keys for these arrays is also possible.
        points : bool, optional, by default None
            True for points (i.e. tuples of len 3 like coords, icoords)
            and False for scaler values (e.g. id). If value is set to
            None, bool value is calculated automatically. Warning:
            this argument must be specified in case that in_data was for
            scaler values in fshape that is (#,..,3) (i.e. not flatten).
            For more information about points automatic calculation,
            check the utility function `utils.ispoints()`.

        Returns
        -------
        ndarray, dict
            array with boundaries removed.

        Raises
        ------
        ValueError
            boundaries are not included or points argument must be
            correctly assigned.
        ValueError
            dtype must be ndarray.

        See Also
        --------
        extract_boundaries: keep only boundary cells from input data.

        ToDo
        ----
        add fmt argument.
        """
        if isinstance(in_data, np.ndarray):
            if points is None:
                points = utils.ispoints(in_data)
            fshape = self.get_fshape(True, points, False)

            if in_data.shape != fshape:
                try:
                    in_data = in_data.reshape(fshape)
                    flatten = True
                except:
                    utils.shape_error(in_data.shape, fshape)
            else:
                flatten = False

            if self.D == 3:
                out_data = in_data[1:-1, 1:-1, 1:-1]
            else:
                if not self.unify:
                    if self.D == 0:
                        out_data = in_data
                    elif self.D == 1:
                        out_data = in_data[1:-1]
                    elif self.D == 2:
                        out_data = in_data[1:-1, 1:-1]
                    else:
                        raise ValueError("Unknown shape.")
                else:
                    fdir = self.get_fdir()
                    if fdir == "-":
                        out_data = in_data
                    elif fdir == "x":
                        out_data = in_data[:, :, 1:-1]
                    elif fdir == "y":
                        out_data = in_data[:, 1:-1, :]
                    elif fdir == "z":
                        out_data = in_data[1:-1, :, :]
                    elif fdir == "xy":
                        out_data = in_data[:, 1:-1, 1:-1]
                    elif fdir == "xz":
                        out_data = in_data[1:-1, :, 1:-1]
                    elif fdir == "yz":
                        out_data = in_data[1:-1, 1:-1, :]
                    else:
                        raise ValueError("Unknown shape.")

            if flatten:
                if not points:
                    out_data = out_data.flatten()
                else:
                    out_data = out_data.reshape((-1, 3))

            return out_data

        elif isinstance(in_data, dict):
            for k, v in in_data.items():
                in_data[k] = self.remove_boundaries(v)
            return in_data
        else:
            raise ValueError("dtype must be ndarray.")

    def extract_boundaries(self, in_data, points=None, fmt="tuple"):
        """Extract boundary cells from ndarrays.

        Parameters
        ----------
        in_data : ndarray
            input array must contain all cells including boundary cell.
        points : bool, optional, by default None
            True for points (i.e. tuples of len 3 like coords, icoords)
            and False for scaler values (e.g. id). If value is set to
            None, bool value is calculated automatically. Warning:
            this argument must be specified in case that in_data was for
            scaler values in fshape that is (#,..,3) (i.e. not flatten).
            For more information about points automatic calculation,
            check the utility function `utils.ispoints()`.
        fmt : str, optional, by default "tuple"
            format of output data as str in ['tuple', 'list', 'set',
            'array'].

        Returns
        -------
        ndarray, list, set
            output data based on fmt argument.

        Raises
        ------
        ValueError
            'fmt is unknown' when fmt is not in ['tuple','list','array']
        ValueError
            'dtype must be ndarray' when in_data is not numpy array.

        See Also
        --------
        remove_boundaries: remove boundary cells from input data.

        ToDo
        ----
        The fshape might be checked automatically based on the provided data.
        Confirm the behavior of when self.unify set to True.
        """
        if isinstance(in_data, np.ndarray):
            if points is None:
                points = utils.ispoints(in_data)
            fshape = self.get_fshape(True, points, False)

            if in_data.shape != fshape:
                try:
                    in_data = in_data.reshape(fshape)
                except:
                    utils.shape_error(in_data.shape, fshape)

            if self.D == 3:
                out_data = np.concatenate(
                    [
                        in_data[:, [0, -1], :].flatten(),
                        in_data[:, 1:-1, [0, -1]].flatten(),
                        in_data[[0, -1], 1:-1, 1:-1].flatten(),
                    ]
                )
            else:
                if not self.unify:
                    if self.D == 0:
                        out_data = in_data
                    elif self.D == 1:
                        out_data = in_data[[0, -1]].flatten()
                    elif self.D == 2:
                        out_data = np.concatenate(
                            [
                                in_data[[0, -1], :].flatten(),
                                in_data[1:-1, [0, -1]].flatten(),
                            ]
                        )
                    else:
                        raise ValueError("unknown shape.")
                else:
                    fdir = self.get_fdir()
                    if fdir == "-":
                        out_data = in_data
                    elif fdir == "x":
                        out_data = in_data[:, :, [0, -1]]
                    elif fdir == "y":
                        out_data = in_data[:, [0, -1], :]
                    elif fdir == "z":
                        out_data = in_data[[0, -1], :, :]
                    elif fdir == "xy":
                        out_data = np.concatenate(
                            [
                                in_data[:, [0, -1], :].flatten(),
                                in_data[:, 1:-1, [0, -1]].flatten(),
                            ]
                        )
                    elif fdir == "xz":
                        out_data = np.concatenate(
                            [
                                in_data[[0, -1], :, :].flatten(),
                                in_data[1:-1, :, [0, -1]].flatten(),
                            ]
                        )
                    elif fdir == "yz":
                        out_data = np.concatenate(
                            [
                                in_data[[0, -1], :, :].flatten(),
                                in_data[1:-1, [0, -1], :].flatten(),
                            ]
                        )
                    else:
                        raise ValueError("unknown shape.")

            if not points:
                out_data = np.sort(out_data.flatten())
            else:
                out_data = out_data.reshape((-1, 3))

            return utils.reformat(out_data, fmt)
        else:
            raise ValueError("dtype must be ndarray.")

    @_lru_cache(maxsize=2)
    def get_boundaries(self, by="id", fmt="tuple"):
        """Return all boundary cells by id or coords.

        Parameters
        ----------
        by : str, optional, by default 'id'
            output boundaries as 'id' or 'coords'. Other undefined str
            values will raise ValueError.
        fmt : str, optional, by default "tuple"
            format of output data as str in ['tuple', 'list', 'set',
            'array']. When option 'array' is used, utils.isin() must be
            used to check if a tuple of 3 is in the array. For a better
            performance, use 'set' to check if an item is in or not and
            use tuples to iterate through items.

        Returns
        -------
        ndarray, list, set
            boundaries by id or coords based on fmt argument.

        Raises
        ------
        ValueError
            by argument must be either 'id' or 'coords'.
        """
        if by == "id":
            cells_id = self.get_cells_id(True, True, "array")
            self.boundaries_id = self.extract_boundaries(cells_id, False, fmt)
            return self.boundaries_id
        elif by == "coords":
            cells_coords = self.get_cells_coords(True, True, "array")
            self.boundaries_coords = self.extract_boundaries(cells_coords, True, fmt)
            return self.boundaries_coords
        else:
            raise ValueError("'by' argument must be either 'id' or 'coords'.")

    @_lru_cache(maxsize=None)
    def get_cell_boundaries(self, id=None, coords=None, fmt="dict"):
        """Return b

        Parameters
        ----------
        id : _type_, optional
            _description_, by default None
        coords : _type_, optional
            _description_, by default None
        fmt : str, optional
            _description_, by default "dict"

        Returns
        -------
        _type_
            _description_

        Raises
        ------
        ValueError
            _description_
        ValueError
            _description_
        """
        if id is not None:
            boundaries = self.get_boundaries("id", "set")
            cells_id_b = self.get_cells_id(True, False)
            assert id in cells_id_b, f"cell id is out of range {cells_id_b}."
            cell_neighbors = self.get_cell_neighbors(id=id, boundary=True, fmt=fmt)
        elif coords is not None:
            boundaries = self.get_boundaries("coords", "set")
            cells_coords = self.get_cells_coords(True, False)
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
                set(cell_neighbors["x"]).intersection(boundaries)
            )
            cell_boundaries["y"] = list(
                set(cell_neighbors["y"]).intersection(boundaries)
            )
            cell_boundaries["z"] = list(
                set(cell_neighbors["z"]).intersection(boundaries)
            )
            return cell_boundaries
        elif fmt == "list":
            return list(set(cell_neighbors).intersection(boundaries))
        else:
            raise ValueError("format argument must be either 'dict' or 'list'.")

    # -------------------------------------------------------------------------
    # Properties
    # -------------------------------------------------------------------------

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
        # if not hasattr(self, "is_homogeneous"):
        #     self.get_is_homogeneous()

    # set_properties = set_props

    def set_phi(self, phi, id=None, coords=None):
        """ """
        if id is None and coords is None:
            self.get_ones(True, False)
            self.phi = self.ones * phi
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
            self.get_ones(True, False)
            self.kx = self.ones * kx
            # self.permeability_x = self.kx
            s = "all cells"
        else:
            if id is not None:
                kx_f = self.kx.flatten()
                kx_f[id] = kx
                fshape = self.get_fshape(True, False, False)
                self.kx = kx_f.reshape(fshape)
            if coords is not None:
                icoords = self.get_cell_icoords(coords)
                self.kx[icoords] = kx
            s = "cell " + str(coords)
            self.get_is_homogeneous()
        print(f"- Permeability at x-direction (kx) is set to {kx} for {s}.")

    # set_permeability_x = set_kx

    def set_ky(self, ky, id=None, coords=None):
        if id is None and coords is None:
            self.get_ones(True, False)
            self.ky = self.ones * ky
            # self.permeability_y = self.ky
            s = "all cells"
        else:
            if id is not None:
                ky_f = self.ky.flatten()
                ky_f[id] = ky
                fshape = self.get_fshape(True, False, False)
                self.ky = ky_f.reshape(fshape)
            if coords is not None:
                icoords = self.get_cell_icoords(coords)
                self.ky[icoords] = ky
            s = "cell " + str(coords)
            self.get_is_homogeneous()
        print(f"- Permeability at y-direction (ky) is set to {ky} for {s}.")

    # set_permeability_y = set_ky

    def set_kz(self, kz, id=None, coords=None):
        if id is None and coords is None:
            self.get_ones(True, False)
            self.kz = self.ones * kz
            # self.permeability_z = self.kz
            s = "all cells"
        else:
            if id is not None:
                kz_f = self.kz.flatten()
                kz_f[id] = kz
                fshape = self.get_fshape(True, False, False)
                self.kz = kz_f.reshape(fshape)
            if coords is not None:
                icoords = self.get_cell_icoords(coords)
                self.kz[icoords] = kz
            s = "cell " + str(coords)
            self.get_is_homogeneous()
        print(f"- Permeability at z-direction (kz) is set to {kz} for {s}.")

    # set_permeability_z = set_kz

    @_lru_cache(maxsize=6)
    def get_k(self, dir, boundary=True):
        self.get_fdir()
        assert dir in self.fdir, f"k{dir} is not in flow direction"
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
            self.get_ones(True, False)
            self.z = self.ones * z
            # self.tops = self.z
            s = "all cells"
        else:
            if id is not None:
                z_f = self.z.flatten()
                z_f[id] = z
                fshape = self.get_fshape(True, False, False)
                self.z = z_f.reshape(fshape)
            if coords is not None:
                icoords = self.get_cell_icoords(coords)
                self.z[icoords] = z
            s = "cell " + str(coords)
        print(f"- Tops (z) is set to {z} for {s}.")

    set_tops = set_z

    def get_is_homogeneous(self):
        self.get_fdir()
        lst = []

        phi = self.remove_boundaries(self.phi, False).flatten()
        if np.all(phi == self.phi[0]):
            is_homogeneous_phi = True
        else:
            is_homogeneous_phi = True
        lst.append(is_homogeneous_phi)

        if "x" in self.fdir:
            kx = self.remove_boundaries(self.kx, False).flatten()
            if np.all(kx == self.kx[0]):
                is_homogeneous_x = True
            else:
                is_homogeneous_x = False
            lst.append(is_homogeneous_x)

        if "y" in self.fdir:
            ky = self.remove_boundaries(self.ky, False).flatten()
            if np.all(ky == self.ky[0]):
                is_homogeneous_y = True
            else:
                is_homogeneous_y = False
            lst.append(is_homogeneous_y)

        if "z" in self.fdir:
            kz = self.remove_boundaries(self.kz, False).flatten()
            if np.all(kz == self.kz[0]):
                is_homogeneous_z = True
            else:
                is_homogeneous_z = False
            lst.append(is_homogeneous_z)

        if all(lst):
            self.is_homogeneous = True
        else:
            self.is_homogeneous = False
        print(f"- Flag is_homogeneous is set to {self.is_homogeneous}.")
        return self.is_homogeneous

    # -------------------------------------------------------------------------
    # Properties - G
    # -------------------------------------------------------------------------

    def get_G(self, dir):
        self.get_fdir()
        if dir in self.fdir:
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
        else:
            print(f"- G{dir} is not in fdir of {self.fdir}.")

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
            return dx / (area * k)
        elif self.D == 1:
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

    # -------------------------------------------------------------------------
    # Pyvista
    # -------------------------------------------------------------------------

    @_lru_cache(maxsize=2)
    def get_pyvista_grid(self, boundary=True):
        """
        https://docs.pyvista.org/api/core/_autosummary/pyvista.ExplicitStructuredGrid.html
        """
        if boundary:
            dims = np.array((self.nx_b, self.ny_b, self.nz_b)) + 1
        else:
            dims = np.array((self.nx, self.ny, self.nz)) + 1

        corners = self.get_corners(boundary)
        pyvista_grid = pv.ExplicitStructuredGrid(dims, corners)

        if self.verbose:
            s = "with boundary" if boundary else "without boundary"
            print(f"- Pyvista grid (pyvista_grid) {s} was created.")

        return pyvista_grid

    @_lru_cache(maxsize=2)
    def get_corners(self, boundary=True):
        """

        (Reference: https://docs.pyvista.org/examples/00-load/create-explicit-structured-grid.html)
        """

        if "x" in self.fdir:
            xcorn = np.insert(self.dxx.cumsum(), 0, 0)
        else:
            xcorn = np.arange(0, (self.nx + 1) * self.dxx[0], self.dxx[0])

        if "y" in self.fdir:
            ycorn = np.insert(self.dyy.cumsum(), 0, 0)
        else:
            ycorn = np.arange(0, (self.ny + 1) * self.dyy[0], self.dyy[0])

        if "z" in self.fdir:
            zcorn = np.insert(self.dzz.cumsum(), 0, 0)
        else:
            zcorn = np.arange(0, (self.nz + 1) * self.dzz[0], self.dzz[0])

        # Boundary:
        if boundary:
            ix = 2 if "x" in self.fdir else 0
            iy = 2 if "y" in self.fdir else 0
            iz = 2 if "z" in self.fdir else 0
        else:
            ix = 0
            iy = 0
            iz = 0
            if "x" in self.fdir:
                xcorn = xcorn[1:-1]
            if "y" in self.fdir:
                ycorn = ycorn[1:-1]
            if "z" in self.fdir:
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

        if self.verbose:
            s = "with boundary" if boundary else "without boundary"
            print(f"- Grid corners {s} were calculated.")
            print(
                "    - xcorn shape:",
                xcorn.shape,
                "- ycorn shape:",
                ycorn.shape,
                "- zcorn shape:",
                zcorn.shape,
            )

        # Combine corners:
        corners = np.stack((xcorn, ycorn, zcorn))
        corners = corners.transpose()

        return corners

    # -------------------------------------------------------------------------
    # Dimensions
    # -------------------------------------------------------------------------

    def get_cells_dims(self, dx, dy, dz, verbose=True):
        self.get_n_max(True)
        self.get_fshape(True, False, False)
        m_list = []
        if "x" in self.fdir:
            self.dxx = np.ones(self.nx_b, dtype="int") * dx
            m_list.append(self.dxx)
        else:
            self.dxx = np.ones(self.n_max, dtype="int") * dx
            m_list.append(dx)
        if "y" in self.fdir:
            self.dyy = np.ones(self.ny_b, dtype="int") * dy
            m_list.append(self.dyy)
        else:
            self.dyy = np.ones(self.n_max, dtype="int") * dy
            m_list.append(dy)
        if "z" in self.fdir:
            self.dzz = np.ones(self.nz_b, dtype="int") * dz
            m_list.append(self.dzz)
        else:
            self.dzz = np.ones(self.n_max, dtype="int") * dz
            m_list.append(dz)

        self.dx, self.dy, self.dz = np.meshgrid(*m_list, copy=False)
        fshape = self.get_fshape(True, False, False)
        self.dx = np.transpose(self.dx, axes=(0, 2, 1)).reshape(fshape)
        self.dy = np.transpose(self.dy, axes=(2, 0, 1)).reshape(fshape)
        self.dz = np.transpose(self.dz, axes=(2, 1, 0)).reshape(fshape)

        if verbose:
            print(f"- Cells dims (dx, dy, dz) were computed.")

    @_lru_cache(maxsize=5)
    def get_cells_d(self, dir, boundary=True, fshape=True):
        if dir == "x":
            cells_d = self.dx
        elif dir == "y":
            cells_d = self.dy
        elif dir == "z":
            cells_d = self.dz

        if fshape:
            shape = self.get_fshape(boundary, False, False)
            cells_d = cells_d.reshape(shape)
        if not boundary:
            cells_d = self.remove_boundaries(cells_d)
        return cells_d

    @_lru_cache(maxsize=2)
    def get_cells_dx(self, boundary=True, fshape=True):
        return self.get_cells_d("x", boundary, fshape)

    @_lru_cache(maxsize=2)
    def get_cells_dy(self, boundary=True, fshape=True):
        return self.get_cells_d("y", boundary, fshape)

    @_lru_cache(maxsize=2)
    def get_cells_dz(self, boundary=True, fshape=True):
        return self.get_cells_d("z", boundary, fshape)

    @_lru_cache(maxsize=None)
    def get_cell_d(self, dir, id=None, coords=None):
        cells_d = self.get_cells_d(dir=dir, boundary=True, fshape=True)
        if id is not None:
            return cells_d.flatten()[id]
        elif coords is not None:
            icoords = self.get_cell_icoords(coords)
            return cells_d[icoords]
        else:
            raise ValueError("at least id or coords argument must be defined.")

    @_lru_cache(maxsize=None)
    def get_cell_dx(self, id=None, coords=None):
        return self.get_cell_d("x", id, coords)

    @_lru_cache(maxsize=None)
    def get_cell_dy(self, id=None, coords=None):
        return self.get_cell_d("y", id, coords)

    @_lru_cache(maxsize=None)
    def get_cell_dz(self, id=None, coords=None):
        return self.get_cell_d("z", id, coords)

    # -------------------------------------------------------------------------
    # Helpers
    # -------------------------------------------------------------------------

    def get_verbose_str(self, boundary, fshape):
        s1 = "with fshape" if fshape else "without fshape"
        s2 = "with boundary" if boundary else "without boundary"
        return s1, s2

    def add_boundary(self, d):
        if isinstance(d, [list, tuple, np.array]):
            self.get_n_max(True)
            if len(d) == self.n_max - 2:
                d = np.concatenate((d[0], d, d[-1]))
        return d

    # -------------------------------------------------------------------------
    # Centers
    # -------------------------------------------------------------------------

    @_lru_cache(maxsize=4)
    def get_cells_center(self, boundary=True, fshape=False):
        pyvista_grid = self.get_pyvista_grid(True)
        self.cells_center = pyvista_grid.cell_centers().points

        if fshape:
            shape = self.get_fshape(boundary, True, False)
            self.cells_center = self.cells_center.reshape(shape)

        if not boundary:
            self.cells_center = self.remove_boundaries(self.cells_center, True)

        if self.verbose:
            s1, s2 = self.get_verbose_str(boundary, fshape)
            print(f"- Cells center (cells_center) was computed ({s1} - {s2}).")

        return self.cells_center

    # -------------------------------------------------------------------------
    # Area
    # -------------------------------------------------------------------------

    @_lru_cache(maxsize=4)
    def get_cells_area_x(self, boundary=True, fshape=True):
        self.area_x = self.dy.flatten() * self.dz.flatten()
        if fshape:
            shape = self.get_fshape(boundary, False, False)
            self.area_x = self.area_x.reshape(shape)
        if not boundary:
            self.area_x = self.remove_boundaries(self.area_x)
        return self.area_x

    @_lru_cache(maxsize=4)
    def get_cells_area_y(self, boundary=True, fshape=True):
        self.area_y = self.dx.flatten() * self.dz.flatten()
        if fshape:
            shape = self.get_fshape(boundary, False, False)
            self.area_y = self.area_y.reshape(shape)
        if not boundary:
            self.area_y = self.remove_boundaries(self.area_y)
        return self.area_y

    @_lru_cache(maxsize=4)
    def get_cells_area_z(self, boundary=True, fshape=True):
        self.area_z = self.dx.flatten() * self.dy.flatten()
        if fshape:
            shape = self.get_fshape(boundary, False, False)
            self.area_z = self.area_z.reshape(shape)
        if not boundary:
            self.area_z = self.remove_boundaries(self.area_z)
        return self.area_z

    @_lru_cache(maxsize=8)
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

    @_lru_cache(maxsize=None)
    def get_cell_area(self, dir, id=None, coords=None):
        if id is not None:
            cells_area = self.get_cells_area(dir, True, False)
            return cells_area[id]
        elif coords is not None:
            cells_area = self.get_cells_area(dir, True, True)
            return cells_area[coords[2], coords[1], coords[0]]
        else:
            raise ValueError("At least id or coords argument must be defined.")

    @_lru_cache(maxsize=None)
    def get_cell_area_x(self, id=None, coords=None):
        return self.get_cell_area("x", id, coords)

    @_lru_cache(maxsize=None)
    def get_cell_area_y(self, id=None, coords=None):
        return self.get_cell_area("y", id, coords)

    @_lru_cache(maxsize=None)
    def get_cell_area_z(self, id=None, coords=None):
        return self.get_cell_area("z", id, coords)

    # -------------------------------------------------------------------------
    # Volume
    # -------------------------------------------------------------------------

    @_lru_cache(maxsize=2)
    def get_volume(self, boundary=True):
        pyvista_grid = self.get_pyvista_grid(boundary)
        self.V = pyvista_grid.volume
        return self.V

    @_lru_cache(maxsize=6)
    def get_cells_volume(self, boundary=True, fshape=False, pyvista=False):

        if pyvista:
            pyvista_grid = self.get_pyvista_grid(True)
            self.cells_V = pyvista_grid.compute_cell_sizes()["Volume"]
            self.cells_V = self.cells_V.round(2)
        else:
            self.cells_V = self.dx.flatten() * self.dy.flatten() * self.dz.flatten()

        if fshape:
            shape = self.get_fshape(boundary, False, False)
            self.cells_V = self.cells_V.reshape(shape)
        if not boundary:
            self.cells_V = self.remove_boundaries(self.cells_V)
        print("- Cells volumes (cells_V) was computed.")
        return self.cells_V

    @_lru_cache(maxsize=None)
    def get_cell_volume(self, id=None, coords=None):
        if id is not None:
            cells_volume = self.get_cells_volume(True, False)
            return cells_volume.flatten()[id]
        elif coords is not None:
            cells_volume = self.get_cells_volume(True, True)
            return cells_volume[coords[2], coords[1], coords[0]]
        else:
            raise ValueError("at least id or coords argument must be defined.")

    def show(
        self,
        label=None,  # 'coords' or 'id',
        boundary=False,
        corners=False,
    ):
        """
        - centers_label: str ('coords', 'id')
        """
        self.get_n_max(True)
        pyvista_grid = self.get_pyvista_grid(boundary)

        if self.n_max > 12:
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
                labels = self.get_cells_coords(boundary, False, "tuple")
            elif label == "id":
                labels = self.get_cells_id(boundary, False, "tuple")
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
        title = "{}D model by {} (flow at {}-direction {})".format(
            self.D, label, self.fdir, s
        )

        pl.add_title(title, font="courier", color="white", font_size=8)
        pl.add_camera_orientation_widget()
        pl.enable_fly_to_right_click()
        pl.show_axes()
        if self.D == 1:
            if self.fdir == "y":
                fdir = "yz"
            else:
                fdir = "xz"
        elif self.D == 2:
            fdir = self.fdir
        else:
            fdir = "xz"
        pl.camera_position = fdir
        pl.set_background("black", top="gray")
        pl.show(title="openresim 3D show", full_screen=True)

    # get_flow_shape = get_fshape
    # self.flow_shape = self.fshape
    # get_dimension = get_D
    # self.dimension = self.D


if __name__ == "__main__":

    dx = 11  # [11,21,31,41]
    dy = 12  # [12,22,32,42]
    dz = 13  # [13,23,33,43]
    grid = CartGrid(
        nx=3,
        ny=3,
        nz=1,
        dx=dx,
        dy=dy,
        dz=dz,
        kx=270,
        ky=10,
        kz=20,
        phi=0.27,
    )

    print("Test 1:")
    cells = grid.get_cells_coords(True, False)
    print(cells)
    # print(grid.extract_boundaries(cells, True, "tuple"))
    # print(grid.extract_boundaries(cells, True, "set"))
    # print(grid.extract_boundaries(cells, True, "array"))
    b1 = grid.extract_boundaries(cells, True, "array")
    b2 = np.array([(0, 1, 0), (4, 1, 0)])
    print(utils.intersection(b1, b2, "list"))
    # print(grid.get_cell_boundaries(coords=(3, 3, 0), fmt="list"))

    # grid.show(boundary=True, label="coords")
    # boundaries = grid.get_cell_boundaries(coords=coords, fmt="list")
    # print("Neighbors:", boundaries)

    # grid.set_kx(kx=10, coords=(1,1,1))
    # print(grid.kx)
    # print(grid.get_cells_id())
    # print(grid.kx.shape)
    # print(grid.cells_volume_b.shape)
    # grid.show(boundary=True, label="id", corners=False)
    # grid.show(boundary=False, label="id", corners=False)
    # grid.show(boundary=True, label="dy")
    # grid.show(boundary=True, label="dz")
    # grid.show(boundary=True, label="area_x")
    # grid.show(boundary=True, label="area_y")
    # grid.show(boundary=True, label="area_z")

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
