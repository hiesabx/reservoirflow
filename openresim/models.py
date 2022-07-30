"""
Model classes to create reservoir simulation models.

This module contains all model classes used to create a reservoir 
simulation model in combination with with a Fluid class and Grid class.
"""
import time
from openresim.base import Base
from openresim import grids, fluids, wells, plots
from openresim.utils import _lru_cache
import numpy as np
import sympy as sym
import scipy.sparse as ss
import scipy.sparse.linalg as ssl
import matplotlib.pyplot as plt
from tqdm import tqdm


class Model(Base):
    """Model class used to create a reservoir simulation model.

    Model class represents the fluid flow process in a reservoir
    due to pressure change cause by production or injection wells.

    Parameters
    ----------
    Base : class
        Base class with universal settings.

    Returns
    -------
    Model
        Model object.
    """

    name = "Model"

    def __init__(
        self,
        grid: grids.Grid,
        fluid: fluids.Fluid,
        well: wells.Well = None,
        pi: int = None,
        dt: int = 1,
        dtype: str = "double",
        unit="field",
        verbose=True,
    ):
        """Reservoir simulation Model class.

        Parameters
        ----------
        grid : grids.Grid
            Grid module.
        fluid : fluids.SinglePhaseFluid
            Fluid module.
        well : wells.Well, optional, by default None
            Well module.
        pi : int, optional, by default None
            Initial reservoir pressure.
        dt : int, optional, by default 1
            Time duration for each time step.
        dtype : str or `np.dtype`, optional, by default 'double'
            data type used in all arrays. Numpy dtype such as
            `np.single` or `np.double` can be used.
        unit : str ('field', 'metric'), optional, by default 'field'
            units used in input and output. Parameters can be defined as
            `unit='field'` (default) or `unit='metric'`. `units`
            attribute can be accessed from this class using
            (`Model.units`).
        verbose : bool, optional, by default False
            print information for debugging.
        """
        super().__init__(unit, dtype, verbose)
        self.grid = grid
        self.fluid = fluid
        assert self.dtype == grid.dtype, "grid dtype is not compatible."
        assert self.dtype == fluid.dtype, "Fluid dtype is not compatible."

        self.cells_terms = {}
        self.dt = dt
        self.A = None
        self.nsteps = 1
        self.tstep = 0

        self.__initialize__(pi, well)
        self.__calc_comp()
        self.__calc_T()
        self.__calc_RHS()

    # -------------------------------------------------------------------------
    # Basic:
    # -------------------------------------------------------------------------

    def __initialize__(self, pi, well):
        """Initialize reservoir pressures, rates, and wells.

        Parameters
        ----------
        pi : int, float
            initial reservoir pressure.
        well : Well
            well class.
        """
        ones = self.grid.get_ones(True, False, False)[np.newaxis]
        self.pressures = ones * np.nan
        self.rates = self.grid.get_zeros(True, False, False)[np.newaxis]

        self.pi = pi
        if pi is not None:
            cells_id = self.grid.get_cells_id(False, False, "list")
            self.pressures[0][cells_id] = pi

        self.wells = {}
        if well is not None:
            self.set_well(well)

        if self.verbose:
            print("[info] the model was initialized.")

    # -------------------------------------------------------------------------
    # Properties:
    # -------------------------------------------------------------------------

    def __calc_comp(self):
        """Calculates total compressibility."""
        if self.fluid.comp_type == self.grid.comp_type == "incompressible":
            self.set_comp(0)
        else:
            self.set_comp(self.fluid.comp + self.grid.comp)

        if self.verbose:
            print("[info] model compressibility (comp) was calculated.")

    def __calc_T(self):
        """Calculates transmissibility at every flow direction (fdir)."""
        self.T = {}
        for dir in self.grid.get_fdir():
            self.T[dir] = self.get_T(dir, False)

        if self.verbose:
            print("[info] transmissibility in all directions was computed.")

    def get_T(self, dir="x", fshape=False):
        """Returns transmissibility (T) at all cells' boundaries.

        Parameters
        ----------
        dir : str
            direction as string in ['x', 'y', 'z'].
        fshape : bool, optional, by default False
            values in flow shape (True) or flatten (False). In this
            case, fshape is for cells' boundaries.

        Returns
        -------
        ndarray
            array of transmissibility at all cells' boundaries.
        """
        return self.grid.get_cells_G(dir, fshape) / (self.fluid.mu * self.fluid.B)

    @_lru_cache(maxsize=None)
    def get_cell_T(self, id=None, coords=None, boundary=False):
        """Returns transmissibility (T) at all cell faces.

        Parameters
        ----------
        id : int, iterable of int, by default None
            cell id based on natural order as int. For multiple cells,
            list of int [id,id,..] or tuple of int (id,id,...).
            NotFullyImplemented.
        coords : iterable of int, iterable of tuples of int, by default
            None cell coordinates (i,j,k) as a tuple of int. For
            multiple cells, tuple of tuples of int as
            ((i,j,k),(i,j,k),..). NotFullyImplemented.
        boundary : bool, optional, by default True
            values in flow shape (True) or flatten (False).

        Returns
        -------
        ndarray
            array of G based on dir argument.

        ToDo
        ----
        - for now use only with id
        """

        cell_G = self.grid.get_cell_G(id, coords, boundary)
        cell_T = {}
        for id_n in cell_G:
            cell_T[id_n] = cell_G[id_n] / (self.fluid.mu * self.fluid.B)

        return cell_T

    # -------------------------------------------------------------------------
    # Wells:
    # -------------------------------------------------------------------------

    def __calc_well_G(self, id=None):
        """Calculates well Geometry factor (G).

        Parameters
        ----------
        id : int, optional, by default None
            cell id based on natural order as int.

        Returns
        -------
        float
            well geometry factor G.

        ToDo
        ----
        - use k and d based on well direction.
        """
        fdir = self.grid.get_fdir()
        if "x" == fdir:
            k_H = self.grid.k["x"][id]
        elif "xy" == fdir:
            k_H = (self.grid.k["x"][id] * self.grid.k["y"][id]) ** 0.5
        else:
            raise ValueError(f"k for fdir='{fdir}' is not defined.")
        G_n = (
            2
            * np.pi
            * self.factors["transmissibility conversion"]
            * k_H
            * self.grid.d["z"][id]
        )

        G_d = np.log(self.wells[id]["r_eq"] / self.wells[id]["r"] * 12)

        if "s" in self.wells[id].keys():
            G_d += self.wells[id]["s"]

        return G_n / G_d

    def __calc_well_r_eq(self, id):
        """Calculates well equivalent radius (r_eq).

        Parameters
        ----------
        id : int, optional, by default None
            cell id based on natural order as int.

        Returns
        -------
        float
            well equivalent radius (r_eq).

        ToDo
        ----
        - use k and d based on well direction.
        """
        D = self.grid.get_D()
        fdir = self.grid.get_fdir()
        if "x" == fdir:
            d = self.grid.d["x"][id] ** 2 + self.grid.d["x"][id] ** 2
            return 0.14 * d**0.5
        if "y" == fdir:
            d = self.grid.d["y"][id] ** 2 + self.grid.d["y"][id] ** 2
            return 0.14 * d**0.5
        elif "xy" == fdir:
            kx_ky = self.grid.k["x"][id] / self.grid.k["y"][id]
            ky_kx = self.grid.k["y"][id] / self.grid.k["x"][id]
            return (
                0.28
                * (
                    ky_kx**0.5 * self.grid.d["x"][id] ** 2
                    + kx_ky**0.5 * self.grid.d["y"][id] ** 2
                )
                ** 0.5
                / (ky_kx**0.25 + kx_ky**0.25)
            )
        else:
            raise ValueError(f"k for fdir='{fdir}' is not defined.")

    def set_well(self, well=None, id=None, q=None, pwf=None, r=None, s=None):
        """Set a well in a specific cell

        Parameters
        ----------
        well : Well class, optional, by default None
            well information. If this class was used as input, all other
            arguments will be ignored except id will be used instead of
            well.id.
        id : int, optional, by default None
            well location using cell id based on natural order as int.
            This value is given a higher priority over well.id.
        q : int, float, optional, by default None
            well rate as positive for injection or negative for
            production
        pwf : int, float, optional, by default None
            bottom hole flowing pressure (BHFP).
        r : int, float, optional, by default None
            well radius.
        s : int, float, optional, by default None
            well skin factor

        ToDo
        ----
        - Change production to positive and injection to negative.
        """

        if well is not None:
            if id is None:
                id = well.id
            self.wells[id] = vars(well)
        else:
            assert id is not None, "id must be defined"
            if not id in self.wells.keys():
                self.wells[id] = {}
            if q is not None:
                self.wells[id]["q"] = q
            if pwf is not None:
                self.wells[id]["pwf"] = pwf
            if r is not None:
                self.wells[id]["r"] = r
            if s is not None:
                self.wells[id]["s"] = s
        if "q" in self.wells[id]:
            self.rates[0][id] = self.wells[id]["q"]
        self.wells[id]["r_eq"] = self.__calc_well_r_eq(id)
        self.wells[id]["G"] = self.__calc_well_G(id)

        if self.verbose:
            print(f"[info] a well in cell {id} was set.")

    # -------------------------------------------------------------------------
    # Boundaries:
    # -------------------------------------------------------------------------

    def set_boundary(self, id: int, cond: str, v: float):
        """Set a boundary condition in a cell.

        Parameters
        ----------
        id : int, optional, by default None
            well location using cell id based on natural order as int.
        cond : str
            boundary constant condition. Three conditions are possible:
            (2) Constant rate: str in ['rate', 'q'],
            (1) Constant pressure: str in ['pressure', 'press', 'p'],
            (3) Constant pressure gradient: str in ['gradient', 'grad',
            'g'].
        v : int, float
            constant value to specify the condition in cond argument.

        ToDo
        ----
        - d is taken at x direction for gradient.
        """

        if cond in ["rate", "q"]:
            self.rates[0][id] = v
        elif cond in ["pressure", "press", "p"]:
            self.pressures[0][id] = v
        elif cond in ["gradient", "grad", "g"]:
            # n = self.grid.get_cell_neighbors(id, None, False, "tuple")
            # print("neighbors:", n)
            # T = self.get_cell_T(id, None, False)
            # print("trans:", T)
            # print("q:", self.rates)
            # print("T:", self.T)
            # print("d:", self.grid.d)
            # self.rates[0][id] = self.T["x"][id] * self.grid.d["x"][id] * v
            self.rates[0][id] = np.nan
        else:
            raise ValueError(f"cond argument {cond} is unknown.")

        if self.verbose:
            print(f"[info] boundary cond in cell {id} was set to constant {cond}.")

    def set_boundaries(self, bdict: dict):
        """Set boundary conditions using a dictionary.

        Parameters
        ----------
        bdict : dict
            boundary condition dictionary where keys correspond to cells
            id and values are a tuple of cond and value (e.g.
            {0: ("pressure", 4000), 1: ("rate", 0)})
        """
        self.bdict = bdict
        boundaries = self.grid.get_boundaries("id", "set")
        for id in bdict.keys():
            assert id in boundaries, f"cell {id} is not a boundary cell."
            cond, v = bdict[id]
            self.set_boundary(id, cond, v)

    def set_all_boundaries(self, cond, v):
        """Set the same boundary condition in all boundaries.

        Parameters
        ----------
        cond : str
            boundary constant condition. Three conditions are possible:
            (2) Constant rate: str in ['rate', 'q'],
            (1) Constant pressure: str in ['pressure', 'press', 'p'],
            (3) Constant pressure gradient: str in ['gradient', 'grad',
            'g'].
        v : int, float
            constant value to specify the condition in cond argument.
        """
        boundaries = self.grid.get_boundaries("id", "set")
        for id in boundaries:
            self.set_boundary(id, cond, v)

    # -------------------------------------------------------------------------
    # Flow Equations:
    # -------------------------------------------------------------------------

    def __calc_RHS(self):
        """Calculates flow equation for RHS.

        ToDo
        ----
        - make sure RHS is suitable in case of floats.

        Backup
        ------
        - in case of incompressible:
            n = self.grid.get_n_cells(True)
            self.RHS =  np.zeros(n, dtype=self.dtype)
            or
            self.RHS = 0
        """
        if self.comp_type == "incompressible":
            n = self.grid.get_n_cells(True)
            self.RHS = np.zeros(n, dtype=self.dtype)
            # or self.RHS = 0
        elif self.comp_type == "compressible":
            RHS_n = self.grid.V * self.grid.phi * self.comp
            RHS_d = self.factors["volume conversion"] * self.fluid.B * self.dt
            self.RHS = RHS_n / RHS_d
        else:
            raise ValueError("compressibility type is unknown.")

    def __calc_n_terms(self, id, n_id, p_id, n_T):
        """Calculates cell flow equation with a neighbor cell.

        This function derives flow terms between a specific cell (id)
        and another neighbor cell (n_id).

        Parameters
        ----------
        id : int
            cell id based on natural order as int.
        n_id : int
            neighbor cell id based on natural order as int.
        p_id : Symbol
            pressure symbol at cell id.
        n_T : float
            transmissibility with neighbor cell.

        Returns
        -------
        _type_
            _description_

        Backup
        ------
        - calc without exec:
            exec(f"dp_ = p{n_id} - p{id}")
            dp = locals()["dp_"]
            trans = self.T[dir][min(n_id, id)]
            h = self.grid.z[n_id] - self.grid.z[id]
            acc = self.fluid.g * h
            n_term = trans * (dp - acc)
        - T with matrix:
            self.T[dir][min(n_id,id)]
        """
        exec(f"p{n_id} = sym.Symbol('p{n_id}')")
        exec(
            f"n_term = n_T * ((p{n_id} - p_id) - "
            + "(self.fluid.g * (self.grid.z[n_id] - self.grid.z[id])))"
        )

        return locals()["n_term"]

    def __calc_b_terms(self, id, b_id, p_id, b_T):
        """Calculates cell flow equation with a boundary cell.

        This function derives flow terms between a specific cell (id)
        and another boundary cell (b_id).

        Parameters
        ----------
        id : int
            cell id based on natural order as int.
        b_id : int
            boundary cell id based on natural order as int.
        p_id : Symbol
            pressure symbol at cell id.
        b_T : float
            transmissibility with boundary cell.

        Returns
        -------
        _type_
            _description_

        Backup
        ------
        - exec with boundary pressure:
        exec(f"p{b_id}=sym.Symbol('p{b_id}')")
        if not np.isnan(self.pressures[self.tstep][b_id]):
            exec(
                f"b_term = b_T * 2 * ((p{b_id} - p_id)"
                + " - (self.fluid.g * (self.grid.z[b_id]-self.grid.z[id])))"
            )
            exec(f"b_term = b_term.subs(p{b_id}, self.pressures[self.tstep][b_id])")
        else:
            exec(f"b_term = self.rates[self.tstep][b_id]")

        return locals()["b_term"]

        - T with matrix:
            self.T[dir][min(b_id,id)]
        """
        p_id_b = self.pressures[self.tstep][b_id]

        if not np.isnan(p_id_b):
            b_term = (
                b_T
                * 2
                * (
                    (p_id_b - p_id)
                    - (self.fluid.g * (self.grid.z[b_id] - self.grid.z[id]))
                )
            )
        else:
            if b_id in self.bdict and self.bdict[b_id][0] == "gradient":
                v = self.bdict[b_id][1]
                # ToDo: 0 or self.tstep?
                self.rates[self.tstep][b_id] = b_T * self.grid.d["x"][id] * v
            b_term = self.rates[self.tstep][b_id]

        return b_term

    def __calc_w_terms(self, id, p_id):
        """Calculates cell flow equation for the well.

        This function derives flow terms between a specific cell (id)
        and its well (if exists).

        Parameters
        ----------
        id : int
            cell id based on natural order as int.
        p_id : Symbol
            cell pressure symbol at cell id.

        Returns
        -------
        _type_
            _description_

        Backup
        ------
        - exec:
            exec(
                f"w_term = - self.wells[id]['G']"
                + f"/ (self.fluid.B*self.fluid.mu)"
                + f"* (p_id - self.wells[id]['pwf'])"
            )
            return locals()["w_term"]

        """
        if "q" in self.wells[id]:
            return self.wells[id]["q"]
        else:
            return (
                -self.wells[id]["G"]
                / (self.fluid.B * self.fluid.mu)
                * (p_id - self.wells[id]["pwf"])
            )

    def __calc_a_term(self, id, p_id):
        """Calculates cell flow equation for the accumulation term.

        Parameters
        ----------
        id : _type_
            _description_

        Returns
        -------
        _type_
            _description_

        Raises
        ------
        Exception
            _description_

        ToDo
        ----
        - consider unifying RHS or if cond.
        """
        if self.comp_type == "incompressible":
            return 0
        else:
            try:
                return self.RHS[id] * (p_id - self.pressures[self.tstep][id])
            except:
                raise Exception("Initial pressure (pi) must be specified")

    # @_lru_cache(maxsize=None)
    def get_cell_eq(self, id):
        """Return cell equation.

        Parameters
        ----------
        id : int, optional
            cell id based on natural order as int.

        Returns
        -------
        tuple
            cell equation as a tuple of (lhs, rhs).

        ToDo
        ----
        - n_term naming.

        Backup
        ------
        - constant pressure:
            # exec(f"p{i}=sym.Symbol('p{i}')")
            # ToDo: keep pressure constant at specific cell (requires A adjust)
            # if not np.isnan(self.pressures[self.tstep][i]):
            #     exec(f"p{i} = {self.pressures[self.tstep][i]}")
        - n_term to use pressure values:
            # To Do: keep pressure constant at specific cell (requires A adjust)
            # if not np.isnan(self.pressures[self.tstep][neighbor]):
            #     exec(f"p{neighbor} = {self.pressures[self.tstep][neighbor]}")
        - n_term in one calc.
        exec(
            f"n_term = self.T[dir][min(neighbor,id)] * ((p{n_id} - p{id})
            - (self.fluid.g * (self.grid.z[neighbor] - self.grid.z[id])))"
        )
        - dict sort:
            cell_eq_lhs = dict(
                sorted(
                    cell_eq.lhs.as_coefficients_dict().items(),
                    key=lambda x: int(str(x[0])[1:]),
                )
            )
        """

        cells_id = self.grid.get_cells_id(False, False, "set")
        assert id in cells_id, f"id is out of range {cells_id}."
        exec(f"p{id}=sym.Symbol('p{id}')")
        p_id = locals()[f"p{id}"]

        if id not in self.cells_terms:
            neighbors = self.grid.get_cell_neighbors(id=id, boundary=False, fmt="tuple")
            boundaries = self.grid.get_cell_boundaries(id=id, fmt="tuple")
            terms = []
            T = self.get_cell_T(id, None, True)
            if self.verbose:
                print(f"[info] cell id: {id}")
                print(f"[info]    - Neighbors: {neighbors}")
                print(f"[info]    - Boundaries: {boundaries}")
            for n_id in neighbors:
                n_terms = self.__calc_n_terms(id, n_id, p_id, T[n_id])
                terms.append(n_terms)
                if self.verbose:
                    print(f"[info] Neighbor terms: {n_terms}")
            for b_id in boundaries:
                b_terms = self.__calc_b_terms(id, b_id, p_id, T[b_id])
                terms.append(b_terms)
                if self.verbose:
                    print(f"[info] Boundary terms: {b_terms}")
            if id in self.wells:
                w_terms = self.__calc_w_terms(id, p_id)
                terms.append(w_terms)
                if self.verbose:
                    print(f"[info] Well terms: {w_terms}")
            if self.verbose:
                print("[info] terms:", terms)

            self.cells_terms[id] = terms
        else:
            terms = self.cells_terms[id]

        a_term = self.__calc_a_term(id, p_id)
        if self.verbose:
            print("[info] accumulation terms:", a_term)

        cell_eq = sym.Eq(sum(terms), a_term)

        if (
            cell_eq.lhs.as_coefficients_dict()[1] != 0
            or cell_eq.rhs.as_coefficients_dict()[1] != 0
        ):
            cell_eq = cell_eq.simplify()

        if self.verbose:
            print("[info] flow equation:", cell_eq)

        cell_eq_lhs = dict(cell_eq.lhs.as_coefficients_dict())

        return cell_eq_lhs, cell_eq.rhs

    # @lru_cache(maxsize=None)
    def get_cells_eq(self):
        """Return flow equations for all internal cells.

        Backup
        ------
        - internal cells id:
            for i in self.grid.order[self.grid.i_blocks.astype("bool")]:
        """
        cells_eq = {}
        for id in self.grid.get_cells_id(False, False, "tuple"):
            cells_eq[id] = self.get_cell_eq(id)
            if self.verbose:
                print(f"[info] cell id: {id}")
                print(f"[info]      - lhs: {cells_eq[id][0]}")
                print(f"[info]      - rhs: {cells_eq[id][1]}")

        return cells_eq

    # -------------------------------------------------------------------------
    # Coefficient Matrix:
    # -------------------------------------------------------------------------

    # @lru_cache(maxsize=None)
    def __update_matrix(self, id, sparse=True):
        """
        Update coefficient matrix (A) and result vector (d).
        Note: arrays are passed by reference.
        """
        i_lhs, i_rhs = self.get_cell_eq(id)
        i_lhs = list(i_lhs.values())
        self.d[id - 1] = np.array(i_rhs).astype(self.dtype)

        if self.tstep == 0:
            if id == 1:
                self.A[id - 1, id - 1 : id + 1] = i_lhs
            elif id == self.grid.nx:
                self.A[id - 1, id - 2 : id] = i_lhs
            else:
                self.A[id - 1, id - 2 : id + 1] = i_lhs

    # @lru_cache(maxsize=None)
    def get_matrix(self, sparse=True, verbose=False):
        """Create coefficient matrix (A) and result vector (d)."""
        if self.grid.D == 1:

            # Construct d vector:
            if all(self.RHS == 0):
                self.d = ss.lil_matrix((self.grid.nx, 1), dtype=self.dtype)
            else:
                try:
                    self.d = ss.lil_matrix(
                        (-self.RHS[1:-1] * self.pressures[self.tstep][1:-1]).reshape(
                            -1, 1
                        )
                    )
                except:
                    raise Exception("Initial pressure (pi) must be specified")
            if not sparse:
                self.d = self.d.toarray()  # ss.lil_matrix(self.d, dtype=self.dtype)

            # Construct A matrix:
            if self.tstep == 0:
                self.A = ss.diags(
                    [
                        -self.T[1:] - self.T[:-1] - self.RHS[1:-1],
                        self.T[1:-1],  # East trans for interior blocks
                        self.T[1:-1],
                    ],  # West trans for interior blocks
                    [0, 1, -1],
                    shape=(self.grid.nx, self.grid.nx),
                    format="lil",  # “dia”, “csr”, “csc”, “lil”
                    dtype=self.dtype,
                )
                self.A[0, 0] = -self.T[1] - self.RHS[1]
                self.A[-1, -1] = -self.T[-2] - self.RHS[-1]
                if not sparse:
                    self.A = self.A.toarray()

            # Update matrix if there is pressure or flow in 'west' boundary:
            if (
                not np.isnan(self.pressures[self.tstep][0])
                or self.rates[self.tstep][0] != 0
            ):
                self.__update_matrix(1, sparse, verbose)

            # Update matrix if there is pressure or flow in 'east' boundary:
            if (
                not np.isnan(self.pressures[self.tstep][-1])
                or self.rates[self.tstep][-1] != 0
            ):
                # at last grid: self.grid.nx or -2
                self.__update_matrix(self.grid.nx, sparse, verbose)

            # Update matrix in wells i_blocks:
            for i in self.wells.keys():
                self.__update_matrix(i, sparse, verbose)

            if verbose:
                if sparse:
                    print("- A:\n", self.A.toarray())
                    print("- d:\n", self.d.toarray())
                else:
                    print("- A:\n", self.A)
                    print("- d:\n", self.d)

            return self.A, self.d

    def get_d(self, sparse=False):
        if self.comp_type == "incompressible":
            n_cells = self.grid.get_n_cells(False)
            self.d = ss.lil_matrix((n_cells, 1), dtype=self.dtype)
        else:
            cells_id = self.grid.get_cells_id(False, False, "array")
            pressures = self.pressures[self.tstep][cells_id]
            RHS = self.RHS[cells_id]
            try:
                self.d = ss.lil_matrix((-RHS * pressures).reshape(-1, 1))
            except:
                raise Exception("Initial pressure (pi) must be specified")

        if not sparse:
            self.d = self.d.toarray()

        if self.verbose:
            print("[info] - d:\n", self.d)

        return self.d

    def init_matrices(self, sparse=False):

        if self.tstep == 0:
            n = self.grid.get_n_cells(False)
            if sparse:
                self.d = ss.lil_matrix((n, 1), dtype=self.dtype)
                self.A = ss.lil_matrix((n, n), dtype=self.dtype)
            else:
                self.d = np.zeros((n, 1), dtype=self.dtype)
                self.A = np.zeros((n, n), dtype=self.dtype)

        cells_eq = self.get_cells_eq()
        cells_id = self.grid.get_cells_id(False, False, "tuple")
        for i, id in enumerate(cells_id):
            cell_lhs, cell_rhs = cells_eq[id]  # or self.get_cell_eq(id)
            ids = [cells_id.index(int(str(s)[1:])) for s in cell_lhs.keys()]
            self.d[i] = np.array(cell_rhs).astype(self.dtype)
            self.A[i, ids] = np.array(list(cell_lhs.values())).astype(self.dtype)

            if self.verbose:
                print(f"[info] cell id: {id}")
                print(f"[info]      - ids: {ids}")
                print(f"[info]      - lhs: {cell_lhs}")
                print(f"[info]      - rhs: {cell_rhs}")

        if self.verbose:
            print("[info] - A:\n", self.A)
            print("[info] - d:\n", self.d)

        return self.A, self.d

    # -------------------------------------------------------------------------
    # Numerical Solution:
    # -------------------------------------------------------------------------

    # @lru_cache(maxsize=None)
    def solve(self, sparse=True, check_MB=True, update=True, verbose=False):

        self.init_matrices(sparse)

        if sparse:
            pressures = ssl.spsolve(self.A.tocsc(), self.d)
        else:
            pressures = np.linalg.solve(self.A, self.d).flatten()
            # same as: np.dot(np.linalg.inv(self.A), self.d)

        if self.grid.D == 1:
            # Update pressures:
            if update:
                self.tstep += 1
                self.pressures = np.vstack([self.pressures, self.pressures[-1]])
                # self.pressures[self.tstep] = self.pressures[self.tstep-1].copy()
                self.pressures[self.tstep][1:-1] = pressures
                self.rates = np.vstack([self.rates, self.rates[-1]])
                # self.rates[self.tstep] = self.rates[self.tstep-1].copy()

                # Update wells:
                for i in self.wells.keys():
                    if "q" in self.wells[i]:
                        self.wells[i]["pwf"] = self.pressures[self.tstep][i] + (
                            self.wells[i]["q"]
                            * self.fluid.B
                            * self.fluid.mu
                            / self.wells[i]["G"]
                        )
                    if "pwf" in self.wells[i]:
                        self.wells[i]["q"] = (
                            -self.wells[i]["G"]
                            / (self.fluid.B * self.fluid.mu)
                            * (self.pressures[self.tstep][i] - self.wells[i]["pwf"])
                        )
                        self.rates[self.tstep][i] = self.wells[i]["q"]

                # Update boundaries:
                for boundary in self.grid.boundaries:
                    i = 1 if boundary == 0 else boundary - 1
                    if not np.isnan(self.pressures[self.tstep][boundary]):
                        self.rates[self.tstep][boundary] = (
                            self.T[min(i, boundary)]
                            * 2
                            * (
                                (
                                    self.pressures[self.tstep][boundary]
                                    - self.pressures[self.tstep][i]
                                )
                                - (
                                    self.fluid.g
                                    * (self.grid.z[boundary] - self.grid.z[i])
                                )
                            )
                        )

                if check_MB:
                    self.check_MB(verbose)

        if verbose:
            print("- Pressures:\n", self.pressures[self.tstep])
            print("- rates:\n", self.rates[self.tstep])

        return pressures

    def run(self, nsteps=10, sparse=True, check_MB=True, verbose=False):
        self.nsteps += nsteps
        start_time = time.time()
        for _ in tqdm(
            range(1, nsteps + 1), unit="steps", colour="green", position=0, leave=True
        ):
            self.solve(sparse, check_MB, update=True, verbose=verbose)
        duration = round(time.time() - start_time, 2)
        print(f"Simulation run of {nsteps} steps is finished in {duration} seconds.")

    # -------------------------------------------------------------------------
    # Material Balance:
    # -------------------------------------------------------------------------

    def check_MB(self, verbose=False, error_threshold=0.1):
        """Material Balance Check"""
        if verbose:
            print(f"Error in step {self.tstep}")
        if self.comp_type == "incompressible":
            self.error = self.rates[self.tstep].sum()  # must add up to 0
            if verbose:
                print(f"    - Error: {self.error}")
        elif self.comp_type == "compressible":
            # Check MB error over a time step:
            self.incremental_error = (
                self.RHS[1:-1]
                * (
                    self.pressures[self.tstep][1:-1]
                    - self.pressures[self.tstep - 1][1:-1]
                )
            ).sum() / self.rates[self.tstep].sum()
            # Check MB error from initial state to current time step: (less accurate)
            self.cumulative_error = (
                self.RHS[1:-1]
                * self.dt
                * (self.pressures[self.tstep][1:-1] - self.pressures[0][1:-1])
            ).sum() / (self.dt * self.tstep * self.rates.sum())
            self.error = abs(self.incremental_error - 1)
            if verbose:
                print(f"    - Incremental Error: {self.incremental_error}")
                print(f"    -  Cumulative Error: {self.cumulative_error}")
                print(
                    f"    -       Total Error: {self.incremental_error+self.cumulative_error}"
                )

        assert (
            abs(self.error) < error_threshold
        ), f"""
        Material balance error ({self.error}) higher than the allowed error ({error_threshold})."""

    # -------------------------------------------------------------------------
    # Visualization:
    # -------------------------------------------------------------------------

    def plot(self, property: str = "pressures", i: int = None, tstep: int = None):

        if tstep is None:
            tstep = self.tstep

        if i is not None:
            exec(f"plt.plot(self.{property}[:, i].flatten())")
            plt.xlabel("Days")
        elif tstep is not None:
            exec(f"plt.plot(self.{property}[tstep, :].flatten())")
            plt.xlabel("Grid (i)")
            plt.xticks(ticks=range(0, self.grid.nx + 2))
        plt.grid()
        plt.show()

    def plot_grid(self, property: str = "pressures", tstep: int = None):
        if tstep is None:
            tstep = self.tstep
        exec(f"plt.imshow(self.{property}[tstep][1:-1][np.newaxis, :])")
        plt.colorbar(label=f"{property.capitalize()} ({self.units[property[:-1]]})")
        plt.title(f"{property.capitalize()} Distribution")
        plt.yticks([])
        plt.xlabel("Grid (i)")
        plt.xticks(ticks=range(0, 4), labels=range(1, 5))
        plt.show()

    def show_grid(
        self, property: str, show_centers=True, show_boundary=False, show_bounds=False
    ):
        plots.show_grid(self, property, show_centers, show_boundary, show_bounds)

    def copy(self):
        """Copy model (under development)

        Returns:
            _type_: _description_
        """
        # https://stackoverflow.com/questions/48338847/how-to-copy-a-python-class-instance-if-deepcopy-does-not-work
        copy_model = Model(
            grid=self.grid,
            fluid=self.fluid,
            pi=self.pi,
            dt=self.dt,
            dtype=self.dtype,
            unit=self.unit,
        )
        # for w in self.wells:
        #     well = wells.Well(self.wells[w])
        #     copy_model.set_well(well)
        # copy_model.set_boundaries(self.b_dict)
        return copy_model

    # -------------------------------------------------------------------------
    # Synonyms:
    # -------------------------------------------------------------------------

    def allow_synonyms(self):
        self.set_transmissibility = self.set_trans
        self.transmissibility = self.T
        self.set_properties = self.__calc_comp

    # -------------------------------------------------------------------------
    # End
    # -------------------------------------------------------------------------


if __name__ == "__main__":
    grid = grids.Cartesian(
        nx=4,
        ny=1,
        nz=1,
        dx=300,
        dy=350,
        dz=40,
        phi=0.27,
        kx=270,
        ky=200,
        kz=200,
        comp=1 * 10**-6,
        dtype="double",
    )
    fluid = fluids.SinglePhase(mu=0.5, B=1, rho=50, comp=1 * 10**-5, dtype="double")

    model = Model(grid, fluid, pi=4000, dtype="double")
    model.set_well(id=5, q=-600, s=1.5, r=3.5)
    model.set_well(id=6, pwf=1000, s=1.5, r=3.5)
    model.set_boundaries({0: ("pressure", 4000), 1: ("rate", 0)})
    # print(model.RHS)
    # print(model.wells)
    # print(grid.get_cell_G(id=5, boundary=True))
    # print(model.get_cell_T(5, None, True))
    # print(model.get_cell_eq(5))
    # print(model.get_cell_eq(6))
    # model.get_cells_eq()
    print(grid.phi)
    print(model.get_d())
    # model.run(nsteps=6, sparse=False, check_MB=True, verbose=False)
    # print(grid.get_Gx(fshape=True))
