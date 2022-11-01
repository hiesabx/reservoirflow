"""
Model classes to create reservoir simulation models.

This module contains all model classes used to create a reservoir 
simulation model in combination with a Fluid class and Grid class.
"""
import time
from openresim.base import Base
from openresim import grids, fluids, wells, plots, profme, utils
from openresim.utils import _lru_cache
import numpy as np
import sympy as sym
import scipy.linalg as sl
import scipy.sparse as ss
import scipy.sparse.linalg as ssl
import tensorflow as tf
import matplotlib.pyplot as plt
from tqdm import tqdm
import warnings
import pandas as pd
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
import pyvista as pv
from datetime import date

# from numba import jit


class Model(Base):
    """Model class used to create a reservoir simulation model.

    Model class represents the fluid flow process in a reservoir
    due to pressure change cause by boundary conditions or by
    (production or injection) wells.

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
        start_date: date = None,
        dtype: str = "double",
        unit="field",
        verbose=False,
    ):
        """Reservoir simulation Model class.

        Parameters
        ----------
        grid : grids.Grid
            Grid object .
        fluid : fluids.Fluid
            Fluid object.
        well : wells.Well, optional, by default None
            Well object.
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
        assert self.dtype == fluid.dtype, "fluid dtype is not compatible."

        self.cells_terms = {}
        self.dt = dt
        self.nsteps = 0
        self.tstep = 0
        self.ctime = 0

        self.__initialize__(pi, start_date, well)
        self.__calc_comp()
        # self.__calc_dir_T()
        self.__calc_RHS()
        self.bdict = {}
        # self.set_boundaries({})
        # self.set_boundaries({0: ("rate", 0)})

    # -------------------------------------------------------------------------
    # Basic:
    # -------------------------------------------------------------------------

    def __initialize__(self, pi, start_date, well):
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
            self.pressures[0, self.grid.cells_id] = pi

        if start_date is None:
            self.start_date = date.today()
        else:
            self.start_date = start_date

        self.wells = {}
        self.w_pressures = defaultdict(list)
        if well is not None:
            self.set_well(well)

        self.n = self.grid.get_n(False)
        self.cells_i = self.grid.get_cells_i(False)
        self.cells_id = self.grid.get_cells_id(False, False, "list")
        self.cells_i_dict = dict(zip(self.cells_id, self.cells_i))
        self.boundaries_id = self.grid.get_boundaries("id", "array")

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

    def __calc_dir_T(self):
        """Calculates transmissibility at every flow direction (fdir)."""
        self.T = self.get_cells_T_vect(False, True)

    # def get_cells_T(self, dir=None, boundary=False, fshape=False):
    #     """Returns transmissibility (T) for all cells.

    #     Parameters
    #     ----------
    #     dir : str
    #         direction as string in ['x', 'y', 'z']. If None or 'all',
    #         transmissibility in all grid flow directions (fdir) will be
    #         included in a dictionary.
    #     fshape : bool, optional, by default False
    #         values in flow shape (True) or flatten (False). In this
    #         case, fshape is for cells' boundaries.

    #     Returns
    #     -------
    #     ndarray
    #         array of transmissibility at all cells' boundaries.
    #     """
    #     if dir in ["x", "y", "z"]:
    #         if self.verbose:
    #             s0 = dir + " direction"
    #         T = self.grid.get_cells_G(dir, boundary, fshape) / (
    #             self.fluid.mu * self.fluid.B
    #         )
    #     elif dir in [None, "all", "-"]:
    #         if self.verbose:
    #             s0 = "all directions"
    #         T = {}
    #         for dir in self.grid.get_fdir():
    #             T[dir] = self.get_cells_T(dir, boundary, fshape)
    #     else:
    #         raise ValueError("dir argument is given a wrong value.")

    #     if self.verbose:
    #         s1, s2 = utils.get_verbose_str(boundary, fshape)
    #         print(f"[info] transmissibility (T) in {s0} was computed ({s1} - {s2}).")

    #     return T

    @_lru_cache(maxsize=None)
    def get_cell_T(self, id=None, coords=None, boundary=False):
        """Returns transmissibility (T) at all cell faces.

        Parameters
        ----------
        id : int, iterable of int, by default None
            cell id based on natural order as int.
        coords : iterable of int, by default None
            cell coordinates (i,j,k) as a tuple of int.
        boundary : bool, optional, by default True
            values with boundary (True) or without boundary (False).

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

    def get_cells_T_loop(self, boundary=False, sparse=False):
        """_summary_

        Parameters
        ----------
        boundary : bool, optional
            _description_, by default True
        sparse : bool, optional
            _description_, by default False

        Returns
        -------
        _type_
            _description_
        """
        n = self.grid.get_n(boundary)
        if sparse:
            T_array = ss.lil_matrix((n, n), dtype=self.dtype)
        else:
            T_array = np.zeros((n, n), dtype=self.dtype)

        if boundary:
            for id in self.grid.cells_id:
                T = self.get_cell_T(id, None, False)
                for id_n in T.keys():
                    T_array[id, id_n] = T[id_n]
        else:
            for id in self.grid.cells_id:
                i = self.cells_i_dict[id]
                T = self.get_cell_T(id, None, False)
                cells_i_n = [self.cells_i_dict[x] for x in T.keys()]
                for id_n, i_n in zip(T.keys(), cells_i_n):
                    T_array[i, i_n] = T[id_n]
        return T_array

    def get_cells_T_vect(self, boundary=False, sparse=False):
        return self.grid.get_cells_G(boundary, sparse) / (self.fluid.mu * self.fluid.B)

    def get_cells_T_diag(self, boundary=False, diag_n=1):
        if diag_n == 3:
            diag, _ = self.grid.get_cells_G_diag_3(boundary)
        elif diag_n == 2:
            diag, _ = self.grid.get_cells_G_diag_2(boundary)
            if self.grid.D > 2:
                self.grid.get_cells_G_diag_3(boundary, diag)
        elif diag_n == 1:
            diag = self.grid.get_cells_G_diag_1(boundary)
            if self.grid.D > 1:
                self.grid.get_cells_G_diag_2(boundary, diag)
        return diag / (self.fluid.mu * self.fluid.B)

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
        if fdir == "x":
            k_H = self.grid.k["x"][id]
        elif fdir == "xy":
            k_H = (self.grid.k["x"][id] * self.grid.k["y"][id]) ** 0.5
        elif fdir == "xyz":
            # print(f"[warning] __calc_well_G at {fdir} has to be verified.")
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
        fdir = self.grid.get_fdir()
        if fdir in ["x", "y"]:
            d = self.grid.d["x"][id] ** 2 + self.grid.d["y"][id] ** 2
            return 0.14 * d**0.5
        elif fdir == "xy":
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
        elif fdir == "xyz":
            # print(f"[warning] __calc_well_r_eq at {fdir} has to be verified.")
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
            bottom hole flowing pressure (BHFP). If was not defined,
            None value will be set to zero.
        r : int, float, optional, by default None
            well radius.
        s : int, float, optional, by default None
            well skin factor

        ToDo
        ----
        - Change production to positive and injection to negative.
        """
        assert id in self.grid.cells_id, "a well must be placed within the reservoir"
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
                self.wells[id]["q_sp"] = q
                self.wells[id]["constrain"] = "q"
            if pwf is not None:
                self.wells[id]["pwf"] = pwf
                self.wells[id]["pwf_sp"] = pwf
                if "q" not in self.wells[id].keys():
                    self.wells[id]["constrain"] = "pwf"
                self.w_pressures[id].append(self.pressures[self.tstep, id])
            if "constrain" not in self.wells[id].keys():
                self.wells[id]["constrain"] = None
            if r is not None:
                self.wells[id]["r"] = r
            if s is not None:
                self.wells[id]["s"] = s

        self.wells[id]["r_eq"] = self.__calc_well_r_eq(id)
        self.wells[id]["G"] = self.__calc_well_G(id)
        if "pwf" not in self.wells[id].keys():
            self.wells[id]["pwf"] = 0
            self.wells[id]["pwf_sp"] = 0
            self.w_pressures[id].append(self.pressures[self.tstep, id])

        if self.verbose:
            print(f"[info] a well in cell {id} was set.")

    # -------------------------------------------------------------------------
    # Boundaries:
    # -------------------------------------------------------------------------

    def set_boundary(self, id_b: int, cond: str, v: float):
        """Set a boundary condition in a cell.

        Parameters
        ----------
        id_b : int, optional, by default None
            boundary cell id based on natural order as int.
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
            self.rates[self.tstep, id_b] = v
        elif cond in ["pressure", "press", "p"]:
            self.pressures[self.tstep, id_b] = v
        elif cond in ["gradient", "grad", "g"]:
            ((id, T),) = self.get_cell_T(id_b, None, False).items()
            n_dict = self.grid.get_cell_neighbors(id_b, None, False, "dict")
            dir = [dir for dir in n_dict if id in n_dict[dir]][0]
            self.rates[self.tstep, id_b] = T * self.grid.d[dir][id] * v
        else:
            raise ValueError(f"cond argument {cond} is unknown.")

        if self.verbose:
            print(f"[info] boundary in cell {id_b} was set to constant {cond}.")

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
        for id in bdict:
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
        boundaries = self.grid.get_boundaries("id", "tuple")
        for id in boundaries:
            self.set_boundary(id, cond, v)

    # -------------------------------------------------------------------------
    # Flow Equations:
    # -------------------------------------------------------------------------

    @_lru_cache(maxsize=1)
    def __calc_RHS(self):
        """Calculates flow equation for RHS.

        ToDo
        ----
        - make sure RHS is suitable in case of floats.
        """
        if self.comp_type == "incompressible":
            n = self.grid.get_n(True)
            self.RHS = np.zeros(n, dtype=self.dtype)
        elif self.comp_type == "compressible":
            RHS_n = self.grid.V * self.grid.phi * self.comp
            RHS_d = self.factors["volume conversion"] * self.fluid.B * self.dt
            self.RHS = RHS_n / RHS_d
        else:
            raise ValueError("compressibility type is unknown.")
        return self.RHS

    @_lru_cache(maxsize=None)
    def __calc_n_terms(self, id, id_n, p, T):
        """Calculates cell flow equation with a neighbor cell.

        This function derives flow terms between a specific cell (id)
        and another neighbor cell (n_id).

        Parameters
        ----------
        id : int
            cell id based on natural order as int.
        id_n : int
            neighbor cell id based on natural order as int.
        p : Symbol
            pressure symbol at cell id.
        T : float
            transmissibility with the neighbor cell.

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
        - exec implementation:
            exec(f"p{n_id} = sym.Symbol('p{n_id}')")
            exec(
                f"n_term = n_T * ((p{n_id} - p_id) - "
                + "(self.fluid.g * (self.grid.z[n_id] - self.grid.z[id])))"
            )
            return locals()["n_term"]

        """
        p_n = eval(f"sym.Symbol('p{id_n}')")
        dz = self.grid.z[id_n] - self.grid.z[id]
        return T * ((p_n - p) - (self.fluid.g * dz))

    @_lru_cache(maxsize=None)
    def __calc_b_terms(self, id, id_b, p, T):
        """Calculates cell flow equation with a boundary cell.

        This function derives flow terms between a specific cell (id)
        and another boundary cell (b_id).

        Parameters
        ----------
        id : int
            cell id based on natural order as int.
        id_b : int
            boundary cell id based on natural order as int.
        p : Symbol
            pressure symbol at cell id.
        T : float
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
        p_b = self.pressures[self.tstep, id_b]

        if not np.isnan(p_b):
            dz = self.grid.z[id_b] - self.grid.z[id]
            b_term = T * 2 * ((p_b - p) - (self.fluid.g * dz))
        else:
            b_term = self.rates[self.tstep, id_b]

        return b_term

    def __calc_w_terms(self, id, p):
        """Calculates cell flow equation for the well.

        This function derives flow terms between a specific cell (id)
        and its well (if exists).

        Parameters
        ----------
        id : int
            cell id based on natural order as int.
        p : Symbol or value
            cell pressure symbol or value at cell id.

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
        if "q" in self.wells[id] and self.wells[id]["constrain"] == "q":
            return self.wells[id]["q"]
        else:
            return (
                -self.wells[id]["G"]
                / (self.fluid.B * self.fluid.mu)
                * (p - self.wells[id]["pwf"])
            )

    def __calc_a_term(self, id, p):
        """Calculates cell flow equation for the accumulation term.

        Parameters
        ----------
        id : _type_
            cell id based on natural order as int.

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
                return self.RHS[id] * (p - self.pressures[self.tstep, id])
            except:
                raise Exception("Initial pressure (pi) must be specified")

    def __simplify_eq(self, cell_eq):
        if (
            cell_eq.lhs.as_coefficients_dict()[1] != 0
            or cell_eq.rhs.as_coefficients_dict()[1] != 0
        ):
            cell_eq = cell_eq.simplify()
        return cell_eq

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
        """
        p = eval(f"sym.Symbol('p{id}')")

        if not id in self.cells_terms:
            assert id in self.grid.cells_id, f"id is out of range {self.grid.cells_id}."
            neighbors = self.grid.get_cell_neighbors(id=id, boundary=False, fmt="array")
            boundaries = self.grid.get_cell_boundaries(id=id, fmt="array")
            terms = {"f_terms": [], "a_term": 0}
            T = self.get_cell_T(id, None, True)
            if self.verbose:
                print(f"[info] cell id: {id}")
                print(f"[info]    - Neighbors: {neighbors}")
                print(f"[info]    - Boundaries: {boundaries}")

            for id_n in neighbors:
                n_terms = self.__calc_n_terms(id, id_n, p, T[id_n])
                terms["f_terms"].append(n_terms)
                if self.verbose:
                    print(f"[info] Neighbor terms: {n_terms}")

            for id_b in boundaries:
                b_terms = self.__calc_b_terms(id, id_b, p, T[id_b])
                terms["f_terms"].append(b_terms)
                if self.verbose:
                    print(f"[info] Boundary terms: {b_terms}")

            if id in self.wells.keys():
                w_terms = self.__calc_w_terms(id, p)
                terms["f_terms"].append(w_terms)
                if self.verbose:
                    print(f"[info] Well terms: {w_terms}")

            terms["a_term"] = self.__calc_a_term(id, p)
            if self.verbose:
                print("[info] Accumulation term:", terms["a_term"])

            self.cells_terms[id] = terms
            if self.verbose:
                print("[info] terms:", terms)

        else:
            terms = self.cells_terms[id]
            if id in self.wells.keys() and self.wells[id]["constrain"] == "pwf":
                w_terms = self.__calc_w_terms(id, p)
                if self.verbose:
                    print(f"[info] Well terms (updated): {w_terms}")
                terms["f_terms"][-1] = w_terms
            if self.comp_type == "compressible":
                terms["a_term"] = self.__calc_a_term(id, p)

        cell_eq = sym.Eq(sum(terms["f_terms"]), terms["a_term"])
        cell_eq = self.__simplify_eq(cell_eq)
        cell_eq_lhs = cell_eq.lhs.as_coefficients_dict()

        if self.verbose:
            print(f"[info] Flow equation {id}:", cell_eq)

        return cell_eq_lhs, cell_eq.rhs

    def get_cells_eq(self, threading=False):
        """Return flow equations for all internal cells."""
        cells_eq = {}
        n_threads = self.n // 2
        if threading:
            with ThreadPoolExecutor(n_threads) as executor:
                # with ProcessPoolExecutor(2) as executor:
                equations = executor.map(self.get_cell_eq, self.grid.cells_id)
                for id, eq in zip(self.grid.cells_id, equations):
                    cells_eq[id] = eq
        else:
            for id in self.grid.cells_id:
                cells_eq[id] = self.get_cell_eq(id)
                if self.verbose:
                    print(f"[info] cell id: {id}")
                    print(f"[info]      - lhs: {cells_eq[id][0]}")
                    print(f"[info]      - rhs: {cells_eq[id][1]}")

        return cells_eq

    # -------------------------------------------------------------------------
    # Matrices:
    # -------------------------------------------------------------------------

    def __update_matrices(self, id):
        """Update flow equations' matrices (A, d).

        Parameters
        ----------
        id : int
            cell id based on natural order as int.

        Backup
        ------
        - arrays for lhs and rhs:
            self.d[i] = np.array(cell_rhs).astype(self.dtype)
            self.A[i, ids] = np.array(list(cell_lhs.values())).astype(self.dtype)
        - finding cell i:
            ids = [self.cells_id.index(int(str(s)[1:])) for s in cell_lhs.keys()]
        """
        cell_lhs, cell_rhs = self.cells_eq[id]
        ids = [self.cells_i_dict[int(str(s)[1:])] for s in cell_lhs.keys()]
        self.d[self.cells_i_dict[id]] = cell_rhs
        self.A[self.cells_i_dict[id], ids] = list(cell_lhs.values())
        if self.verbose:
            print(f"[info] cell id: {id}")
            print(f"[info]      - ids: {ids}")
            print(f"[info]      - lhs: {cell_lhs}")
            print(f"[info]      - rhs: {cell_rhs}")

    def init_matrices(self, sparse=False, threading=False):
        """Initialize flow equations' matrices (A, d).

        Parameters
        ----------
        sparse : bool, optional, by default False
            use sparse matrices instead of dense matrices.
        threading : bool, optional, by default False
            use multiple threads for concurrence workers. The maximum
            number of threads are set to the half number of grids.

        Returns
        -------
        _type_
            _description_

        ToDo
        ----
        - Update only required cells.
        """
        self.cells_eq = self.get_cells_eq(threading)

        if self.tstep == 0:
            if sparse:
                self.d = ss.lil_matrix((self.n, 1), dtype=self.dtype)
                self.A = ss.lil_matrix((self.n, self.n), dtype=self.dtype)
            else:
                self.d = np.zeros((self.n, 1), dtype=self.dtype)
                self.A = np.zeros((self.n, self.n), dtype=self.dtype)

        if threading:
            with ThreadPoolExecutor(self.n) as executor:
                # with ProcessPoolExecutor(2) as executor:
                executor.map(self.__update_matrices, self.cells_id)
        else:
            for id in self.cells_id:
                self.__update_matrices(id)

        if self.verbose:
            print("[info] - A:\n", self.A)
            print("[info] - d:\n", self.d)

        return self.A, self.d

    # -------------------------------------------------------------------------
    # Matrices: vectorized
    # -------------------------------------------------------------------------

    def init_A(self, sparse=False):
        """Initialize A matrix.

        A matrix initialization is needed on in timestep zero.

        Parameters
        ----------
        sparse : bool, optional, by default False
            _description_

        Returns
        -------
        ndarray
            matrix A is initialized in place and can be accessed by
            self.A_.
        """
        # T = self.get_cells_T_array(True).toarray()
        # self.A_ = T[:, self.cells_id][self.cells_id]
        # self.A_[self.cells_i, self.cells_i] = (
        #     -self.A_[self.cells_i, :].sum(axis=1) - self.RHS[self.cells_id]
        # )
        # if sparse:
        #     self.A_ = ss.lil_matrix(self.A_, dtype=self.dtype)

        # return self.A_
        # self.A_ = self.get_cells_T_array(False, True).toarray()
        self.A_ = self.get_cells_T_vect(False, sparse)
        v1 = -self.A_[self.cells_i, :].sum(axis=1).flatten()
        v2 = self.RHS[self.cells_id].flatten()
        v3 = v1 - v2
        self.A_[self.cells_i, self.cells_i] = v3
        # if sparse:
        #     self.A_ = ss.lil_matrix(self.A_, dtype=self.dtype)
        return self.A_

    def __init_d(self, sparse):
        """Initialize d vector.

        d vector initialization is needed in every timestep especially
        if the system if compressible (i.g. changes in pressure affect
        vector d values).

        Parameters
        ----------
        sparse : bool, optional, by default False
            _description_

        Returns
        -------
        None
            vector d is initialized in place and can be accessed by
            self.d_.
        Raises
        ------
        Exception
            in case the initial reservoir pressure was not defined.
        """
        if sparse:
            self.d_ = ss.lil_matrix((self.n, 1), dtype=self.dtype)
        else:
            self.d_ = np.zeros((self.n, 1), dtype=self.dtype)

        if self.comp_type == "compressible":
            try:
                self.d_[:] = (
                    -self.RHS[self.grid.cells_id]
                    * self.pressures[self.tstep, self.grid.cells_id]
                ).reshape(-1, 1)
            except:
                raise Exception("Initial pressure (pi) must be specified")

    def __update_z(self):
        """_summary_



        ToDo
        ----
        - T for different geometries is still not ready.
        """
        # all 1D in x direction.
        z = self.grid.z[self.grid.cells_id]
        if not np.all(z == z[0]):
            z_l = np.append(z[1:], np.nan)
            z_u = np.append(np.nan, z[:-1])
            dz_l = self.fluid.g * np.nan_to_num(z_l - z)
            dz_u = self.fluid.g * np.nan_to_num(z_u - z)
            # T = self.T["x"][self.grid.cells_id]
            # T = np.diag(self.get_cells_T(True, False), 1)[self.cells_i]
            T = self.get_cells_T_diag(True, 1)[self.cells_id]
            v = T * dz_l + T * dz_u
            self.d_ += v.reshape(-1, 1)

    def get_matrices_vectorized(self, sparse=False, threading=False):
        """_summary_

        Parameters
        ----------
        sparse : bool, optional
            _description_, by default False
        threading : bool, optional
            _description_, by default False

        Returns
        -------
        _type_
            _description_
        """
        update_z = False
        if self.tstep == 0:
            self.resolve = defaultdict(lambda: False)
            self.init_A(sparse)
            bdict_keys = [
                id_b for id_b in self.bdict.keys() if self.bdict[id_b][0] == "pressure"
            ]
            self.bdict_v = {}
            for id_b in bdict_keys:
                ((id, T),) = self.get_cell_T(id_b, None, False).items()
                p = eval(f"sym.Symbol('p{id}')")
                b_term = self.__calc_b_terms(id, id_b, p, T)
                v0, v1 = b_term.as_coefficients_dict().values()
                self.bdict_v[id_b] = (v0, v1, id)
                self.A_[self.cells_i_dict[id], self.cells_i_dict[id]] += v1
                # self.A_[self.cells_i_dict[id], self.cells_i_dict[id]] -= T * 2

        self.__init_d(sparse)

        for id in self.wells.keys():
            if self.wells[id]["constrain"] == "q":
                w_term = self.__calc_w_terms(id, self.pressures[self.tstep, id])
                self.d_[self.cells_i_dict[id], 0] -= w_term
                update_z = True
            elif self.wells[id]["constrain"] == "pwf":
                p = eval(f"sym.Symbol('p{id}')")
                w_term = self.__calc_w_terms(id, p)
                v = w_term.as_coefficients_dict().values()
                if len(v) == 1:
                    ((v0),) = v
                    v1 = 0
                elif len(v) == 2:
                    v0, v1 = v
                else:
                    raise ValueError("unknown length")
                self.d_[self.cells_i_dict[id], 0] -= v0
                if not self.resolve[id]:
                    self.A_[self.cells_i_dict[id], self.cells_i_dict[id]] += v1
                    self.resolve[id] = True
            else:
                pass  # no constrain

        for id_b in self.bdict.keys():
            id = self.grid.get_cell_neighbors(id_b, None, False, "list")[0]
            if self.bdict[id_b][0] == "pressure":
                self.d_[self.cells_i_dict[id], 0] -= self.bdict_v[id_b][0]
            else:  # elif self.bdict[id_b][0] in ["gradient", "rate"]:
                self.d_[self.cells_i_dict[id], 0] -= self.rates[self.tstep, id_b]

        if update_z:
            self.__update_z()

        return self.A_, self.d_

    def get_d(self, sparse=False):
        if self.comp_type == "incompressible":
            self.d = ss.lil_matrix((self.n, 1), dtype=self.dtype)
        else:
            pressures = self.pressures[self.tstep][self.grid.cells_id]
            RHS = self.RHS[self.grid.cells_id]
            try:
                self.d = ss.lil_matrix((-RHS * pressures).reshape(-1, 1))
            except:
                raise Exception("Initial pressure (pi) must be specified")

        if not sparse:
            self.d = self.d.toarray()

        if self.verbose:
            print("[info] - d:\n", self.d)

        return self.d

    # -------------------------------------------------------------------------
    # Numerical Solution:
    # -------------------------------------------------------------------------

    def __update_wells(self):
        """_summary_


        Backup
        ------
        - well q calc:
            self.__calc_w_terms(
                    id, self.pressures[self.tstep][id]
                )
            or
            self.wells[id]["q"] = (
                -self.wells[id]["G"]
                / (self.fluid.B * self.fluid.mu)
                * (self.pressures[self.tstep][id] - self.wells[id]["pwf"])
            )
        - all calc original:
            if "q" in self.wells[id]:
                self.wells[id]["pwf"] = self.pressures[self.tstep][id] + (
                    self.wells[id]["q"]
                    * self.fluid.B
                    * self.fluid.mu
                    / self.wells[id]["G"]
                )
            self.w_pressures[id].append(self.wells[id]["pwf"])
            if "pwf" in self.wells[id]:
                self.wells[id]["q"] = (
                    -self.wells[id]["G"]
                    / (self.fluid.B * self.fluid.mu)
                    * (self.pressures[self.tstep][id] - self.wells[id]["pwf"])
                )
                self.rates[self.tstep][id] = self.wells[id]["q"]
        """
        resolve = False
        tstep_w_pressures = {}
        for id in self.wells.keys():
            if "q_sp" in self.wells[id]:
                pwf_est = self.pressures[self.tstep, id] + (
                    self.wells[id]["q_sp"]
                    * self.fluid.B
                    * self.fluid.mu
                    / self.wells[id]["G"]
                )
            else:
                pwf_est = self.wells[id]["pwf"]

            if pwf_est > self.wells[id]["pwf_sp"]:
                self.wells[id]["constrain"] = "q"
            else:
                if (
                    pwf_est < self.wells[id]["pwf_sp"]
                    and self.wells[id]["q"] == self.wells[id]["q_sp"]
                ):
                    resolve = True

                self.wells[id]["constrain"] = "pwf"
                pwf_est = self.wells[id]["pwf_sp"]

            self.wells[id]["pwf"] = pwf_est
            q_est = self.__calc_w_terms(id, self.pressures[self.tstep, id])
            self.wells[id]["q"] = self.rates[self.tstep, id] = q_est

            if resolve:
                return True
            else:
                tstep_w_pressures[id] = pwf_est

        for id in self.wells.keys():
            self.w_pressures[id].append(tstep_w_pressures[id])

        return False

    def __update_boundaries(self):
        for id_b in self.bdict.keys():
            ((id_n, T),) = self.get_cell_T(id_b, None, False).items()
            p_n = self.pressures[self.tstep, id_n]
            b_terms = self.__calc_b_terms(id_n, id_b, p_n, T)
            self.rates[self.tstep, id_b] = b_terms

    def __print_arrays(self, sparse):
        if sparse:
            A, d = self.A.toarray(), self.d.toarray()
            A_, d_ = self.A_.toarray(), self.d_.toarray()
        else:
            A, d = self.A, self.d
            A_, d_ = self.A_, self.d_
        print("step:", self.tstep)
        print(np.concatenate([A, A_, A - A_], axis=0))
        print(np.concatenate([d, d_, d - d_], axis=1))
        print()

    def solve(
        self,
        sparse=True,
        threading=False,
        vectorize=True,
        check_MB=True,
        update=True,
        print_arrays=False,
        isolver="cgs",
    ):
        """Solve a single simulation tstep.

        Parameters
        ----------
        sparse : bool, optional, by default True
            _description_
        threading : bool, optional, by default True
            _description_
        vectorize : bool, optional, by default True
            _description_
        check_MB : bool, optional, by default True
            _description_
        update : bool, optional, by default True
            _description_
        print_arrays : bool, optional, by default False
            _description_
        isolver : str, optional, by default "cgs"
            iterative solver for sparse matrices. Available solvers are
            ["bicg", "bicgstab", "cg", "cgs", "gmres", "lgmres",
            "minres", "qmr", "gcrotmk", "tfqmr"].
            If None, direct solver is used. Only relevant when argument
            sparse=True. Option "cgs" is recommended to increase
            performance while option "minres" is not recommended due to
            high MB error. For more information check [1][2].

        Backup
        ------
        - Direct solutions can also be obtained using matrix dot product
        (usually slower) as following:
            pressures = np.dot(np.linalg.inv(A), d).flatten()

        References
        ----------
        [1] https://docs.scipy.org/doc/scipy/reference/sparse.linalg.html#solving-linear-problems
        [2] https://scipy-lectures.org/advanced/scipy_sparse/solvers.html#iterative-solvers
        """
        if print_arrays:
            A, d = self.init_matrices(sparse, threading)  #  has to be first
            self.get_matrices_vectorized(sparse, threading)
            self.__print_arrays(sparse)
        else:
            if vectorize:
                A, d = self.get_matrices_vectorized(sparse, threading)
            else:
                A, d = self.init_matrices(sparse, threading)

        if sparse:
            if isolver:
                if isolver == "bicg":
                    solver = ssl.bicg
                elif isolver == "bicgstab":
                    solver = ssl.bicgstab
                elif isolver == "cg":
                    solver = ssl.cg
                elif isolver == "cgs":
                    solver = ssl.cgs
                elif isolver == "gmres":
                    solver = ssl.gmres
                elif isolver == "lgmres":
                    solver = ssl.lgmres
                elif isolver == "minres":
                    solver = ssl.minres
                    warnings.warn("option isolver='minres' is not recommended.")
                elif isolver == "qmr":
                    solver = ssl.qmr
                elif isolver == "gcrotmk":
                    solver = ssl.gcrotmk
                # elif isolver == "tfqmr":
                # solver = ssl.tfqmr
                else:
                    raise ValueError("isolver is unknown.")
                pressures, exit_code = solver(A.tocsc(), d.todense(), atol=0)
                assert exit_code == 0, "unsuccessful convergence"
            else:
                pressures = ssl.spsolve(A.tocsc(), d.todense(), use_umfpack=True)
        else:
            pressures = sl.solve(A, d).flatten()

        if update:
            self.tstep += 1
            self.pressures = np.vstack([self.pressures, self.pressures[-1]])
            self.pressures[self.tstep, self.grid.cells_id] = pressures
            self.rates = np.vstack([self.rates, self.rates[-1]])
            self.__update_boundaries()
            resolve = self.__update_wells()
            if resolve:
                self.rates = self.rates[: self.tstep]
                self.pressures = self.pressures[: self.tstep]
                self.tstep -= 1
                self.solve(sparse, threading, vectorize, False, True, False)
                if self.verbose:
                    print(f"[info] Time step {self.tstep} was resolved.")

            if check_MB:
                self.check_MB()

        if self.verbose:
            print("[info] Pressures:\n", self.pressures[self.tstep])
            print("[info] rates:\n", self.rates[self.tstep])

    def run(
        self,
        nsteps=10,
        sparse=True,
        threading=True,
        vectorize=True,
        check_MB=True,
        print_arrays=False,
        isolver="cgs",
    ):
        """Perform a simulation run for nsteps.

        Parameters
        ----------
        nsteps : int, optional, by default 10
            _description_
        sparse : bool, optional, by default True
            _description_
        threading : bool, optional, by default True
            _description_
        check_MB : bool, optional, by default True
            _description_
        isolver : str, optional, by default "cgs"
            iterative solver for sparse matrices. Available solvers are
            ["bicg", "bicgstab", "cg", "cgs", "gmres", "lgmres",
            "minres", "qmr", "gcrotmk", "tfqmr"].
            If None, direct solver is used. Only relevant when argument
            sparse=True. Option "cgs" is recommended to increase
            performance while option "minres" is not recommended due to
            high MB error. For more information check [1][2].

        References
        ----------
        [1] https://docs.scipy.org/doc/scipy/reference/sparse.linalg.html#solving-linear-problems
        [2] https://scipy-lectures.org/advanced/scipy_sparse/solvers.html#iterative-solvers
        """
        start_time = time.time()
        self.nsteps += nsteps
        self.run_ctime = 0
        if self.verbose:
            self.verbose = False
            verbose_restore = True
        else:
            verbose_restore = False
        print(f"[info] Simulation run started: {nsteps} timesteps.")
        for _ in tqdm(
            range(1, nsteps + 1),
            unit="steps",
            colour="green",
            position=0,
            leave=True,
        ):
            self.solve(
                sparse,
                threading,
                vectorize,
                check_MB,
                True,
                print_arrays,
                isolver,
            )

        self.run_ctime = round(time.time() - start_time, 2)
        self.ctime += self.run_ctime
        print(
            f"[info] Simulation run of {nsteps} steps",
            f"finished in {self.run_ctime} seconds.",
        )
        if check_MB:
            print(f"\n[info] Material Balance Error: {self.error}.")

        if verbose_restore:
            self.verbose = True

    # -------------------------------------------------------------------------
    # Material Balance:
    # -------------------------------------------------------------------------

    def check_MB(self, verbose=False, error_threshold=0.1):
        """Material Balance Check

        Parameters
        ----------
        verbose : bool, optional, by default False
            _description_
        error_threshold : float, optional, by default 0.1
            _description_
        """
        if verbose:
            print(f"[info] Error in step {self.tstep}")

        if self.comp_type == "incompressible":
            # rates must add up to 0:
            self.error = self.rates[self.tstep].sum()
            if verbose:
                print(f"[info]    - Error: {self.error}")
        elif self.comp_type == "compressible":
            # error over a timestep:
            self.error = (
                self.RHS[self.grid.cells_id]
                * (
                    self.pressures[self.tstep, self.grid.cells_id]
                    - self.pressures[self.tstep - 1, self.grid.cells_id]
                )
            ).sum() / self.rates[self.tstep].sum()
            # error from initial timestep to current timestep: (less accurate)
            self.cumulative_error = (
                self.RHS[self.grid.cells_id]
                * self.dt
                * (
                    self.pressures[self.tstep, self.grid.cells_id]
                    - self.pressures[0, self.grid.cells_id]
                )
            ).sum() / (self.dt * self.tstep * self.rates.sum())
            self.error = abs(self.error - 1)
            if self.verbose:
                print(f"[info]    - Incremental Error: {self.error}")
                print(f"[info]    -  Cumulative Error: {self.cumulative_error}")
                print(
                    f"[info]    -       Total Error: {self.error+self.cumulative_error}"
                )

        if abs(self.error) > error_threshold:
            warnings.warn("High material balance error.")
            print(
                f"[warning] Material balance error ({self.error}) higher",
                f" than the allowed error ({error_threshold}).",
            )

    # -------------------------------------------------------------------------
    # Data:
    # -------------------------------------------------------------------------

    def __concat(self, data, df):
        if df is not None:
            df = pd.concat([df, data], axis=1)
            return df
        return data

    def __add_time(self, units, df=None):
        if units:
            time_str = f" [{self.units['time']}]"
        else:
            time_str = ""
        time = np.arange(0, (self.tstep + 1) * self.dt, self.dt)
        data = pd.Series(time, name="Time" + time_str)
        return self.__concat(data, df)

    def __add_date(self, units, df=None):
        if units:
            date_str = f" [d.m.y]"
        else:
            date_str = ""
        date_series = pd.date_range(
            start=self.start_date,
            periods=self.tstep + 1,
            freq=str(self.dt) + "D",
        ).strftime("%d.%m.%Y")
        data = pd.Series(date_series, name="Date" + date_str)
        return self.__concat(data, df)

    def __add_cells_rate(self, units, boundary, df=None):
        if units:
            rate_str = f" [{self.units['rate']}]"
        else:
            rate_str = ""
        cells_id = self.grid.get_cells_id(boundary, False, "array")
        cells = [id for id in cells_id if id not in self.wells.keys()]
        labels = [f"q{str(id)}" + rate_str for id in cells]
        array = self.rates[:, cells]
        data = pd.DataFrame(array, columns=labels)
        return self.__concat(data, df)

    def __add_cells_pressures(self, units, boundary, df=None):
        if units:
            press_str = f" [{self.units['pressure']}]"
        else:
            press_str = ""
        cells_id = self.grid.get_cells_id(boundary, False, "array")
        labels = [f"P{str(id)}" + press_str for id in cells_id]
        array = self.pressures[:, cells_id]
        data = pd.DataFrame(array, columns=labels)
        return self.__concat(data, df)

    def __add_wells_rate(self, units, boundary, df=None):
        if units:
            rate_str = f" [{self.units['rate']}]"
        else:
            rate_str = ""
        labels = [f"Q{str(id)}" + rate_str for id in self.wells.keys()]
        array = self.rates[:, list(self.wells.keys())]
        data = pd.DataFrame(array, columns=labels)
        return self.__concat(data, df)

    def __add_wells_pressures(self, units, boundary, df=None):
        if units:
            press_str = f" [{self.units['pressure']}]"
        else:
            press_str = ""
        labels = [f"Pwf{str(id)}" + press_str for id in self.w_pressures]
        data = pd.DataFrame(self.w_pressures)
        data.columns = labels
        return self.__concat(data, df)

    def get_dataframe(
        self,
        boundary=True,
        units=True,
        columns=["time", "date", "wells"],
        save=True,
        drop_nan=True,
        drop_zero=True,
    ):
        """Returns simulation data as a dataframe.

        Parameters
        ----------
        boundary : bool, optional, by default True
            values with boundary (True) or without boundary (False).
            It is only relevant when cells columns are selected.
        units : bool, optional, by default True
            column names with units (True) or without units (False).
        columns : list, optional, by default ["time", "date", "wells"]
            selected columns to be added to the dataframe. The following
            options are available:
            "time": for time steps as specified in dt.
            "date": for dates as specified in dt and start_date.
            "q", "rates": for all rates including cells and wells.
            "p", "pressures": for all pressures including cells and wells.
            "cells": for all cells rates and pressures.
            "wells": for all wells rates and pressures.
            "cells_rate": for all cells rates.
            "cells_pressure": for all cells pressures.
            "wells_rate": for all wells rates.
            "wells_pressure": for all wells pressures.
        save : bool, optional, by default True
            save output as a csv file.
        drop_nan : bool, optional, by default True
            drop columns which contain only nan values.
        drop_zero : bool, optional, by default True
            drop columns which contain only zero values.

        Returns
        -------
        DataFrame
            simulation data as a dataframe.
        """

        col_dict = {
            "time": ["t", "time"],
            "date": ["d", "date"],
            "cells_rate": ["q", "rates", "cells", "cells_rate"],
            "cells_pressure": ["p", "pressures", "cells", "cells_pressure"],
            "wells_rate": ["q", "rates", "wells", "wells_rate"],
            "wells_pressure": ["p", "pressures", "wells", "wells_pressure"],
        }

        df = None
        for c in columns:
            if c.lower() in col_dict["time"]:
                df = self.__add_time(units, df)
            if c.lower() in col_dict["date"]:
                df = self.__add_date(units, df)
            if c.lower() in col_dict["cells_rate"]:
                df = self.__add_cells_rate(units, boundary, df)
            if c.lower() in col_dict["cells_pressure"]:
                df = self.__add_cells_pressures(units, boundary, df)
            if c.lower() in col_dict["wells_rate"]:
                df = self.__add_wells_rate(units, boundary, df)
            if c.lower() in col_dict["wells_pressure"]:
                df = self.__add_wells_pressures(units, boundary, df)

        if drop_nan:
            df = df.dropna(axis=1, how="all")

        if drop_zero:
            df = df.loc[:, (df != 0).any(axis=0)]

        df.index.name = "steps"
        if save:
            df.to_csv("model_data.csv")
            print("[info] Model data was successfully saved.")

        return df

    # -------------------------------------------------------------------------
    # Visualization:
    # -------------------------------------------------------------------------

    def plot(self, prop: str = "pressures", id: int = None, tstep: int = None):

        if tstep is None:
            tstep = self.tstep

        if id is not None:
            exec(f"plt.plot(self.{prop}[:, id].flatten())")
            plt.xlabel("Days")
        elif tstep is not None:
            exec(f"plt.plot(self.{prop}[tstep, :].flatten())")
            plt.xlabel("Grid (id)")
            plt.xticks(ticks=range(0, self.grid.nx + 2))
        plt.grid()
        plt.show()

    def plot_grid(self, property: str = "pressures", tstep: int = None):
        if tstep is None:
            tstep = self.tstep
        cells_id = self.grid.get_cells_id(False, False, "list")
        exec(f"plt.imshow(self.{property}[tstep][cells_id][np.newaxis, :])")
        plt.colorbar(label=f"{property.capitalize()} ({self.units[property[:-1]]})")
        plt.title(f"{property.capitalize()} Distribution")
        plt.yticks([])
        plt.xlabel("Grid (i)")
        plt.xticks(ticks=range(0, 4), labels=range(1, 5))
        plt.show()

    # def __show_wells(self, pl):
    #     for w in model.wells:
    #         # x = model.grid.dx[1:w+1].sum() + model.grid.dx[w]//2
    #         # y = model.grid.dy[w]//2
    #         # z = 100
    #         height = model.grid.dz[w] * 10
    #         # well_cell_i = w if boundary else w - 1
    #         well_cell_center = list(
    #             model.grid.get_pyvista_grid(True).extract_cells(w).GetCenter()
    #         )
    #         well_cell_center[2] = height // 2
    #         well = pv.Cylinder(
    #             center=well_cell_center,
    #             height=height,
    #             radius=model.wells[w]["r"],
    #             direction=(0, 0, 1),
    #         )
    #         pl.add_mesh(well)

    #     return pl

    def show(self, property: str, centers=False, boundary=False, bounds=False):
        plots.show(self, property, centers, boundary, bounds)

    # def get_gif(self, prop, boundary=False, wells=True):
    #     if prop in ['p', 'press', 'pressure', 'pressures']:
    #         values = self.pressures
    #     elif prop in ['q', 'Q', 'rete', 'rates']:
    #         values = self.rates

    #     grid = model.grid.get_pyvista_grid(boundary)
    #     grid.cell_data[prop] = values[1]

    #     pl = pv.Plotter(notebook=False, off_screen=True)

    #     if wells:
    #         pl = self.__show_wells(pl)

    #     pl.add_mesh(
    #         grid,
    #         clim=limits,
    #         # style='wireframe',
    #         show_edges=True,
    #         opacity=0.7,
    #         lighting=True,
    #         ambient=0.2,
    #         n_colors=5,
    #         colormap="Blues",
    #         label=property,
    #         categories=True,
    #         # nan_color='gray',
    #         nan_opacity=0.7,
    #         # use_transparency=True,
    #         scalars=values[-1],  # or str 'pressures'
    #         scalar_bar_args=cbar_opt,
    #         show_scalar_bar=True,
    #         # annotations=annotations,
    #     )
    #     pl.open_gif("images/grid.gif")

    #     pts = grid.points.copy()
    #     for step in range(model.nsteps):
    #         pl.update_coordinates(pts, render=False)
    #         pl.update_scalars(values[step], render=False)
    #         pl.render()
    #         pl.write_frame()
    #         # time.sleep(2)

    #     pl.close()

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

    # -------------------------------------------------------------------------
    # End
    # -------------------------------------------------------------------------


if __name__ == "__main__":

    def create_model_example_7_1():
        grid = grids.Cartesian(
            nx=4, ny=1, nz=1, dx=300, dy=350, dz=40, phi=0.27, kx=270, dtype="double"
        )
        fluid = fluids.SinglePhase(mu=0.5, B=1, dtype="double")
        model = Model(grid, fluid, dtype="double", verbose=False)
        model.set_well(id=4, q=-600, s=1.5, r=3.5)
        model.set_boundaries({0: ("pressure", 4000), 5: ("rate", 0)})
        return model

    def create_model():
        grid = grids.Cartesian(
            nx=100,
            ny=100,
            nz=2,
            dx=300,
            dy=350,
            dz=20,
            phi=0.27,
            kx=1,
            ky=1,
            kz=0.1,
            comp=1 * 10**-6,
            dtype="double",
        )
        fluid = fluids.SinglePhase(
            mu=0.5, B=1, rho=50, comp=1 * 10**-5, dtype="double"
        )
        model = Model(
            grid, fluid, pi=6000, dt=5, start_date="10.10.2018", dtype="double"
        )

        cells_id = model.grid.get_cells_id(False, True)[-1].flatten()
        wells = np.random.choice(cells_id, 6, False)
        for id in wells:
            model.set_well(id=id, q=-300, pwf=100, s=1.5, r=3.5)

        wells = np.random.choice(cells_id, 6, False)
        for id in wells:
            model.set_well(id=id, q=100, s=0, r=3.5)

        return model

    model = create_model()
    # model.run(10, isolver="cgs")
    model.run(10, isolver=None)
    model.show("pressures")
