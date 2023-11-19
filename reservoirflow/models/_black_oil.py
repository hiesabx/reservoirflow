import time
import warnings
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from datetime import date

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.linalg as sl
import scipy.sparse as ss
import scipy.sparse.linalg as ssl
import sympy as sym
from tqdm import tqdm

import reservoirflow as rf
from reservoirflow import fluids, grids, scalers, utils, wells
from reservoirflow.models._model import _Model
from reservoirflow.utils.helpers import _lru_cache

# if __name__ == "__main__":
#     from _model import _Model
# else:
#     from ._model import _Model


class BlackOil(_Model):
    """Black oil model class.

    Returns
    -------
    Model
        BlackOil model object.
    """

    name = "BlackOil Model"

    def __init__(
        self,
        grid: grids.RegularCartesian,
        fluid: fluids.SinglePhase,
        well=None,  # wells.Well = None,
        pi: int = None,
        dt: int = 1,
        start_date: date = None,
        dtype: str = "double",
        unit="field",
        verbose=False,
    ):
        """Create a BlackOil reservoir simulation model.

        Parameters
        ----------
        grid : rf.grids.Grid
            Grid object.
        fluid : rf.fluids.Fluid
            Fluid object.
        well : rf.wells.Well, optional
            Well object.
        pi : int, optional
            Initial reservoir pressure.
        dt : int, optional
            Time step duration.
        start_date : date, optional
            Start date of the simulation run. If None, today's date is
            used.
        dtype : str or `np.dtype`, optional
            data type used in all arrays. Numpy dtype such as
            `np.single` or `np.double` can be used.
        unit : str ('field', 'metric', 'lab'), optional
            units used in input and output. Parameters can be defined as
            `unit='field'` (default), `unit='metric'`, or `unit='lab'`.
            `units`attribute can be accessed from this class using
            (`Model.units`).
        verbose : bool, optional
            print information for debugging.

        Notes
        -----
        .. note::
            Units are defined based on `unit` argument, for more
            details, check
            `Units & Factors </user_guide/units_factors/units_factors.html>`_.
            For definitions, check
            `Glossary </user_guide/glossary/glossary.html>`_.
        """
        super().__init__(unit, dtype, verbose)
        self.grid = grid
        self.fluid = fluid
        assert self.dtype == grid.dtype, "grid dtype is not compatible."
        assert self.dtype == fluid.dtype, "fluid dtype is not compatible."

        self.cells_terms = {}
        self.dt = dt
        self.nsteps = 1
        self.tstep = 0
        self.ctime = 0

        self.__initialize__(pi, start_date, well)
        self.__calc_comp()
        self.__calc_RHS()
        self.bdict = {}
        self.bdict_update = []

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

        Notes
        -----
        Initialization with initial pressure (Pi) :
            - setting boundaries with Pi have wrong effect on
            __calc_b_terms() method where pressure is taken instead of
            taking rate specified at the boundary (implementation 1).
            >>> self.pressures[0, :] = pi
            - while implementation 2 solves this issue, specifying
            initial pressure at boundaries will be carried (copied)
            in the following steps as if this is a constant boundary
            cond which is misleading.
        """
        ones = self.grid.get_ones(True, False, False)[np.newaxis]
        self.pressures = ones * np.nan
        self.rates = self.grid.get_zeros(True, False, False)[np.newaxis]
        self.ds = self.grid.get_zeros(False, False, False)[np.newaxis]
        n = self.grid.get_n(False)
        self.As = np.zeros((1, n * n), dtype=self.dtype)

        self.pi = pi
        if pi is not None:
            self.pressures[0, self.grid.cells_id] = pi
        else:
            warnings.warn("Initial reservoir pressure is not defined.")
            print(f"[warning] Pi is by default set to {self.pi}.")

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
        self.cells_id = self.grid.get_cells_id(False, False, "array")
        self.cells_i_dict = dict(zip(self.cells_id, self.cells_i))
        self.boundaries_id = self.grid.get_boundaries("id", "array")

        self.scalers_dict = {
            "time": ["MinMaxScaler", (0, 1)],
            "space": ["MinMaxScaler", (-1, 1)],
            "pressure": ["MinMaxScaler", (-1, 1)],
            "rate": [None, None],
        }
        self.set_scalers(self.scalers_dict)

        if self.verbose:
            print("[info] the model was initialized.")

    # -------------------------------------------------------------------------
    # Properties:
    # -------------------------------------------------------------------------

    def get_shape(self, boundary: bool = True) -> tuple:
        """Solution shape.

        Parameters
        ----------
        boundary : bool, optional
            with grid boundary (True) or without grid boundary (False).

        Returns
        -------
        tuple
            tuple as (number of time steps, number of girds)
        """
        return (self.nsteps, self.grid.get_n(boundary))

    def __calc_comp(self):
        """Calculates total compressibility."""
        if self.fluid.comp_type == self.grid.comp_type == "incompressible":
            self.set_comp(0)
        else:
            self.set_comp(self.fluid.comp + self.grid.comp)

        if self.verbose:
            print("[info] model compressibility (comp) was calculated.")

    @_lru_cache(maxsize=None)
    def get_cell_trans(
        self,
        cell_id=None,
        cell_coords=None,
        boundary: bool = False,
    ):
        """Returns transmissibility (T) at all cell faces.

        Parameters
        ----------
        id : int, iterable of int
            cell id based on natural order as int.
        coords : iterable of int
            cell coordinates (i,j,k) as a tuple of int.
        boundary : bool, optional
            include boundary cells.

        Returns
        -------
        ndarray
            array of G based on dir argument.

        """
        # ToDo
        # ----
        # - for now use only with id

        cell_G = self.grid.get_cell_G(cell_id, cell_coords, boundary)
        muB = self.fluid.mu * self.fluid.B
        return {k: v / muB for k, v in cell_G.items()}

    def get_cells_trans(
        self,
        boundary: bool = False,
        sparse: bool = False,
        vectorize: bool = True,
    ):
        """_summary_

        Parameters
        ----------
        boundary : bool, optional
            _description_
        sparse : bool, optional
            _description_

        Returns
        -------
        _type_
            _description_
        """
        if vectorize:
            return self.grid.get_cells_G(boundary, sparse) / (
                self.fluid.mu * self.fluid.B
            )

        n = self.grid.get_n(boundary)
        if sparse:
            T_array = ss.lil_matrix((n, n), dtype=self.dtype)
        else:
            T_array = np.zeros((n, n), dtype=self.dtype)

        if boundary:
            for cell_id in self.grid.cells_id:
                T = self.get_cell_trans(cell_id, None, False)
                for cell_n_id in T.keys():
                    T_array[cell_id, cell_n_id] = T[cell_n_id]
        else:
            for cell_id in self.grid.cells_id:
                cell_i = self.cells_i_dict[cell_id]
                T = self.get_cell_trans(cell_id, None, False)
                cells_n_i = [self.cells_i_dict[x] for x in T.keys()]
                for cell_n_id, cell_n_i in zip(T.keys(), cells_n_i):
                    T_array[cell_i, cell_n_i] = T[cell_n_id]
        return T_array

    def get_cells_trans_diag(self, boundary: bool = False, diag_n=1):
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

    def __calc_well_G(self, cell_id=None):
        """Calculates well Geometry factor (G).

        Parameters
        ----------
        id : int, optional
            cell id based on natural order as int.

        Returns
        -------
        float
            well geometry factor G.

        """
        # ToDo
        # ----
        # - use k and d based on well direction.
        fdir = self.grid.get_fdir()
        if fdir == "x":
            k_H = self.grid.k["x"][cell_id]
        elif fdir == "xy":
            k_H = (self.grid.k["x"][cell_id] * self.grid.k["y"][cell_id]) ** 0.5
        elif fdir == "xyz":
            # print(f"[warning] __calc_well_G at {fdir} has to be verified.")
            k_H = (self.grid.k["x"][cell_id] * self.grid.k["y"][cell_id]) ** 0.5
        else:
            raise ValueError(f"k for fdir='{fdir}' is not defined.")
        G_n = (
            2
            * np.pi
            * self.factors["transmissibility conversion"]
            * k_H
            * self.grid.d["z"][cell_id]
        )

        G_d = np.log(self.wells[cell_id]["r_eq"] / self.wells[cell_id]["r"] * 12)

        if "s" in self.wells[cell_id].keys():
            G_d += self.wells[cell_id]["s"]

        return G_n / G_d

    def __calc_well_r_eq(self, cell_id):
        """Calculates well equivalent radius (r_eq).

        Parameters
        ----------
        id : int, optional
            cell id based on natural order as int.

        Returns
        -------
        float
            well equivalent radius (r_eq).

        """
        # ToDo
        # ----
        # - use k and d based on well direction.
        fdir = self.grid.get_fdir()
        if fdir in ["x", "y"]:
            d = self.grid.d["x"][cell_id] ** 2 + self.grid.d["y"][cell_id] ** 2
            return 0.14 * d**0.5
        elif fdir == "xy":
            kx_ky = self.grid.k["x"][cell_id] / self.grid.k["y"][cell_id]
            ky_kx = self.grid.k["y"][cell_id] / self.grid.k["x"][cell_id]
            return (
                0.28
                * (
                    ky_kx**0.5 * self.grid.d["x"][cell_id] ** 2
                    + kx_ky**0.5 * self.grid.d["y"][cell_id] ** 2
                )
                ** 0.5
                / (ky_kx**0.25 + kx_ky**0.25)
            )
        elif fdir == "xyz":
            # print(f"[warning] __calc_well_r_eq at {fdir} has to be verified.")
            kx_ky = self.grid.k["x"][cell_id] / self.grid.k["y"][cell_id]
            ky_kx = self.grid.k["y"][cell_id] / self.grid.k["x"][cell_id]
            return (
                0.28
                * (
                    ky_kx**0.5 * self.grid.d["x"][cell_id] ** 2
                    + kx_ky**0.5 * self.grid.d["y"][cell_id] ** 2
                )
                ** 0.5
                / (ky_kx**0.25 + kx_ky**0.25)
            )
        else:
            raise ValueError(f"k for fdir='{fdir}' is not defined.")

    def set_well(self, well=None, cell_id=None, q=None, pwf=None, r=None, s=None):
        """Set a well in a specific cell

        Parameters
        ----------
        well : Well class, optional
            well information. If this class was used as input, all other
            arguments will be ignored except id will be used instead of
            well.id.
        id : int, optional
            well location using cell id based on natural order as int.
            This value is given a higher priority over well.id.
        q : int, float, optional
            well rate as positive for injection or negative for
            production
        pwf : int, float, optional
            bottom hole flowing pressure (BHFP). If was not defined,
            None value will be set to zero.
        r : int, float, optional
            well radius.
        s : int, float, optional
            well skin factor

        """
        # ToDo
        # ----
        # - Change production to positive and injection to negative.
        if well is not None:
            if cell_id is None:
                cell_id = well.cell_id
            assert (
                cell_id in self.grid.cells_id
            ), "a well must be placed within the reservoir"
            self.wells[cell_id] = vars(well)
        else:
            assert cell_id is not None, "id must be defined"
            assert (
                cell_id in self.grid.cells_id
            ), "a well must be placed within the reservoir"
            if cell_id not in self.wells:
                self.wells[cell_id] = {}
            if q is not None:
                self.wells[cell_id]["q"] = q
                self.wells[cell_id]["q_sp"] = q
                self.wells[cell_id]["constrain"] = "q"
            if pwf is not None:
                self.wells[cell_id]["pwf"] = pwf
                self.wells[cell_id]["pwf_sp"] = pwf
                if "q" not in self.wells[cell_id].keys():
                    self.wells[cell_id]["constrain"] = "pwf"
                self.w_pressures[cell_id].append(self.pressures[self.tstep, cell_id])
            if "constrain" not in self.wells[cell_id].keys():
                self.wells[cell_id]["constrain"] = None
            if r is not None:
                self.wells[cell_id]["r"] = r
            if s is not None:
                self.wells[cell_id]["s"] = s

        self.wells[cell_id]["r_eq"] = self.__calc_well_r_eq(cell_id)
        self.wells[cell_id]["G"] = self.__calc_well_G(cell_id)
        if "pwf" not in self.wells[cell_id].keys():
            self.wells[cell_id]["pwf"] = 0
            self.wells[cell_id]["pwf_sp"] = 0
            self.w_pressures[cell_id].append(self.pressures[self.tstep, cell_id])

        if self.verbose:
            print(f"[info] a well in cell {cell_id} was set.")

    # -------------------------------------------------------------------------
    # Boundaries:
    # -------------------------------------------------------------------------

    def set_boundary(self, cell_b_id: int, cond: str, v: float):
        """Set a boundary condition in a cell.

        Parameters
        ----------
        id_b : int, optional
            boundary cell id based on natural order as int.
        cond : str
            boundary constant condition. Three conditions are possible:
            (2) Constant rate: str in ['rate', 'q'],
            (1) Constant pressure: str in ['pressure', 'press', 'p'],
            (3) Constant pressure gradient: str in ['gradient', 'grad',
            'g'].
        v : int, float
            constant value to specify the condition in cond argument.

        """
        # ToDo
        # ----
        # - d is taken at x direction for gradient.
        cond = cond.lower()
        if cond in ["rate", "q"]:
            self.rates[self.tstep, cell_b_id] = v
        elif cond in ["pressure", "press", "p"]:
            self.pressures[self.tstep, cell_b_id] = v
        elif cond in ["gradient", "grad", "g"]:
            ((cell_id, T),) = self.get_cell_trans(cell_b_id, None, False).items()
            cell_n = self.grid.get_cell_neighbors(cell_b_id, None, False, "dict")
            dir = [dir for dir in cell_n if cell_id in cell_n[dir]][0]
            self.rates[self.tstep, cell_b_id] = T * self.grid.d[dir][cell_id] * v
        else:
            raise ValueError(f"cond argument {cond} is unknown.")

        if self.verbose:
            print(f"[info] boundary in cell {cell_b_id} was set to constant {cond}.")

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
        for cell_id in bdict:
            assert cell_id in boundaries, f"cell {cell_id} is not a boundary cell."
            cond, v = bdict[cell_id]
            self.set_boundary(cell_id, cond, v)

        self.bdict_update = [
            id_b for id_b in self.bdict.keys() if self.bdict[id_b][0] == "pressure"
        ]

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
        for cell_id in boundaries:
            self.set_boundary(cell_id, cond, v)

    # -------------------------------------------------------------------------
    # Flow Equations:
    # -------------------------------------------------------------------------

    @_lru_cache(maxsize=1)
    def __calc_RHS(self):
        """Calculates flow equation for RHS."""
        # ToDo
        # ----
        # - make sure RHS is suitable in case of floats.
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
    def __calc_n_term(
        self,
        cell_id,
        cell_n_id,
        cell_p,
        trans,
    ) -> float:
        """Calculates neighbor flow term.

        This function calculates the neighbor flow term between
        a specific cell (cell_id) and its neighbor cell (cell_n_id).

        Parameters
        ----------
        cell_id : int
            cell id based on natural order as int.
        cell_n_id : int
            neighbor cell id based on natural order as int.
        cell_p : Symbol
            cell pressure at cell_id.
        trans : float
            transmissibility between cell_id and cell_b_id.

        Returns
        -------
        float
            neighbor flow term (n_term).
        """
        cell_n_p = eval(f"sym.Symbol('p{cell_n_id}')")
        dz = self.grid.z[cell_n_id] - self.grid.z[cell_id]
        return trans * ((cell_n_p - cell_p) - (self.fluid.g * dz))

    @_lru_cache(maxsize=None)
    def __calc_b_term(
        self,
        cell_id,
        cell_b_id,
        cell_p,
        trans,
    ) -> float:
        """Calculates boundary flow term.

        This function calculates the boundary flow term between a
        specific cell (cell_id) and its boundary cell (cell_b_id).

        Parameters
        ----------
        cell_id : int
            cell id based on natural order as int.
        cell_b_id : int
            boundary cell id based on natural order as int.
        cell_p : Symbol
            pressure symbol at cell id.
        trans : float
            transmissibility between cell_id and cell_b_id.

        Returns
        -------
        float
            boundary flow term (b_term).
        """

        # implementation 1:
        # problamatic in case initial pressure is set at boundaries.
        # cell_b_p = self.pressures[self.tstep, cell_b_id]
        # if not np.isnan(cell_b_p):
        #     dz = self.grid.z[cell_b_id] - self.grid.z[cell_id]
        #     b_term = trans * 2 * ((cell_b_p - cell_p) - (self.fluid.g * dz))
        # else:
        #     b_term = self.rates[self.tstep, cell_b_id]

        # implementation 2:
        if cell_b_id in self.bdict:
            cond, v = self.bdict[cell_b_id]
            if cond.lower() in ["pressure", "press", "p"]:
                dz = self.grid.z[cell_b_id] - self.grid.z[cell_id]
                return trans * 2 * ((v - cell_p) - (self.fluid.g * dz))
            else:  # elif cond in ["rate", "q", "gradient", "grad", "g"]:
                return v
        else:
            return 0.0

    def __calc_w_term(
        self,
        cell_id,
        cell_p,
    ) -> float:
        """Calculates well flow term.

        This function calculates the well flow term between a
        specific cell (cell_id) and its well (if exists).

        Parameters
        ----------
        cell_id : int
            cell id based on natural order as int.
        cell_p : Symbol or value
            cell pressure symbol or value at cell id.

        Returns
        -------
        float
            well flow term (w_term).
        """
        if "q" in self.wells[cell_id] and self.wells[cell_id]["constrain"] == "q":
            return self.wells[cell_id]["q"]
        else:
            return (
                -self.wells[cell_id]["G"]
                / (self.fluid.B * self.fluid.mu)
                * (cell_p - self.wells[cell_id]["pwf"])
            )

    def __calc_a_term(
        self,
        cell_id,
        cell_p,
    ):
        """Calculates accumulation term.

        Parameters
        ----------
        cell_id : int
            cell id based on natural order as int.

        Returns
        -------
        float
            accumulation term (a_term).

        Raises
        ------
        ValueError
            Initial pressure was not defined.
        """
        # ToDo
        # ----
        # - consider unifying RHS or if cond.
        if self.comp_type == "incompressible":
            return 0.0
        else:
            try:
                return self.RHS[cell_id] * (
                    cell_p - self.pressures[self.tstep, cell_id]
                )
            except:
                raise ValueError("Initial pressure (pi) must be specified")

    def __simplify_eq(self, cell_eq):
        if (
            cell_eq.lhs.as_coefficients_dict()[1] != 0
            or cell_eq.rhs.as_coefficients_dict()[1] != 0
        ):
            cell_eq = cell_eq.simplify()
        return cell_eq

    def get_cell_eq(self, cell_id):
        """Return cell equation.

        Parameters
        ----------
        id : int, optional
            cell id based on natural order as int.

        Returns
        -------
        tuple
            cell equation as a tuple of (lhs, rhs).

        """
        # ToDo
        # ----
        # - n_term naming.

        # Backup
        # ------
        # - constant pressure:
        #     # exec(f"p{i}=sym.Symbol('p{i}')")
        #     # ToDo: keep pressure constant at specific cell (requires A adjust)
        #     # if not np.isnan(self.pressures[self.tstep][i]):
        #     #     exec(f"p{i} = {self.pressures[self.tstep][i]}")
        # - n_term to use pressure values:
        #     # To Do: keep pressure constant at specific cell (requires A adjust)
        #     # if not np.isnan(self.pressures[self.tstep][neighbor]):
        #     #     exec(f"p{neighbor} = {self.pressures[self.tstep][neighbor]}")
        # - n_term in one calc.
        # exec(
        #     f"n_term = self.T[dir][min(neighbor,id)] * ((p{n_id} - p{id})
        #     - (self.fluid.g * (self.grid.z[neighbor] - self.grid.z[id])))"
        # )
        cell_p = eval(f"sym.Symbol('p{cell_id}')")

        if cell_id not in self.cells_terms:
            assert (
                cell_id in self.grid.cells_id
            ), f"id is out of range {self.grid.cells_id}."
            neighbors = self.grid.get_cell_neighbors(
                id=cell_id, boundary=False, fmt="array"
            )
            boundaries = self.grid.get_cell_boundaries(id=cell_id, fmt="array")
            # f_terms: flow terms, a_terms: accumulation terms
            terms = {"f_terms": [], "a_term": 0}
            T = self.get_cell_trans(cell_id, None, True)
            if self.verbose:
                print(f"[info] cell id: {cell_id}")
                print(f"[info]    - Neighbors: {neighbors}")
                print(f"[info]    - Boundaries: {boundaries}")

            for id_n in neighbors:
                n_term = self.__calc_n_term(cell_id, id_n, cell_p, T[id_n])
                terms["f_terms"].append(n_term)
                if self.verbose:
                    print(f"[info] Neighbor terms: {n_term}")

            for cell_b_id in boundaries:
                b_term = self.__calc_b_term(cell_id, cell_b_id, cell_p, T[cell_b_id])
                terms["f_terms"].append(b_term)
                if self.verbose:
                    print(f"[info] Boundary terms: {b_term}")

            if cell_id in self.wells.keys():
                w_term = self.__calc_w_term(cell_id, cell_p)
                terms["f_terms"].append(w_term)
                if self.verbose:
                    print(f"[info] Well terms: {w_term}")

            terms["a_term"] = self.__calc_a_term(cell_id, cell_p)
            if self.verbose:
                print("[info] Accumulation term:", terms["a_term"])

            self.cells_terms[cell_id] = terms
            if self.verbose:
                print("[info] terms:", terms)

        else:
            terms = self.cells_terms[cell_id]
            if (
                cell_id in self.wells.keys()
                and self.wells[cell_id]["constrain"] == "pwf"
            ):
                w_term = self.__calc_w_term(cell_id, cell_p)
                if self.verbose:
                    print(f"[info] Well terms (updated): {w_term}")
                terms["f_terms"][-1] = w_term
            if self.comp_type == "compressible":
                terms["a_term"] = self.__calc_a_term(cell_id, cell_p)

        cell_eq = sym.Eq(sum(terms["f_terms"]), terms["a_term"])
        cell_eq = self.__simplify_eq(cell_eq)
        cell_eq_lhs = cell_eq.lhs.as_coefficients_dict()

        if self.verbose:
            print(f"[info] Flow equation {cell_id}:", cell_eq)

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

    def __update_matrices(self, cell_id):
        """Update flow equations' matrices (A, d).

        Parameters
        ----------
        id : int
            cell id based on natural order as int.

        Notes
        -----
        - arrays for lhs and rhs:
            self.d[i] = np.array(cell_rhs).astype(self.dtype)
            self.A[i, ids] = np.array(list(cell_lhs.values())).astype(self.dtype)
        - finding cell i:
            ids = [self.cells_id.index(int(str(s)[1:])) for s in cell_lhs.keys()]
        """
        cell_lhs, cell_rhs = self.cells_eq[cell_id]
        ids = [self.cells_i_dict[int(str(s)[1:])] for s in cell_lhs.keys()]
        self.d[self.cells_i_dict[cell_id]] = cell_rhs
        self.A[self.cells_i_dict[cell_id], ids] = list(cell_lhs.values())
        if self.verbose:
            print(f"[info] cell id: {cell_id}")
            print(f"[info]      - ids: {ids}")
            print(f"[info]      - lhs: {cell_lhs}")
            print(f"[info]      - rhs: {cell_rhs}")

    def init_matrices(self, sparse=False, threading=False):
        """Initialize flow equations' matrices (A, d).

        Parameters
        ----------
        sparse : bool, optional
            use sparse matrices instead of dense matrices.
        threading : bool, optional
            use multiple threads for concurrence workers. The maximum
            number of threads are set to the half number of grids.

        Returns
        -------
        _type_
            _description_

        """
        # ToDo
        # ----
        # - Update only required cells.
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
            for cell_id in self.cells_id:
                self.__update_matrices(cell_id)

        if self.verbose:
            print("[info] - A:\n", self.A)
            print("[info] - d:\n", self.d)

        return self.A, self.d

    # -------------------------------------------------------------------------
    # Matrices: vectorized
    # -------------------------------------------------------------------------

    def init_A(self, sparse: bool = False):
        """Initialize ceofficient matrix (`A`).

        For a system of linear equations `Au=d`, `A` is the
        ceofficient matrix (known), `d` is the constant vector (known),
        and `u` is the variable vector (unknown e.g., pressure).

        This function initialize the ceofficient matrix (`A`) which is
        needed only at initial timestep (i.e., `timestep=0`).

        Parameters
        ----------
        sparse : bool, optional
            use sparse matrices instead of dense matrices.

        Returns
        -------
        ndarray
            ceofficient A is initialized in place.
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
        self.A_ = self.get_cells_trans(False, sparse, True)
        v1 = -self.A_[self.cells_i, :].sum(axis=1).flatten()
        v2 = self.RHS[self.cells_id].flatten()
        v3 = v1 - v2
        self.A_[self.cells_i, self.cells_i] = v3
        if sparse:
            self.A_ = ss.lil_matrix(self.A_, dtype=self.dtype)
        return self.A_

    def init_d(self, sparse):
        """Initialize constant vector (`d`).

        For a system of linear equations `Au=d`, `A` is the
        ceofficient matrix (known), `d` is the constant vector (known),
        and `u` is the variable vector (unknown e.g., pressure).

        This function initialize the constant vector (`d`) which is
        needed at every timestep in case of a compressible system. In
        case of an incompressible system, a constant zero vector is
        used.

        Parameters
        ----------
        sparse : bool, optional
            _description_

        Returns
        -------
        ndarray
            vector d is initialized in place and can be accessed by
            `self.d_`.

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

        return self.d_

    def __update_z(self):
        """_summary_"""
        # ToDo
        # ----
        # - T for different geometries is still not ready.

        # all 1D in x direction.
        z = self.grid.z[self.grid.cells_id]
        if not np.all(z == z[0]):
            z_l = np.append(z[1:], np.nan)
            z_u = np.append(np.nan, z[:-1])
            dz_l = self.fluid.g * np.nan_to_num(z_l - z)
            dz_u = self.fluid.g * np.nan_to_num(z_u - z)
            # T = self.T["x"][self.grid.cells_id]
            # T = np.diag(self.get_cells_T(True, False), 1)[self.cells_i]
            T = self.get_cells_trans_diag(True, 1)[self.cells_id]
            v = T * dz_l + T * dz_u
            self.d_ += v.reshape(-1, 1)

    def get_matrices_vectorized(self, sparse=False, threading=False):
        """_summary_

        Parameters
        ----------
        sparse : bool, optional
            _description_
        threading : bool, optional
            _description_

        Returns
        -------
        _type_
            _description_
        """
        update_z = False
        if self.tstep == 0:
            self.resolve = defaultdict(lambda: False)
            self.init_A(sparse)
            self.bdict_v = {}
            for id_b in self.bdict_update:
                ((id, T),) = self.get_cell_trans(id_b, None, False).items()
                p = eval(f"sym.Symbol('p{id}')")
                b_term = self.__calc_b_term(id, id_b, p, T)
                v0, v1 = b_term.as_coefficients_dict().values()
                self.bdict_v[id_b] = (v0, v1, id)
                self.A_[self.cells_i_dict[id], self.cells_i_dict[id]] += v1
                # self.A_[self.cells_i_dict[id], self.cells_i_dict[id]] -= T * 2

        self.init_d(sparse)

        for id in self.wells.keys():
            if self.wells[id]["constrain"] == "q":
                w_term = self.__calc_w_term(id, self.pressures[self.tstep, id])
                self.d_[self.cells_i_dict[id], 0] -= w_term
                update_z = True
            elif self.wells[id]["constrain"] == "pwf":
                p = eval(f"sym.Symbol('p{id}')")
                w_term = self.__calc_w_term(id, p)
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
            id = self.grid.get_cell_neighbors(id_b, None, False, "list")
            if len(id) > 0:
                id = id[0]
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
            pressures = self.pressures[self.tstep, self.grid.cells_id]
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


        Notes
        -----
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
            q_est = self.__calc_w_term(id, self.pressures[self.tstep, id])
            self.wells[id]["q"] = self.rates[self.tstep, id] = q_est

            if resolve:
                return True
            else:
                tstep_w_pressures[id] = pwf_est

        for id in self.wells.keys():
            self.w_pressures[id].append(tstep_w_pressures[id])

        return False

    def __update_boundaries(self):
        for id_b in self.bdict_update:
            ((id_n, T),) = self.get_cell_trans(id_b, None, False).items()
            p_n = self.pressures[self.tstep, id_n]
            b_terms = self.__calc_b_term(id_n, id_b, p_n, T)
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
        sparse : bool, optional
            _description_
        threading : bool, optional
            _description_
        vectorize : bool, optional
            _description_
        check_MB : bool, optional
            _description_
        update : bool, optional
            _description_
        print_arrays : bool, optional
            _description_
        isolver : str, optional
            iterative solver for sparse matrices. Available solvers are
            `["bicg", "bicgstab", "cg", "cgs", "gmres", "lgmres",
            "minres", "qmr", "gcrotmk", "tfqmr"]`.
            If None, direct solver is used. Only relevant when argument
            sparse=True. Option "cgs" is recommended to increase
            performance while option "minres" is not recommended due to
            high MB error.

        Notes
        -----
        Direct solutions can also be obtained using matrix dot product
        (usually slower) as following:

        >>> pressures = np.dot(np.linalg.inv(A), d).flatten()
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
            A, d = A.tocsc(), d.todense()
            if isolver:
                solver = utils.solvers.get_isolver(isolver)
                pressures, exit_code = solver(
                    A,
                    d,
                    atol=0,
                    # x0=self.pressures[self.tstep, self.cells_id],
                )
                assert exit_code == 0, "unsuccessful convergence"
            else:
                pressures = ssl.spsolve(A, d, use_umfpack=True)
            A = A.todense()
        else:
            pressures = sl.solve(A, d).flatten()

        if update:
            self.tstep += 1
            self.pressures = np.vstack([self.pressures, self.pressures[-1]])
            self.pressures[self.tstep, self.grid.cells_id] = pressures
            self.rates = np.vstack([self.rates, self.rates[-1]])
            self.As = np.vstack([self.As, A.reshape(1, -1)])
            self.ds = np.vstack([self.ds, d.reshape(1, -1)])
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
        isolver=None,
    ):
        """Perform a simulation run for nsteps.

        Parameters
        ----------
        nsteps : int, optional
            _description_
        sparse : bool, optional
            _description_
        threading : bool, optional
            _description_
        check_MB : bool, optional
            _description_
        isolver : str, optional
            iterative solver for sparse matrices. Available solvers are
            ["bicg", "bicgstab", "cg", "cgs", "gmres", "lgmres",
            "minres", "qmr", "gcrotmk", "tfqmr"].
            If None, direct solver is used. Only relevant when argument
            sparse=True. Direct solver is recommended for more accurate
            calculations. To improve performance, "cgs" is recommended
            to increase performance while option "minres" is not recommended due to
            high MB error. For more information check [1][2].

        References
        ----------
        - scipy: `Solving Linear Problems <https://docs.scipy.org/doc/scipy/reference/sparse.linalg.html#solving-linear-problems>`_.
        - scipy: `Iterative Solvers <https://scipy-lectures.org/advanced/scipy_sparse/solvers.html#iterative-solvers>`_.
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
        pbar = tqdm(
            range(1, nsteps + 1),
            unit="steps",
            colour="green",
            position=0,
            leave=True,
        )

        for step in pbar:
            self.solve(
                sparse,
                threading,
                vectorize,
                check_MB,
                True,
                print_arrays,
                isolver,
            )
            pbar.set_description(f"[step] {step}")

        self.run_ctime = round(time.time() - start_time, 2)
        self.ctime += self.run_ctime
        print(
            f"[info] Simulation run of {nsteps} steps",
            f"finished in {self.run_ctime} seconds.",
        )
        if check_MB:
            print(f"[info] Material Balance Error: {self.error}.")

        if verbose_restore:
            self.verbose = True

    # -------------------------------------------------------------------------
    # Material Balance:
    # -------------------------------------------------------------------------

    def check_MB(self, verbose=False, error_threshold=0.1):
        """Material Balance Check

        Parameters
        ----------
        verbose : bool, optional
            _description_
        error_threshold : float, optional
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
                f"[warning] Material balance error ({self.error})",
                f"in step {self.tstep}",
                f"is higher than the allowed error ({error_threshold}).",
            )

    # -------------------------------------------------------------------------
    # Data:
    # -------------------------------------------------------------------------

    def __concat(self, data, df):
        if df is not None:
            df = pd.concat([df, data], axis=1)
            return df
        return data

    def __add_time(self, units, melt, boundary, scale, df=None):
        if units:
            if scale:
                time_str = " [Scaled]"
            else:
                time_str = f" [{self.units['time']}]"
        else:
            time_str = ""
        time = self.__get_time_vector()
        if scale:
            time = self.time_scaler.transform(time).flatten()
        if melt:
            n_cells = self.grid.get_n(boundary)
            time = np.repeat(time, n_cells)
        data = pd.Series(time, name="Time" + time_str)
        return self.__concat(data, df)

    def __add_date(self, units, melt, boundary, df=None):
        if units:
            date_str = " [d.m.y]"
        else:
            date_str = ""
        date_series = pd.date_range(
            start=self.start_date,
            periods=self.tstep + 1,
            freq=str(self.dt) + "D",
        ).strftime("%d.%m.%Y")
        data = pd.Series(date_series, name="Date" + date_str)
        if melt:
            n_cells = self.grid.get_n(boundary)
            data = data.repeat(n_cells).reset_index(drop=True)
        return self.__concat(data, df)

    def __add_cells_rate(self, units, melt, boundary, scale, df=None):
        if units:
            if scale:
                rate_str = " [Scaled]"
            else:
                rate_str = f" [{self.units['rate']}]"
        else:
            rate_str = ""
        cells_id = self.grid.get_cells_id(boundary, False, "array")
        array = self.rates[:, cells_id]
        if scale:
            array = self.rates_scaler.transform(array)
        if melt:
            labels = ["Q" + rate_str]
            array = array.flatten()
        else:
            labels = [f"Q{str(id)}" + rate_str for id in cells_id]
        data = pd.DataFrame(array, columns=labels)
        return self.__concat(data, df)

    def __add_cells_pressures(self, units, melt, boundary, scale, df=None):
        if units:
            if scale:
                press_str = " [Scaled]"
            else:
                press_str = f" [{self.units['pressure']}]"
        else:
            press_str = ""
        cells_id = self.grid.get_cells_id(boundary, False, "array")
        array = self.pressures[:, cells_id]
        if scale:
            array = self.pressures_scaler.transform(array)
        if melt:
            labels = ["P" + press_str]
            array = array.flatten()
        else:
            labels = [f"P{str(id)}" + press_str for id in cells_id]
        data = pd.DataFrame(array, columns=labels)
        return self.__concat(data, df)

    def __add_wells_rate(self, units, scale, df=None):
        if units:
            if scale:
                rate_str = " [Scaled]"
            else:
                rate_str = f" [{self.units['rate']}]"
        else:
            rate_str = ""
        labels = [f"Qw{str(id)}" + rate_str for id in self.wells.keys()]
        array = self.rates[:, list(self.wells.keys())]
        if scale:
            array = self.rates_scaler.transform(array)
        data = pd.DataFrame(array, columns=labels)
        return self.__concat(data, df)

    def __add_wells_pressures(self, units, scale, df=None):
        if units:
            if scale:
                press_str = " [Scaled]"
            else:
                press_str = f" [{self.units['pressure']}]"
        else:
            press_str = ""
        labels = [f"Pwf{str(id)}" + press_str for id in self.w_pressures]
        data = pd.DataFrame(self.w_pressures)
        if scale:
            data = self.pressures_scaler.transform(data.values)
            data = pd.DataFrame(data)
        data.columns = labels
        return self.__concat(data, df)

    def __add_xyz(self, boundary, melt, scale, df=None):
        if melt:
            cells_center, fdir = self.get_centers(scale, boundary)
            array = np.tile(cells_center.flatten(), self.nsteps).reshape(-1, len(fdir))
            data = pd.DataFrame(array, columns=fdir)
            return self.__concat(data, df)
        return df

    def __get_time_vector(self):
        return np.arange(0, (self.tstep + 1) * self.dt, self.dt)

    def __update_time_scaler(self):
        time = self.__get_time_vector()
        self.time_scaler.fit(time, axis=0)

    def __update_space_scaler(self, boundary):
        cells_center, _ = self.get_centers(False, boundary)
        self.space_scaler.fit(cells_center, axis=0)

    def __update_pressures_scaler(self, boundary):
        config = self.__get_scalers_config(boundary)
        pressures = self.get_df(columns=["cells_pressure"], **config).values
        self.pressures_scaler.fit(pressures, axis=None)

    def __update_rates_scaler(self, boundary):
        """Flow rates scaler

        Parameters
        ----------
        boundary : bool, optional
            include boundary cells.

        Notes
        -----
        Disabled:
            This scaler is disabled since the scaling argument for both
            __add_cells_rate() and __add_wells_rate() is set to False.
        Error:
            When boundary argument is False but no wells, the scaler
            will fail since rows with zeros and nans are removed, no
            rates will remain for scaling (i.e. empty array).
        """
        config = self.__get_scalers_config(boundary)
        rates = self.get_df(columns=["rates"], **config).values
        self.rates_scaler.fit(rates, axis=None)

    def __get_scalers_config(self, boundary):
        return {
            "boundary": boundary,
            "units": False,
            "melt": False,
            "scale": False,
            "save": False,
            "drop_nan": True,
            "drop_zero": True,
        }

    def update_scalers(self, boundary):
        self.__update_time_scaler()  # get_time
        self.__update_space_scaler(boundary)  # get_centers() and __add_xyz()
        self.__update_pressures_scaler(boundary)
        self.__update_rates_scaler(boundary)

    def get_time(self, scale=False):
        time = self.__get_time_vector().reshape(-1, 1)
        if scale:
            self.__update_time_scaler()
            time = self.time_scaler.transform(time)
        return time

    def get_centers(self, scale=False, boundary: bool = True):
        def get_fdir_cols(s):
            if s == "x":
                return 0
            elif s == "y":
                return 1
            elif s == "z":
                return 2
            else:
                return None

        centers = self.grid.get_cells_center(boundary, False, False)
        fdir = list(self.grid.get_fdir())
        fdir_cols = list(map(get_fdir_cols, fdir))
        centers = centers[:, fdir_cols]
        if scale:
            self.__update_space_scaler(True)
            centers = self.space_scaler.transform(centers)
        return centers, fdir

    def get_domain(self, scale, boundary):
        t = self.get_time(scale)
        centers, _ = self.get_centers(scale, boundary)
        return t, centers

    def set_scalers(self, scalers_dict: dict):
        """Set scalers configuration

        To change the scaling settings. Current settings can be shown
        under scalers_dict. By default the following settings are used:

        .. highlight:: python
        .. code-block:: python

            scalers_dict = {
                'time':['MinMaxScaler', (0,1)],
                'space':['MinMaxScaler', (-1,1)],
                'pressure':['MinMaxScaler', (-1,1)],
                'rate':[None,None],
            }

        Note that be default rates are not scaled, time is scaled
        between 0 and 1, while space and pressure are scaled between
        -1 and 1. By default, MinMaxScaler is used for all dimensions.

        Parameters
        ----------
        scalers_dict : dict
            scalers setting as dict in the following
            format: {'time': [scaler_type, scaler_range]} were
            scaler_type can be e.g. 'MinMaxScaler' and scaler_range is a
            tuple of output_range of the scaler e.g. (-1,1). The keys
            must be in ['time', 'space', 'pressure', 'rate'].

        """

        def create_scaler(scaler_type, output_range):
            if scaler_type is None or output_range is None:
                return scalers.Dummy(None), None
            elif scaler_type.lower() in ["minmax", "minmaxscaler"]:
                return scalers.MinMax(output_range=output_range), "MinMax"
            else:
                raise ValueError("scaler type is not defined.")

        col_dict = {
            "time": ["t", "time", "all"],
            "space": ["space", "spatial", "xyz", "x", "y", "z", "all"],
            "pressure": ["p", "pressures", "pressure", "all"],
            "rate": ["q", "rates", "rate", "all"],
        }
        col_vals = sum(col_dict.values(), [])

        for column in scalers_dict.keys():
            scaler_type = scalers_dict[column][0]
            output_range = scalers_dict[column][1]
            if column.lower() in col_dict["time"]:
                self.time_scaler, s_str = create_scaler(scaler_type, output_range)
                column_str = "time"
            elif column.lower() in col_dict["space"]:
                self.space_scaler, s_str = create_scaler(scaler_type, output_range)
                column_str = "space"
            elif column.lower() in col_dict["pressure"]:
                self.pressures_scaler, s_str = create_scaler(scaler_type, output_range)
                column_str = "pressure"
            elif column.lower() in col_dict["rate"]:
                self.rates_scaler, s_str = create_scaler(scaler_type, output_range)
                column_str = "rate"
            else:  # if column.lower() not in col_vals:
                raise ValueError(f"column {column} does not exist.")

            if scaler_type is None or output_range is None:
                self.scalers_dict[column_str] = [None, None]
            else:
                self.scalers_dict[column_str] = [s_str, scalers_dict[column][1]]

    def get_df(
        self,
        columns: list = ["time", "cells", "wells"],
        boundary: bool = True,
        units: bool = False,
        melt: bool = False,
        scale: bool = False,
        save: bool = False,
        drop_nan: bool = True,
        drop_zero: bool = True,
    ):
        """Returns simulation data as a dataframe.

        Parameters
        ----------
        columns : list[str], optional
            selected columns to be added to the dataframe. The following
            options are available:
            "time": for time steps as specified in dt.
            "date": for dates as specified in dt and start_date.
            "q", "rates": for all (cells and wells) rates.
            "p", "pressures": for all (cells and wells) pressures.
            "cells": for all cells rates and pressures.
            "wells": for all wells rates and pressures.
            "cells_rate": for all cells rates (including wells' cells).
            "cells_pressure": for all cells pressures.
            "wells_rate": for all wells rates.
            "wells_pressure": for all wells pressures.
        boundary : bool, optional
            include boundary cells.
            It is only relevant when cells columns are selected.
        units : bool, optional
            column names with units (True) or without units (False).
        melt : bool, optional
            to melt columns of the same property to one column. By
            default, cells id, xyz (based on grid fdir), step columns
            are included while wells columns (wells_rate,
            wells_pressure) are ignored.
        scale : bool, optional
            scale time, space (x, y, z), rates, and pressures. To change
            the scaling settings use `set_scalers() </api/reservoirflow.models.html#reservoirflow.models.BlackOil.set_scalers>`_.
            Current settings can be shown under scalers_dict.
            By default:

            .. highlight:: python
            .. code-block:: python

                scalers_dict = {
                    'time':['MinMaxScaler', (0,1)],
                    'space':['MinMaxScaler', (-1,1)],
                    'pressure':['MinMaxScaler', (-1,1)],
                    'rate':[None,None]
                }

        save : bool, optional
            save output as a csv file.
        drop_nan : bool, optional
            drop columns which contain only nan values if melt is False.
            if melt is True, drop rows which contain any nan values.
        drop_zero : bool, optional
            drop columns which contain only zero values. This argument
            is ignored if melt is True.

        Returns
        -------
        DataFrame :
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
        col_vals = sum(col_dict.values(), [])

        df = pd.DataFrame()

        if scale:
            self.update_scalers(True)

        if melt:
            n_cells = self.grid.get_n(boundary)
            cells_id = self.grid.get_cells_id(boundary, False, "array")
            df["id"] = np.tile(cells_id, self.nsteps)
            df["Step"] = np.repeat(np.arange(self.nsteps), n_cells)
            df = self.__add_xyz(boundary, melt, scale, df)

        for c in columns:
            if c.lower() not in col_vals:
                raise ValueError(f"column {c} does not exist.")
            if c.lower() in col_dict["time"]:
                df = self.__add_time(units, melt, boundary, scale, df)
            if c.lower() in col_dict["date"]:
                df = self.__add_date(units, melt, boundary, df)
            if c.lower() in col_dict["cells_rate"]:
                df = self.__add_cells_rate(units, melt, boundary, scale, df)
            if c.lower() in col_dict["cells_pressure"]:
                df = self.__add_cells_pressures(units, melt, boundary, scale, df)
            if c.lower() in col_dict["wells_rate"] and not melt:
                df = self.__add_wells_rate(units, scale, df)
            if c.lower() in col_dict["wells_pressure"] and not melt:
                df = self.__add_wells_pressures(units, scale, df)

        if melt:
            if drop_nan:
                df = df.dropna(axis=0, how="any")
                df.reset_index(drop=True, inplace=True)
        else:
            if drop_nan:
                df = df.dropna(axis=1, how="all")
            if drop_zero:
                df = df.loc[:, (df != 0).any(axis=0)]
            df.index.name = "Step"
        df.index = df.index.astype(int, copy=False)

        if save:
            df.to_csv("model_data.csv")
            print("[info] Model data was successfully saved.")

        return df

    # -------------------------------------------------------------------------
    # Visualization:
    # -------------------------------------------------------------------------

    show = utils.plots.show_model
    save_gif = utils.plots.save_gif

    # -------------------------------------------------------------------------
    # Plotting:
    # -------------------------------------------------------------------------

    def plot(self, prop: str = "pressures", id: int = None, tstep: int = None):
        """Show values in a cartesian plot.

        Parameters
        ----------
        prop : str, optional
            property name from ["rates", "pressures"].
        id : int, optional
            cell id. If None, all cells are selected.
        tstep : int, optional
            time step. If None, the last time step is selected.
        """
        if tstep is None:
            tstep = self.tstep

        if id is not None:
            exec(f"plt.plot(self.{prop}[:, id].flatten())")
            plt.xlabel("Days")
        elif tstep is not None:
            exec(f"plt.plot(self.{prop}[tstep, :].flatten())")
            plt.xlabel("Grid (id)")
            boundary = True
            x_skip = 5
            x_ticks = np.arange(0, self.grid.get_n(boundary), x_skip)
            x_labels = self.grid.get_cells_id(boundary, fshape=False, fmt="array")
            x_labels = x_labels[::x_skip]
            plt.xticks(ticks=x_ticks, labels=x_labels)
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

    def copy(self):
        """Copy model (under development)

        Returns:
            _type_: _description_
        """
        # https://stackoverflow.com/questions/48338847/how-to-copy-a-python-class-instance-if-deepcopy-does-not-work
        copy_model = _Model(
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
        grid = grids.RegularCartesian(
            nx=4, ny=1, nz=1, dx=300, dy=350, dz=40, phi=0.27, kx=270, dtype="double"
        )
        fluid = fluids.SinglePhase(mu=0.5, B=1, dtype="double")
        model = BlackOil(grid, fluid, dtype="double", verbose=False)
        model.set_well(cell_id=4, q=-600, s=1.5, r=3.5)
        model.set_boundaries({0: ("pressure", 4000), 5: ("rate", 0)})
        return model

    def create_model():
        grid = grids.RegularCartesian(
            nx=10,
            ny=10,
            nz=3,
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
        model = BlackOil(
            grid, fluid, pi=6000, dt=5, start_date="10.10.2018", dtype="double"
        )

        cells_id = model.grid.get_cells_id(False, True)[-1].flatten()
        wells = np.random.choice(cells_id, 6, False)
        for id in wells:
            model.set_well(cell_id=id, q=-300, pwf=100, s=1.5, r=3.5)

        wells = np.random.choice(cells_id, 6, False)
        for id in wells:
            model.set_well(cell_id=id, q=100, s=0, r=3.5)

        return model

    model = create_model()
    model.run(10, sparse=True, vectorize=True, isolver="cgs")
    # model.run(10, isolver=None)
    model.show("pressures")
    model = create_model()
    model.run(10, sparse=False, vectorize=True, isolver=None)
    # model.run(10, isolver=None)
    model.show("pressures")
    model = create_model()
    model.run(10, sparse=True, vectorize=False, isolver="cgs")
    # model.run(10, isolver=None)
    model.show("pressures")
    model = create_model()
    model.run(10, sparse=False, vectorize=False, isolver=None)
    # model.run(10, isolver=None)
    model.show("pressures")
