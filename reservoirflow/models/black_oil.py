"""
BlackOil
========
"""

import warnings
from collections import defaultdict
from datetime import date

# import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as ss

# import reservoirflow as rf
from reservoirflow.grids import Grid
from reservoirflow.fluids import Fluid
from reservoirflow import scalers, utils, plots
from reservoirflow.models.model import Model
from reservoirflow.utils.helpers import _lru_cache


class BlackOil(Model):
    """BlackOil model class.

    Returns
    -------
    Model
        Model object.
    """

    name = "BlackOil"

    def __init__(
        self,
        grid: Grid,
        fluid: Fluid,
        well=None,  # wells.Well = None,
        pi: int = None,
        dt: int = 1,
        start_date: date = None,
        unit: str = "field",
        dtype: str = "double",
        verbose: bool = False,
    ):
        """Create BlackOil Model.

        Parameters
        ----------
        grid : Grid
            Grid object.
        fluid : Fluid
            Fluid object.
        well : Well, optional
            Well object. Wells can be added latter using ``set_well()``
            method.
        pi : int, optional
            Initial reservoir pressure.
        dt : int, optional
            Time step duration.
        start_date : date, optional
            Start date of the simulation run. If None, today's date is
            used.
        unit : str ('field', 'metric', 'lab'), optional
            unit used in input and output. Both `units` and `factors`
            attributes will be updated based on the selected `unit` and
            can be accessed directly from this class.
        dtype : str or `np.dtype`, optional
            data type used in all arrays. Numpy dtype such as
            `np.single` or `np.double` can be used.
        verbose : bool, optional
            print information for debugging.

        Notes
        -----
        .. note::
            Both attributes units and factors are defined based on `unit`
            argument, for more details, check
            `Units & Factors </user_guide/units_factors/units_factors.html>`_.
            For definitions, check
            `Glossary </user_guide/glossary/glossary.html>`_.

        .. attention::
            deciding sparsity must be at the initialization step and not
            in ``solve()`` or ``run()`` methods because this will require
            that the ceofficient matrix must be reinitialized each time
            sparsity is changed which requires tracking this change.
        """
        super().__init__(unit, dtype, verbose)
        self.grid = grid
        self.fluid = fluid
        assert self.dtype == grid.dtype, "grid dtype is not compatible."
        assert self.dtype == fluid.dtype, "fluid dtype is not compatible."

        self.cells_terms = {}
        # newtest:
        self.dt = dt
        # self.nsteps = 1
        # self.tstep = 0
        # self.ctime = 0

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
        start_date : date, optional
            Start date of the simulation run. If None, today's date is
            used.
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
        self.n = self.grid.get_n(False)
        self.cells_i = self.grid.get_cells_i(False)
        self.cells_id = self.grid.get_cells_id(False, False, "array")
        self.cells_i_dict = dict(zip(self.cells_id, self.cells_i))
        self.boundaries_id = self.grid.get_boundaries("id", "array")

        # newtest
        # ones = self.grid.get_ones(True, False, False)[np.newaxis]
        # self.pressures = ones * np.nan
        # self.rates = self.grid.get_zeros(True, False, False)[np.newaxis]
        self.init_pressures, self.init_rates = self.__get_arrays()

        self.pi = pi
        if pi is not None:
            # newtest
            # self.pressures[0, self.grid.cells_id] = pi
            self.init_pressures[0, self.grid.cells_id] = pi
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

        self.scalers_dict = {
            "time": ["MinMax", (0, 1)],
            "space": ["MinMax", (-1, 1)],
            "pressure": ["MinMax", (-1, 1)],
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
            include boundary cells.

        Returns
        -------
        tuple
            tuple as (number of time steps, number of girds)
        """
        assert self.solution is not None, "Model was not compiled."
        return (self.solution.nsteps, self.grid.get_n(boundary))

    def __get_arrays(self) -> tuple:
        ones = self.grid.get_ones(True, False, False)[np.newaxis]
        pressures = ones * np.nan
        rates = self.grid.get_zeros(True, False, False)[np.newaxis]
        return pressures, rates

    def get_init_arrays(self) -> tuple:
        """Initialization arrays.

        Returns
        -------
        tuple
            tuple as (pressures_array, rates_array)
        """
        return self.init_pressures.copy(), self.init_rates.copy()

    def __get_tstep(self):
        if self.solution is None:
            return 0
        else:
            return self.solution.tstep

    def __get_pressure(self, tstep, cell_id):
        if tstep == 0:
            return self.init_pressures[tstep, cell_id]
        else:
            assert self.solution is not None, "Model was not compiled."
            if tstep is None:
                return self.solution.pressures[:, cell_id]
            else:
                return self.solution.pressures[tstep, cell_id]

    def __set_pressure(self, tstep, cell_id, value):
        if tstep == 0:
            self.init_pressures[tstep, cell_id] = value
        else:
            assert self.solution is not None, "Model was not compiled."
            if tstep is None:
                self.solution.pressures[:, cell_id] = value
            else:
                self.solution.pressures[tstep, cell_id] = value

    def __get_rate(self, tstep, cell_id):
        if tstep == 0:
            return self.init_rates[tstep, cell_id]
        else:
            assert self.solution is not None, "Model was not compiled."
            if tstep is None:
                return self.solution.rates[:, cell_id]
            else:
                return self.solution.rates[tstep, cell_id]

    def __set_rate(self, tstep, cell_id, value):
        if tstep == 0:
            self.init_rates[tstep, cell_id] = value
        else:
            assert self.solution is not None, "Model was not compiled."
            if tstep is None:
                self.solution.rates[:, cell_id] = value
            else:
                self.solution.rates[tstep, cell_id] = value

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
            include boundary cells.
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
                # newtest
                # self.w_pressures[cell_id].append(self.pressures[self.get_tstep(), cell_id])
                self.w_pressures[cell_id].append(
                    self.__get_pressure(self.__get_tstep(), cell_id)
                )
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
            # newtest
            # self.w_pressures[cell_id].append(self.pressures[self.get_tstep(), cell_id])
            self.w_pressures[cell_id].append(
                self.__get_pressure(self.__get_tstep(), cell_id)
            )

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
            # newtest
            # self.rates[self.get_tstep(), cell_b_id] = v
            self.__set_rate(self.__get_tstep(), cell_b_id, v)
        elif cond in ["pressure", "press", "p"]:
            # newtest
            # self.pressures[self.get_tstep(), cell_b_id] = v
            self.__set_pressure(self.__get_tstep(), cell_b_id, v)
        elif cond in ["gradient", "grad", "g"]:
            ((cell_id, T),) = self.get_cell_trans(cell_b_id, None, False).items()
            cell_n = self.grid.get_cell_neighbors(cell_b_id, None, False, "dict")
            dir = [dir for dir in cell_n if cell_id in cell_n[dir]][0]
            # newtest
            # self.rates[self.get_tstep(), cell_b_id] = T * self.grid.d[dir][cell_id] * v
            self.__set_rate(
                self.__get_tstep(), cell_b_id, T * self.grid.d[dir][cell_id] * v
            )
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

    @_lru_cache(maxsize=1)
    def get_alpha(
        self,
        method="mean",
    ):
        """Calculates alpha parameter.
        This parameter is used to calculate alpha constant.

        Note
        ----
        - This method is used only for 1D in x-direction flow.

        Parameters
        ----------
        method : str, optional, default "mean"
            method to calculate alpha parameter. The following methods
            are available:
            - first value: `"first"`.
            - last value: `"last"`.
            - mean value: `"mean"`.
            - vector or array: `None`.

        Returns
        -------
        float
            alpha parameter.
        """
        fdir = self.grid.get_fdir()
        if fdir == "x":
            k = self.grid.k["x"]
        else:
            raise ValueError(f"k for fdir='{fdir}' is not defined.")

        lhs_f = (self.factors["transmissibility conversion"] * k) / (
            self.fluid.mu * self.fluid.B
        )
        rhs_f = (self.grid.phi * self.comp) / (
            self.factors["volume conversion"] * self.fluid.B
        )
        alpha = lhs_f / rhs_f

        if method == "mean":
            self.alpha = np.mean(alpha)
        elif method == "first":  # or np.all(alpha == alpha[0]):
            self.alpha = alpha[0]
        elif method == "last":
            self.alpha = alpha[-1]
        elif method in [None, "vector", "array"]:
            pass
        else:  #
            raise ValueError(f"k for fdir='{fdir}' is not defined.")

        return self.alpha

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
            if scale and self.scalers_dict["time"][0] is not None:
                time_str = " [scaled]"
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
            periods=self.__get_tstep() + 1,
            freq=str(self.dt) + "D",
        ).strftime("%d.%m.%Y")
        data = pd.Series(date_series, name="Date" + date_str)
        if melt:
            n_cells = self.grid.get_n(boundary)
            data = data.repeat(n_cells).reset_index(drop=True)
        return self.__concat(data, df)

    def __add_cells_rate(self, units, melt, boundary, scale, df=None):
        if units:
            if scale and self.scalers_dict["rate"][0] is not None:
                rate_str = " [scaled]"
            else:
                rate_str = f" [{self.units['rate']}]"
        else:
            rate_str = ""
        cells_id = self.grid.get_cells_id(boundary, False, "array")
        # newtest
        # array = self.rates[:, cells_id]
        array = self.__get_rate(None, cells_id)
        if scale:
            array = self.rate_scaler.transform(array)
        if melt:
            labels = ["Q" + rate_str]
            array = array.flatten()
        else:
            labels = [f"Q{str(id)}" + rate_str for id in cells_id]
        data = pd.DataFrame(array, columns=labels)
        return self.__concat(data, df)

    def __add_cells_pressures(self, units, melt, boundary, scale, df=None):
        if units:
            if scale and self.scalers_dict["pressure"][0] is not None:
                press_str = " [scaled]"
            else:
                press_str = f" [{self.units['pressure']}]"
        else:
            press_str = ""
        cells_id = self.grid.get_cells_id(boundary, False, "array")
        # newtest
        # array = self.pressures[:, cells_id]
        array = self.__get_pressure(None, cells_id)
        if scale:
            array = self.pressure_scaler.transform(array)
        if melt:
            labels = ["P" + press_str]
            array = array.flatten()
        else:
            labels = [f"P{str(id)}" + press_str for id in cells_id]
        data = pd.DataFrame(array, columns=labels)
        return self.__concat(data, df)

    def __add_wells_rate(self, units, scale, df=None):
        if units:
            if scale and self.scalers_dict["rate"][0] is not None:
                rate_str = " [scaled]"
            else:
                rate_str = f" [{self.units['rate']}]"
        else:
            rate_str = ""
        labels = [f"Qw{str(id)}" + rate_str for id in self.wells.keys()]
        # newtest
        # array = self.rates[:, list(self.wells.keys())]
        array = self.__get_rate(None, list(self.wells.keys()))
        if scale:
            array = self.rate_scaler.transform(array)
        data = pd.DataFrame(array, columns=labels)
        return self.__concat(data, df)

    def __add_wells_pressures(self, units, scale, df=None):
        if units:
            if scale and self.scalers_dict["pressure"][0] is not None:
                press_str = " [scaled]"
            else:
                press_str = f" [{self.units['pressure']}]"
        else:
            press_str = ""
        labels = [f"Pwf{str(id)}" + press_str for id in self.w_pressures]
        data = pd.DataFrame(self.w_pressures)
        if scale:
            data = self.pressure_scaler.transform(data.values)
            data = pd.DataFrame(data)
        data.columns = labels
        return self.__concat(data, df)

    def __add_xyz(self, boundary, melt, scale, df=None):
        if melt:
            cells_center, fdir = self.get_centers(scale, boundary)
            array = np.tile(cells_center.flatten(), self.solution.nsteps).reshape(
                -1, len(fdir)
            )
            data = pd.DataFrame(array, columns=fdir)
            return self.__concat(data, df)
        return df

    def __get_time_vector(self):
        return np.arange(0, (self.__get_tstep() + 1) * self.dt, self.dt)

    def __update_time_scaler(self):
        time = self.__get_time_vector()
        self.time_scaler.fit(time, axis=0)

    def __update_space_scaler(self, boundary):
        cells_center, _ = self.get_centers(False, boundary)
        self.space_scaler.fit(cells_center, axis=0)

    def __update_pressure_scaler(self, boundary):
        config = self.__get_scalers_config(boundary)
        pressures = self.get_df(columns=["cells_pressure"], **config).values
        if self.solution and self.solution.__class__.__name__ in ["D1P1"]:
            if self.verbose:
                print(
                    "[warning] To avoid scaling based on unstable solutions, "
                    "pressures scaler is updated based on the final timestep "
                    "for D1P1 solution."
                )
            self.pressure_scaler.fit(pressures[-1], axis=None)
        else:
            self.pressure_scaler.fit(pressures, axis=None)

    def __update_rate_scaler(self, boundary):
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
        self.rate_scaler.fit(rates, axis=None)

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
        self.__update_pressure_scaler(boundary)
        self.__update_rate_scaler(boundary)

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
        """Set scalers configuration.

        To change the scaling settings. Current settings can be shown
        under ``scalers_dict``. By default the following settings are
        used:

        .. code-block:: python

            scalers_dict = {
                'time':['MinMax', (0,1)],
                'space':['MinMax', (-1,1)],
                'pressure':['MinMax', (-1,1)],
                'rate':[None,None],
            }

        Note that by default rates are not scaled, time is scaled
        between 0 and 1, while space and pressure are scaled between
        -1 and 1. By default, MinMax scaler is used for all dimensions.

        Parameters
        ----------
        scalers_dict : dict
            scalers setting as dict in the following
            format: {'time': [scaler_type, scaler_range]} were:
            * scaler_type is a string with the scaler name from the
            ``scalers`` module (e.g. 'MinMax').
            * scaler_range is a tuple of the output_range of the scaler
            e.g. (-1,1).
            The keys must be in ['time', 'space', 'pressure', 'rate'].
        """

        def create_scaler(scaler_type, output_range):
            if scaler_type is None or output_range is None:
                return scalers.Dummy(None), None
            elif scaler_type.lower() in ["minmax", "minmaxscaler"]:
                return scalers.MinMax(output_range=output_range), "MinMax"
            elif scaler_type.lower() in ["logminmax", "logminmaxscaler"]:
                return scalers.LogMinMax(output_range=output_range), "LogMinMax"
            else:
                raise ValueError("scaler type is not defined.")

        col_dict = {
            "time": ["t", "time", "all"],
            "space": ["space", "spatial", "xyz", "x", "y", "z", "all"],
            "pressure": ["p", "pressures", "pressure", "all"],
            "rate": ["q", "rates", "rate", "all"],
        }
        # col_vals = sum(col_dict.values(), [])

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
                self.pressure_scaler, s_str = create_scaler(scaler_type, output_range)
                column_str = "pressure"
            elif column.lower() in col_dict["rate"]:
                self.rate_scaler, s_str = create_scaler(scaler_type, output_range)
                column_str = "rate"
            else:  # if column.lower() not in col_vals:
                raise ValueError(f"column {column} does not exist.")

            if scaler_type is None or output_range is None:
                self.scalers_dict[column_str] = [None, None]
            else:
                self.scalers_dict[column_str] = [s_str, scalers_dict[column][1]]

    def get_scalers(
        self,
        name: str = None,
    ):
        """Get scalers configuration.

        Parameters
        ----------
        name : str, optional
            name of the scaler to be returned. If None or "all" or "*",
            all scalers will be returned. Otherwise, the scaler with the
            specified name will be returned.

        Returns
        -------
        dict
            scalers configuration as a dictionary.
        """
        # self.time_scaler
        # self.space_scaler
        # self.rate_scaler
        # self.pressure_scaler
        if name in self.scalers_dict:
            return self.scalers_dict[name]
        elif name is None or name in ["all", "*"]:
            return self.scalers_dict
        else:
            raise ValueError(
                f"name must be defined. Use: {list(self.scalers_dict.keys)}."
            )

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

                - ``"time"``: for time steps as specified in dt.
                - ``"date"``: for dates as specified in dt and start_date.
                - ``"q"``, ``"rates"``: for all (cells and wells) rates.
                - ``"p"``, ``"pressures"``: for all (cells and wells) pressures.
                - ``"cells"``: for all cells rates and pressures.
                - ``"wells"``: for all wells rates and pressures.
                - ``"cells_rate"``: for all cells rates (including wells' cells).
                - ``"cells_pressure"``: for all cells pressures.
                - ``"wells_rate"``: for all wells rates.
                - ``"wells_pressure"``: for all wells pressures.

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

            .. code-block:: python

                scalers_dict = {
                    'time':['MinMax', (0,1)],
                    'space':['MinMax', (-1,1)],
                    'pressure':['MinMax', (-1,1)],
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
            df["id"] = np.tile(cells_id, self.solution.nsteps)
            df["Step"] = np.repeat(np.arange(self.solution.nsteps), n_cells)
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
        df.index = df.index.astype("int32", copy=False)

        if save:
            df.to_csv("model_data.csv")
            print("[info] Model data was successfully saved.")

        return df

    # -------------------------------------------------------------------------
    # Visualization:
    # -------------------------------------------------------------------------

    show = utils.pyvista.show_model
    save_gif = utils.pyvista.save_gif

    # -------------------------------------------------------------------------
    # Plotting:
    # -------------------------------------------------------------------------

    def get_values(
        self,
        scale=False,
        boundary=True,
    ):

        df = self.get_df(
            columns=["time", "cells_pressure"],
            boundary=boundary,
            units=False,
            melt=True,
            scale=scale,
            drop_zero=False,
            drop_nan=False,
        )
        shape = self.get_shape(boundary)
        fdir = list(self.grid.get_fdir())
        ncols = len(fdir) + 2  # +2 for time and P
        values = df[["Time", *fdir, "P"]].values.reshape(*shape, ncols)
        X, Y = values[:, :, :2], values[:, :, 2]
        return X, Y

    def plot(
        self,
        type: str = "line",
        scale: bool = False,
        boundary: bool = True,
        nrows: int = 3,
        ncols: int = 3,
        solution: str = None,
    ):
        """Plot solution values.

        Parameters
        ----------
        type : str, optional
            type of the plot in ['line', 'contour']. Currently only 1D
            plots are supported. Default is "line".
        scale : bool, optional
            scale values using scalers. Default is False.
        boundary : bool, optional
            include boundary cells in the plot. Default is True.
        nrows : int, optional
            number of rows in the plot. Default is 3.
        ncols : int, optional
            number of columns in the plot. Default is 3.
        solution : str, optional
            name of the solution to plot. If "all" or "*", all solutions
            will be plotted. If None, the current solution will be
            plotted. Default is None.

        Returns
        -------
        plotter : Plot1D, Plot2D, or Plot3D
            plotter object with the plotted solution values.
        """
        if self.grid.D == 1:
            if type == "line":
                plotter = plots.Plot1D(
                    nrows=nrows,
                    ncols=ncols,
                    verbose=self.verbose,
                )
            elif type == "contour":
                plotter = plots.Contour1D(
                    verbose=self.verbose,
                )
            else:
                raise ValueError(f"Plot type {type} is not supported for 1D grid.")
        elif self.grid.D == 2:
            raise NotImplementedError("2D plotting is not implemented yet.")
        elif self.grid.D == 3:
            raise NotImplementedError("3D plotting is not implemented yet.")
        else:
            raise ValueError(f"Grid dimension {self.grid.D} is not supported.")

        if solution in ["all", "*"]:
            for solution_name in self.solutions:
                self.set_solution(solution_name)
                X, Y = self.get_values(scale, boundary)
                plotter.add(
                    x=X,
                    y=Y,
                    name=solution_name,
                )
        else:
            if solution is not None:
                self.set_solution(solution)
            X, Y = self.get_values(scale, boundary)
            plotter.add(
                x=X,
                y=Y,
                name=self.solution.name,
            )

        return plotter.plot()

    plot_solution = plot

    def plot_solutions(
        self,
        scale=True,
        boundary=True,
        nrows: int = 3,
        ncols: int = 3,
    ):
        """Plot all solutions in a grid.

        Parameters
        ----------
        scale : bool, optional
            scale values using scalers. Default is True.
        boundary : bool, optional
            include boundary cells in the plot. Default is True.
        nrows : int, optional
            number of rows in the plot. Default is 3.
        ncols : int, optional
            number of columns in the plot. Default is 3.

        Returns
        -------
        plotter : Plot1D, Plot2D, or Plot3D
            plotter object with the plotted solution values.
        """
        return self.plot(
            scale=scale,
            boundary=boundary,
            nrows=nrows,
            ncols=ncols,
            solution="all",
        )

    def contour(self):

        return None

    # def plot(self, prop: str = "pressures", id: int = None, tstep: int = None):
    #     """Show values in a cartesian plot.

    #     Parameters
    #     ----------
    #     prop : str, optional
    #         property name from ["rates", "pressures"].
    #     id : int, optional
    #         cell id. If None, all cells are selected.
    #     tstep : int, optional
    #         time step. If None, the last time step is selected.
    #     """
    #     if tstep is None:
    #         tstep = self.__get_tstep()

    #     if id is not None:
    #         exec(f"plt.plot(self.{prop}[:, id].flatten())")
    #         plt.xlabel("Days")
    #     elif tstep is not None:
    #         exec(f"plt.plot(self.{prop}[tstep, :].flatten())")
    #         plt.xlabel("Grid (id)")
    #         boundary = True
    #         x_skip = 5
    #         x_ticks = np.arange(0, self.grid.get_n(boundary), x_skip)
    #         x_labels = self.grid.get_cells_id(boundary, fshape=False, fmt="array")
    #         x_labels = x_labels[::x_skip]
    #         plt.xticks(ticks=x_ticks, labels=x_labels)
    #     plt.grid()
    #     plt.show()

    # def plot_grid(self, property: str = "pressures", tstep: int = None):
    #     if tstep is None:
    #         tstep = self.__get_tstep()
    #     cells_id = self.grid.get_cells_id(False, False, "list")
    #     exec(f"plt.imshow(self.{property}[tstep][cells_id][np.newaxis, :])")
    #     plt.colorbar(label=f"{property.capitalize()} ({self.units[property[:-1]]})")
    #     plt.title(f"{property.capitalize()} Distribution")
    #     plt.yticks([])
    #     plt.xlabel("Grid (i)")
    #     plt.xticks(ticks=range(0, 4), labels=range(1, 5))
    #     plt.show()

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
        """Allow full descriptions.

        This function maps functions as following:

        .. code-block:: python

            self.set_transmissibility = self.set_trans
            self.transmissibility = self.T

        """
        self.set_transmissibility = self.set_trans
        self.transmissibility = self.T

    # -------------------------------------------------------------------------
    # End
    # -------------------------------------------------------------------------


if __name__ == "__main__":

    import reservoirflow as rf

    def create_model_example_7_1():
        grid = rf.grids.RegularCartesian(
            nx=4, ny=1, nz=1, dx=300, dy=350, dz=40, phi=0.27, kx=270, dtype="double"
        )
        fluid = rf.fluids.SinglePhase(mu=0.5, B=1, dtype="double")
        model = BlackOil(grid, fluid, dtype="double", pi=4000, verbose=False)
        model.set_well(cell_id=4, q=-600, s=1.5, r=3.5)
        model.set_boundaries({0: ("pressure", 4000), 5: ("rate", 0)})
        return model

    def create_model():
        grid = rf.grids.RegularCartesian(
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
        fluid = rf.fluids.SinglePhase(
            mu=0.5,
            B=1,
            rho=50,
            comp=1 * 10**-5,
            dtype="double",
        )
        model = BlackOil(
            grid,
            fluid,
            pi=6000,
            dt=5,
            start_date="10.10.2018",
            dtype="double",
        )

        cells_id = model.grid.get_cells_id(False, True)[-1].flatten()
        wells = np.random.choice(cells_id, 6, False)
        for id in wells:
            model.set_well(cell_id=id, q=-300, pwf=100, s=1.5, r=3.5)

        wells = np.random.choice(cells_id, 6, False)
        for id in wells:
            model.set_well(cell_id=id, q=100, s=0, r=3.5)

        return model

    sparse = True

    model = create_model()
    model.compile(stype="numerical", method="FDM", sparse=sparse)
    model.run(nsteps=10, vectorize=True, isolver="cgs")
    model.show("pressures")

    model = create_model()
    model.compile(stype="numerical", method="FDM", sparse=sparse)
    model.run(nsteps=10, vectorize=True, isolver=None)
    model.show("pressures")
