"""
FDM
---

Finite Difference Method (FDM) class.
"""
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

# from ..solution import Solution
import reservoirflow as rf
from reservoirflow.solutions.solution import Solution
from reservoirflow.utils.helpers import _lru_cache


class FDM(Solution):
    """FDM solution class.

    FDM is a Finite-Difference-Method.

    Returns
    -------
    Solution
        Solution object.
    """

    name = "FDM"

    # -------------------------------------------------------------------------
    # Flow Equations: symbolic
    # -------------------------------------------------------------------------

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
        dz = self.model.grid.z[cell_n_id] - self.model.grid.z[cell_id]
        return trans * ((cell_n_p - cell_p) - (self.model.fluid.g * dz))

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
        # cell_b_p = self.model.pressures[self.model.tstep, cell_b_id]
        # if not np.isnan(cell_b_p):
        #     dz = self.model.grid.z[cell_b_id] - self.model.grid.z[cell_id]
        #     b_term = trans * 2 * ((cell_b_p - cell_p) - (self.model.fluid.g * dz))
        # else:
        #     b_term = self.model.rates[self.model.tstep, cell_b_id]

        # implementation 2:
        if cell_b_id in self.model.bdict:
            cond, v = self.model.bdict[cell_b_id]
            if cond.lower() in ["pressure", "press", "p"]:
                dz = self.model.grid.z[cell_b_id] - self.model.grid.z[cell_id]
                return trans * 2 * ((v - cell_p) - (self.model.fluid.g * dz))
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
        if (
            "q" in self.model.wells[cell_id]
            and self.model.wells[cell_id]["constrain"] == "q"
        ):
            return self.model.wells[cell_id]["q"]
        else:
            return (
                -self.model.wells[cell_id]["G"]
                / (self.model.fluid.B * self.model.fluid.mu)
                * (cell_p - self.model.wells[cell_id]["pwf"])
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
        if self.model.comp_type == "incompressible":
            return 0.0
        else:
            try:
                return self.model.RHS[cell_id] * (
                    cell_p - self.model.pressures[self.model.tstep, cell_id]
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
        #     # if not np.isnan(self.model.pressures[self.model.tstep][i]):
        #     #     exec(f"p{i} = {self.model.pressures[self.model.tstep][i]}")
        # - n_term to use pressure values:
        #     # To Do: keep pressure constant at specific cell (requires A adjust)
        #     # if not np.isnan(self.model.pressures[self.model.tstep][neighbor]):
        #     #     exec(f"p{neighbor} = {self.model.pressures[self.model.tstep][neighbor]}")
        # - n_term in one calc.
        # exec(
        #     f"n_term = self.T[dir][min(neighbor,id)] * ((p{n_id} - p{id})
        #     - (self.model.fluid.g * (self.model.grid.z[neighbor] - self.model.grid.z[id])))"
        # )
        cell_p = eval(f"sym.Symbol('p{cell_id}')")

        if cell_id not in self.model.cells_terms:
            assert (
                cell_id in self.model.grid.cells_id
            ), f"id is out of range {self.model.grid.cells_id}."
            neighbors = self.model.grid.get_cell_neighbors(
                id=cell_id, boundary=False, fmt="array"
            )
            boundaries = self.model.grid.get_cell_boundaries(id=cell_id, fmt="array")
            # f_terms: flow terms, a_terms: accumulation terms
            terms = {"f_terms": [], "a_term": 0}
            T = self.model.get_cell_trans(cell_id, None, True)
            if self.model.verbose:
                print(f"[info] cell id: {cell_id}")
                print(f"[info]    - Neighbors: {neighbors}")
                print(f"[info]    - Boundaries: {boundaries}")

            for id_n in neighbors:
                n_term = self.__calc_n_term(cell_id, id_n, cell_p, T[id_n])
                terms["f_terms"].append(n_term)
                if self.model.verbose:
                    print(f"[info] Neighbor terms: {n_term}")

            for cell_b_id in boundaries:
                b_term = self.__calc_b_term(cell_id, cell_b_id, cell_p, T[cell_b_id])
                terms["f_terms"].append(b_term)
                if self.model.verbose:
                    print(f"[info] Boundary terms: {b_term}")

            if cell_id in self.model.wells.keys():
                w_term = self.__calc_w_term(cell_id, cell_p)
                terms["f_terms"].append(w_term)
                if self.model.verbose:
                    print(f"[info] Well terms: {w_term}")

            terms["a_term"] = self.__calc_a_term(cell_id, cell_p)
            if self.model.verbose:
                print("[info] Accumulation term:", terms["a_term"])

            self.model.cells_terms[cell_id] = terms
            if self.model.verbose:
                print("[info] terms:", terms)

        else:
            terms = self.model.cells_terms[cell_id]
            if (
                cell_id in self.model.wells.keys()
                and self.model.wells[cell_id]["constrain"] == "pwf"
            ):
                w_term = self.__calc_w_term(cell_id, cell_p)
                if self.model.verbose:
                    print(f"[info] Well terms (updated): {w_term}")
                terms["f_terms"][-1] = w_term
            if self.model.comp_type == "compressible":
                terms["a_term"] = self.__calc_a_term(cell_id, cell_p)

        cell_eq = sym.Eq(sum(terms["f_terms"]), terms["a_term"])
        cell_eq = self.__simplify_eq(cell_eq)
        cell_eq_lhs = cell_eq.lhs.as_coefficients_dict()

        if self.model.verbose:
            print(f"[info] Flow equation {cell_id}:", cell_eq)

        return cell_eq_lhs, cell_eq.rhs

    def get_cells_eq(self, threading=False):
        """Return flow equations for all internal cells."""
        cells_eq = {}
        n_threads = self.model.n // 2
        if threading:
            with ThreadPoolExecutor(n_threads) as executor:
                # with ProcessPoolExecutor(2) as executor:
                equations = executor.map(self.get_cell_eq, self.model.grid.cells_id)
                for id, eq in zip(self.model.grid.cells_id, equations):
                    cells_eq[id] = eq
        else:
            for id in self.model.grid.cells_id:
                cells_eq[id] = self.get_cell_eq(id)
                if self.model.verbose:
                    print(f"[info] cell id: {id}")
                    print(f"[info]      - lhs: {cells_eq[id][0]}")
                    print(f"[info]      - rhs: {cells_eq[id][1]}")

        return cells_eq

    # -------------------------------------------------------------------------
    # Matrices: symbolic
    # -------------------------------------------------------------------------

    def __update_matrices_symb(self, cell_id):
        """Update flow equations' matrices (A, d).

        Parameters
        ----------
        id : int
            cell id based on natural order as int.

        Notes
        -----
        - arrays for lhs and rhs:
            self.d[i] = np.array(cell_rhs).astype(self.model.dtype)
            self.A[i, ids] = np.array(list(cell_lhs.values())).astype(self.model.dtype)
        - finding cell i:
            ids = [self.model.cells_id.index(int(str(s)[1:])) for s in cell_lhs.keys()]
        """
        cell_lhs, cell_rhs = self.cells_eq[cell_id]
        ids = [self.model.cells_i_dict[int(str(s)[1:])] for s in cell_lhs.keys()]
        self.d[self.model.cells_i_dict[cell_id]] = cell_rhs
        self.A[self.model.cells_i_dict[cell_id], ids] = list(cell_lhs.values())
        if self.model.verbose:
            print(f"[info] cell id: {cell_id}")
            print(f"[info]      - ids: {ids}")
            print(f"[info]      - lhs: {cell_lhs}")
            print(f"[info]      - rhs: {cell_rhs}")

    def get_matrices_symb(self, threading=False):
        """Initialize flow equations' matrices (A, d).

        Parameters
        ----------
        sparse : bool, optional
            use sparse matrices instead of dense matrices.
        threading : bool, optional
            use multiple threads for concurrence workers. The maximum
            number of threads are set to the half number of cells.

        Returns
        -------
        _type_
            _description_

        """
        # ToDo
        # ----
        # - Update only required cells.
        self.cells_eq = self.get_cells_eq(threading)

        if self.model.tstep == 0 or not hasattr(self, "A") or not hasattr(self, "d"):
            # second and third conditions allow switching vectorize
            # True/False after timestep=0.
            if self.model.sparse:
                self.d = ss.lil_matrix((self.model.n, 1), dtype=self.model.dtype)
                self.A = ss.lil_matrix(
                    (self.model.n, self.model.n), dtype=self.model.dtype
                )
            else:
                self.d = np.zeros((self.model.n, 1), dtype=self.model.dtype)
                self.A = np.zeros((self.model.n, self.model.n), dtype=self.model.dtype)

        if threading:
            with ThreadPoolExecutor(self.model.n) as executor:
                # with ProcessPoolExecutor(2) as executor:
                executor.map(self.__update_matrices_symb, self.model.cells_id)
        else:
            for cell_id in self.model.cells_id:
                self.__update_matrices_symb(cell_id)

        if self.model.verbose:
            print("[info] - A:\n", self.A)
            print("[info] - d:\n", self.d)

        return self.A, self.d

    # -------------------------------------------------------------------------
    # Matrices: vectorized
    # -------------------------------------------------------------------------

    def __init_A(self):
        """Initialize ceofficient matrix (`A`).

        For a system of linear equations `Au=d`, `A` is the
        ceofficient matrix (known), `d` is the constant vector (known),
        and `u` is the variable vector (unknown e.g., pressure).

        This function initialize the ceofficient matrix (`A`) which is
        needed only at initial timestep (i.e., `timestep=0`).

        Returns
        -------
        ndarray
            ceofficient A is initialized in place.
        """
        # T = self.get_cells_T_array(True).toarray()
        # self.A_ = T[:, self.model.cells_id][self.model.cells_id]
        # self.A_[self.model.cells_i, self.model.cells_i] = (
        #     -self.A_[self.model.cells_i, :].sum(axis=1) - self.model.RHS[self.model.cells_id]
        # )
        # if sparse:
        #     self.A_ = ss.lil_matrix(self.A_, dtype=self.model.dtype)

        # return self.A_
        # self.A_ = self.get_cells_T_array(False, True).toarray()
        self.A_ = self.model.get_cells_trans(False, self.model.sparse, True)
        v1 = -self.A_[self.model.cells_i, :].sum(axis=1).flatten()
        v2 = self.model.RHS[self.model.cells_id].flatten()
        v3 = v1 - v2
        self.A_[self.model.cells_i, self.model.cells_i] = v3
        if self.model.sparse:
            self.A_ = ss.lil_matrix(self.A_, dtype=self.model.dtype)
        return self.A_

    def __init_d(self):
        """Initialize constant vector (`d`).

        For a system of linear equations `Au=d`, `A` is the
        ceofficient matrix (known), `d` is the constant vector (known),
        and `u` is the variable vector (unknown e.g., pressure).

        This function initialize the constant vector (`d`) which is
        needed at every timestep in case of a compressible system. In
        case of an incompressible system, a constant zero vector is
        used.

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
        if self.model.sparse:
            self.d_ = ss.lil_matrix((self.model.n, 1), dtype=self.model.dtype)
        else:
            self.d_ = np.zeros((self.model.n, 1), dtype=self.model.dtype)

        if self.model.comp_type == "compressible":
            try:
                self.d_[:] = (
                    -self.model.RHS[self.model.grid.cells_id]
                    * self.model.pressures[self.model.tstep, self.model.grid.cells_id]
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
        z = self.model.grid.z[self.model.grid.cells_id]
        if not np.all(z == z[0]):
            z_l = np.append(z[1:], np.nan)
            z_u = np.append(np.nan, z[:-1])
            dz_l = self.model.fluid.g * np.nan_to_num(z_l - z)
            dz_u = self.model.fluid.g * np.nan_to_num(z_u - z)
            # T = self.T["x"][self.model.grid.cells_id]
            # T = np.diag(self.get_cells_T(True, False), 1)[self.model.cells_i]
            T = self.model.get_cells_trans_diag(True, 1)[self.model.cells_id]
            v = T * dz_l + T * dz_u
            self.d_ += v.reshape(-1, 1)

    def get_matrices_vect(self, threading=False):
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
        if self.model.tstep == 0 or not hasattr(self, "A_"):
            # second condition allow switching vectorize True/False after
            # timestep=0.
            self.resolve = defaultdict(lambda: False)
            self.__init_A()
            self.bdict_v = {}
            for id_b in self.model.bdict_update:
                ((id, T),) = self.model.get_cell_trans(id_b, None, False).items()
                p = eval(f"sym.Symbol('p{id}')")
                b_term = self.__calc_b_term(id, id_b, p, T)
                v0, v1 = b_term.as_coefficients_dict().values()
                self.bdict_v[id_b] = (v0, v1, id)
                self.A_[self.model.cells_i_dict[id], self.model.cells_i_dict[id]] += v1
                # self.A_[self.model.cells_i_dict[id], self.model.cells_i_dict[id]] -= T * 2

        self.__init_d()

        for id in self.model.wells.keys():
            if self.model.wells[id]["constrain"] == "q":
                w_term = self.__calc_w_term(
                    id, self.model.pressures[self.model.tstep, id]
                )
                self.d_[self.model.cells_i_dict[id], 0] -= w_term
                update_z = True
            elif self.model.wells[id]["constrain"] == "pwf":
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
                self.d_[self.model.cells_i_dict[id], 0] -= v0
                if not self.resolve[id]:
                    self.A_[
                        self.model.cells_i_dict[id], self.model.cells_i_dict[id]
                    ] += v1
                    self.resolve[id] = True
            else:
                pass  # no constrain

        for id_b in self.model.bdict.keys():
            id = self.model.grid.get_cell_neighbors(id_b, None, False, "list")
            if len(id) > 0:
                id = id[0]
                if self.model.bdict[id_b][0] == "pressure":
                    self.d_[self.model.cells_i_dict[id], 0] -= self.bdict_v[id_b][0]
                else:  # elif self.model.bdict[id_b][0] in ["gradient", "rate"]:
                    self.d_[self.model.cells_i_dict[id], 0] -= self.model.rates[
                        self.model.tstep, id_b
                    ]

        if update_z:
            self.__update_z()

        return self.A_, self.d_

    # -------------------------------------------------------------------------
    # Numerical Solution:
    # -------------------------------------------------------------------------

    def __update_wells(self):
        """_summary_


        Notes
        -----
        - well q calc:
            self.__calc_w_terms(
                    id, self.model.pressures[self.model.tstep][id]
                )
            or
            self.model.wells[id]["q"] = (
                -self.model.wells[id]["G"]
                / (self.model.fluid.B * self.model.fluid.mu)
                * (self.model.pressures[self.model.tstep][id] - self.model.wells[id]["pwf"])
            )
        - all calc original:
            if "q" in self.model.wells[id]:
                self.model.wells[id]["pwf"] = self.model.pressures[self.model.tstep][id] + (
                    self.model.wells[id]["q"]
                    * self.model.fluid.B
                    * self.model.fluid.mu
                    / self.model.wells[id]["G"]
                )
            self.model.w_pressures[id].append(self.model.wells[id]["pwf"])
            if "pwf" in self.model.wells[id]:
                self.model.wells[id]["q"] = (
                    -self.model.wells[id]["G"]
                    / (self.model.fluid.B * self.model.fluid.mu)
                    * (self.model.pressures[self.model.tstep][id] - self.model.wells[id]["pwf"])
                )
                self.model.rates[self.model.tstep][id] = self.model.wells[id]["q"]
        """
        resolve = False
        tstep_w_pressures = {}
        for id in self.model.wells.keys():
            if "q_sp" in self.model.wells[id]:
                pwf_est = self.model.pressures[self.model.tstep, id] + (
                    self.model.wells[id]["q_sp"]
                    * self.model.fluid.B
                    * self.model.fluid.mu
                    / self.model.wells[id]["G"]
                )
            else:
                pwf_est = self.model.wells[id]["pwf"]

            if pwf_est > self.model.wells[id]["pwf_sp"]:
                self.model.wells[id]["constrain"] = "q"
            else:
                if (
                    pwf_est < self.model.wells[id]["pwf_sp"]
                    and self.model.wells[id]["q"] == self.model.wells[id]["q_sp"]
                ):
                    resolve = True

                self.model.wells[id]["constrain"] = "pwf"
                pwf_est = self.model.wells[id]["pwf_sp"]

            self.model.wells[id]["pwf"] = pwf_est
            q_est = self.__calc_w_term(id, self.model.pressures[self.model.tstep, id])
            self.model.wells[id]["q"] = self.model.rates[self.model.tstep, id] = q_est

            if resolve:
                return True
            else:
                tstep_w_pressures[id] = pwf_est

        for id in self.model.wells.keys():
            self.model.w_pressures[id].append(tstep_w_pressures[id])

        return False

    def __update_boundaries(self):
        for id_b in self.model.bdict_update:
            ((id_n, T),) = self.model.get_cell_trans(id_b, None, False).items()
            p_n = self.model.pressures[self.model.tstep, id_n]
            b_terms = self.__calc_b_term(id_n, id_b, p_n, T)
            self.model.rates[self.model.tstep, id_b] = b_terms

    def __print_arrays(self, sparse):
        if sparse:
            A, d = self.A.toarray(), self.d.toarray()
            A_, d_ = self.A_.toarray(), self.d_.toarray()
        else:
            A, d = self.A, self.d
            A_, d_ = self.A_, self.d_
        print("step:", self.model.tstep)
        print(np.concatenate([A, A_, abs(A) - abs(A_)], axis=0))
        print(np.concatenate([d, d_, abs(d) - abs(d_)], axis=1))
        print()

    def solve(
        self,
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
        # sparse = self.model.sparse
        if print_arrays:
            A, d = self.get_matrices_symb(threading)  #  has to be first
            self.get_matrices_vect(threading)
            print(f"self.A : {type(self.A)}")
            print(f"self.d : {type(self.d)}")
            print(f"self.A_: {type(self.A_)}")
            print(f"self.d_: {type(self.d_)}")
            self.__print_arrays(self.model.sparse)
        else:
            if vectorize:
                A, d = self.get_matrices_vect(threading)
            else:
                A, d = self.get_matrices_symb(threading)

        if self.model.sparse:
            A, d = A.tocsc(), d.toarray()
            if isolver:
                solver = rf.solutions.numerical.solvers.get_isolver(isolver)
                pressures, exit_code = solver(
                    A,
                    d,
                    atol=0,
                    # x0=self.model.pressures[self.model.tstep, self.model.cells_id],
                )
                assert exit_code == 0, "unsuccessful convergence"
            else:
                pressures = ssl.spsolve(A, d, use_umfpack=True)
            A = A.toarray()
        else:
            pressures = sl.solve(A, d).flatten()

        if update:
            self.model.tstep += 1
            self.model.pressures = np.vstack(
                [self.model.pressures, self.model.pressures[-1]]
            )
            self.model.pressures[self.model.tstep, self.model.grid.cells_id] = pressures
            self.model.rates = np.vstack([self.model.rates, self.model.rates[-1]])
            self.model.As = np.vstack([self.model.As, A.reshape(1, -1)])
            self.model.ds = np.vstack([self.model.ds, d.reshape(1, -1)])
            self.__update_boundaries()
            resolve = self.__update_wells()
            if resolve:
                self.model.rates = self.model.rates[: self.model.tstep]
                self.model.pressures = self.model.pressures[: self.model.tstep]
                self.model.tstep -= 1
                self.solve(threading, vectorize, False, True, False)
                if self.model.verbose:
                    print(f"[info] Time step {self.model.tstep} was resolved.")

            if check_MB:
                self.check_MB()

        if self.model.verbose:
            print("[info] Pressures:\n", self.model.pressures[self.model.tstep])
            print("[info] rates:\n", self.model.rates[self.model.tstep])

    def run(
        self,
        nsteps=10,
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
        - SciPy: `Solving Linear Problems <https://docs.scipy.org/doc/scipy/reference/sparse.linalg.html#solving-linear-problems>`_.
        - SciPy: `Iterative Solvers <https://scipy-lectures.org/advanced/scipy_sparse/solvers.html#iterative-solvers>`_.
        """
        start_time = time.time()
        self.model.nsteps += nsteps
        self.run_ctime = 0
        if self.model.verbose:
            self.model.verbose = False
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
                # sparse,
                threading,
                vectorize,
                check_MB,
                True,
                print_arrays,
                isolver,
            )
            pbar.set_description(f"[step] {step}")

        self.run_ctime = round(time.time() - start_time, 2)
        self.model.ctime += self.run_ctime
        print(
            f"[info] Simulation run of {nsteps} steps",
            f"finished in {self.run_ctime} seconds.",
        )
        if check_MB:
            print(f"[info] Material Balance Error: {self.model.error}.")

        if verbose_restore:
            self.model.verbose = True

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
            print(f"[info] Error in step {self.model.tstep}")

        if self.model.comp_type == "incompressible":
            # rates must add up to 0:
            self.model.error = self.model.rates[self.model.tstep].sum()
            if verbose:
                print(f"[info]    - Error: {self.model.error}")
        elif self.model.comp_type == "compressible":
            # error over a timestep:
            self.model.error = (
                self.model.RHS[self.model.grid.cells_id]
                * (
                    self.model.pressures[self.model.tstep, self.model.grid.cells_id]
                    - self.model.pressures[
                        self.model.tstep - 1, self.model.grid.cells_id
                    ]
                )
            ).sum() / self.model.rates[self.model.tstep].sum()
            # error from initial timestep to current timestep: (less accurate)
            self.model.cumulative_error = (
                self.model.RHS[self.model.grid.cells_id]
                * self.model.dt
                * (
                    self.model.pressures[self.model.tstep, self.model.grid.cells_id]
                    - self.model.pressures[0, self.model.grid.cells_id]
                )
            ).sum() / (self.model.dt * self.model.tstep * self.model.rates.sum())
            self.model.error = abs(self.model.error - 1)
            if self.model.verbose:
                print(f"[info]    - Incremental Error: {self.model.error}")
                print(f"[info]    -  Cumulative Error: {self.model.cumulative_error}")
                print(
                    f"[info]    -       Total Error: {self.model.error+self.model.cumulative_error}"
                )

        if abs(self.model.error) > error_threshold:
            warnings.warn("High material balance error.")
            print(
                f"[warning] Material balance error ({self.model.error})",
                f"in step {self.model.tstep}",
                f"is higher than the allowed error ({error_threshold}).",
            )
