"""
Model classes to create reservoir simulation models.

This module contains all model classes used to create a reservoir 
simulation model in combination with with a Fluid class and Grid class.
"""
import time
from openresim.base import Base
from openresim import grids, fluids, wells, plots
import numpy as np
import sympy as sym
import scipy.sparse as ss
import scipy.sparse.linalg as ssl
import matplotlib.pyplot as plt
from tqdm import tqdm


class Model(Base):
    """Model class to create a reservoir simulation model.

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

        self.flow_equations_terms_dict = {}
        self.dt = dt
        self.A = None
        self.nsteps = 1
        self.tstep = 0

        self.__initialize__(pi, well)
        self.__calc_comp()
        self.__calc_trans()
        self.__calc_RHS()

    # -------------------------------------------------------------------------
    # Basic:
    # -------------------------------------------------------------------------

    def __initialize__(self, pi, well):
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

    # -------------------------------------------------------------------------
    # Properties:
    # -------------------------------------------------------------------------

    def __calc_comp(self):
        """ """
        if self.fluid.comp_type == self.grid.comp_type == "incompressible":
            self.set_comp(0)
        else:
            self.set_comp(self.fluid.comp + self.grid.comp)

    def __calc_trans(self):
        self.trans = {}
        for dir in self.grid.get_fdir():
            self.trans[dir] = self.get_trans(dir, False)

    def get_trans(self, dir="x", fshape=False):
        """Returns transmissibility.

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
            array of transmissibility based on dir argument.
        """

        G = self.grid.get_G(dir, fshape)
        trans = G / (self.fluid.mu * self.fluid.B)

        return trans

    # -------------------------------------------------------------------------
    # Wells:
    # -------------------------------------------------------------------------

    def __calc_well_G(self, id=None):
        G_n = (
            2
            * np.pi
            * self.factors["transmissibility conversion"]
            * self.grid.k["x"][id]
            * self.grid.d["z"][id]
        )

        G_d = np.log(self.wells[id]["r_eq"] / self.wells[id]["r"] * 12)

        if "s" in self.wells[id].keys():
            G_d += self.wells[id]["s"]

        return G_n / G_d

    def __calc_well_r_eq(self, id):
        return 0.14 * (self.grid.dx[id] ** 2 + self.grid.dy[id] ** 2) ** 0.5

    def set_well(self, id=None, well=None, q=None, pwf=None, r=None, s=None):
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

    # -------------------------------------------------------------------------
    # Boundaries:
    # -------------------------------------------------------------------------

    def set_boundary(self, i: int, cond: str, v: float):
        """Set model boundary condition

        Args:
            i (int): boundary grid index
                - can be obtained from Grid class using get_boundaries() method.
            cond (str): boundary condition
                - contant pressure: 'pressure', 'press', 'p'
                - constant rate: 'rate', 'q'
                - constant pressure gradient: 'pressure gradient', 'press grad', 'gradient', 'grad', 'pg', 'g'
            v (float): value
        """
        if cond in ["rate", "q"]:
            self.rates[0][i] = v
        if cond in ["pressure gradient", "press grad", "gradient", "grad", "pg", "g"]:
            self.rates[0][i] = self.trans[i] * self.grid.dx[i] * v
        if cond in ["pressure", "press", "p"]:
            self.pressures[0][i] = v

    def set_boundaries(self, b_dict: dict):
        """

        Args:
            b_dict (dict): _description_
        """
        self.b_dict = b_dict
        for i in b_dict.keys():
            [(cond, v)] = b_dict[i].items()
            self.set_boundary(i, cond, v)

    # -------------------------------------------------------------------------
    # Flow Equations:
    # -------------------------------------------------------------------------

    def __calc_RHS(self):
        if self.comp_type == "incompressible":
            n = self.grid.get_n_cells(True)
            self.RHS = 0  # np.zeros(n, dtype=self.dtype)
        elif self.comp_type == "compressible":
            RHS_n = self.grid.V * self.grid.phi * self.comp
            RHS_d = self.factors["volume conversion"] * self.fluid.B * self.dt
            self.RHS = RHS_n / RHS_d

    def get_i_flow_equation(self, i, verbose=False):

        assert i > 0 and i <= self.grid.nx, "grid index is out of range."
        exec(f"p{i}=sym.Symbol('p{i}')")

        if i not in self.flow_equations_terms_dict:

            i_neighbors = self.grid.get_neighbors(i)
            i_boundaries = self.grid.get_boundaries(i)

            if verbose:
                print(f"i: {i} - Neighbors: {i_neighbors} - Boundaries: {i_boundaries}")

            # exec(f"p{i}=sym.Symbol('p{i}')")
            # ToDo: keep pressure constant at specific cell (requires A adjust)
            # if not np.isnan(self.pressures[self.tstep][i]):
            #     exec(f"p{i} = {self.pressures[self.tstep][i]}")
            terms = []

            # 1. Flow from grid neighbors:
            for neighbor in i_neighbors:
                exec(f"p{neighbor} = sym.Symbol('p{neighbor}')")
                # To Do: keep pressure constant at specific cell (requires A adjust)
                # if not np.isnan(self.pressures[self.tstep][neighbor]):
                #     exec(f"p{neighbor} = {self.pressures[self.tstep][neighbor]}")
                exec(
                    f"""n_term = self.trans[{min(neighbor,i)}] * ((p{neighbor} - p{i})
                                - (self.fluid.g * (self.grid.z[{neighbor}] - self.grid.z[{i}])))
                """
                )
                terms.append(locals()["n_term"])

            # 2. Flow from grid boundaries:
            for boundary in i_boundaries:
                exec(f"p{boundary}=sym.Symbol('p{boundary}')")
                if not np.isnan(self.pressures[self.tstep][boundary]):
                    exec(
                        f"""b_term = self.trans[{min(boundary,i)}] * 2 * ((p{boundary} - p{i}) 
                                    - (self.fluid.g * (self.grid.z[{boundary}]-self.grid.z[{i}])))
                    """
                    )
                    exec(
                        f"b_term = b_term.subs(p{boundary}, {self.pressures[self.tstep][boundary]})"
                    )
                else:
                    exec(f"b_term = self.rates[self.tstep][{boundary}]")
                terms.append(locals()["b_term"])

            # 3. Flow from grid well:
            if i in self.wells:
                if "q" in self.wells[i]:
                    terms.append(self.wells[i]["q"])
                else:
                    exec(
                        f"""w_term = - self.wells[{i}]['G'] / (self.fluid.B*self.fluid.mu) * (p{i} - self.wells[{i}]['pwf'])"""
                    )
                    terms.append(locals()["w_term"])

            if verbose:
                print("terms:", terms)

            self.flow_equations_terms_dict[i] = terms

        else:
            terms = self.flow_equations_terms_dict[i]

        # 4. Accumulation term:
        if self.RHS[i] == 0:
            a_term = 0
        else:
            try:
                exec(
                    f"accumulation = self.RHS[{i}] * (p{i} - {self.pressures[self.tstep][i]})"
                )
            except:
                raise Exception("Initial pressure (pi) must be specified")
            a_term = locals()["accumulation"]

        # 4. Overall grid flow equation:
        i_flow_equation = sym.Eq(sum(terms), a_term)
        if (
            i_flow_equation.lhs.as_coefficients_dict()[1] != 0
            or i_flow_equation.rhs.as_coefficients_dict()[1] != 0
        ):
            i_flow_equation = i_flow_equation.simplify()

        # 5. Find lhs and rhs:
        i_lhs = dict(
            sorted(
                i_flow_equation.lhs.as_coefficients_dict().items(),
                key=lambda x: str(x[0]),
            )
        )
        i_rhs = i_flow_equation.rhs
        return i_lhs, i_rhs

    # @lru_cache(maxsize=None)
    def get_flow_equations(self, verbose=False):
        for i in self.grid.order[self.grid.i_blocks.astype("bool")]:
            i_lhs, i_rhs = self.get_i_flow_equation(i, verbose)
            print(f"Grid {i}: {i_lhs}, {i_rhs}")

    # -------------------------------------------------------------------------
    # Coefficient Matrix:
    # -------------------------------------------------------------------------

    # @lru_cache(maxsize=None)
    def update_matrix(self, i, sparse=True, verbose=False):
        """
        Update coefficient matrix (A) and result vector (d).
        Note: arrays are passed by reference.
        """
        i_lhs, i_rhs = self.get_i_flow_equation(i, verbose)
        i_lhs = list(i_lhs.values())
        self.d[i - 1] = np.array(i_rhs).astype(self.dtype)

        if self.tstep == 0:
            if i == 1:
                self.A[i - 1, i - 1 : i + 1] = i_lhs
            elif i == self.grid.nx:
                self.A[i - 1, i - 2 : i] = i_lhs
            else:
                self.A[i - 1, i - 2 : i + 1] = i_lhs

    def __get_matrix(self, sparse=False, verbose=False):
        """Create coefficient matrix (A) and result vector (d).

        Very slow and used only for testing purposes.
        """
        if sparse:
            self.d = ss.lil_matrix((self.grid.nx, 1), dtype=self.dtype)
            self.A = ss.lil_matrix((self.grid.nx, self.grid.nx), dtype=self.dtype)
        else:
            self.d = np.zeros((self.grid.nx, 1), dtype=self.dtype)
            self.A = np.zeros((self.grid.nx, self.grid.nx), dtype=self.dtype)

        for i in range(1, self.grid.nx + 1):
            i_lhs, i_rhs = self.get_i_flow_equation(i, verbose)
            self.d[i - 1] = np.array(i_rhs).astype(self.dtype)
            i_lhs = np.array(list(i_lhs.values())).astype(self.dtype)
            if i == 1:
                self.A[i - 1, i - 1 : i + 1] = i_lhs
            elif i == self.grid.nx:
                self.A[i - 1, i - 2 : i] = i_lhs
            else:
                self.A[i - 1, i - 2 : i + 1] = i_lhs

        if verbose:
            print("- A:\n", self.A)
            print("- d:\n", self.d)

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
                        -self.trans[1:] - self.trans[:-1] - self.RHS[1:-1],
                        self.trans[1:-1],  # East trans for interior blocks
                        self.trans[1:-1],
                    ],  # West trans for interior blocks
                    [0, 1, -1],
                    shape=(self.grid.nx, self.grid.nx),
                    format="lil",  # “dia”, “csr”, “csc”, “lil”
                    dtype=self.dtype,
                )
                self.A[0, 0] = -self.trans[1] - self.RHS[1]
                self.A[-1, -1] = -self.trans[-2] - self.RHS[-1]
                if not sparse:
                    self.A = self.A.toarray()

            # Update matrix if there is pressure or flow in 'west' boundary:
            if (
                not np.isnan(self.pressures[self.tstep][0])
                or self.rates[self.tstep][0] != 0
            ):
                self.update_matrix(1, sparse, verbose)

            # Update matrix if there is pressure or flow in 'east' boundary:
            if (
                not np.isnan(self.pressures[self.tstep][-1])
                or self.rates[self.tstep][-1] != 0
            ):
                # at last grid: self.grid.nx or -2
                self.update_matrix(self.grid.nx, sparse, verbose)

            # Update matrix in wells i_blocks:
            for i in self.wells.keys():
                self.update_matrix(i, sparse, verbose)

            if verbose:
                if sparse:
                    print("- A:\n", self.A.toarray())
                    print("- d:\n", self.d.toarray())
                else:
                    print("- A:\n", self.A)
                    print("- d:\n", self.d)

            return self.A, self.d

    # -------------------------------------------------------------------------
    # Numerical Solution:
    # -------------------------------------------------------------------------

    # @lru_cache(maxsize=None)
    def solve(self, sparse=True, check_MB=True, update=True, verbose=False):

        self.get_matrix(sparse, verbose)
        # self.__get_matrix(sparse, verbose)

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
                            self.trans[min(i, boundary)]
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
        self.transmissibility = self.trans
        self.set_properties = self.__calc_comp

    # -------------------------------------------------------------------------
    # End
    # -------------------------------------------------------------------------


if __name__ == "__main__":
    grid = grids.Cartesian(
        nx=4,
        ny=2,
        nz=2,
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
    fluid = fluids.SinglePhaseFluid(
        mu=0.5, B=1, rho=50, comp=1 * 10**-5, dtype="double"
    )

    model = Model(grid, fluid, pi=4000, dtype="double")
    model.set_well(id=4, q=-600, s=1.5, r=3.5)
    model.set_boundaries({0: {"pressure": 4000}, -1: {"rate": 0}})
    print(model.RHS)
    # model.run(nsteps=6, sparse=False, check_MB=True, verbose=False)
    # print(grid.get_Gx(fshape=True))
    # print(model.pressures)
    # print(model.rates)
