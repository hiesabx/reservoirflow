from abc import ABC, abstractmethod
import warnings

# import reservoirflow as rf
# from reservoirflow.models import Model


class Solution(ABC):
    """Abstract solution class.

    .. attention::

        This is an abstract class and can't be instantiated. This class
        is only used as a parent for other classes of ``solutions``
        submodules.

    Returns
    -------
    Solution
        Solution object.
    """

    def __init__(
        self,
        model,  #: Model,  #: rf.models.Model,
        sparse,
    ):
        """Construct solution object.

        Parameters
        ----------
        model : Model
            a model object from ``models`` module.
        sparse : bool
            using sparse computing for a better performance.
        """
        self.model = model
        self.sparse = sparse
        self.pressures, self.rates = self.model.get_init_arrays()
        self.nsteps = 1
        self.tstep = 0
        self.ctime = 0
        self.tstep_error = 0  # timestep error
        self.ctime_error = 0  # comulative error

    @property
    def total_error(self):
        return self.tstep_error + self.ctime_error

    @abstractmethod
    def solve(self):
        """Solve a single timestep.

        .. attention::
            This is an abstract method.
        """

    @abstractmethod
    def run(self):
        """Solve multiple timesteps.

        .. attention::
            This is an abstract method.
        """

    # -------------------------------------------------------------------------
    # Material Balance:
    # -------------------------------------------------------------------------

    def check_MB(self, tstep, verbose=False, error_threshold=0.1):
        """Material Balance Check

        Parameters
        ----------
        verbose : bool, optional
            print values.
        error_threshold : float, optional
            maximum limit of allowed error.
        """
        if verbose:
            print(f"[info] Error in step {tstep}")

        if self.model.comp_type == "incompressible":
            # rates must add up to 0:
            self.tstep_error = self.rates[tstep].sum()
            if verbose:
                print(f"[info]    - Error: {self.tstep_error}")
        elif self.model.comp_type == "compressible":
            # error over a timestep:
            self.tstep_error = (
                self.model.RHS[self.model.grid.cells_id]
                * (
                    self.pressures[tstep, self.model.grid.cells_id]
                    - self.pressures[tstep - 1, self.model.grid.cells_id]
                )
            ).sum() / self.rates[tstep].sum()
            # error from initial timestep to current timestep: (less accurate)
            self.ctime_error = (
                self.model.RHS[self.model.grid.cells_id]
                * self.model.dt
                * (
                    self.pressures[tstep, self.model.grid.cells_id]
                    - self.pressures[0, self.model.grid.cells_id]
                )
            ).sum() / (self.model.dt * tstep * self.rates.sum())
            self.tstep_error = abs(self.tstep_error - 1)
            if verbose:
                print(f"[info]    - Incremental Error: {self.tstep_error}")
                print(f"[info]    -  Cumulative Error: {self.ctime_error}")
                print(f"[info]    -       Total Error: {self.total_error}")

        if verbose and abs(self.tstep_error) > error_threshold:
            warnings.warn("High material balance error.")
            print(
                f"[warn] Material balance error ({self.tstep_error})",
                f"in step {tstep}",
                f"is higher than the allowed error ({error_threshold}).",
            )

    def check_MB_nsteps(self, verbose=False, error_threshold=0.1):
        for tstep in range(1, self.nsteps):
            self.check_MB(
                tstep=tstep,
                verbose=verbose,
                error_threshold=error_threshold,
            )

    def __str__(self):
        return (
            f"{self.name}"
            # + f"-{self.stype}"
            # + f"-{self.method}"
            # + f"-{self.sparse}"
        )

    def __repr__(self):
        return (
            "rf.solutions.Solution("
            + f"name='{self.name}'"
            # + f"model='{self.model.name}'"
            # + f", sparse={self.sparse}"
            # + f", name='{self.name}'"
            # + f", solution='{self.model.solution.name}'"
            + ")"
        )
