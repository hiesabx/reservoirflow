import time

import numpy as np
from tqdm import tqdm

from reservoirflow import scalers
from reservoirflow.solutions.solution import Solution

# from reservoirflow.utils.profme import cProfiler


class D1P1(Solution):
    """D1P1 solution class.

    D1P1 is a 1-Dimension-1-Phase.

    .. caution::
        This class is not available.

    Returns
    -------
    Solution
        Solution object.
    """

    name = "D1P1"

    def calc_solution(
        self,
        N=1000,
        scale=False,
        output_range=[-1, 1],
        clean=True,
    ):

        # Independent variables: t, x
        t, x = self.model.get_domain(scale=False, boundary=True)
        L = x.max() - x.min()
        xD = (x - x.min()) / L
        tD = self.model.get_alpha() * t / (L**2)
        tD_values, xD_values = np.meshgrid(tD, xD, indexing="ij")

        # Dependent variable: p
        p = self.pressures
        input_range = [0, 1]
        input_scaler = scalers.MinMax(input_range).fit(p, axis=None)
        pDi = input_scaler.transform(self.model.pi)
        pD0 = input_scaler.transform(p[0, 0])
        pDn = input_scaler.transform(p[0, -1])

        # Analytical solution:
        N_range = np.arange(1, N + 1)
        pDi0 = pDi - pD0
        pDin = pDn - pDi
        xDpi = np.pi * xD_values
        tDpi = -np.pi**2 * tD_values
        pDsum = np.zeros_like(xD_values, dtype="double")
        for n in N_range:
            pDsum += (
                (pDi0 / n + pDin * ((-1) ** n) / n)
                * np.sin(n * xDpi)
                * np.exp((n**2) * tDpi)
            )
        pD = pD0 + (pDn - pD0) * xD_values + 2 / np.pi * pDsum

        # Remove values out of range:
        if clean:
            pD[pD < input_range[0]] = input_range[0]
            pD[pD > input_range[1]] = input_range[1]

        shape = self.model.get_shape(True)
        t, x = self.model.get_domain(scale=scale, boundary=True)
        t_values, x_values = np.meshgrid(t, x, indexing="ij")
        X = np.stack((t_values, x_values), axis=-1).reshape(*shape, 2)

        if scale:
            output_scaler = scalers.MinMax(output_range, input_range)
            p = output_scaler.transform(pD)
            # xD_values = output_scaler.transform(xD_values)
            # X = np.stack((tD_values, xD_values), axis=-1).reshape(*shape, 2)
        else:
            p = input_scaler.inverse_transform(pD)
            # t_values, x_values = np.meshgrid(t, x, indexing='ij')
            # X = np.stack((t_values, x_values), axis=-1).reshape(*shape, 2)

        return X, p

    def solve(self):
        raise NotImplementedError

    def run(
        self,
        nsteps=10,
        N=100,
        clean=False,
        # threading=True,
        vectorize=True,
        # check_MB=True,
        # print_arrays=False,
        # isolver=None,
    ):
        """Run the simulation for a given number of steps.

        Parameters
        ----------
        nsteps : int, optional
            Number of steps to run the simulation, by default 10.
        N : int, optional
            Number of terms in the analytical solution, by default 100.
        clean : bool, optional
            If True, remove values out of range in the analytical solution,
            by default False.
        """
        start_time = time.time()
        self.tstep += nsteps
        self.nsteps += nsteps
        self.N = N
        self.run_ctime = 0
        if self.model.verbose:
            self.model.verbose = False
            verbose_restore = True
        else:
            verbose_restore = False

        print(f"[info] Simulation run started: {nsteps} timesteps.")

        N_range = np.arange(1, self.N + 1)

        # scale = False
        # output_range = [-1, 1]
        # clean = False

        # for step in pbar:
        #     pbar.set_description(f"[step] {step}")
        # self.solve(
        #     threading,
        #     vectorize,
        #     check_MB,
        #     True,
        #     print_arrays,
        #     isolver,
        # )
        # self.pressures = self.calc_solution(
        #     N=1000,
        #     scale=False,
        #     output_range=[-1, 1],
        #     clean=True,
        #     )

        # Independent variables: t, x
        t, x = self.model.get_domain(scale=False, boundary=True)
        L = x.max() - x.min()
        xD = (x - x.min()) / L
        tD = self.model.get_alpha() * t / (L**2)
        tD_values, xD_values = np.meshgrid(tD, xD, indexing="ij")

        # Dependent variable: p
        p = self.pressures
        input_range = [0, 1]  # Analytical solution domain
        input_scaler = scalers.MinMax(output_range=input_range).fit(p, axis=None)
        pDi = input_scaler.transform(self.model.pi)
        pD0 = input_scaler.transform(p[0, 0])
        pDn = input_scaler.transform(p[0, -1])

        # Analytical solution:
        pDi0 = pDi - pD0
        pDin = pDn - pDi
        xDpi = np.pi * xD_values
        tDpi = -np.pi**2 * tD_values
        pDsum = np.zeros_like(xD_values, dtype="double")

        for n in tqdm(
            N_range,
            unit="steps",
            colour="green",
            position=0,
            leave=True,
            desc="[step]",
        ):
            # pbar.set_description(f"[step] {n}")
            pDsum += (
                (pDi0 / n + pDin * ((-1) ** n) / n)
                * np.sin(n * xDpi)
                * np.exp((n**2) * tDpi)
            )
        pD = pD0 + (pDn - pD0) * xD_values + 2 / np.pi * pDsum

        # Remove values out of range:
        if clean:
            pD[pD < input_range[0]] = input_range[0]
            pD[pD > input_range[1]] = input_range[1]
            self.pressures = np.vstack(
                [
                    self.pressures[0, :],
                    input_scaler.inverse_transform(pD[1:, :]),
                ]
            )
        else:
            self.pressures = input_scaler.inverse_transform(pD)

        # shape = self.model.get_shape(True)
        # t, x = self.model.get_domain(scale=scale, boundary=True)
        # t_values, x_values = np.meshgrid(t, x, indexing="ij")
        # X = np.stack((t_values, x_values), axis=-1).reshape(*shape, 2)

        # if scale:
        #     output_scaler = scalers.MinMax(output_range, input_range)
        #     p = output_scaler.transform(pD)
        #     # xD_values = output_scaler.transform(xD_values)
        #     # X = np.stack((tD_values, xD_values), axis=-1).reshape(*shape, 2)
        # else:
        #     p = input_scaler.inverse_transform(pD)
        #     # t_values, x_values = np.meshgrid(t, x, indexing='ij')
        #     # X = np.stack((t_values, x_values), axis=-1).reshape(*shape, 2)

        self.run_ctime = round(time.time() - start_time, 2)
        self.ctime = self.run_ctime
        print(
            f"[info] Simulation run of {nsteps} steps",
            f"finished in {self.run_ctime} seconds.",
        )
        # if check_MB:
        # print(f"[info] Material Balance Error: {self.error}.")

        if verbose_restore:
            self.model.verbose = True
