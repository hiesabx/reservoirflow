import time
from collections import defaultdict
from tqdm import tqdm

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm

import torch
import torch.autograd as autograd
import torch.nn as nn

# from reservoirflow.solutions.solution import Solution
from reservoirflow.models.model import Model
from reservoirflow.solutions.neurical import Network
from reservoirflow import scalers
from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import ProcessPoolExecutor


class PINN(Network, nn.Module):
    """PINN solution class.

    PINN is a Physics-Informed-Neural-Network.

    .. caution::
        This class is not available.

    Returns
    -------
    Solution
        Solution object.
    """

    name = "PINN"

    def __init__(
        self,
        model: Model,  #: rf.models.Model,
        sparse: bool = True,
    ):
        """Create Network Solution.

        Parameters
        ----------
        model : Model
            Model object.
        sparse : bool, optional, default: True
            using sparse computing for a better performance.
        """
        super().__init__(model, sparse)
        nn.Module.__init__(self)
        self.dtype = torch.float32
        # self.dtype = torch.float32
        # self.X_train = self.as_tensor(X_train, False)
        # self.Y_train = self.as_tensor(Y_train, False)
        # self.X_test = self.as_tensor(X_test, False)
        # self.Y_test = self.as_tensor(Y_test, False)
        # self.X_phy = self.as_tensor(X_phy, True)
        # self.rates_V_RHS = self.as_tensor(rates_V_RHS, False)
        # self.LHS_RHS = self.as_tensor(LHS_RHS, False)
        # self.rates_V = self.as_tensor(rates_V, False)
        # self.LHSf = self.as_tensor(LHSf, False)
        # self.RHSf = self.as_tensor(RHSf, False)
        # self.Y_phy = torch.zeros(self.X_phy.shape[0], dtype=self.dtype)

        # X, Y = model.get_values(boundary=True, scale=True)
        # self.add(X, Y, "train")

        self.neurons = [2, *[64] * 6, 1]
        lr = 0.01

        self.epochs = 0
        self.epoch = 0

        self.activation = nn.Tanh
        # self.activation = nn.GELU
        self.loss_func = nn.MSELoss
        self.network = self.get_network(self.neurons)
        self.lr = lr

    def assert_data(self):
        assert len(self.X_train) == len(self.Y_train), "input length does not match."
        assert self.X_train.shape[1] == self.neurons[0], "input shape does not match."
        assert (
            self.Y_train.shape[1] == self.neurons[-1] + 1
        ), "output shape does not match."

        assert len(self.X_test) == len(self.Y_test), "input length does not match."
        assert self.X_test.shape[1] == self.neurons[0], "input shape does not match."
        assert (
            self.Y_test.shape[1] == self.neurons[-1] + 1
        ), "output shape does not match."

        assert len(self.X_phy) == len(self.rates_V_RHS), "input length does not match."
        assert len(self.X_phy) == len(self.LHS_RHS), "input length does not match."
        assert self.X_phy.shape[1] == self.neurons[0], "input shape does not match."

    def as_tensor(self, a, grad=False):
        if isinstance(a, pd.DataFrame):
            a = a.values
        return torch.tensor(a, dtype=self.dtype, requires_grad=grad)

    def init_optimizer(self, lr, opt):
        if lr is None:
            lr = self.lr
        if opt.lower() == "adam":
            return torch.optim.Adam(self.parameters(), lr=lr, weight_decay=1e-6)
        elif opt.lower() == "lbfgs":
            return torch.optim.LBFGS(
                self.parameters(),
                lr,
                max_iter=100,
                max_eval=None,
                tolerance_grad=1e-11,
                tolerance_change=1e-11,
                history_size=100,
                line_search_fn="strong_wolfe",
            )
        else:
            raise ValueError("Optimizer is unknown. Use ['adam', 'lbfgs']")

    def get_network(self, neurons):
        bias = True
        self.N_layers = len(neurons)
        self.N_range = range(self.N_layers - 1)
        self.linears = [
            nn.Linear(
                neurons[i],
                neurons[i + 1],
                dtype=self.dtype,
                bias=bias,
            )
            for i in self.N_range
        ]
        self.layers = []
        for i in self.N_range:
            nn.init.xavier_normal_(
                self.linears[i].weight,
                gain=1.0,
            )
            if bias:
                nn.init.zeros_(self.linears[i].bias)
            if i < self.N_layers - 2 and self.activation is not None:
                self.layers.append([self.linears[i], self.activation()])
            else:
                self.layers.append([self.linears[i]])
        return nn.Sequential(*sum(self.layers, []))

    def forward(self, x):
        if not torch.is_tensor(x):
            x = self.as_tensor(x, False)
        with torch.no_grad():
            return self.network(x)

    predict = forward

    def loss_train(self, x, y, reduction="mean"):
        return self.loss_func(reduction=reduction)(self.network(x).flatten(), y[:, 0])

    def loss_phy(self, x):
        p = self.network(x)
        Dp = self.gradient(p, x)
        DDp = self.gradient(Dp[:, 1], x)
        Y_phy_h = (
            self.LHS_RHS * DDp[:, 1].flatten() - Dp[:, 0].flatten() + self.rates_V_RHS
        )  # perfect
        return torch.mean(Y_phy_h**2)
        # return self.loss_func(reduction="mean")(Y_phy_h, self.Y_phy)

    def loss(self):
        """_summary_

        Important Notes
        ---------------
        Loss should be divided into two components as following:
        1. training loss: includes both initial and boundary condition.
        Using separate loss instead for initial and boundary conditions
        is not working.
        2. physics loss: includes the residual of a continuous
        dimensionless PDE. Note that at least one of the derivative
        terms should be separated from any factors (preferably the lower
        order term e.g. Dpt). Factors and derivative must be flatten()
        to have a better performance.

        Returns
        -------
        _type_
            _description_
        """
        loss_train = self.loss_train(self.X_train, self.Y_train)
        loss_phy = 1e-1 * self.loss_phy(self.X_phy)
        return loss_train + loss_phy

    def gradient(self, y_h, x):
        """
        y_h: output (predicted)
        x: input (requires_grad=True)
        Dy: output derivative with respect to input x

        Remarks
        -------
        - create_graph is set to True allowing to compute higher
        order derivative products.
        """
        return autograd.grad(
            y_h,
            x,
            torch.ones_like(y_h),
            retain_graph=True,
            create_graph=True,
        )[0]

    def closure(self):
        self.optimizer.zero_grad()
        loss = self.loss()
        loss.backward()
        # loss.backward(retain_graph=True)
        return loss

    def fit(self, freq=None, epochs=5000, lr=None, opt="adam"):
        run_epoch = self.epochs
        self.epochs += epochs
        self.run_ctime = 0
        if freq is None:
            if epochs < 1:
                raise ValueError("epochs must be higher than 0.")
            elif epochs == 1:
                freq = 1
            elif epochs < 11:
                freq = epochs // 2
            elif epochs < 101:
                freq = epochs // 10
            else:
                freq = epochs // 100

        self.optimizer = self.init_optimizer(lr=lr, opt=opt)

        if not hasattr(self, "results"):
            self.run_id = 1
            self.results = defaultdict(list)
        else:
            self.run_id += 1

        pbar = tqdm(
            range(1, epochs + 1),
            unit="epochs",
            colour="green",
            position=0,
            leave=True,
        )

        init_start_time = time.time()
        for epoch in pbar:
            step_start_time = time.time()
            self.epoch = run_epoch + epoch
            self.optimizer.step(self.closure)
            with torch.no_grad():
                loss_value = (
                    self.loss_train(self.X_test, self.Y_test).detach().float().numpy()
                )
                pbar.set_description(f"[epoch] {epoch} - [loss] {loss_value}")
                if epoch % freq == 0:
                    self.results["run_id"].append(self.run_id)
                    self.results["epoch"].append(self.epoch)
                    self.results["loss"].append(loss_value)
                    self.run_ctime = round(time.time() - step_start_time, 2)
                    self.results["time"].append(self.run_ctime)
                if max(loss_value, 0) == 0.0:
                    print(f"early stop at loss {loss_value}")
                    break

        self.run_ctime = round((time.time() - init_start_time) / 60, 3)
        self.ctime += self.run_ctime
        print(
            f"[info] Pinn fitting of {epochs} epochs",
            f"finished in {self.run_ctime} minutes",
            f"(loss: {loss_value}).",
        )

        plt.plot(self.results["epoch"], self.results["loss"], label="loss")
        plt.ylabel("MSE")
        plt.xlabel("epoch")
        plt.show()

        return pd.DataFrame(self.results)

    def solve(self):
        raise NotImplementedError

    def run(
        self,
        nsteps=10,
        N=100,
        clean=False,
        threading=True,
        vectorize=True,
        check_MB=True,
        # print_arrays=False,
        # isolver=None,
    ):
        """Run the simulation for a given number of steps.

        Parameters
        ----------
        nsteps : int, optional
            Number of steps to run the simulation, by default 10.
        N : int, optional
            Number of epochs in the neurical solution, by default 100.
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

        # Independent variables: t, x
        alpha = self.model.get_alpha(method="mean")
        t, x = self.model.get_domain(scale=False, boundary=True)
        L = x.max() - x.min()
        xD = (x - x.min()) / L
        # tD = alpha * t / (L**2)
        tD = -np.pi**2 * alpha * t / (L**2)
        tD_values, xD_values = np.meshgrid(
            tD,
            xD,
            # sparse=self.sparse,
            indexing="ij",
        )

        # Dependent variable: p
        p = self.pressures
        input_range = [0, 1]  # Analytical solution domain
        input_scaler = scalers.MinMax(output_range=input_range).fit(
            p,
            axis=None,
        )

        # Neurical solution:
        progress = tqdm(
            np.arange(1, self.N + 1),
            unit="steps",
            colour="green",
            position=0,
            leave=True,
            desc="[step]",
        )
        if threading:
            with ThreadPoolExecutor(self.model.n_threads) as executor:
                executor.map(self.__update_pressures_sum, progress)
        else:
            for n in progress:
                self.__update_pressures_sum(n)

        p_pred = pD0

        # Remove values out of range:
        if clean:
            p_pred[p_pred < input_range[0]] = input_range[0]
            p_pred[p_pred > input_range[1]] = input_range[1]

            self.pressures = np.vstack(
                [
                    self.pressures[0, :],
                    input_scaler.inverse_transform(p_pred[1:, :]),
                ]
            )
        else:
            self.pressures = input_scaler.inverse_transform(p_pred)

        self.rates = np.repeat(self.rates, repeats=nsteps + 1, axis=0)
        self.model.update_boundaries_rates_nsteps()

        self.run_ctime = round(time.time() - start_time, 2)
        self.ctime = self.run_ctime
        print(
            f"[info] Simulation run of {nsteps} steps",
            f"finished in {self.run_ctime} seconds.",
        )
        # if check_MB:
        if check_MB:
            self.check_MB_nsteps(self.model.verbose)
            print(f"[info] Material Balance Error: {self.tstep_error}.")

        if verbose_restore:
            self.model.verbose = True
