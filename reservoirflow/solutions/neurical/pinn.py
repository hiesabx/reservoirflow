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
from torch.autograd.functional import hessian

# from reservoirflow.solutions.solution import Solution
from reservoirflow.models.model import Model
from reservoirflow.solutions.neurical import Network
from reservoirflow import scalers
from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import ProcessPoolExecutor


class RelativeL2Loss(nn.Module):
    def __init__(self, epsilon=1e-6, reduction="mean"):
        super(RelativeL2Loss, self).__init__()
        self.epsilon = epsilon
        self.reduction = reduction

    def forward(self, y_h, y):
        # Flatten along all but batch dimension
        diff = y_h - y
        num = torch.sum(diff**2, dim=tuple(range(1, diff.dim())))
        denom = torch.sum(y**2, dim=tuple(range(1, y.dim()))) + self.epsilon
        rel_loss = num / denom
        if self.reduction == "mean":
            return torch.mean(rel_loss)
        elif self.reduction == "sum":
            return torch.sum(rel_loss)
        else:
            return rel_loss


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

        # Neural network parameters:
        self.epochs = 0
        self.epoch = 0
        self.neurons = [2, *[64] * 6, 1]
        self.bias = True
        # self.activation = nn.Tanh  # nn.Tanh, nn.Sigmoid, nn.GELU
        self.loss_func = nn.MSELoss()  # nn.MSELoss, nn.L1Loss
        # self.loss_func = RelativeL2Loss()  # nn.MSELoss, nn.L1Loss
        self.network = self.get_network()
        # self.network = torch.compile(self.network)
        self.w_y = 0.95  # weight for the training loss
        self.w_r = 1 - self.w_y  # weight for the physics loss

        # Early stopping:
        self.patience = 10
        self.min_delta = 1e-4
        self.counter = 0
        self.min_validation_loss = float("inf")

    def early_stop(self, validation_loss):
        if validation_loss < self.min_validation_loss:
            self.min_validation_loss = validation_loss
            self.counter = 0
        elif validation_loss > (self.min_validation_loss + self.min_delta):
            self.counter += 1
            if self.counter >= self.patience:
                return True
        return False

    def assert_data(self):
        assert len(self.X_data) == len(self.Y_data), "input length does not match."
        assert self.X_data.shape[1] == self.neurons[0], "input shape does not match."
        assert self.Y_data.shape[1] == self.neurons[-1], "output shape does not match."
        assert len(self.X_test) == len(self.Y_test), "input length does not match."
        assert self.X_test.shape[1] == self.neurons[0], "input shape does not match."
        assert self.Y_test.shape[1] == self.neurons[-1], "output shape does not match."
        # assert len(self.X_r) == len(self.F_px_pt), "input length does not match."
        assert self.X_r.shape[1] == self.neurons[0], "input shape does not match."

    def as_tensor(self, a, grad=False):
        if isinstance(a, pd.DataFrame):
            a = a.values
        return torch.tensor(a, dtype=self.dtype, requires_grad=grad)

    def init_optimizer(
        self,
        opt="adam",
        lr=None,
        weight_decay=None,
    ):
        if opt.lower() == "adam":
            return torch.optim.Adam(
                self.parameters(),
                lr=lr,
                weight_decay=weight_decay,
            )
        elif opt.lower() == "adamw":
            return torch.optim.AdamW(
                self.parameters(),
                lr=lr,
                weight_decay=weight_decay,
            )
        elif opt.lower() == "sgd":
            return torch.optim.SGD(
                self.parameters(),
                lr=lr,
                momentum=0.9,
                weight_decay=weight_decay,
            )
        elif opt.lower() == "adagrad":
            return torch.optim.Adagrad(
                self.parameters(),
                lr=lr,
                weight_decay=weight_decay,
            )
        elif opt.lower() == "rmsprop":
            return torch.optim.RMSprop(
                self.parameters(),
                lr=lr,
                alpha=0.99,
                eps=1e-8,
                weight_decay=weight_decay,
            )
        elif opt.lower() == "adamax":
            return torch.optim.Adamax(
                self.parameters(),
                lr=lr,
                betas=(0.9, 0.999),
                eps=1e-8,
                weight_decay=weight_decay,
            )
        elif opt.lower() == "adadelta":
            return torch.optim.Adadelta(
                self.parameters(),
                lr=lr,
                rho=0.9,
                eps=1e-6,
                weight_decay=weight_decay,
            )
        elif opt.lower() == "lbfgs":
            return torch.optim.LBFGS(
                self.parameters(),
                lr=lr,
                max_iter=3,
                max_eval=4,
                history_size=10,
                tolerance_grad=1e-7,
                tolerance_change=1e-7,  # 1 * np.finfo(float).eps
                line_search_fn="strong_wolfe",
            )
        else:
            raise ValueError("Optimizer is unknown. Use ['adam', 'lbfgs']")

    def get_network(self):
        self.N_layers = len(self.neurons)
        self.N_range = range(self.N_layers - 1)
        self.linears = [
            nn.Linear(
                self.neurons[i],
                self.neurons[i + 1],
                dtype=self.dtype,
                bias=self.bias,
            )
            for i in self.N_range
        ]
        activation = nn.Tanh()  # nn.GELU(), nn.SiLU()
        n_act = self.N_layers - 2
        activations = [
            *[activation] * n_act,
            None,
        ]
        self.layers = []
        for i in self.N_range:
            nn.init.xavier_normal_(
                self.linears[i].weight,
                gain=1.0,
            )
            if self.bias:
                nn.init.zeros_(self.linears[i].bias)
            # if i < self.N_layers - 3 and self.activation is not None:
            # self.layers.append([self.linears[i], self.activation()])
            # else:
            # self.layers.append([self.linears[i]])
            if activations[i] is not None:
                self.layers.append([self.linears[i], activations[i]])
            else:
                self.layers.append([self.linears[i]])
        return nn.Sequential(*sum(self.layers, []))

    def forward(self, x):
        if not torch.is_tensor(x):
            x = self.as_tensor(x, False)
        with torch.no_grad():
            return self.network(x)

    predict = forward

    def loss_y(self, x, y):
        return self.loss_func(self.network(x), y)

    def loss_y_value(self, x, y):
        return self.loss_y(x, y).detach().float().numpy()

    def loss_r(self, x):
        """Physics loss function.

        This function calculates the residual of a PDE as:
            y_phy_h = LHS_RHS * d2p_dx2 - dp_dt
        Where:
            - x is the input in form of (time, x)
            - p is the pressure
            - dp_dt is the first temporal derivative
            - d2p_dx2 is the second spatial derivative.
            - LHS_RHS is the ratio of the left-hand side factor and right-hand side factor of the PDE.

        Parameters
        ----------
        x : torch.Tensor
            Input tensor of shape (N, 2), where N is the number of samples.

        Returns
        -------
        torch.Tensor
            Residual of the PDE, averaged over all samples.
        """

        # prediction: p(t, x)
        p = self.network(x)
        # # rates:
        # # q = self.update_rates()
        # # first derivative: dp/dt, dp/dx
        dp = self.gradient(p, x)
        # # first temporal derivative: dp/dt
        dp_dt = dp[:, 0]  # .flatten()
        # # second derivative: d2p/dt2, d2p/dx2
        d2p = self.gradient(dp[:, 1], x)
        # # second spacial derivative: d2p/dx2
        d2p_dx2 = d2p[:, 1]  # .flatten()
        # residual of the PDE:
        # r = self.F_px_pt * d2p_dx2 - dp_dt  # + self.F_q_pt * q
        # return the mean of the squared residual:
        # return torch.mean(r**2)
        return self.loss_func(self.F_px_pt * d2p_dx2, dp_dt)

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
            outputs=y_h,
            inputs=x,
            grad_outputs=torch.ones_like(y_h),
            retain_graph=True,
            create_graph=True,  # allows to compute higher order derivatives
        )[0]

    def loss(self):
        """Calculate the total loss.

        Important Notes
        ---------------
        Loss should be divided into two components as following:
        1. training loss: includes both initial and boundary condition.
        Using separate loss instead for initial and boundary conditions
        is not working.
        2. physics loss: includes the residual of a continuous
        dimensionless PDE. Note that at least one of the derivative
        terms should be separated from any factors (preferably the lower
        order term e.g. dp_dt). Factors and derivative must be flatten()
        to have a better performance.

        Returns
        -------
        torch.Tensor
            Total loss, which is the sum of training loss and physics loss.
        """
        # Training loss:
        loss_d = self.loss_y(self.X_data, self.Y_data)
        # Physics loss:
        loss_r = self.loss_r(self.X_r)
        # Total loss:
        return self.w_d * loss_d + self.w_r * loss_r

    def closure(self):
        self.optimizer.zero_grad()
        # Data loss:
        l_d = self.loss_y(self.X_data, self.Y_data)
        # Residual loss:
        l_r = self.loss_r(self.X_r)
        # Total loss:
        l_t = self.w_d * l_d + self.w_r * l_r
        l_t.backward()  # loss.backward(retain_graph=True)
        # Update losses:
        self.l_d = l_d.detach().float().numpy()
        self.l_r = l_r.detach().float().numpy()
        self.l_t = l_t.detach().float().numpy()
        return l_t

    def fit(
        self,
        freq=None,
        epochs=5000,
        opt="adam",
        lr=0.01,
        weight_decay=None,
    ):
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

        self.optimizer = self.init_optimizer(
            opt=opt,
            lr=lr,
            weight_decay=weight_decay,
        )

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
                loss_test = self.loss_y_value(self.X_test, self.Y_test)
                pbar.set_description(
                    (
                        f"[epoch] {epoch} - "
                        f"[loss total] {self.l_t:.6f} - "
                        f"[loss test] {loss_test:.6f} - "
                        f"[loss data] {self.l_d:.6f} - "
                        f"[loss residual] {self.l_r:.6f}"
                    )
                )
                if epoch % freq == 0:
                    self.results["run_id"].append(self.run_id)
                    self.results["epoch"].append(self.epoch)
                    self.results["loss (total)"].append(self.l_t)
                    self.results["loss (test)"].append(loss_test)
                    self.results["loss (data)"].append(self.l_d)
                    self.results["loss (residual)"].append(self.l_r)
                    self.results["loss_d + loss_r"].append(self.l_d + self.l_r)
                    self.run_ctime = round(time.time() - step_start_time, 2)
                    self.results["time"].append(self.run_ctime)
                # if self.early_stop(loss_y_test):
                #     break

        # Reset early stopping counter:
        self.counter = 0

        self.run_ctime = round((time.time() - init_start_time) / 60, 3)
        self.ctime += self.run_ctime
        print(
            f"[info] Pinn fitting of {epochs} epochs",
            f"finished in {self.run_ctime} minutes - ",
            f"loss (test): {loss_test:.6f}.",
        )

    def plot(self):
        plt.plot(
            self.results["epoch"],
            self.results["loss (total)"],
            label="loss (total: self.w_y * loss_d + self.w_r * loss_r)",
        )
        plt.plot(
            self.results["epoch"],
            self.results["loss_d + loss_r"],
            label="loss (total: loss_d + loss_r)",
        )
        plt.plot(
            self.results["epoch"],
            self.results["loss (test)"],
            label="loss (test)",
        )
        plt.plot(
            self.results["epoch"],
            self.results["loss (data)"],
            label="loss (data)",
        )
        plt.plot(
            self.results["epoch"],
            self.results["loss (residual)"],
            label="loss (residual)",
        )
        plt.title(f"PINN {self.name} - {self.model.name} - {self.run_id}")
        plt.legend()
        plt.grid()
        plt.ylabel("Loss")
        plt.xlabel("epoch")
        plt.show()

        return pd.DataFrame(self.results)

    def solve(self):
        raise NotImplementedError

    def update_pressures(
        self,
        clean=True,
    ):
        p = self.network(self.X)
        output_shape = self.model.get_shape(boundary=True)
        p = p.detach().numpy().reshape(output_shape)
        if clean:
            Vmin = self.model.pressure_scaler.Vmin
            Vmax = self.model.pressure_scaler.Vmax
            p[p < Vmin] = Vmin
            p[p > Vmax] = Vmax

            self.pressures = np.vstack(
                [
                    self.pressures[0, :],
                    self.model.pressure_scaler.inverse_transform(p[1:, :]),
                ]
            )
        else:
            self.pressures = self.model.pressure_scaler.inverse_transform(p)

    def update_rates(self):
        self.update_pressures(clean=True)
        self.model.update_boundaries_rates_nsteps()
        q = self.model.get_df(
            columns=["cells_rate"],
            boundary=True,
            scale=False,
            units=False,
            melt=True,
            drop_zero=False,
            drop_nan=True,
        )["Q"].values
        return self.as_tensor(q, False)

    def run(
        self,
        nsteps=10,
        N=100,
        clean=False,
        reference=None,
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
        self.model.update_shapes()
        self.N = N
        self.run_ctime = 0
        if self.model.verbose:
            self.model.verbose = False
            verbose_restore = True
        else:
            verbose_restore = False

        print(f"[info] Simulation run started: {nsteps} timesteps.")

        # Scale the model:
        scale = True

        # Training Data: initial and boundary conditions
        X, Y = self.model.get_data(
            times_id=None,
            cells_id=None,
            scale=scale,
            drop_nan=True,
            shuffle=False,
        )
        self.X_data = self.as_tensor(X, False)
        self.Y_data = self.as_tensor(Y, False)
        n_d = self.X_data.shape[0]

        # Testing Data: based on a reference solution if available
        if reference is not None:
            self.model.set_solution(reference)
            X, Y = self.model.get_data(
                times_id=slice(1, -1, 1),
                cells_id=None,
                scale=scale,
                drop_nan=False,
                shuffle=False,
            )
            self.model.set_solution(self.name)
        self.X_test = self.as_tensor(X, False)
        self.Y_test = self.as_tensor(Y, False)

        # Full Domain as X(t, x):
        X = self.model.get_X(
            times_id=None,  # None for all times
            cells_id=None,  # None for all cells
            scale=scale,
            shuffle=False,
        )
        self.X = self.as_tensor(X, False)

        # Residual Domain as X(t, x):
        X_r = self.model.get_X(
            times_id=slice(0, -1, 1),
            cells_id=self.model.grid.get_cells_id(boundary=True),
            scale=scale,
            shuffle=False,
        )
        self.X_r = self.as_tensor(X_r, True)
        n_r = self.X_r.shape[0]

        n_t = n_d + n_r
        self.w_r = n_d / n_t
        self.w_d = n_r / n_t
        print("w_r", self.w_r)
        print("w_d", self.w_d)

        # Differential Factors of the PDE:
        F_px, F_pt, F_q = self.model.get_factors(
            boundary=True,
            scale=scale,
            method="mean",
        )
        # F_px *= 1.0404  # diff_term
        # F_px *= 1 + self.w_r
        self.F_px_pt = F_px / F_pt

        # if initial_r:
        #     nsteps_r = self.nsteps
        # else:
        #     nsteps_r = self.nsteps - 1

        # F_px_pt = F_px / F_pt
        # F_px_pt = np.tile(F_px_pt, nsteps_r)
        # self.F_px_pt = self.as_tensor(F_px_pt, False)

        # F_pt_px = F_pt / F_px
        # F_pt_px = np.tile(F_pt_px, nsteps_r)
        # self.F_pt_px = self.as_tensor(F_pt_px, False)

        # F_px = np.tile(F_px, nsteps_r)
        # self.F_px = self.as_tensor(F_px, False)
        # F_pt = np.tile(F_pt, nsteps_r)
        # self.F_pt = self.as_tensor(F_pt, False)

        # Update Rates: update_rates (did not work as expected)
        # F_q_pt = F_q / F_pt
        # F_q_pt = np.tile(F_q_pt, nsteps_r)
        # self.F_q_pt = self.as_tensor(F_q_pt, False)

        # Check if the data is correct:
        self.assert_data()

        # opts: ["adam", "adamw", "sgd", "adagrad", "rmsprop", "adamax", "adadelta", "lbfgs"]
        self.fit(epochs=100, opt="adam", lr=0.01, weight_decay=1e-6)
        # self.fit(epochs=300, opt="adam", lr=0.001, weight_decay=1e-6)
        self.fit(epochs=100, opt="lbfgs", lr=0.01)
        # self.fit(epochs=400, opt="adam", lr=0.0001, weight_decay=1e-6)
        # self.fit(epochs=10, opt="lbfgs", lr=0.01)
        # self.fit(epochs=400, opt="adam", lr=0.0001)
        # self.fit(epochs=100, lr=0.0001, opt="adam")
        # self.fit(epochs=400, lr=0.0001, opt="adam")
        # self.fit(epochs=1000, lr=0.00001, opt="adam")
        # self.fit(epochs=1000, lr=0.00001, opt="adam")

        # Prediction:
        self.update_pressures(clean=clean)

        # Update boundaries rates with the new pressures:
        # self.model.update_rates_shape()
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
