"""
2D Line Plot
============
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

from reservoirflow.plots.plot import Plot


class Plot2D(Plot):
    """1D Line Plot class for visualizing 2D data.

    Parameters
    ----------
    Plot : _type_
        _description_
    """

    name = "2D Plot"

    def __init__(
        self,
        verbose=False,
        error=False,
        tsteps=None,
    ):
        """_summary_

        Parameters
        ----------
        dtype : str, optional
            _description_, by default "double"
        unit : str, optional
            _description_, by default "field"
        verbose : bool, optional
            _description_, by default False
        """
        super().__init__(verbose, error)
        self.tsteps = tsteps
        self.nrows = 3
        self.ncols = 3

        if tsteps is not None:
            assert (
                len(tsteps) == self.nrows * self.ncols
            ), "tsteps are not compatible with the number of rows and columns."

    def __set_axis_labels(self, axs):
        for i, ax in enumerate(axs.ravel()):
            if i >= 6:
                ax.set_xlabel("x")
            if i in [0, 3, 6]:
                ax.set_ylabel("p")

    def __add_legends(self, fig, N):
        # https://stackoverflow.com/a/59393045/11549398
        labels_handles = {
            label: handle
            for ax in fig.axes
            for handle, label in zip(*ax.get_legend_handles_labels())
        }

        fig.legend(
            labels_handles.values(),
            labels_handles.keys(),
            loc="upper center",
            bbox_to_anchor=(0.5, 0.03),
            ncol=6,  # ncol=N,
        )

    def __get_lims_ticks(self, x, Y, ylims):

        xmin = x.min()
        xmax = x.max()
        xmax_ = round(xmax, 1)
        xstep = (xmax_ - xmin) / 4
        xmin_ = xmin - xstep * 0.25
        xlim = (xmin_, xmax * 1.1)
        xticks = np.linspace(xmin, xmax_, 5)

        if ylims is not None:
            ymin = ylims[0]
            ymax = ylims[1]
        else:
            ymin = Y.min()
            ymax = Y.max()
        # n = 0 if ymax < 10 else -3
        ymax_ = ymax  # round(ymax, n)
        ystep = (ymax_ - ymin) / 4
        ymin_ = ymin - ystep * 0.25
        ylim = (ymin_, ymax_)
        if ymax_ > ymax:
            ystep = ystep if ymax < 2 else 0
        m = 5 if ystep == 0 else 6
        yticks = np.linspace(ymin, ymax_ + ystep, m)

        return (xlim, xticks), (ylim, yticks)

    def plot(
        self,
        ylims,
        # id=None,
    ):

        # fig = plt.figure()
        # ax = fig.add_subplot(projection="3d")
        # arr = generate_arr(rsm)
        # # Initial points:
        # ax.scatter(
        #     arr[0, :, 0], arr[0, :, 1], arr[0, :, 2], s=10, c="g", label="Initial"
        # )
        # # Boundary points:
        # ids_b = np.append(rsm.boundaries_id, list(rsm.wells.keys()), axis=0).astype(int)
        # ax.scatter(
        #     arr[:, ids_b, 0],
        #     arr[:, ids_b, 1],
        #     arr[:, ids_b, 2],
        #     s=10,
        #     c="r",
        #     label="Boundary",
        # )
        # # FDM:
        # ax.plot_wireframe(
        #     arr[:, :, 0],
        #     arr[:, :, 1],
        #     arr[:, :, 2],
        #     rstride=5,
        #     cstride=5,
        #     alpha=0.3,
        #     color="k",
        #     label="linear interp",
        # )
        # # PINN:
        # Y_test_h = (
        #     rsm_pinn.predict(arr[:, :, :2]).reshape(arr[:, :, 0].shape).detach().numpy()
        # )
        # ax.plot_wireframe(
        #     arr[:, :, 0],
        #     arr[:, :, 1],
        #     Y_test_h,
        #     rstride=5,
        #     cstride=5,
        #     alpha=0.3,
        #     color="m",
        #     label="linear interp",
        # )
        # ax.set_xlabel("t")
        # ax.set_ylabel("x")
        # ax.set_zlabel("p")
        # plt.legend()
        # plt.show()

        fig, axs = plt.subplots(
            self.nrows,
            self.ncols,
            figsize=(10, 6),
            sharey=True,
            sharex=True,
            subplot_kw=dict(projection="3d"),
        )

        plt.subplots_adjust(hspace=0.3, wspace=0.2)
        alpha = 1
        # tsteps = [0, 1, 5, 10, 40, 50, 60, 90, 100]

        # assert len(tsteps) == self.nrows * self.ncols, "tsteps are not compatible"

        # self.__set_axis_labels(axs)
        for i, ax in enumerate(axs.ravel()):
            if i >= 6:
                ax.set_xlabel("x")
                ax.set_ylabel("y")
            if i in [0, 3, 6]:
                ax.set_zlabel("p")

        # if id is None:
        #     names = list(self.Data.keys())
        #     N = len(names)
        # else:
        #     if isinstance(id, list):
        #         lst = list(self.Data.keys())
        #         names = [lst[i] for i in id]
        #         N = len(names)
        #     else:
        #         names = [list(self.Data.keys())[id]]
        #         N = 1

        names = list(self.Data.keys())
        N = len(names)
        for n in range(N):
            name = names[n]
            lstyle = "-" if "ana" in name.lower() else "--"
            X, Y = self.Data[name]
            print(X.shape, Y.shape)
            t = X[:, :, :, :, 0]
            x = X[0, :, :, :, 1].squeeze()
            y = X[0, :, :, :, 2].squeeze()
            print(
                x.shape,
                y.shape,
                Y[0].squeeze().shape,
            )
            if self.tsteps is None:
                pass
                # tsteps = get_sampled_array_combined_distribution(
                #     t,
                #     self.nrows * self.ncols,
                # )
            else:
                tsteps = self.tsteps

            assert len(tsteps) == self.nrows * self.ncols, "tsteps are not compatible"

            for i, ax in enumerate(axs.ravel()):
                tstep = tsteps[i]

                # ax.plot_wireframe(
                ax.plot_wireframe(
                    x,
                    y,
                    Y[tstep].squeeze(),
                    rstride=5,
                    cstride=5,
                    alpha=0.3,
                    # color='m',
                    cmap=cm.coolwarm,
                    label=name,
                )
                # ax.plot(
                #     x,
                #     Y[tstep, :],
                #     label=name,
                #     linestyle=lstyle,
                #     alpha=alpha,
                # )
                ax.tick_params(top=True, right=True)
                # ax.set_title(f"t={t[tstep]}")
                ax.grid(True)
                if i == 0 and n == 0:
                    (xlim, xticks), (ylim, yticks) = self.__get_lims_ticks(
                        x,
                        Y,
                        ylims,
                    )
                    ax.set_xlim(xlim)
                    ax.set_ylim(ylim)
                    ax.set_xticks(xticks)
                    ax.set_yticks(yticks)

        self.__add_legends(fig, N)
        return plt.show()

    # -------------------------------------------------------------------------
    # End
    # -------------------------------------------------------------------------


if __name__ == "__main__":
    plot = Plot2D()
    print(plot)
