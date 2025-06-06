"""
1D Plot
=======
"""

import numpy as np
import matplotlib.pyplot as plt

from reservoirflow.plots.plot import Plot


class Contour1D(Plot):
    """1D Contour Plot class for visualizing 1D data.

    Parameters
    ----------
    Plot : _type_
        _description_
    """

    name = "1D Contour Plot"

    def __init__(
        self,
        # nrows=1,
        # ncols=1,
        verbose=False,
    ):
        """Create 1D Plot object.

        Parameters
        ----------
        nrows : int
            number of rows in the plot.
        ncols : int
            number of columns in the plot.
        verbose : bool, optional
            print verbose output, by default False
        """
        super().__init__(verbose)
        # self.nrows = nrows
        # self.ncols = ncols

    def plot(
        self,
    ):
        """Draw pressures.

        Parameters
        ----------
        X : ndarray
            2D-array where rows are time steps and columns are cells centers.
        Y : ndarray
            _description_
        ax : _type_, optional
            _description_, by default None
        fig_type : str, optional
            string from a list ['contourf', 'contour', 'pcolor'], by default
            'contourf'
        title : str, optional
            string, by default 'FDM'
        cbar_num : int, optional
            number of contour colors, by default 11

        Raises
        ------
        ValueError
            _description_
        """

        name = list(self.Data.keys())[0]
        ax = None
        fig_type = "contourf"
        cbar_num = 11
        Y_range = None

        if ax is None:
            fig, ax = plt.subplots(figsize=(6, 3))

        X = self.Data[name][0]
        Y = self.Data[name][1]
        t = X[:, :, 0]
        x = X[:, :, 1]
        if Y_range is None:
            ymin = Y.min()
            ymax = Y.max()
        else:
            ymin = Y_range[0]
            ymax = Y_range[1]

        # xmin = x.min()
        # xmax = x.max()

        if ymax >= 2:
            ylim = (x.min(), x.max())
            yticks = np.linspace(ylim[0], ylim[1], 5)
            # cbar_ticks = np.linspace(round(ymin, -2), round(ymax, -2), cbar_num)
            cbar_ticks = np.linspace(ymin, round(ymax, -2), cbar_num)
        else:
            ylim = (x.min(), x.max())
            yticks = np.linspace(ylim[0], ylim[1], 5)
            cbar_ticks = np.linspace(ymin, ymax, cbar_num)

        # config = dict(vmin=ymin, vmax=ymax, cmap="coolwarm",)
        config = dict(
            vmin=cbar_ticks.min(),
            vmax=cbar_ticks.max(),
            cmap="coolwarm",
        )
        if fig_type != "pcolor":
            config.update(levels=cbar_ticks, extend="both")

        if fig_type == "contour":
            subplot = ax.contour(t, x, Y, **config)
        elif fig_type == "contourf":
            subplot = ax.contourf(t, x, Y, **config)
        elif fig_type == "pcolor":
            subplot = ax.pcolor(t, x, Y, **config)
        else:
            raise ValueError("fig_type is not known.")
        ax.set_title(name)
        ax.set_xlabel("t")
        # ax.set_xlim(t.min(), t.max())
        # ax.set_xticks(xticks)
        ax.set_ylabel("x")
        ax.set_ylim(*ylim)
        ax.set_yticks(yticks)
        ax.tick_params(top=True, right=True)

        plt.colorbar(subplot, ax=ax, ticks=cbar_ticks)

        return plt.show()

    def compare_draw_pressures(
        self,
        X,
        Y,
        name,
        X_h,
        Y_h,
        name_h,
        fig_type="contourf",
        error_type="abs",
        Y_range=None,
    ):
        """_summary_

        Parameters
        ----------
        X : ndarray
            2D-array where rows are time steps and columns are cells centers.
        Y : ndarray
            output based on FDM.
        Y_h : ndarray
            output based on PINNs.
        fig_type : str, optional
            string from a list ['contourf', 'contour', 'pcolor'], by default
            'contourf'.

        Raises
        ------
        ValueError
            fig_type is not known.
        """
        fig, axs = plt.subplots(1, 3, figsize=(18, 3))
        plt.subplots_adjust(hspace=0.5)

        if error_type.lower() in ["d", "dif", "diff", "difference"]:
            error = Y - Y_h
            error_str = f"Error = ({name}-{name_h})"
        elif error_type.lower() in ["a", "abs", "absolute"]:
            error = np.abs(Y - Y_h)
            error_str = f"Error = ABS({name}-{name_h})"
        elif error_type.lower() in ["r", "rel", "relative"]:
            error = np.abs(Y - Y_h) / np.abs(Y)
            error_str = f"Error = ABS({name}-{name_h})/{name}"
        else:
            raise ValueError("Unknown error_type argument.")

        names = [name, name_h, error_str]
        X_list = [X, X_h, X]
        Y_list = [Y, Y_h, error]
        fig_types = [fig_type, fig_type, fig_type]
        Y_range_list = [Y_range, Y_range, None]

        for i in range(3):
            self.draw_pressures(
                X_list[i],
                Y_list[i],
                names[i],
                axs[i],
                fig_types[i],
                11,
                Y_range_list[i],
            )

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

    def __get_lims_ticks(self, x, Y):

        xmin = x.min()
        xmax = x.max()
        xmax_ = round(xmax, 1)
        xstep = (xmax_ - xmin) / 4
        xmin_ = xmin - xstep * 0.25
        xlim = (xmin_, xmax * 1.1)
        xticks = np.linspace(xmin, xmax_, 5)

        ymin = Y.min()
        ymax = Y.max()
        n = 0 if ymax < 2 else -3
        ymax_ = round(ymax, n)
        ystep = (ymax_ - ymin) / 4
        ymin_ = ymin - ystep * 0.25
        ylim = (ymin_, ymax_)
        if ymax_ > ymax:
            ystep = ystep if ymax < 2 else 0
        m = 5 if ystep == 0 else 6
        yticks = np.linspace(ymin, ymax_ + ystep, m)

        return (xlim, xticks), (ylim, yticks)

    # def plot(self, id=None):
    #     fig, axs = plt.subplots(
    #         self.nrows, self.ncols, figsize=(10, 6), sharey=True, sharex=True
    #     )
    #     plt.subplots_adjust(hspace=0.3, wspace=0.2)
    #     alpha = 1
    #     tsteps = [0, 1, 5, 10, 40, 50, 60, 90, 100]
    #     assert len(tsteps) == self.nrows * self.ncols, "tsteps are not compatible"
    #     self.__set_axis_labels(axs)

    #     if id is None:
    #         names = list(self.Data.keys())
    #         N = len(names)
    #     else:
    #         if isinstance(id, list):
    #             lst = list(self.Data.keys())
    #             names = [lst[i] for i in id]
    #             N = len(names)
    #         else:
    #             names = [list(self.Data.keys())[id]]
    #             N = 1

    #     for n in range(N):
    #         name = names[n]
    #         lstyle = "-" if "ana" in name.lower() else "--"
    #         X = self.Data[name][0]
    #         Y = self.Data[name][1]
    #         t = X[:, 0, 0]
    #         x = X[0, :, 1]

    #         for i, ax in enumerate(axs.ravel()):
    #             tstep = tsteps[i]
    #             ax.plot(x, Y[tstep, :], label=name, linestyle=lstyle, alpha=alpha)
    #             ax.tick_params(top=True, right=True)
    #             ax.set_title(f"t={t[tstep]}")
    #             ax.grid(True)
    #             if i == 0 and n == 0:
    #                 (xlim, xticks), (ylim, yticks) = self.__get_lims_ticks(x, Y)
    #                 ax.set_xlim(xlim)
    #                 ax.set_ylim(ylim)
    #                 ax.set_xticks(xticks)
    #                 ax.set_yticks(yticks)

    #     self.__add_legends(fig, N)
    #     plt.show()

    # -------------------------------------------------------------------------
    # End
    # -------------------------------------------------------------------------


if __name__ == "__main__":
    contour = Contour1D()
    print(contour)
