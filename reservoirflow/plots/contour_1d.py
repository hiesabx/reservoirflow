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
        error=False,
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
        super().__init__(verbose, error)
        # self.nrows = nrows
        # self.ncols = ncols

    def plot_case(
        self,
        name,
        ax=None,
        fig_type="contourf",
        cbar_num=11,
        Y_range=None,
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

        if ax is None:
            fig, ax = plt.subplots(figsize=(6, 3))

        X, Y = self.Data[name]
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

        ylim = (x.min(), x.max())
        yticks = np.linspace(ylim[0], ylim[1], 5)
        cbar_ticks = np.linspace(ymin, ymax, cbar_num)
        # if ymax >= 2:
        #     # cbar_ticks = np.linspace(round(ymin, -2), round(ymax, -2), cbar_num)
        #     cbar_ticks = np.linspace(ymin, round(ymax, -1), cbar_num)
        # else:
        #     cbar_ticks = np.linspace(ymin, round(ymax, 1), cbar_num)

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

    def plot(
        self,
        ylims,
        # X,
        # Y,
        # name,
        # X_h,
        # Y_h,
        # name_h,
        # fig_type="contourf",
        # error_type="abs",
        # Y_range=None,
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

        # if error_type.lower() in ["d", "dif", "diff", "difference"]:
        #     error = Y - Y_h
        #     error_str = f"Error = ({name}-{name_h})"
        # elif error_type.lower() in ["a", "abs", "absolute"]:
        #     error = np.abs(Y - Y_h)
        #     error_str = f"Error = ABS({name}-{name_h})"
        # elif error_type.lower() in ["r", "rel", "relative"]:
        #     error = np.abs(Y - Y_h) / np.abs(Y)
        #     error_str = f"Error = ABS({name}-{name_h})/{name}"
        # else:
        #     raise ValueError("Unknown error_type argument.")

        # names = [name, name_h, error_str]
        # X_list = [X, X_h, X]
        # Y_list = [Y, Y_h, error]
        # fig_types = [fig_type, fig_type, fig_type]
        # Y_range_list = [Y_range, Y_range, None]

        names = list(self.Data.keys())
        N = len(names)

        if N == 1:
            fig, axs = plt.subplots(1, N, figsize=(6, 3))
            axs = [axs]
        else:
            if self.error and N == 2:
                N += 1
                X1, Y1 = self.Data[names[0]]
                X2, Y2 = self.Data[names[1]]
                Y3 = np.abs(Y1 - Y2)  # / np.abs(Y1)
                self.add(X1, Y3, f"Absolute Error (sum={Y3.sum().round(3)})")
                names = list(self.Data.keys())
            fig, axs = plt.subplots(1, N, figsize=(N * 6, N + 1))
        plt.subplots_adjust(hspace=0.5)

        if ylims is not None:
            Y_range = ylims
        else:
            X, Y = self.Data[names[0]]
            Y_range = (Y.min(), Y.max())
        for n in range(N):
            if "error" in names[n].lower():
                Y_range = None
            self.plot_case(
                names[n],
                axs[n],
                Y_range=Y_range,
            )

        return plt.show()

    # -------------------------------------------------------------------------
    # End
    # -------------------------------------------------------------------------


if __name__ == "__main__":
    contour = Contour1D()
    print(contour)
