"""
1D Plot
=======
"""

import numpy as np
import matplotlib.pyplot as plt

from reservoirflow.plots.plot import Plot


class Plot1D(Plot):
    """1D Line Plot class for visualizing 1D data.

    Parameters
    ----------
    Plot : _type_
        _description_
    """

    name = "1D Line Plot"

    def __init__(
        self,
        nrows=1,
        ncols=1,
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
        self.nrows = nrows
        self.ncols = ncols

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

    def plot(self, id=None):
        fig, axs = plt.subplots(
            self.nrows,
            self.ncols,
            figsize=(10, 6),
            sharey=True,
            sharex=True,
        )
        plt.subplots_adjust(hspace=0.3, wspace=0.2)
        alpha = 1
        tsteps = [0, 1, 5, 10, 40, 50, 60, 90, 100]

        assert len(tsteps) == self.nrows * self.ncols, "tsteps are not compatible"

        self.__set_axis_labels(axs)

        if id is None:
            names = list(self.Data.keys())
            N = len(names)
        else:
            if isinstance(id, list):
                lst = list(self.Data.keys())
                names = [lst[i] for i in id]
                N = len(names)
            else:
                names = [list(self.Data.keys())[id]]
                N = 1

        for n in range(N):
            name = names[n]
            lstyle = "-" if "ana" in name.lower() else "--"
            X = self.Data[name][0]
            Y = self.Data[name][1]
            t = X[:, 0, 0]
            x = X[0, :, 1]

            for i, ax in enumerate(axs.ravel()):
                tstep = tsteps[i]
                ax.plot(x, Y[tstep, :], label=name, linestyle=lstyle, alpha=alpha)
                ax.tick_params(top=True, right=True)
                ax.set_title(f"t={t[tstep]}")
                ax.grid(True)
                if i == 0 and n == 0:
                    (xlim, xticks), (ylim, yticks) = self.__get_lims_ticks(x, Y)
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
    plot = Plot1D()
    print(plot)
