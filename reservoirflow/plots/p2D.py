from reservoirflow.plots.plot import Plot


class p2D(Plot):

    name = "2D Plot"

    def __init__(
        self,
        dtype="double",
        unit="field",
        verbose=False,
    ):
        super().__init__()

    # -------------------------------------------------------------------------
    # End
    # -------------------------------------------------------------------------


if __name__ == "__main__":
    plot = p2D()
    print(plot)