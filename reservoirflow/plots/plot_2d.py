"""
2D
==
"""

from reservoirflow.plots.plot import Plot


class Plot2D(Plot):
    """_summary_

    Parameters
    ----------
    Plot : _type_
        _description_
    """

    name = "2D Plot"

    def __init__(
        self,
        dtype="double",
        unit="field",
        verbose=False,
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
        super().__init__()

    # -------------------------------------------------------------------------
    # End
    # -------------------------------------------------------------------------


if __name__ == "__main__":
    plot = Plot2D()
    print(plot)
