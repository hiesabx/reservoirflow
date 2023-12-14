"""
2D
==
"""
from reservoirflow.plots._plot import _Plot


class p2D(_Plot):
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
    plot = p2D()
    print(plot)
