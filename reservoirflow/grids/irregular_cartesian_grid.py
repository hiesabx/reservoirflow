from .grid import Grid


class IrregularCartesianGrid(Grid):
    """IrregularCartesian grid class.

    .. caution::
        This class is not available.

    Returns
    -------
    Grid
        Grid object.
    """

    name = "IrregularCartesianGrid"

    def __init__(self, **kwargs):
        raise NotImplementedError("This class is not implemented.")
