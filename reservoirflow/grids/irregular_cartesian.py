from reservoirflow.grids.grid import Grid


class IrregularCartesian(Grid):
    """IrregularCartesian grid class.

    .. caution::
        This class is not available.

    Returns
    -------
    Grid
        Grid object.
    """

    name = "IrregularCartesian"

    def __init__(self, **kwargs):
        raise NotImplementedError("This class is not implemented.")
