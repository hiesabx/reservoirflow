from reservoirflow.grids.grid import Grid


class Radial(Grid):
    """Radial grid class.

    .. caution::
        This class is not available.

    Returns
    -------
    Grid
        Grid object.
    """

    name = "Radial"

    def __init__(self, **kwargs):
        raise NotImplementedError("This class is not implemented.")
