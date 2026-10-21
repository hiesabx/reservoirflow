from .grid import Grid


class RadialGrid(Grid):
    """Radial grid class.

    .. caution::
        This class is not available.

    Returns
    -------
    Grid
        Grid object.
    """

    name = "RadialGrid"

    def __init__(self, **kwargs):
        raise NotImplementedError("This class is not implemented.")
