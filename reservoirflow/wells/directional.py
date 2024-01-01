from reservoirflow.wells.well import Well


class Directional(Well):
    """Directional well class.

    .. caution::
        This class is not available.

    Returns
    -------
    Well
        Well object.
    """

    name = "Directional"

    def __init__(self, **kwargs):
        raise NotImplementedError("This class is not implemented.")
