from reservoirflow.wells.well import Well


class MultiCell(Well):
    """MultiCell well class.

    .. caution::
        This class is not available.

    Returns
    -------
    Well
        Well object.
    """

    name = "MultiCell"

    def __init__(self, **kwargs):
        raise NotImplementedError("This class is not implemented.")
