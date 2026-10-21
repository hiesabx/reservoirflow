from reservoirflow.wells.well import Well


class MultiCellWell(Well):
    """MultiCellWell well class.

    .. caution::
        This class is not available.

    Returns
    -------
    Well
        Well object.
    """

    name = "MultiCellWell"

    def __init__(self, **kwargs):
        raise NotImplementedError("This class is not implemented.")
