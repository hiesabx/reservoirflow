from reservoirflow.wells.well import Well


class DirectionalWell(Well):
    """DirectionalWell well class.

    .. caution::
        This class is not available.

    Returns
    -------
    Well
        Well object.
    """

    name = "DirectionalWell"

    def __init__(self, **kwargs):
        raise NotImplementedError("This class is not implemented.")
