"""
MultiPhase
==========
"""

from reservoirflow.fluids.fluid import Fluid


class MultiPhase(Fluid):
    """MultiPhase fluid class.

    .. caution::
        This class is not available.

    Returns
    -------
    Fluid
        Fluid object.
    """

    name = "MultiPhase"

    def __init__(self, **kwargs):
        raise NotImplementedError("This class is not implemented.")


if __name__ == "__main__":
    pass
