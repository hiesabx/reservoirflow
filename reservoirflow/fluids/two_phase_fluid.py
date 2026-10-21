"""
TwoPhase
========
"""

from .fluid import Fluid


class TwoPhaseFluid(Fluid):
    """TwoPhase fluid class.

    .. caution::
        This class is not available.

    Returns
    -------
    Fluid
        Fluid object.
    """

    name = "TwoPhaseFluid"

    def __init__(self, **kwargs):
        raise NotImplementedError("This class is not implemented.")


if __name__ == "__main__":
    pass
