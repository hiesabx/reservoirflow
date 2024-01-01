from reservoirflow.solutions.solution import Solution


class PINN(Solution):
    """PINN solution class.

    PINN is a Physics-Informed-Neural-Network.

    .. caution::
        This class is not available.

    Returns
    -------
    Solution
        Solution object.
    """

    name = "PINN"

    def __init__(self, **kwargs):
        raise NotImplementedError("This class is not implemented.")

    def solve(self):
        raise NotImplementedError

    def run(self):
        raise NotImplementedError
