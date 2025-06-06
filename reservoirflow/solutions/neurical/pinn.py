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

    def fit(self):
        raise NotImplementedError

    def solve(self):
        raise NotImplementedError

    def run(self):
        raise NotImplementedError
