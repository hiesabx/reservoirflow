from reservoirflow.solutions.solution import Solution


class D2P1Solution(Solution):
    """D2P1Solution solution class.

    D2P1Solution is a 2-Dimension-1-Phase.

    .. caution::
        This class is not available.

    Returns
    -------
    Solution
        Solution object.
    """

    name = "D2P1Solution"

    def __init__(self, **kwargs):
        raise NotImplementedError("This class is not implemented.")

    def solve(self):
        raise NotImplementedError

    def run(self):
        raise NotImplementedError
