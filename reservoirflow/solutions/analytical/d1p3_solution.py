from reservoirflow.solutions.solution import Solution


class D1P3Solution(Solution):
    """D1P3Solution solution class.

    D1P3Solution is a 1-Dimension-3-Phase.

    .. caution::
        This class is not available.

    Returns
    -------
    Solution
        Solution object.
    """

    name = "D1P3Solution"

    def __init__(self, **kwargs):
        raise NotImplementedError("This class is not implemented.")

    def solve(self):
        raise NotImplementedError

    def run(self):
        raise NotImplementedError
