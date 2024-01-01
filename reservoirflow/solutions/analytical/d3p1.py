from reservoirflow.solutions.solution import Solution


class D3P1(Solution):
    """D3P1 solution class.

    D3P1 is a 3-Dimension-1-Phase.

    .. caution::
        This class is not available.

    Returns
    -------
    Solution
        Solution object.
    """

    name = "D3P1"

    def __init__(self, **kwargs):
        raise NotImplementedError("This class is not implemented.")

    def solve(self):
        raise NotImplementedError

    def run(self):
        raise NotImplementedError
