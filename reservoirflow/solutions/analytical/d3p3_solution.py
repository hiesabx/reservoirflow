from reservoirflow.solutions.solution import Solution


class D3P3Solution(Solution):
    """D3P3Solution solution class.

    D3P3Solution is a 3-Dimension-3-Phase.

    .. caution::
        This class is not available.

    Returns
    -------
    Solution
        Solution object.
    """

    name = "D3P3Solution"

    def __init__(self, **kwargs):
        raise NotImplementedError("This class is not implemented.")

    def solve(self):
        raise NotImplementedError

    def run(self):
        raise NotImplementedError
