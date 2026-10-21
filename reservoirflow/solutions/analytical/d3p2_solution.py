from reservoirflow.solutions.solution import Solution


class D3P2Solution(Solution):
    """D3P2Solution solution class.

    D3P2Solution is a 3-Dimension-2-Phase.

    .. caution::
        This class is not available.

    Returns
    -------
    Solution
        Solution object.
    """

    name = "D3P2Solution"

    def __init__(self, **kwargs):
        raise NotImplementedError("This class is not implemented.")

    def solve(self):
        raise NotImplementedError

    def run(self):
        raise NotImplementedError
