from reservoirflow.solutions.solution import Solution


class D2P3(Solution):
    """D2P3 solution class.

    D2P3 is a 2-Dimension-3-Phase.

    .. caution::
        This class is not available.

    Returns
    -------
    Solution
        Solution object.
    """

    name = "D2P3"

    def __init__(self, **kwargs):
        raise NotImplementedError("This class is not implemented.")

    def solve(self):
        raise NotImplementedError

    def run(self):
        raise NotImplementedError
