from reservoirflow.solutions.solution import Solution


class FVM(Solution):
    """FVM solution class.

    FVM is a Finite-Volume-Method.

    .. caution::
        This class is not available.

    Returns
    -------
    Solution
        Solution object.
    """

    name = "FVM"

    def __init__(self, **kwargs):
        raise NotImplementedError("This class is not implemented.")

    def solve(self):
        raise NotImplementedError

    def run(self):
        raise NotImplementedError
