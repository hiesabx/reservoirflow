from reservoirflow.solutions.solution import Solution


class FemSolution(Solution):
    """FemSolution solution class.

    FemSolution is a Finite-Element-Method.

    .. caution::
        This class is not available.

    Returns
    -------
    Solution
        Solution object.
    """

    name = "FemSolution"

    def __init__(self, **kwargs):
        raise NotImplementedError("This class is not implemented.")

    def solve(self):
        raise NotImplementedError

    def run(self):
        raise NotImplementedError
