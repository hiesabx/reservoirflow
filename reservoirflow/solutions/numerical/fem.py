from reservoirflow.solutions.solution import Solution


class FEM(Solution):
    """FEM solution class.

    FEM is a Finite-Element-Method.

    .. caution::
        This class is not available.

    Returns
    -------
    Solution
        Solution object.
    """

    name = "FEM"

    def __init__(self, **kwargs):
        raise NotImplementedError("This class is not implemented.")

    def solve(self):
        raise NotImplementedError

    def run(self):
        raise NotImplementedError
