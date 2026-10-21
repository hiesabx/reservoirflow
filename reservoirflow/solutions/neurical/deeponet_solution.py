# from reservoirflow.solutions.solution import Solution
from .network_solution import NetworkSolution


class DeepONetSolution(NetworkSolution):
    """DeepONetSolution solution class.

    DeepONetSolution is a Deep-Operator-NetworkSolution.

    .. caution::
        This class is not available.

    Returns
    -------
    Solution
        Solution object.
    """

    name = "DeepONetSolution"

    def __init__(self, **kwargs):
        raise NotImplementedError("This class is not implemented.")

    def solve(self):
        raise NotImplementedError

    def run(self):
        raise NotImplementedError
