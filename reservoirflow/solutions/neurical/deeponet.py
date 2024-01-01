from reservoirflow.solutions.solution import Solution


class DeepONet(Solution):
    """DeepONet solution class.

    DeepONet is a Deep-Operator-Network.

    .. caution::
        This class is not available.

    Returns
    -------
    Solution
        Solution object.
    """

    name = "DeepONet"

    def __init__(self, **kwargs):
        raise NotImplementedError("This class is not implemented.")

    def solve(self):
        raise NotImplementedError

    def run(self):
        raise NotImplementedError
