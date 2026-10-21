from reservoirflow.models.model import Model
from reservoirflow.solutions.solution import Solution


class NetworkSolution(Solution):
    """PinnSolution solution class.

    PinnSolution is a Physics-Informed-Neural-NetworkSolution.

    .. caution::
        This class is not available.

    Returns
    -------
    Solution
        Solution object.
    """

    name = "PinnSolution"

    def __init__(
        self,
        model: Model,  #: rf.models.Model,
        sparse: bool = True,
    ):
        """Create NetworkSolution Solution.

        Parameters
        ----------
        model : Model
            Model object.
        sparse : bool, optional, default: True
            using sparse computing for a better performance.
        """
        super().__init__(model, sparse)
        self.Data = {}

    def add(self, x, y, name):
        if self.verbose:
            if name in self.Data.keys():
                print(f"[Info] Solution: {name} was updated.")
            else:
                print(f"[Info] Solution: {name} was added.")
        self.Data[name] = [x, y]
        return self

    def train(self):
        raise NotImplementedError

    def solve(self):
        raise NotImplementedError

    def run(self):
        raise NotImplementedError
