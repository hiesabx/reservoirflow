from reservoirflow.solutions.solution import Solution
from reservoirflow.models.model import Model


class Network(Solution):
    """PINN solution class.

    PINN is a Physics-Informed-Neural-Network.

    .. caution::
        This class is not available.

    Returns
    -------
    Solution
        Solution object.
    """

    name = "PINN"

    def __init__(
        self,
        model: Model,  #: rf.models.Model,
        sparse: bool = True,
    ):
        """Create Network Solution.

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
