from abc import ABC, abstractmethod
import reservoirflow as rf

class Solution(ABC):
    """Abstract solution class.

    .. attention::

        This is an abstract class and can't be instantiated. This class
        is only used as a parent for other classes of ``solutions``
        submodules.

    Returns
    -------
    Solution
        Solution object.
    """

    def __init__(
        self,
        model,#: rf.models.Model,
        sparse,
    ):
        """Construct solution object.

        Parameters
        ----------
        model : Model
            a model object from ``models`` module.
        sparse : bool
            using sparse computing for a better performance.
        """
        self.model = model
        self.sparse = sparse
        self.pressures, self.rates = self.model.get_init_arrays()
        
        self.nsteps = 1
        self.tstep = 0
        self.ctime = 0

    @abstractmethod
    def solve(self):
        """Solve a single timestep.

        .. attention::
            This is an abstract method.
        """

    @abstractmethod
    def run(self):
        """Solve multiple timesteps.

        .. attention::
            This is an abstract method.
        """
