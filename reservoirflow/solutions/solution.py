from abc import ABC, abstractmethod


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
        model,
    ):
        """Construct solution object.

        Parameters
        ----------
        model : Model
            a model object from ``models`` module.
        """
        self.model = model

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
