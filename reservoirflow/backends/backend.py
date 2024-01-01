from abc import ABC, abstractmethod


class Backend(ABC):
    """Abstract backend class.

    .. attention::

        This is an abstract class and can't be instantiated. This class
        is only used as a parent for other classes of ``backends``
        module.

    Returns
    -------
    Backend
        Backend object.
    """

    def __init__(
        self,
        name,
    ):
        """Construct backend object.

        Parameters
        ----------
        name : str
            backend name.
        """
        self.name = name

    @abstractmethod
    def transpose(self):
        """Transpose."""

    @abstractmethod
    def ones(self):
        """Ones."""


if __name__ == "__main__":
    backend = Backend("NumPy")
