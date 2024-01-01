from abc import ABC, abstractmethod

from reservoirflow.base import Base


class Well(ABC, Base):
    """Abstract well class.

    .. attention::

        This is an abstract class and can't be instantiated. This class
        is only used as a parent for other classes of ``wells`` module.

    Returns
    -------
    Well
        Well object.
    """

    name = "Well"

    def __init__(self):
        pass


if __name__ == "__main__":
    well = Well()
    print(well)
