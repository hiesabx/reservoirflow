from abc import ABC, abstractmethod

# from reservoirflow._base import _Base


class Scaler(ABC):
    """Abstract scaler class.

    .. attention::

        This is an abstract class and can't be instantiated. This class
        is only used as a parent for other classes of ``scalers`` module.

    Returns
    -------
    Scaler
        Scaler object.
    """

    name: str = "Scaler"
    """Class name."""

    def __init__(
        self,
        output_range: tuple | None = None,
        input_range: tuple | None = None,
    ):
        """Construct scaler object.

        Parameters
        ----------
        output_range : tuple, optional
            output range used in the transformation.
        input_range : tuple, optional
            input range used in the transformation.
        """

    @abstractmethod
    def set_input_range(self, input_range: tuple):
        """Set input range for the scaler.

        Parameters
        ----------
        input_range : tuple
            input range used in the transformation.
        """

    @abstractmethod
    def set_output_range(self, output_range: tuple):
        """Set output range for the scaler.

        Parameters
        ----------
        output_range : tuple
            output range used in the transformation.
        """

    @abstractmethod
    def fit(self, v, axis: int = 0):
        """Fit scaler with input values.

        Parameters
        ----------
        v : array
            values before transformation to transform.
        axis : int, by default 0
            use ``axis=0`` for vertical (i.e., across rows) operations,
            and use ``axis=1`` for horizontal (i.e., across columns)
            operations. For a table with multiple features as columns,
            using ``axis=0`` is desired and the length of the output is
            equal to the number of features. For a vector of values,
            either vertical or horizontal, ``axis=0`` is required.
        """

    @abstractmethod
    def transform(self, v):
        """Transform input based on output range.

        Parameters
        ----------
        v : array
            values before transformation to transform.
        """

    # scale = transform

    @abstractmethod
    def inverse_transform(self, vbar):
        """Transform input back to the original (input) range.

        Parameters
        ----------
        vbar : array
            values after transformation to inverse.
        """

    # descale = inverse_transform

    @abstractmethod
    def fit_transform(self, v, axis=0):
        """Fit scaler and transform input based on output range.

        Parameters
        ----------
        v : array
            values before transformation.
        axis : int, by default 0
            use ``axis=0`` for vertical (i.e., across rows) operations,
            and use ``axis=1`` for horizontal (i.e., across columns)
            operations. For a table with multiple features as columns,
            using ``axis=0`` is desired and the length of the output is
            equal to the number of features.
        """


if __name__ == "__main__":
    # dtype = "double"
    # unit = "field"
    # verbose = False
    # scaler = Scaler(unit, dtype, verbose)
    # print(scaler)
    scaler = Scaler()
