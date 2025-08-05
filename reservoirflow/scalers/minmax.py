"""
MinMax
------
"""

import numpy as np
from reservoirflow.scalers.scaler import Scaler


class MinMax(Scaler):
    """MinMax scaler class.

    This scaler is used to scale input data based on
    ``output_range=(min,max)``. If ``input_range`` is set to ``None``
    instead of ``input_range=(min,max)``, then ``input_range`` is
    inferred based on input data.

    .. hint::

        Using ``input_range=(min,max)`` is useful in some cases to match
        the scaling with other solutions when ``input_range`` can't be
        inferred from input data (e.g. unstable solution).

    .. note::

        Note that if the input array has multiple feature each with its
        own range (not unified), then using ``input_range=None`` is
        required to infer ``input_range`` for each feature.

    Returns
    -------
    Scaler
        Scaler object.
    """

    name = "MinMax"

    def __init__(
        self,
        output_range: tuple,
        input_range: tuple | None = None,
    ):
        self.Vmin = output_range[0]
        self.Vmax = output_range[1]
        if input_range is not None:
            self.vmin = input_range[0]
            self.vmax = input_range[1]
        else:
            self.vmin = None
            self.vmax = None

    def set_input_range(self, input_range: tuple):
        self.vmin = input_range[0]
        self.vmax = input_range[1]
        return self

    @property
    def input_range(self):
        if self.vmin is None or self.vmax is None:
            raise ValueError("Input range is not set.")
        return (self.vmin, self.vmax)

    @property
    def output_range(self):
        return (self.Vmin, self.Vmax)

    def get_input_range(self):
        return (self.vmin, self.vmax)

    def set_output_range(self, output_range: tuple):
        self.Vmin = output_range[0]
        self.Vmax = output_range[1]
        return self

    def get_output_range(self):
        return (self.Vmin, self.Vmax)

    def get_factors(self):
        return (self.vmax - self.vmin) / (self.Vmax - self.Vmin)

    def fit(self, v, axis=0):
        if len(v.shape) > 2 and axis == 0:
            msg = (
                "axis=0 is not allowed with input len(shape) > 2. "
                + "Use axis=None instead. "
                + "Note that in this case overall min and max are used for scaling."
            )
            raise ValueError(msg)
        self.vmin = np.nanmin(v, axis=axis)  # v.min(axis=axis)
        self.vmax = np.nanmax(v, axis=axis)  # v.nanmax(axis=axis)
        if (self.vmin.ndim == 0 and self.vmin == self.vmax) or (
            self.vmin.ndim > 0 and all(self.vmin == self.vmax)
        ):
            self.vmax = 1 + self.vmin
        return self

    def transform(self, v):
        self.__check_vmin_vmax__()
        vbar = (self.Vmax - self.Vmin) * (v - self.vmin) / (
            self.vmax - self.vmin
        ) + self.Vmin
        return vbar  #: transformed input values.

    def inverse_transform(self, vbar):
        self.__check_vmin_vmax__()
        v = (self.vmax - self.vmin) * (vbar - self.Vmin) / (
            self.Vmax - self.Vmin
        ) + self.vmin
        return v  #: inverse_transformed values (back to original).


if __name__ == "__main__":
    scaler = MinMax(output_range=(0, 1))
    print(scaler)
