"""
MinMax
------
"""

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

    def set_output_range(self, output_range: tuple):
        self.Vmin = output_range[0]
        self.Vmax = output_range[1]
        return self

    def fit(self, v, axis=0):
        if len(v.shape) > 2 and axis == 0:
            msg = (
                "axis=0 is not allowed with input len(shape) > 2. "
                + "Use axis=None instead. "
                + "Note that in this case overall min and max are used for scaling."
            )
            raise ValueError(msg)
        self.vmin = v.min(axis=axis)  #: input minimum value.
        self.vmax = v.max(axis=axis)  #: input maximum value.
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

    def fit_transform(self, v, axis=0):
        self.fit(v, axis)
        return self.transform(v)

    def __check_vmin_vmax__(self):
        if self.vmin is None or self.vmax is None:
            msg = (
                "input_range=[vmin,vmax] is not defined. "
                + "Use fit (or fit_transform) or define input_range"
                + " in initialization."
            )
            raise ValueError(msg)


if __name__ == "__main__":
    scaler = MinMax(output_range=(0, 1))
    print(scaler)
