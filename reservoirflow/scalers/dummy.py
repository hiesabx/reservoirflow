"""
Dummy
=====
"""

from reservoirflow.scalers.scaler import Scaler


class Dummy(Scaler):
    """Dummy scaler class.

    This scaler mimics the scalers' behavior but does not apply any
    transformations (i.e. returns data without scaling).

    Returns
    -------
    Scaler
        Scaler object.
    """

    name = "Dummy"

    def set_input_range(self, input_range):
        return self

    def set_output_range(self, output_range):
        return self

    def fit(self, v, axis=0):
        return self

    def transform(self, v):
        return v

    def inverse_transform(self, vbar):
        return vbar

    def fit_transform(self, v, axis=0):
        self.fit(v, axis)
        return self.transform(v)


if __name__ == "__main__":
    scaler = Dummy(output_range=(0, 1))
    print(scaler)
