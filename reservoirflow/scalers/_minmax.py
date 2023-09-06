import numpy as np

from reservoirflow.scalers._scaler import _Scaler


class MinMax(_Scaler):
    name = "Min Max Scaler"

    def __init__(self, output_range, input_range=None):
        self.Vmin = output_range[0]
        self.Vmax = output_range[1]
        if input_range is not None:
            self.vmin = input_range[0]
            self.vmax = input_range[1]
        else:
            self.vmin = None
            self.vmax = None

    def set_output_range(self, output_range):
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
        self.vmin = v.min(axis=axis)
        self.vmax = v.max(axis=axis)
        return self

    def transform(self, v):
        self.__check_vmin_vmax__()
        vbar = (self.Vmax - self.Vmin) * (v - self.vmin) / (
            self.vmax - self.vmin
        ) + self.Vmin
        return vbar

    scale = transform

    def inverse_transform(self, vbar):
        self.__check_vmin_vmax__()
        v = (self.vmax - self.vmin) * (vbar - self.Vmin) / (
            self.Vmax - self.Vmin
        ) + self.vmin
        return v

    descale = inverse_transform

    def fit_transform(self, v, axis=0):
        self.fit(v, axis)
        return self.transform(v)

    def __check_vmin_vmax__(self):
        if self.vmin is None or self.vmax is None:
            msg = (
                "input_range=[vmin,vmax] is not defined.\n"
                + "Use fit (or fit_transform) or define input_range in initialization."
            )
            raise ValueError(msg)


if __name__ == "__main__":
    scaler = MinMax(output_range=(0, 1))
    print(scaler)
