from reservoirflow.scalers._scaler import _Scaler


class Dummy(_Scaler):
    """Dummy Scaler

    Returns data without scaling.
    """

    name = "Dummy Scaler"

    def __init__(self, output_range: tuple = None, input_range: tuple = None):
        """Dummy Scaler

        Parameters
        ----------
        output_range : tuple, optional
            output range
        input_range : tuple, optional
            input range
        """
        pass

    def set_output_range(self, output_range):
        return self

    def fit(self, v, axis=0):
        return self

    def transform(self, v):
        return v

    scale = transform

    def inverse_transform(self, vbar):
        return vbar

    descale = inverse_transform

    def fit_transform(self, v, axis=0):
        self.fit(v, axis)
        return self.transform(v)


if __name__ == "__main__":
    scaler = Dummy(output_range=(0, 1))
    print(scaler)
