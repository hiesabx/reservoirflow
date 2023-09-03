from reservoirflow.scalers._scaler import Scaler


class Dummy(Scaler):
    name = "Dummy Scaler"

    def __init__(self, output_range=None, input_range=None):
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
