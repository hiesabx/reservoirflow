from reservoirflow._base import _Base


class _Scaler(_Base):
    name = "Scaler"

    def __init__(self, unit, dtype, verbose):
        super().__init__(unit, dtype, verbose)


if __name__ == "__main__":
    dtype = "double"
    unit = "field"
    verbose = False
    scaler = _Scaler(unit, dtype, verbose)
    print(scaler)
