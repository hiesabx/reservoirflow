import numpy as np
from reservoirflow.base import Base


class Scaler(Base):
    name = "Scaler"

    def __init__(self, unit, dtype, verbose):
        super().__init__(unit, dtype, verbose)


if __name__ == "__main__":
    dtype = "double"
    unit = "field"
    verbose = False
    scaler = Scaler(unit, dtype, verbose)
    print(scaler)
