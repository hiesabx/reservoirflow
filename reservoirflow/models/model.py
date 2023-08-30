from reservoirflow.base import Base


class Model(Base):
    name = "Model"

    def __init__(self, unit, dtype, verbose):
        super().__init__(unit, dtype, verbose)


if __name__ == "__main__":
    dtype = "double"
    unit = "field"
    verbose = False
    model = Model(unit, dtype, verbose)
    print(model)
