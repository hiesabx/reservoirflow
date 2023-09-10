from reservoirflow._base import _Base


class Well(_Base):
    name = "Well"

    def __init__(self):
        pass


if __name__ == "__main__":
    well = Well()
    print(well)
