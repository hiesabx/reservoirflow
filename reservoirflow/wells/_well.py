from reservoirflow._base import Base


class Well(Base):
    name = "Well"

    def __init__(self):
        pass


if __name__ == "__main__":
    well = Well()
    print(well)
