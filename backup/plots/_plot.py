from reservoirflow.base import Base


class _Plot(Base):
    name = "Plot"

    def __init__(self):
        pass


if __name__ == "__main__":
    plot = _Plot()
    print(plot)
