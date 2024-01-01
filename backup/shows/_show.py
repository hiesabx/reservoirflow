from reservoirflow.base import Base


class _Show(Base):
    name = "Visual"

    def __init__(self):
        pass


if __name__ == "__main__":
    visual = Visual()
    print(visual)
