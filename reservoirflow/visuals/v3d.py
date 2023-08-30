from reservoirflow.visuals.visual import Visual


class v3D(Visual):
    name = "3D Visual"

    def __init__(
        self,
        dtype="double",
        unit="field",
        verbose=False,
    ):
        super().__init__()

    # -------------------------------------------------------------------------
    # End
    # -------------------------------------------------------------------------


if __name__ == "__main__":
    visual = v3D()
    print(visual)
