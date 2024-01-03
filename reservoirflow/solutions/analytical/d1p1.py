from reservoirflow.solutions.solution import Solution


class D1P1(Solution):
    """D1P1 solution class.

    D1P1 is a 1-Dimension-1-Phase.

    .. caution::
        This class is not available.

    Returns
    -------
    Solution
        Solution object.
    """

    name = "D1P1"

    def __init__(self, **kwargs):
        raise NotImplementedError("This class is not implemented.")


def get_alpha(self, method="mean"):
    alpha_n = (
        self.model.factors["transmissibility conversion"] * self.model.grid.kx
    ) / (self.model.fluid.mu * self.model.fluid.B)
    alpha_d = (self.model.grid.phi * self.model.comp) / (
        self.model.factors["volume conversion"] * self.model.fluid.B
    )
    alpha = alpha_n / alpha_d

    if method in ["first", None]:  # or np.all(alpha == alpha[0]):
        alpha = alpha[0]
    elif method in ["average", "avg", "mean"]:
        alpha = alpha.mean()
    elif method in ["vector", "array"]:
        pass
    else:  #
        pass

    # print(f"{alpha=}")

    return alpha

    def solve(self):
        raise NotImplementedError

    def run(self):
        raise NotImplementedError
