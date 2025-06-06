from reservoirflow.base import Base


class Plot(Base):
    name = "Plot"

    def __init__(self, verbose):
        self.verbose = verbose
        self.Data = {}

    def add(self, x, y, name):
        if self.report:
            if name in self.Data.keys():
                print(f"[Info] Solution: {name} was updated.")
            else:
                print(f"[Info] Solution: {name} was added.")
        self.Data[name] = [x, y]
        return self

    # def add_solution(self, solution, name=None):
    #     """Add solution data to the plot.

    #     Parameters
    #     ----------
    #     solution : tuple
    #         A tuple containing (x, y) data.
    #     name : str, optional
    #         Name of the solution, by default None.
    #     """
    #     if name is None:
    #         name = solution.name
    #     x, y = solution
    #     self.add(x, y, name)
    #     return self


if __name__ == "__main__":
    plot = Plot()
    print(plot)
