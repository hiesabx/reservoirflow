STYPES = {
    "numerical": ["n", "num", "numerical"],
    "analytical": ["a", "ana", "analytical"],
    "neurical": ["nn", "neu", "neurical"],
}

METHODS = {
    "numerical": {
        "fdm": ["fdm", "finite difference method", "finite-difference-method"],
        "fvm": ["fvm", "finite volume method", "finite-volume-method"],
        "fem": ["fem", "finite element method", "finite-element-method"],
    },
    "analytical": {
        "1d1p": ["1d1p", "1-dimension-1-phase", "1dimension1phase"],
        "1d2p": ["1d2p", "1-dimension-2-phase", "1dimension2phase"],
        "1d3p": ["1d3p", "1-dimension-3-phase", "1dimension3phase"],
        "2d1p": ["2d1p", "2-dimension-1-phase", "2dimension1phase"],
        "2d2p": ["2d1p", "2-dimension-2-phase", "2dimension2phase"],
        "2d3p": ["2d1p", "2-dimension-3-phase", "2dimension3phase"],
        "3d1p": ["3d1p", "3-dimension-1-phase", "3dimension1phase"],
        "3d2p": ["3d1p", "3-dimension-2-phase", "3dimension2phase"],
        "3d3p": ["3d1p", "3-dimension-3-phase", "3dimension3phase"],
    },
    "neurical": {
        "pinn": [
            "pinn",
            "physics informed neural network",
            "physics-informed-neural-network",
        ],
        "deeponet": [
            "deeponet",
            "deep operator network",
            "deep-operator-network",
        ],
    },
}

MODES = {
    "vectorized": ["v", "vect", "vectorize", "vectorized"],
    "symbolized": ["s", "symb", "symbolize", "symbolized"],
    "availability": {
        "vectorized": [
            "numerical",
            "analytical",
            "neurical",
        ],
        "symbolized": [
            "numerical",
        ],
    },
}

SOLVERS = {
    "direct": ["d", "direct"],
    "iterative": ["i", "iterative"],
    "neurical": ["n", "neurical"],
    "availability": {
        "direct": [
            "numerical",
            "analytical",
            "neurical",
        ],
        "iterative": [
            "numerical",
        ],
        "neurical": [
            "numerical",
        ],
    },
}


class Compiler:
    """Compiler class.

    This class is used to compile a model from ``models`` module.

    Returns
    -------
    Compiler
        Compiler object.
    """

    def __init__(
        self,
        model,
        stype: str,
        method: str,
        # mode: str,
        # solver: str,
    ):
        """Construct compiler object.

        Parameters
        ----------
        model : Model
            a model object from ``models`` module.
        stype : str
            solution type in ['numerical', 'analytical', 'neurical'].
        method : str
            solution method as following:
            - 'numerical' methods: ['FDM', 'FVM', 'FEM'].
            - 'analytical' methods: ['1D1P', '1D2P', etc.].
            - 'neurical' methods: ['PINN', 'DeepONet'].
        """
        # mode : str
        #     solution mode in ['vectorized', 'symbolized'].
        # solver : str
        #     solution solver in ['direct', 'iterative', 'neurical'].

        self.model = model
        self.__set_stype(stype)
        self.__set_method(method)
        # self.__set_mode(mode)
        # self.__set_solver(solver)
        self.__add_solution()

    def __add_solution(self):
        if self.stype == "numerical":
            if self.method == "FDM":
                from reservoirflow.solutions.numerical.fdm import FDM

                self.model.solution = FDM(self.model)
                print("[info] FDM was assigned as model.solution.")
            elif self.method == "FVM":
                from reservoirflow.solutions.numerical.fvm import FVM

                self.model.solution = FVM(self.model)
                print("[info] FVM was assigned as model.solution.")
            elif self.method == "FEM":
                from reservoirflow.solutions.numerical.fem import FEM

                self.model.solution = FEM(self.model)
                print("[info] FEM was assigned as model.solution.")
            else:
                print("[INFO] Numerical EquationSystem was not constructed.")
                raise ValueError("Not ready.")
        elif self.stype == "analytical":
            raise ValueError("Not ready.")
        elif self.stype == "neurical":
            if self.method == "PINN":
                raise ValueError("Not ready.")
            elif self.method == "DeepONet":
                raise ValueError("Not ready.")
            else:
                print("[INFO] Numerical EquationSystem was not constructed.")
                raise ValueError("Not ready.")
        else:
            print("[INFO] EquationSystem was not constructed.")
            raise ValueError("Not ready.")

    def __set_stype(self, stype):
        if stype.lower() in STYPES["numerical"]:
            self.stype = "numerical"
        elif stype.lower() in STYPES["analytical"]:
            self.stype = "analytical"
        elif stype.lower() in STYPES["neurical"]:
            self.stype = "neurical"
        else:
            raise ValueError(
                "Unknown value in stype argument. "
                + f"Value must be in {list(STYPES.keys())}."
            )

    def __set_method(self, method):
        if self.stype == "numerical":
            if method.lower() in METHODS["numerical"]["fdm"]:
                self.method = "FDM"
            elif method.lower() in METHODS["numerical"]["fvm"]:
                self.method = "FVM"
            elif method.lower() in METHODS["numerical"]["fem"]:
                self.method = "FEM"
            else:
                raise ValueError(
                    f"Unknown value in method argument for stype={self.stype}. "
                    + f"Value must be in {list(METHODS['numerical'].keys())}."
                )
        elif self.stype == "analytical":
            self.method = "nan"
        elif self.stype == "neurical":
            if method.lower() in METHODS["neurical"]["pinn"]:
                self.method = "PINN"
            elif method.lower() in METHODS["neurical"]["deeponet"]:
                self.method = "DeepONet"
            else:
                raise ValueError(
                    f"Unknown value in method argument for stype={self.stype}. "
                    + f"Value must be in {list(METHODS['neurical'].keys())}."
                )
        else:
            raise ValueError(
                "Unkown value in stype argument. "
                + f"Value must be in {list(STYPES.keys())}."
            )

    def __set_mode(self, mode):
        if mode.lower() in MODES["vectorized"]:
            self.mode = "vectorized"
        elif mode.lower() in MODES["symbolized"]:
            self.mode = "symbolized"
        else:
            raise ValueError(
                "Unknown value in mode argument. "
                + f"Value must be in {list(MODES.keys())[:-1]}."
            )
        if self.stype not in MODES["availability"][self.mode]:
            raise ValueError(
                f"Selected mode={self.mode} is not available "
                + f"when stype={self.stype} and method={self.method}. "
                + f"Available option are {MODES['availability']}"
            )

    def __set_solver(self, solver):
        if solver.lower() in SOLVERS["direct"]:
            self.solver = "direct"
        elif solver.lower() in SOLVERS["iterative"]:
            self.solver = "iterative"
        elif solver.lower() in SOLVERS["neurical"]:
            self.solver = "neurical"
        else:
            raise ValueError(
                "Unknown value in solver argument. "
                + f"Value must be in {list(SOLVERS.keys())[:-1]}."
            )
        if self.stype not in SOLVERS["availability"][self.solver]:
            raise ValueError(
                f"Selected solver={self.solver} is not available "
                + f"when stype={self.stype} and method={self.method}. "
                + f"Available option are {SOLVERS['availability']}"
            )

    def __repr__(self):
        return (
            f"Compiler(model='{self.model.name}', stype='{self.stype}', "
            + f"method='{self.method}', solution='{self.model.solution.name}')"
        )
        # return (
        #     f"Compiler(model='{self.model.name}', stype='{self.stype}', "
        #     + f"method='{self.method}', mode='{self.mode}', "
        #     + f"solver='{self.solver}', solution='{self.model.solution.name}')"
        # )


if __name__ == "__main__":
    # compiler = Compiler(stype="numerical", method="fdm", mode="s", solver="n")
    compiler = Compiler(stype="numerical", method="fdm")
    print(compiler)
    compiler.fit()
    compiler = Compiler(stype="neurical", method="pinn")
    print(compiler)
    compiler.fit()
