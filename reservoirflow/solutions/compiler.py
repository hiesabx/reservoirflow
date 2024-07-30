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
        "d1p1": ["d1p1", "1d1p", "1-dimension-1-phase", "1dimension1phase"],
        "d1p2": ["d1p2", "1d2p", "1-dimension-2-phase", "1dimension2phase"],
        "d1p3": ["d1p3", "1d3p", "1-dimension-3-phase", "1dimension3phase"],
        "d2p1": ["d2p1", "2d1p", "2-dimension-1-phase", "2dimension1phase"],
        "d2p2": ["d2p2", "2d1p", "2-dimension-2-phase", "2dimension2phase"],
        "d2p3": ["d2p3", "2d1p", "2-dimension-3-phase", "2dimension3phase"],
        "d3p1": ["d3p1", "3d1p", "3-dimension-1-phase", "3dimension1phase"],
        "d3p2": ["d3p2", "3d1p", "3-dimension-2-phase", "3dimension2phase"],
        "d3p3": ["d3p3", "3d1p", "3-dimension-3-phase", "3dimension3phase"],
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
        sparse: bool=True,
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
            
            - numerical methods: ``['FDM', 'FVM', 'FEM']``.
            - analytical methods: ``['1D1P', '1D2P', etc.]``.
            - neurical methods: ``['PINN', 'DeepONet']``.
            
        sparse : bool, optional, default: True
            using sparse computing for a better performance.
        """
        self.model = model
        
        if isinstance(sparse, bool):
            self.sparse = sparse
        else:
            raise ValueError("Argument sparse must be bool.")
        
        self.__set_stype(stype)
        self.__set_method(method)
        self.__add_solution(sparse)

    def __add_solution(self, sparse):
        if self.stype == "numerical":
            if self.method == "FDM":
                from reservoirflow.solutions.numerical.fdm import FDM

                self.model.solution = FDM(self.model, sparse)
            elif self.method == "FVM":
                from reservoirflow.solutions.numerical.fvm import FVM

                self.model.solution = FVM(self.model, sparse)
            elif self.method == "FEM":
                from reservoirflow.solutions.numerical.fem import FEM

                self.model.solution = FEM(self.model, sparse)
            else:
                print("[INFO] Numerical solution could not be assigned.")
                raise ValueError("Not ready.")
        elif self.stype == "analytical":
            if self.method == "D1P1":
                from reservoirflow.solutions.analytical.d1p1 import D1P1

                self.model.solution = D1P1(self.model, sparse)
            elif self.method == "D1P2":
                from reservoirflow.solutions.analytical.d1p2 import D1P2

                self.model.solution = D1P2(self.model, sparse)
            elif self.method == "D1P3":
                from reservoirflow.solutions.analytical.d1p3 import D1P3

                self.model.solution = D1P3(self.model, sparse)
            elif self.method == "D2P1":
                from reservoirflow.solutions.analytical.d2p1 import D2P1

                self.model.solution = D2P1(self.model, sparse)
            elif self.method == "D2P2":
                from reservoirflow.solutions.analytical.d2p2 import D2P2

                self.model.solution = D2P2(self.model, sparse)
            elif self.method == "D2P3":
                from reservoirflow.solutions.analytical.d2p3 import D2P3

                self.model.solution = D2P3(self.model, sparse)
            elif self.method == "D3P1":
                from reservoirflow.solutions.analytical.d3p1 import D3P1

                self.model.solution = D3P1(self.model, sparse)
            elif self.method == "D3P2":
                from reservoirflow.solutions.analytical.d3p2 import D3P2

                self.model.solution = D3P2(self.model, sparse)
            elif self.method == "D3P3":
                from reservoirflow.solutions.analytical.d3p3 import D3P3

                self.model.solution = D3P3(self.model, sparse)
            else:
                print("[INFO] Analytical solution could not be assigned.")
                raise ValueError("Not ready.")
        elif self.stype == "neurical":
            if self.method == "PINN":
                from reservoirflow.solutions.neurical.pinn import PINN

                self.model.solution = PINN(self.model, sparse)
            elif self.method == "DeepONet":
                from reservoirflow.solutions.neurical.deeponet import DeepONet

                self.model.solution = DeepONet(self.model, sparse)
            else:
                print("[INFO] Neurical solution could not be assigned.")
                raise ValueError("Not ready.")
        else:
            print("[INFO] Solution could not be assigned.")
            raise ValueError("Not ready.")
        # self.model.solutions[self.method] = self.model.solution
        print(f"[info] {self.method} was assigned as model.solution.")

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
            if method.lower() in METHODS["analytical"]["d1p1"]:
                self.method = "D1P1"
            elif method.lower() in METHODS["analytical"]["d1p2"]:
                self.method = "D1P2"
            elif method.lower() in METHODS["analytical"]["d1p3"]:
                self.method = "D1P3"
            elif method.lower() in METHODS["analytical"]["d2p1"]:
                self.method = "D2P1"
            elif method.lower() in METHODS["analytical"]["d2p2"]:
                self.method = "D2P2"
            elif method.lower() in METHODS["analytical"]["d2p3"]:
                self.method = "D2P3"
            elif method.lower() in METHODS["analytical"]["d3p1"]:
                self.method = "D3P1"
            elif method.lower() in METHODS["analytical"]["d3p2"]:
                self.method = "D3P2"
            elif method.lower() in METHODS["analytical"]["d3p3"]:
                self.method = "D3P3"
            else:
                raise ValueError(
                    f"Unknown value in method argument for stype={self.stype}. "
                    + f"Value must be in {list(METHODS['analytical'].keys())}."
                )
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

    def __repr__(self):
        return (
              f"Compiler(model='{self.model.name}'"
            + f", stype='{self.stype}'"
            + f", method='{self.method}'"
            + f", sparse={self.sparse}"
            # + f", solution='{self.model.solution.name}'"
            + ")"
        )


if __name__ == "__main__":
    compiler = Compiler(stype="numerical", method="fdm")
    print(compiler)
    compiler.fit()
    compiler = Compiler(stype="neurical", method="pinn")
    print(compiler)
    compiler.fit()
