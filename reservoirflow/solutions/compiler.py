STYPES = {
    "numerical": ["n", "num", "numerical"],
    "analytical": ["a", "ana", "analytical"],
    "neurical": ["nn", "neu", "neurical"],
}

METHODS = {
    "numerical": {
        "FdmSolution": [
            "fdm",
            "fdmsolution",
            "fdm_solution",
            "finite difference method",
            "finite-difference-method",
        ],
        "FvmSolution": [
            "fvm",
            "fvmsolution",
            "fvm_solution",
            "finite volume method",
            "finite-volume-method",
        ],
        "FemSolution": [
            "fem",
            "femsolution",
            "fem_solution",
            "finite element method",
            "finite-element-method",
        ],
    },
    "analytical": {
        "D1P1Solution": [
            "1d1p",
            "d1p1solution",
            "d1p1_solution",
            "1-dimension-1-phase",
            "1dimension1phase",
        ],
        "D1P2Solution": [
            "1d2p",
            "d1p2solution",
            "d1p2_solution",
            "1-dimension-2-phase",
            "1dimension2phase",
        ],
        "D1P3Solution": [
            "1d3p",
            "d1p3solution",
            "d1p3_solution",
            "1-dimension-3-phase",
            "1dimension3phase",
        ],
        "D2P1Solution": [
            "2d1p",
            "d2p1solution",
            "d2p1_solution",
            "2-dimension-1-phase",
            "2dimension1phase",
        ],
        "D2P2Solution": [
            "2d1p",
            "d2p2solution",
            "d2p2_solution",
            "2-dimension-2-phase",
            "2dimension2phase",
        ],
        "D2P3Solution": [
            "2d1p",
            "d2p3solution",
            "d2p3_solution",
            "2-dimension-3-phase",
            "2dimension3phase",
        ],
        "D3P1Solution": [
            "3d1p",
            "d3p1solution",
            "d3p1_solution",
            "3-dimension-1-phase",
            "3dimension1phase",
        ],
        "D3P2Solution": [
            "3d1p",
            "d3p2solution",
            "d3p2_solution",
            "3-dimension-2-phase",
            "3dimension2phase",
        ],
        "D3P3Solution": [
            "3d1p",
            "d3p3solution",
            "d3p3_solution",
            "3-dimension-3-phase",
            "3dimension3phase",
        ],
    },
    "neurical": {
        "PinnSolution": [
            "pinn",
            "pinnsolution",
            "pinn_solution",
            "physics informed neural Network Solution",
            "physics-informed-neural-Network-Solution",
        ],
        "DeepONetSolution": [
            "deeponet",
            "deeponetsolution",
            "deeponet_solution",
            "deep operator Network Solution",
            "deep-operator-Network-Solution",
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
        sparse: bool = True,
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

            - numerical methods: ``['FdmSolution', 'FvmSolution', 'FemSolution']``.
            - analytical methods: ``['1D1P', '1D2P', etc.]``.
            - neurical methods: ``['PinnSolution', 'DeepONetSolution']``.

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
        if self.stype == "Numerical":
            if self.method == "FdmSolution":
                from reservoirflow.solutions.numerical.fdm_solution import FdmSolution

                self.model.solution = FdmSolution(self.model, sparse)
            elif self.method == "FvmSolution":
                from reservoirflow.solutions.numerical.fvm_solution import FvmSolution

                self.model.solution = FvmSolution(self.model, sparse)
            elif self.method == "FemSolution":
                from reservoirflow.solutions.numerical.fem_solution import FemSolution

                self.model.solution = FemSolution(self.model, sparse)
            else:
                print("[INFO] Numerical solution could not be assigned.")
                raise ValueError("Not ready.")
        elif self.stype == "Analytical":
            D = self.model.grid.D
            if "D1" in self.method:
                assert D == 1, "D1 solution is only available for 1D grid."
            elif "D2" in self.method:
                assert D == 2, "D2 solution is only available for 2D grid."
            elif "D3" in self.method:
                assert D == 3, "D3 solution is only available for 3D grid."
            else:
                raise ValueError(
                    "Compiler can not use Grid D={D} for method={self.method}."
                )

            if self.method == "D1P1Solution" and D == 1:
                from reservoirflow.solutions.analytical.d1p1_solution import (
                    D1P1Solution,
                )

                self.model.solution = D1P1Solution(self.model, sparse)
            elif self.method == "D1P2Solution" and D == 1:
                from reservoirflow.solutions.analytical.d1p2_solution import (
                    D1P2Solution,
                )

                self.model.solution = D1P2Solution(self.model, sparse)
            elif self.method == "D1P3Solution" and D == 1:
                from reservoirflow.solutions.analytical.d1p3_solution import (
                    D1P3Solution,
                )

                self.model.solution = D1P3Solution(self.model, sparse)
            elif self.method == "D2P1Solution" and D == 2:
                from reservoirflow.solutions.analytical.d2p1_solution import (
                    D2P1Solution,
                )

                self.model.solution = D2P1Solution(self.model, sparse)
            elif self.method == "D2P2Solution" and D == 2:
                from reservoirflow.solutions.analytical.d2p2_solution import (
                    D2P2Solution,
                )

                self.model.solution = D2P2Solution(self.model, sparse)
            elif self.method == "D2P3Solution" and D == 2:
                from reservoirflow.solutions.analytical.d2p3_solution import (
                    D2P3Solution,
                )

                self.model.solution = D2P3Solution(self.model, sparse)
            elif self.method == "D3P1Solution" and D == 3:
                from reservoirflow.solutions.analytical.d3p1_solution import (
                    D3P1Solution,
                )

                self.model.solution = D3P1Solution(self.model, sparse)
            elif self.method == "D3P2Solution" and D == 3:
                from reservoirflow.solutions.analytical.d3p2_solution import (
                    D3P2Solution,
                )

                self.model.solution = D3P2Solution(self.model, sparse)
            elif self.method == "D3P3Solution" and D == 3:
                from reservoirflow.solutions.analytical.d3p3_solution import (
                    D3P3Solution,
                )

                self.model.solution = D3P3Solution(self.model, sparse)
            else:
                print("[INFO] Analytical solution could not be assigned.")
                raise ValueError("Not ready.")
        elif self.stype == "Neurical":
            if self.method == "PinnSolution":
                from reservoirflow.solutions.neurical.pinn_solution import PinnSolution

                self.model.solution = PinnSolution(self.model, sparse)
            elif self.method == "DeepONetSolution":
                from reservoirflow.solutions.neurical.deeponet_solution import (
                    DeepONetSolution,
                )

                self.model.solution = DeepONetSolution(self.model, sparse)
            else:
                print("[INFO] Neurical solution could not be assigned.")
                raise ValueError("Not ready.")
        else:
            print("[INFO] Solution could not be assigned.")
            raise ValueError("Not ready.")

        self.model.solution.stype = self.stype
        self.model.solution.method = self.method
        self.model.solution.sparse = sparse

        print(f"[info] {self.method} was assigned as model.solution.")

    def __set_stype(self, stype):
        if stype.lower() in STYPES["numerical"]:
            self.stype = "Numerical"
        elif stype.lower() in STYPES["analytical"]:
            self.stype = "Analytical"
        elif stype.lower() in STYPES["neurical"]:
            self.stype = "Neurical"
        else:
            raise ValueError(
                "Unknown value in stype argument. "
                + f"Value must be in {list(STYPES.keys())}."
            )

    def __set_method(self, method):
        if self.stype == "Numerical":
            if method.lower() in METHODS["numerical"]["FdmSolution"]:
                self.method = "FdmSolution"
            elif method.lower() in METHODS["numerical"]["FvmSolution"]:
                self.method = "FvmSolution"
            elif method.lower() in METHODS["numerical"]["FemSolution"]:
                self.method = "FemSolution"
            else:
                raise ValueError(
                    f"Unknown value ('{method}') in method argument for stype={self.stype}. "
                    + f"Value must be in {list(METHODS['numerical'].keys())}."
                )
        elif self.stype == "Analytical":
            if method.lower() in METHODS["analytical"]["D1P1Solution"]:
                self.method = "D1P1Solution"
            elif method.lower() in METHODS["analytical"]["D1P2Solution"]:
                self.method = "D1P2Solution"
            elif method.lower() in METHODS["analytical"]["D1P3Solution"]:
                self.method = "D1P3Solution"
            elif method.lower() in METHODS["analytical"]["D2P1Solution"]:
                self.method = "D2P1Solution"
            elif method.lower() in METHODS["analytical"]["D2P2Solution"]:
                self.method = "D2P2Solution"
            elif method.lower() in METHODS["analytical"]["D2P3Solution"]:
                self.method = "D2P3Solution"
            elif method.lower() in METHODS["analytical"]["D3P1Solution"]:
                self.method = "D3P1Solution"
            elif method.lower() in METHODS["analytical"]["D3P2Solution"]:
                self.method = "D3P2Solution"
            elif method.lower() in METHODS["analytical"]["D3P3Solution"]:
                self.method = "D3P3Solution"
            else:
                raise ValueError(
                    f"Unknown value in method argument for stype={self.stype}. "
                    + f"Value must be in {list(METHODS['analytical'].keys())}."
                )
        elif self.stype == "Neurical":
            if method.lower() in METHODS["neurical"]["PinnSolution"]:
                self.method = "PinnSolution"
            elif method.lower() in METHODS["neurical"]["DeepONetSolution"]:
                self.method = "DeepONetSolution"
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
            "rf.solutions.Compiler("
            + f"model='{self.model.name}'"
            + f", stype='{self.stype}'"
            + f", method='{self.method}'"
            + f", sparse={self.sparse}"
            # + f", solution='{self.model.solution.name}'"
            + ")"
        )

    def __str__(self):
        return (
            f"{self.model.name}"
            + f"-{self.stype}"
            + f"-{self.method}"
            + f"-{self.sparse}"
        )


if __name__ == "__main__":
    compiler = Compiler(stype="numerical", method="FdmSolution")
    print(compiler)
    compiler.fit()
    compiler = Compiler(stype="neurical", method="PinnSolution")
    print(compiler)
    compiler.fit()
