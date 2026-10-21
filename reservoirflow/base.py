import inspect

import numpy as np
from tabulate import tabulate

UNITS = {
    "field": {
        "transmissibility": "stb/(day.psi)",
        "error": "stb/day",
        "pressure": "psia",
        "potential": "psia",
        "time": "days",
        "rate": "stb/day",
        "length": "ft",
        "area": "ft^2",
        "volume": "ft^3",
        "permeability": "md",
        "viscosity": "cp",
        "gas formation volume factor": "bbl/scf",
        "liquid formation volume factor": "bbl/stb",
        "solution gas oil ratio": "scf/stb",
        "phase gravity": "psi/ft",
        "gas flow rate": "scf/day",
        "liquid flow rate": "stb/day",
        "volumetric velocity": "bbl/(day.ft^2)",
        "density": "lbm/ft^3",
        "compressibility": "psi^{-1}",
        "compressibility factor": "dimensionless",
        "temperature": "R",
        "porosity": "fraction",
        "saturation": "fraction",
        "relative permeability": "fraction",
        "angle": "rad",
        "gravitational acceleration": "ft/(sec^2)",
        "transmissibility conversion": "dimensionless",
        "gravity conversion": "dimensionless",
        "volume conversion": "dimensionless",
    },
    "metric": {
        "transmissibility": "m^3/(day.bar)",
        "error": "m^3/day",
        "pressure": "kpa",
        "potential": "kpa",
        "time": "days",
        "rate": "m^3/day",
        "length": "m",
        "area": "m^2",
        "volume": "m^3",
        "permeability": "m^2",
        "viscosity": "mpa.sec",
        "gas formation volume factor": r"m^3/(std\,m^3)",
        "liquid formation volume factor": r"m^3/(std\,m^3)",
        "solution gas oil ratio": r"(std\,m^3)/(std\,m^3)",
        "phase gravity": "kpa/m",
        "gas flow rate": r"(std\,m^3)/day",
        "liquid flow rate": r"(std\,m^3)/day",
        "volumetric velocity": "m/day",
        "density": "kg/m^3",
        "compressibility": "kpa^{-1}",
        "compressibility factor": "dimensionless",
        "temperature": "K",
        "porosity": "fraction",
        "saturation": "fraction",
        "relative permeability": "fraction",
        "angle": "rad",
        "gravitational acceleration": "m/(sec^2)",
        "transmissibility conversion": "dimensionless",
        "gravity conversion": "dimensionless",
        "volume conversion": "dimensionless",
    },
    "lab": {
        "transmissibility": "cm^3/(sec.atm)",
        "error": "cm^3/sec",
        "pressure": "atm",
        "potential": "atm",
        "time": "sec",
        "rate": "cm^3/sec",
        "length": "cm",
        "area": "cm^2",
        "volume": "cm^3",
        "permeability": "darcy",
        "viscosity": "cp",
        "gas formation volume factor": r"cm^3/(std\,cm^3)",
        "liquid formation volume factor": r"cm^3/(std\,cm^3)",
        "solution gas oil ratio": r"(std\,cm^3)/(std\,cm^3)",
        "phase gravity": "atm/cm",
        "gas flow rate": r"(std\,cm^3)/day",
        "liquid flow rate": r"(std\,cm^3)/day",
        "volumetric velocity": "cm/day",
        "density": "g/cm^3",
        "compressibility": "atm^{-1}",
        "compressibility factor": "dimensionless",
        "temperature": "K",
        "porosity": "fraction",
        "saturation": "fraction",
        "relative permeability": "fraction",
        "angle": "rad",
        "gravitational acceleration": "cm/(sec^2)",
        "transmissibility conversion": "dimensionless",
        "gravity conversion": "dimensionless",
        "volume conversion": "dimensionless",
    },
}

FACTORS = {
    "field": {
        "gravitational acceleration": 32.174,  #: {ft}/{sec^2}
        "transmissibility conversion": 0.001127,  #: dimensionless
        "gravity conversion": 0.21584e-3,  #: dimensionless
        "volume conversion": 5.614583,  #: dimensionless
    },
    "metric": {
        "gravitational acceleration": 9.806635,  #: {m}/{sec^2}
        "transmissibility conversion": 0.0864,  #: dimensionless
        "gravity conversion": 0.001,  #: dimensionless
        "volume conversion": 1,  #: dimensionless
    },
    "lab": {
        "gravitational acceleration": 980.6635,  #: {cm}/{sec^2}
        "transmissibility conversion": 1,  #: dimensionless
        "gravity conversion": 0.986923e-6,  #: dimensionless
        "volume conversion": 1,  #: dimensionless
    },
}


NOMENCLATURE = {
    "abbreviation": {
        "transmissibility": "trans",
        "error": "err",
        "pressure": "press",
        "potential": "poten",
        "time": "t",
        "rate": "q",
        "length": "L",
        "area": "A",
        "volume": "V",
        "permeability": "perm",
        "viscosity": "mu",
        "gas formation volume factor": "FVFg",
        "liquid formation volume factor": "FVFl",
        "water formation volume factor": "FVFw",
        "oil formation volume factor": "FVFo",
        "solution gas oil ratio": "GORs",
        "phase gravity": "gamma phase",
        "gas flow rate": "Qg",
        "liquid flow rate": "Ql",
        "water flow rate": "Qw",
        "oil flow rate": "Qo",
        "volumetric velocity": "u",
        "density": "rho",
        "compressibility": "comp",
        "compressibility factor": "z",
        "temperature": "temp",
        "porosity": "phi",
        "saturation": "sat",
        "relative permeability": "relperm",
        "angle": "theta",
        "gravitational acceleration": "g",
        "transmissibility conversion": "beta constant",
        "gravity conversion": "gamma constant",
        "volume conversion": "alpha constant",
    },
    "symbol": {
        "transmissibility": r"\mathbb{T}",
        "error": "e",
        "pressure": "P",
        "potential": r"\Phi",
        "time": "t",
        "rate": "q",
        "length": "L",
        "area": "A",
        "volume": "V",
        "permeability": "k",
        "viscosity": r"\mu",
        "gas formation volume factor": "B_g",
        "liquid formation volume factor": "B_l",
        "water formation volume factor": "B_w",
        "oil formation volume factor": "B_o",
        "solution gas oil ratio": "R_s",
        "phase gravity": r"\gamma_p",
        "gas flow rate": "q_g",
        "liquid flow rate": "q_l",
        "water flow rate": "q_w",
        "oil flow rate": "q_o",
        "volumetric velocity": "u",
        "density": r"\rho",
        "compressibility": "c",
        "compressibility factor": "z",
        "temperature": "T",
        "porosity": r"\phi",
        "saturation": "S",
        "relative permeability": "k_r",
        "angle": r"\theta",
        "gravitational acceleration": "g",
        "transmissibility conversion": r"\beta_c",
        "gravity conversion": r"\gamma_c",
        "volume conversion": r"\alpha_c",
    },
}


class Base:
    """Base Class."""

    name: str = "Base"
    """Returns class name.

    Returns
    -------
    str
        class name.
    """
    unit: str = "field"
    """Returns class unit.

    .. note::
        Both attributes units and factors are defined based on `unit` 
        argument, for more details, check
        `Units & Factors </user_guide/units_factors/units_factors.html>`_.
        For definitions, check
        `Glossary </user_guide/glossary/glossary.html>`_.

    Returns
    -------
    str
        class unit, default 'field'.
    """
    units: dict = UNITS[unit]
    """Returns class units.

    .. note::
        Both attributes units and factors are defined based on `unit` 
        argument, for more details, check
        `Units & Factors </user_guide/units_factors/units_factors.html>`_.
        For definitions, check
        `Glossary </user_guide/glossary/glossary.html>`_.

    Returns
    -------
    dict
        class units, default 'field' units.
    """
    factors: dict = FACTORS[unit]
    """Returns class factors.

    .. note::
        Both attributes units and factors are defined based on `unit` 
        argument, for more details, check
        `Units & Factors </user_guide/units_factors/units_factors.html>`_.
        For definitions, check
        `Glossary </user_guide/glossary/glossary.html>`_.

    Returns
    -------
    dict
        class factors, default 'field' factors.
    """

    def __init__(
        self,
        unit: str = "field",
        dtype="double",
        verbose: bool = True,
    ):
        """Create Base.

        Parameters
        ----------
        unit : str ('field', 'metric', 'lab'), optional
            unit used in input and output. Both `units` and `factors`
            attributes will be updated based on the selected `unit` and
            can be accessed directly from this class.
        dtype : str or `np.dtype`, optional
            data type used in all arrays. Numpy dtype such as
            `np.single` or `np.double` can be used.
        verbose : bool, optional
            print information for debugging.
        """
        self.set_units(unit)
        self.dtype = dtype
        self.verbose = verbose
        self.n_threads = 8

    def set_units(self, unit: str = "field"):
        """Set object units.

        Parameters
        ----------
        unit : str ('field', 'metric', 'lab'), optional
            unit used in input and output. Both `units` and `factors`
            attributes will be updated based on the selected `unit` and
            can be accessed directly from this class.

        Raises
        ------
        ValueError
            Unknown units.
        """
        if unit in UNITS:
            self.unit = unit
            self.units = UNITS[unit]
            self.factors = FACTORS[unit]
        else:
            raise ValueError(f"Unknown units ({unit}).")

    def set_n_threads(self, n: int = 8):
        """Set number of threads.

        Parameters
        ----------
        n : int, optional
            number of threads, by default 8.
        """
        if n > 0:
            self.n_threads = n
        else:
            raise ValueError("Number of threads must be greater than zero.")

    def report(
        self,
        prop: str = None,
        showindex: bool = True,
        ifmt: int = 0,
    ):
        """Print class report.

        Parameters
        ----------
        prop : str, by default None
            class property name. If None, all class properties will be
            printed.
        showindex : bool, by default True
            show table index.
        ifmt : int, by default 0
            integer format based on the following list:

            .. code-block:: python

                tablefmt = [
                            "pipe",
                            "plain",
                            "simple",
                            "fancy_grid",
                            "presto",
                            "orgtbl",
                ]
        """
        props = vars(self)
        tablefmt = [
            "pipe",
            "plain",
            "simple",
            "fancy_grid",
            "presto",
            "orgtbl",
        ]
        ignore_lst = ["units", "factors", "pv_grid", "corners", "centers"]
        if prop == None:
            print(f"{self.name} Information: \n")
            table = tabulate(
                [
                    (str(k), str(v), "-")
                    for k, v in props.items()
                    if k not in ignore_lst
                ],
                headers=["Property", "Value", "Unit"],
                showindex=showindex,
                tablefmt=tablefmt[ifmt],
            )
            print(table)
            print(" ")
        elif prop in props.keys():
            print(f" - {prop}: {props[prop]}")
        else:
            print(
                "Unknown property.",
                "Use prop=None to print all available properties.",
            )

    def __str__(self):
        return f"{self.__class__.__name__}(unit='{self.unit}', dtype='{self.dtype}')"

    def _repr_args(self):
        args = []
        seen = set()
        # walk child -> parent so inherited parameters (unit, dtype, ...) are included
        for cls in type(self).__mro__:
            args.extend(self._repr_args_for_class(cls, seen))
        return args

    def _repr_args_for_class(self, cls, seen):
        init = cls.__dict__.get("__init__")
        if init is None:
            return []

        args = []
        for name, parameter in inspect.signature(init).parameters.items():
            if name == "self" or name in seen:
                continue
            if parameter.kind in (parameter.VAR_POSITIONAL, parameter.VAR_KEYWORD):
                continue
            seen.add(name)
            if not hasattr(self, name):
                continue  # stored under a different name, skip it
            value = getattr(self, name)
            if isinstance(value, np.ndarray):
                value = f"<array{value.shape}>"  # keep repr short
                args.append(f"{name}={value}")
            else:
                args.append(f"{name}={value!r}")
        return args

    def __repr__(self):
        args = self._repr_args()
        return f"{type(self).__name__}({', '.join(args)})"

    # -------------------------------------------------------------------------
    # End
    # -------------------------------------------------------------------------


if __name__ == "__main__":
    b = Base("field", "single", False)
    b.name = "b"
    k = Base("metric", "single", True)
    k.name = "k"
    # print(repr(b))
    b.report()
    k.report()
