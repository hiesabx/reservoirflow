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
        "gravity": "psi/ft",
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
        "gas formation volume factor": "m^3/(std\,m^3)",
        "liquid formation volume factor": "m^3/(std\,m^3)",
        "solution gas oil ratio": "(std\,m^3)/(std\,m^3)",
        "gravity": "kpa/m",
        "gas flow rate": "(std\,m^3)/day",
        "liquid flow rate": "(std\,m^3)/day",
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
        "gas formation volume factor": "cm^3/(std\,cm^3)",
        "liquid formation volume factor": "cm^3/(std\,cm^3)",
        "solution gas oil ratio": "(std\,cm^3)/(std\,cm^3)",
        "gravity": "atm/cm",
        "gas flow rate": "(std\,cm^3)/day",
        "liquid flow rate": "(std\,cm^3)/day",
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


class _Base:
    name = "Base"

    def __init__(
        self,
        unit: str = "field",
        dtype="double",
        verbose: bool = True,
    ):
        self.set_units(unit)
        self.dtype = dtype
        self.verbose = verbose

    def set_units(self, unit: str = "field"):
        if unit in UNITS.keys():
            self.unit = unit
            self.units = UNITS[unit]
            self.factors = FACTORS[unit]
        else:
            raise ValueError(f"The selected unit system ({unit}) is unknown!")

    def report(
        self,
        prop: str = None,
        showindex: bool = True,
        ifmt: int = 0,
    ):
        props = vars(self)
        tablefmt = ["pipe", "plain", "simple", "fancy_grid", "presto", "orgtbl"]
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
        self.report()
        return ""

    def __repr__(self):
        return str(vars(self))

    # -------------------------------------------------------------------------
    # End
    # -------------------------------------------------------------------------


if __name__ == "__main__":
    b = _Base("metric", "single", True)
    b.name = "b"
    print(repr(b))
    k = _Base("metric", "single", True)
    k.name = "k"
    b.report()
    k.report()
