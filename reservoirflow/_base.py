from tabulate import tabulate

UNITS = {
    "field": {
        "transmissibility": "stb/d-psi",
        "error": "stb/d",
        "pressure": "psia",
        "time": "days",
        "rate": "stb/day",
        "length": "ft",
        "area": "ft^2",
        "volume": "ft^3",
    },
    "metric": {
        "transmissibility": "m^3/D-bar",
        "error": "m^3/day",
        "pressure": "kpa",
        "time": "days",
        "rate": "m^3/day",
        "length": "m",
        "area": "m^2",
        "volume": "m^3",
    },
    "lab": {
        "transmissibility": "cm^3/sec-atm",
        "error": "cm^3/sec",
        "pressure": "atm",
        "time": "sec",
        "rate": "cm^3/sec",
        "length": "cm",
        "area": "cm^2",
        "volume": "cm^3",
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

    def __init__(self, unit: str = "field", dtype="double", verbose: bool = True):
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

    def report(self, prop: str = None, showindex: bool = True, ifmt: int = 0):
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
            print("Unknown property! Use prop=None to print all available properties.")

    def __str__(self):
        self.report()
        return ""

    def __repr__(self):  # print(repr(class))
        return str(vars(self))

    # -------------------------------------------------------------------------
    # Synonyms:
    # -------------------------------------------------------------------------

    # def allow_synonyms(self):
    #     self.set_compressibility = self.set_comp
    #     self.compressibility = self.comp
    #     self.compressibility_type = self.comp_type

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
