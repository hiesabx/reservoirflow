from tabulate import tabulate

units_dict = {
    "field": {
        "transmissibility": "STB/D-psi",
        "error": "STB/D",
        "pressure": "Psia",
        "time": "Days",
        "rate": "STB/Day",
    },
    "metric": {"transmissibility": "M3/D-bar", "error": "M3/D"},
}

factors_dict = {
    "field": {
        "transmissibility conversion": 0.001127,
        "gravity conversion": 0.21584 * 10**-3,
        "gravitational acceleration": 32.174,
        "volume conversion": 5.614583,
    },
    "metric": {
        "transmissibility conversion": 0,
        "gravity conversion": 0,
        "gravitational acceleration": 0,
        "volume conversion": 1,
    },
}


class Base:
    name = None
    units_dict = units_dict
    factors_dict = factors_dict

    unit = "field"
    units = units_dict[unit]
    factors = factors_dict[unit]

    def __init__(self, unit="field", dtype="double", verbose=True):
        self.set_units(unit)
        self.dtype = dtype
        self.verbose = verbose

    def set_units(self, unit="field"):
        if unit in self.units_dict.keys():
            self.unit = unit
            self.units = self.units_dict[unit]
            self.factors = self.factors_dict[unit]
        else:
            raise ValueError(f"The selected unit system ({unit}) is unknown!")

    def report(self, prop=None, ifmt=0):
        props = vars(self)
        tablefmt = ["pipe", "plain", "simple", "fancy_grid", "presto", "orgtbl"]
        ignore_lst = ["units", "factors", "pv_grid", "corners", "centers"]
        if prop == None:
            print(f"{self.name} Information: \n")
            table = tabulate(
                [(str(k), str(v)) for k, v in props.items() if k not in ignore_lst],
                headers=["Property", "Value"],
                showindex="always",
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
    b = Base("metric", "single", True)
    b.name = "b"
    print(repr(b))
    k = Base("metric", "single", True)
    k.name = "k"
    b.report()
    k.report()
