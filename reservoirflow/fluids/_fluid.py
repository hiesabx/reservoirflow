from reservoirflow._base import Base


class Fluid(Base):
    name = "Fluid"

    def __init__(self, unit, dtype, verbose):
        super().__init__(unit, dtype, verbose)

    def set_comp(self, comp: float):
        """Fluid compressibility

        Parameters
        ----------
        comp : float
            fluid compressibility.

        Raises
        ------
        ValueError
            Compressibility smaller than zero is not allowed.
        """
        self.comp = comp
        if comp == 0:
            self.comp_type = "incompressible"
        elif comp > 0:
            self.comp_type = "compressible"
        else:
            self.comp_type = None
            raise ValueError("Compressibility smaller than zero is not allowed.")

    # -------------------------------------------------------------------------
    # Synonyms:
    # -------------------------------------------------------------------------

    def allow_synonyms(self):
        self.set_compressibility = self.set_comp
        self.compressibility = self.comp
        self.compressibility_type = self.comp_type

    # -------------------------------------------------------------------------
    # End
    # -------------------------------------------------------------------------


if __name__ == "__main__":
    dtype = "double"
    unit = "field"
    verbose = False
    fluid = Fluid(unit, dtype, verbose)
    print(fluid)
