from reservoirflow._base import Base


class Model(Base):
    """Model class.

    Model class used to create a reservoir simulation model.

    Model class represents the fluid flow process in a reservoir
    due to pressure change cause by boundary conditions or by
    (production or injection) wells.
    """

    name = "Model"

    def __init__(self, unit, dtype, verbose):
        super().__init__(unit, dtype, verbose)

    def set_comp(self, comp: float):
        """Model compressibility

        Parameters
        ----------
        comp : float
            model compressibility which is the total compressibility
            of Grid and Fluid classes.

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
    model = Model(unit, dtype, verbose)
    print(model)
