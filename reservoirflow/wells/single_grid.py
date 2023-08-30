from reservoirflow.wells.well import Well


class SingleGrid(Well):
    """
    Single Grid class.
    """

    name = "Single Grid Well"

    def __init__(self, id=None, q=None, s=None, r=None, constrain=None):
        self.set_props(id, q, s, r, constrain)

    def set_id(self, id):
        self.id = id

    def set_q(self, q):
        self.q = q

    def set_s(self, s):
        self.s = s

    def set_r(self, r):
        self.r = r

    def set_constrain(self, constrain):
        self.constrain = constrain

    def set_props(self, id=None, q=None, s=None, r=None, constrain=None):
        if id:
            self.set_id(id)
        if q:
            self.set_q(q)
        if s:
            self.set_s(s)
        if r:
            self.set_r(r)
        self.set_constrain(constrain)

    # -------------------------------------------------------------------------
    # Synonyms:
    # -------------------------------------------------------------------------

    def allow_synonyms(self):
        self.set_properties = self.set_props
        self.set_rate = self.set_q
        self.rate = self.q
        self.set_skin = self.set_s
        self.skin = self.s
        self.set_radius = self.set_r
        self.radius = self.r

    # -------------------------------------------------------------------------
    # End
    # -------------------------------------------------------------------------


if __name__ == "__main__":
    well = SingleGrid(id=4, q=-600, s=1.5, r=3.5)
    print(well)
