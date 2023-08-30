from reservoirflow.base import Base
import numpy as np


class Fluid(Base):
    
    name = "Fluid"
    
    def __init__(self, unit, dtype, verbose):
        super().__init__(unit, dtype, verbose)
        
        
if __name__ == "__main__":
    dtype="double"
    unit="field"
    verbose=False
    fluid = Fluid(unit, dtype, verbose)
    print(fluid)