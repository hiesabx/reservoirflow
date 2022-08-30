#%%
from openresim import grids, fluids, models, utils
import numpy as np
import pandas as pd


def create_model():
    grid = grids.Cartesian(
        nx=4,
        ny=1,
        nz=1,
        dx=300,
        dy=350,
        dz=40,
        phi=0.27,
        kx=270,
        dtype="double",
    )
    fluid = fluids.SinglePhase(mu=0.5, B=1, dtype="double")
    model = models.Model(grid, fluid, dtype="double", dt=1, verbose=False)
    model.set_well(id=4, q=-600, s=1.5, r=3.5)
    model.set_boundaries({0: ("pressure", 4000), 5: ("rate", 0)})
    return model


@utils.profile("tottime", True, True)
def benchmark(sparse=True, threading=False, check_MB=True, update=True):
    model = create_model()
    model.run(1, sparse, threading, check_MB)


benchmark()
