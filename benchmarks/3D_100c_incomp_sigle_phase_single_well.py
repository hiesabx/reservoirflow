#%%
from openresim import grids, fluids, models
import numpy as np
import pandas as pd
import cProfile
import timeit


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


def benchmark(sparse=True, threading=True, check_MB=True, update=True):
    model = create_model()
    model.solve(sparse, threading, check_MB, update)
    # model.run(3, sparse, threading, check_MB)


results = timeit.repeat(
    "benchmark()",
    "from __main__ import benchmark, create_model",
    number=1,
    repeat=10,
)
print(results)
# #%%
# %timeit benchmark(False, False, True, True)

# #%%
# cProfile.run("benchmark(False, False, False, False)", sort="cumtime")

# #%%
# %timeit benchmark(False, True, True, True)
# #%%
# %timeit benchmark(True, False, True, True)
# #%%
# %timeit benchmark(True, True, True, True)
# # %%
