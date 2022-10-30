#%%
from openresim import grids, fluids, models, profme
import numpy as np
import pandas as pd


def create_model():
    grid = grids.Cartesian(
        nx=100,
        ny=100,
        nz=2,
        dx=300,
        dy=350,
        dz=20,
        phi=0.27,
        kx=1,
        ky=1,
        kz=0.1,
        comp=1 * 10**-6,
        dtype="double",
    )
    fluid = fluids.SinglePhase(mu=0.5, B=1, rho=50, comp=1 * 10**-5, dtype="double")
    model = models.Model(
        grid, fluid, pi=6000, dt=5, start_date="10.10.2018", dtype="double"
    )

    cells_id = model.grid.get_cells_id(False, True)[-1].flatten()
    wells = np.random.choice(cells_id, 6, False)
    for id in wells:
        model.set_well(id=id, q=-300, pwf=100, s=1.5, r=3.5)

    wells = np.random.choice(cells_id, 6, False)
    for id in wells:
        model.set_well(id=id, q=100, s=0, r=3.5)

    return model


isolvers = [
    None,
    "bicg",
    "bicgstab",
    "cg",
    "cgs",
    "gmres",
    "lgmres",
    # "minres",
    "qmr",
    "gcrotmk",
    # "tfqmr",
]
results = {"isolver": [], "error": [], "ctime": []}
for isolver in isolvers:
    model = create_model()
    model.run(10, isolver=isolver)
    results["isolver"].append(isolver)
    results["error"].append(model.error)
    results["ctime"].append(model.ctime)
results = pd.DataFrame(results)
results.sort_values(["error", "ctime"], inplace=True)
print(results)

"""
(nx=100, ny=100, nz=2)
    isolver         error  ctime
0      None  2.981393e-12  14.48
4       cgs  7.852105e-05   6.98 > best option
3        cg  1.074206e-04   6.90
1      bicg  1.152981e-04   7.02
2  bicgstab  2.067467e-04   7.36
9   gcrotmk  6.102345e-04   7.09
6    lgmres  6.294939e-04   7.03
5     gmres  6.756949e-04   6.96
8       qmr  6.790783e-04   7.03
7    minres  3.474882e+02   7.09

(nx=500, ny=500, nz=1)
    isolver         error   ctime
0      None  1.712275e-11  204.03
4       cgs  3.705438e-04  104.20 > best option
2  bicgstab  1.077852e-03  101.92
1      bicg  2.992796e-03   99.31
3        cg  2.992812e-03  100.45
6    lgmres  3.147717e-03  104.89
5     gmres  3.302489e-03  101.69
8       qmr  3.302587e-03  104.15
9   gcrotmk  3.303104e-03  104.72
7    minres  9.061416e+00  105.99
"""
