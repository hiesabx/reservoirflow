# %%
import numpy as np
import pandas as pd

from reservoirflow import fluids, grids, models
from reservoirflow.utils import profme


def create_model():
    grid = grids.RegularCartesian(
        nx=500,
        ny=500,
        nz=4,
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
    model = models.BlackOil(
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
    # None,
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
############################################################
Before solvers module:
############################################################

(nx=100, ny=100, nz=2)
    isolver         error  ctime
0      None  2.981393e-12  14.48
4       cgs  7.852105e-05   6.98 > lowest error
3        cg  1.074206e-04   6.90 > fastest ctime
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
4       cgs  3.705438e-04  104.20 > lowest error
2  bicgstab  1.077852e-03  101.92
1      bicg  2.992796e-03   99.31 > fastest ctime
3        cg  2.992812e-03  100.45
6    lgmres  3.147717e-03  104.89
5     gmres  3.302489e-03  101.69
8       qmr  3.302587e-03  104.15
9   gcrotmk  3.303104e-03  104.72
7    minres  9.061416e+00  105.99

############################################################
After solvers module: 12.11.2022
############################################################

(nx=100, ny=100, nz=2)
    isolver         error  ctime
0      None  2.980949e-12  16.24
4       cgs  8.929843e-05   2.27 > lowest error
3        cg  1.052619e-04   2.19 > fastest ctime
1      bicg  1.070314e-04   2.65
6    lgmres  1.590616e-04   2.46
2  bicgstab  2.298655e-04   2.40
8   gcrotmk  6.567435e-04   2.37
7       qmr  7.232530e-04   2.32
5     gmres  7.392049e-04   2.33

(nx=500, ny=500, nz=1)
    isolver         error   ctime
0      None  9.375833e-12  179.73
4       cgs  3.705177e-04   24.06 > lowest error
2  bicgstab  1.078263e-03   22.77 > fastest ctime
3        cg  2.992928e-03   23.58
1      bicg  2.992964e-03   25.28
8   gcrotmk  3.302225e-03   23.83
7       qmr  3.302482e-03   23.74
5     gmres  3.303131e-03   24.36
6    lgmres  3.303150e-03   24.99

(nx=500, ny=500, nz=4)
    isolver     error   ctime
3       cgs  0.000221  124.22 > lowest error
0      bicg  0.001043  108.34 > fastest ctime
2        cg  0.001043  116.33
1  bicgstab  0.001394  120.88
4     gmres  0.002027  125.36
6       qmr  0.002034  128.71
5    lgmres  0.002043  121.59
7   gcrotmk  0.002050  122.95

############################################################
After using cache in is_homogeneous and is_regular: 20.11.2022
############################################################

(nx=500, ny=500, nz=4)
    isolver     error   ctime
3       cgs  0.000221  111.37 > lowest error
2        cg  0.001043  112.85
0      bicg  0.001043  110.04
1  bicgstab  0.001393  108.95 > fastest ctime
5    lgmres  0.001991  122.98
6       qmr  0.002034  134.56
7   gcrotmk  0.002038  128.41
4     gmres  0.002040  111.76

"""
