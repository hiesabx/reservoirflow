import numpy as np
from scipy import linalg as lin
from openresim import profme, grids, fluids, models

np.set_printoptions(threshold=np.inf, linewidth=100000)


def create_model():
    grid = grids.Cartesian(
        nx=3,
        ny=4,
        nz=7,
        dx=[1, 11, 3, 4, 5],
        dy=[10, 11, 11, 12, 13, 14],
        dz=[2, 3, 2, 3, 4, 4, 3, 4, 5],
        phi=0.27,
        kx=150,
        ky=3,
        kz=10,
        dtype="double",
        unify=False,
    )
    fluid = fluids.SinglePhase(
        mu=3.5,
        B=1,
        rho=50,
        dtype="double",
    )
    model = models.Model(grid, fluid, dtype="double", verbose=False)
    return model


model = create_model()
T1 = model.get_cells_T_loop(False, False)
print(T1)
g = model.grid.get_cells_G_array()
T2 = g / (model.fluid.mu * model.fluid.B)
print(T2)
print(T2 - T1)
print((T2 - T1).sum())
