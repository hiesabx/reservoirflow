from openresim import models, grids, fluids, wells, plots

"""
This is a 1D Model for single-phase incompressible fluid
The workflow of openresim library consists of 8 steps as following:
(default: dtype='double', unit='field')
"""
# Step 1: Define 1D grid (default dtype: 'double')
grid = grids.Cartesian(
    nx=4, ny=1, nz=1, dx=300, dy=350, dz=40, phi=0.27, kx=270, comp=1 * 10**-6
)
# Step 2: Define a fluid (single phase incompressible fluid)
fluid = fluids.SinglePhase(mu=0.5, B=1, rho=50, comp=1 * 10**-5)
# Step 3: Create a model
model = models.Model(grid, fluid, pi=4000)
# Step 4: Add a well
model.set_well(id=4, q=-600, s=1.5, r=3.5)  # well 1 (Producer)
# Step 5: Define boundary conditions
model.set_boundaries(
    {0: {"pressure": 4000}, -1: {"rate": 0}}  # left boundary (constant pressure)
)  # right boundary (constat rate of 0)
# Step 6: Run the model (single time step)
model.solve(sparse=False, check_MB=False, verbose=True)
# print(model)
# # Step 7: Show pressures in 3D grid
model.show(property="pressures", show_boundary=False)
# Step 8: Show report
# model.report()
