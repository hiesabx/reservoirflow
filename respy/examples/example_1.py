from respy import models, grids, fluids, wells, plots
'''
This is a small 1D Model for single phase incompressible fluid
The workflow of respy library consists of 8 steps as following:
'''
# Step 1: define 1D grid
grid = grids.Grid1D(nx=4, ny=1, nz=1, dx=300, dy=350, dz=40, phi=0.27, k=270, dtype='double', unit='field')
#       : change specific grid properties (heterogeneity)
grid.dx[2] = 600 # change dx in grid 2
grid.k[3] = 10 # change permeability in grid 3
# Step 2: define a fluid (single phase incompressible fluid)
fluid = fluids.Fluid(mu=0.5 , B=1, dtype='double', unit='field')
# Step 3: create a model
model = models.Model(grid, fluid, dtype='double')
# Step 4: define well locations
#       : method 1 (directly using model method)
model.set_well(i=4, q=-600, s=1.5, r=3.5) # well 1 (Producer)
#       : method 2 (using wells module)
well_2 = wells.Well(i=1, q=600, s=1.5, r=3.5) # well 2 (Injector)
model.set_well(well=well_2)
# Step 5: define boundary conditions
model.set_boundaries({
     0: {'pressure': 4000}, # left boundary at constant pressure
    -1: {'rate': 0}}) # right boundary at constat rate of 0
# Step 6: run the model (single time step)
model.solve(sparse=False, check_MB=True, verbose=False)
# Step 7: show pressures in 3D grid
model.plot_grid(property='pressures', show_boundary=False)
#       : same but using plots module
plots.plot_grid(model, property='pressures', show_boundary=False)
# Step 8: show report
model.report()