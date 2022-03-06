from openresim import models, grids, fluids, wells, plots
'''
This is a 1D Model for single-phase incompressible fluid
(default: dtype='double', unit='field')
The workflow of openresim library consists of 8 steps as following:
'''
# Step 1: Define 1D grid (default dtype: 'double')
grid = grids.Grid1D(nx=4, ny=1, nz=1, dx=300, dy=350, dz=40, 
                    phi=0.27, k=270)
#       : change specific grid properties to add some heterogeneity
grid.dx[2] = 600    # change dx in grid 2
grid.k[3] = 10      # change permeability in grid 3
# Step 2: Define a fluid (single phase incompressible fluid)
fluid = fluids.SinglePhaseFluid(mu=0.5, B=1)
# Step 3: Create a model
model = models.Model(grid, fluid)
# Step 4: Add wells
#       : method 1 - directly using model method
model.set_well(i=4, q=-600, s=1.5, r=3.5) # well 1 (Producer)
#       : method 2 - using wells module
well_2 = wells.Well(i=1, q=600, s=1.5, r=3.5)   # well 2 (Injector)
model.set_well(well=well_2)
# Step 5: Define boundary conditions
model.set_boundaries({
     0: {'pressure': 4000},     # left boundary (constant pressure)
    -1: {'rate': 0}})           # right boundary (constat rate of 0)
# Step 6: Run the model (single time step)
model.solve(sparse=False, check_MB=True, verbose=False)
# Step 7: Show pressures in 3D grid
#       : method 1 - directly using model method
model.show_grid(property='pressures', show_boundary=False)
#       : method 2 - using wells module
plots.show_grid(model, property='pressures', show_boundary=False, show_bounds=False)
# Step 8: Show report
model.report()
