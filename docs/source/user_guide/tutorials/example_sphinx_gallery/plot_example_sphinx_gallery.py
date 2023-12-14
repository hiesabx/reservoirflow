"""

Example (Sphinx Gallery)
========================

This is a 3D Model for single-phase compressible fluid

Simulation Steps
################

Simulation steps can be done in 8 steps:
  - Define 1D grid
  - Define a fluid
  - Create a model
  - Add a well
  - Define boundary conditions
  - Run the model
  - Show pressures
  - Show DataFrame
"""
# %%

###############################################################################
# Head 1
# ======
#
# This to show that markdown `works` fine here.
# $a=\frac{a}{4}$
#
# __List__
#
# - one
# - two
# - two

# example 1
import reservoirflow as rf

rf.__version__
# sphinx_gallery_thumbnail_path = "user_guide/tutorials/example_sphinx_gallery/_static/example_image.png"

# %%

# Step 1: Define 1D grid (default dtype: 'double')
grid = rf.grids.RegularCartesian(
    nx=4,
    ny=1,
    nz=1,
    dx=300,
    dy=350,
    dz=40,
    phi=0.27,
    kx=270,
    comp=1 * 10**-6,
)

# Step 2: Define a fluid (single phase incompressible fluid)
fluid = rf.fluids.SinglePhase(
    mu=0.5,
    B=1,
    rho=50,
    comp=1 * 10**-5,
)

# Step 3: Create a model
model = rf.models.BlackOil(
    grid,
    fluid,
    pi=4000,
)

# Step 4: Add a well
model.set_well(
    cell_id=4,
    q=-600,
    s=1.5,
    r=3.5,
)  # well 1 (Producer)

# Step 5: Define boundary conditions
model.set_boundaries(
    {
        0: ("pressure", 4000),  # left boundary (constant pressure)
        5: ("rate", 0),  # right boundary (no flow)
    }
)

# Step 6: Run the model (single time step)
model.solve(sparse=False, check_MB=False)


# %%
# show the 3D model

# Step 7: Show pressures in 3D grid
# model.show(property="pressures", boundary=False)

# %%
###############################################################################
# Model DataFrame
# ++++++++++++++++

# Step 8: Show report
df = model.get_df()
df
