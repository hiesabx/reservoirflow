from reservoirflow import fluids, grids, models, wells

grid = grids.RegularCartesian(
    nx=2,
    ny=2,
    nz=1,
    dx=100,
    dy=100,
    dz=10,
    phi=0.3,
    kx=10,
    ky=10,
    kz=10,
    comp=0,
    dtype="double",
)
fluid = fluids.SinglePhase(mu=1.5, B=1, rho=50, comp=1 * 10**-5, dtype="double")
grid.show("id")
model = models.BlackOil(grid, fluid, pi=6000, dt=5, dtype="double")
model.set_well(id=9, q=-400, pwf=1500, s=0, r=3.5)
model.set_boundaries({0: ("rate", 0)})
model.run(nsteps=30, sparse=True)
model.get_df(cells_pressures=True, w_rates=True, w_pressures=True)
