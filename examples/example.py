from openresim import grids, fluids, models, wells

dx = [5,10,20,30,40]

grid = grids.CartGrid(nx=3,ny=1,nz=1, dx=dx, dy=10, dz=5, phi=0.4, k=10)
grid.show(boundary=True, centers='volume')

fluid = fluids.SinglePhaseFluid(mu=10, B=1, rho=64)
w1 = wells.Well(i=(), )
model = models.Model(grid, fluid, )
model.run()


for case in cases:
    case.run(statement, globals=None, locals=None)
    