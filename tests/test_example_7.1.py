from openresim import grids, fluids, models
import numpy as np
import unittest


class TestApp(unittest.TestCase):

    def test_pressures(self):
        desired_pressures = np.array([3989.4367685 , 3968.31030549, 3947.18384248, 3926.05737947])
        np.testing.assert_array_equal(pressures_sparse, pressures_not_sparse)
        np.testing.assert_almost_equal(pressures_sparse, desired_pressures, decimal=5)


if __name__ == '__main__':
    grid = grids.Grid1D(nx=4, ny=1, nz=1, dx=300, dy=350, dz=40, phi=0.27, k=270, dtype='double')
    fluid = fluids.Fluid(mu=0.5 , B=1, dtype='double')
    model = models.Model(grid, fluid, dtype='single')
    model.set_well(i=4, q=-600, s=1.5, r=3.5)
    model.set_boundaries({0: {'pressure': 4000}, -1: {'rate': 0}})
    pressures_sparse = model.solve(sparse=True)
    pressures_not_sparse = model.solve(sparse=False)
    unittest.main()