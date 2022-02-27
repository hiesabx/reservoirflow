from openresim import grids, fluids, models
import numpy as np
import unittest


class TestApp(unittest.TestCase):

    def test_trans(self):
        trans_desired = np.array([28.4004, 28.4004, 28.4004, 28.4004, 28.4004])
        np.testing.assert_array_equal(model.trans, trans_desired)

    def test_fluid(self):
        gravity_desired = 0.347221808
        np.testing.assert_almost_equal(model.fluid.gravity, gravity_desired, decimal=5)

    def test_pressures(self):
        pressures_sparse = model.solve(sparse=True)
        pressures_not_sparse = model.solve(sparse=False)
        pressures_desired = np.array([3978.884697751098, 3936.654093253294, 3894.423488755491, 3852.1928842576876])
        np.testing.assert_array_equal(pressures_sparse, pressures_not_sparse)
        np.testing.assert_almost_equal(pressures_sparse, pressures_desired, decimal=5)

    def test_well(self):
        np.testing.assert_almost_equal(model.wells[4]['q'], -600, decimal=5)
        np.testing.assert_almost_equal(model.wells[4]['r_eq'], 64.53681120105021, decimal=5)
        np.testing.assert_almost_equal(model.wells[4]['G'], 11.08453575337366, decimal=5)
        np.testing.assert_almost_equal(model.wells[4]['pwf'], 3825.1281513112767, decimal=5)

    def test_rates(self):
        rates_desired = np.array([600.0000000000247, 0.0, 0.0, 0.0, -600.0000000000036, 0.0])
        np.testing.assert_almost_equal(model.rates, rates_desired, decimal=5)


if __name__ == '__main__':
    z = [3212.73, 3182.34, 3121.56, 3060.78, 3000, 2969.62]
    grid = grids.Grid1D(nx=4, ny=1, nz=1, dx=300, dy=350, dz=40, z=z, phi=0.27, k=270, dtype='double')
    fluid = fluids.Fluid(mu=0.5, B=1, rho=50, dtype='double')
    model = models.Model(grid, fluid, dtype='double')
    model.set_well(i=4, q=-600, s=1.5, r=3.5)
    model.set_boundaries({0: {'pressure': 4000}, -1: {'rate': 0}})
    model.solve(verbose=False)
    unittest.main()

    