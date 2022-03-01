from openresim import grids, fluids, models
import numpy as np
import unittest


class TestApp(unittest.TestCase):

    def test_trans(self):
        model = create_model()
        trans_desired = np.array([28.4004, 28.4004, 28.4004, 28.4004, 28.4004])
        np.testing.assert_array_equal(model.trans, trans_desired)

    def test_pressures(self):
        model = create_model()
        pressures_sparse = model.solve(sparse=True, update=False, check_MB=False)
        pressures_not_sparse = model.solve(sparse=False, update=False, check_MB=False)
        pressures_desired = np.array([3989.4367685 , 3968.31030549, 3947.18384248, 3926.05737947])
        np.testing.assert_almost_equal(pressures_sparse, pressures_not_sparse, decimal=5)
        np.testing.assert_almost_equal(pressures_sparse, pressures_desired, decimal=5)

    def test_well(self):
        model = create_model()
        model.solve(sparse=True, update=True, check_MB=False)
        np.testing.assert_almost_equal(model.wells[4]['q'], -600, decimal=5)
        np.testing.assert_almost_equal(model.wells[4]['r_eq'], 64.53681120105021, decimal=5)
        np.testing.assert_almost_equal(model.wells[4]['G'], 11.08453575337366, decimal=5)
        np.testing.assert_almost_equal(model.wells[4]['pwf'], 3898.992646527117, decimal=5)

    def test_rates(self):
        model = create_model()
        model.solve(sparse=True, update=True, check_MB=False)
        rates_desired = np.array([600.0000000000169, 0.0, 0.0, 0.0, -600.0000000000036, 0.0])
        np.testing.assert_almost_equal(model.rates, rates_desired, decimal=5)


if __name__ == '__main__':
    def create_model():
        grid = grids.Grid1D(nx=4, ny=1, nz=1, dx=300, dy=350, dz=40, phi=0.27, k=270, dtype='double')
        fluid = fluids.Fluid(mu=0.5 , B=1, dtype='double')
        model = models.Model(grid, fluid, dtype='double')
        model.set_well(i=4, q=-600, s=1.5, r=3.5)
        model.set_boundaries({0: {'pressure': 4000}, -1: {'rate': 0}})
        return model
    unittest.main()

    