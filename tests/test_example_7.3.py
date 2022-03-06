from openresim import grids, fluids, models
import numpy as np
import unittest


class TestApp(unittest.TestCase):

    def test_trans(self):
        trans_desired = np.array([28.4004, 28.4004, 28.4004, 28.4004, 28.4004])
        model = create_model()
        np.testing.assert_array_equal(model.trans, trans_desired)

    def test_pressures(self):
        pressures_desired = np.array([3989.437537513856, 3968.3126125415683, 3947.1876875692815, 3926.0627625969946])
        model = create_model()
        pressures_sparse = model.solve(sparse=True, update=True, check_MB=True)
        model = create_model()
        pressures_not_sparse = model.solve(sparse=False, update=True, check_MB=True)
        np.testing.assert_almost_equal(pressures_sparse, pressures_not_sparse, decimal=5)
        np.testing.assert_almost_equal(pressures_sparse, pressures_desired, decimal=5)

    def test_well(self):
        model = create_model()
        model.solve(sparse=True, update=True, check_MB=True)
        np.testing.assert_almost_equal(model.wells[4]['q'], -599.9563191829004, decimal=5)
        np.testing.assert_almost_equal(model.wells[4]['r_eq'], 64.53681120105021, decimal=5)
        np.testing.assert_almost_equal(model.wells[4]['G'], 11.08453575337366, decimal=5)
        np.testing.assert_almost_equal(model.wells[4]['pwf'], 3899, decimal=5)

    def test_rates(self):
        rates_desired = np.array([599.9563191829617, 0.0, 0.0, 0.0, -599.9563191829004, 0.0])
        model = create_model()
        model.solve(sparse=True, update=True, check_MB=True)
        np.testing.assert_almost_equal(model.rates[1], rates_desired, decimal=5)


if __name__ == '__main__':
    def create_model():
        grid = grids.Grid1D(nx=4, ny=1, nz=1, dx=300, dy=350, dz=40, phi=0.27, k=270, dtype='double')
        fluid = fluids.SinglePhaseFluid(mu=0.5 , B=1, dtype='double')
        model = models.Model(grid, fluid, dtype='double')
        model.set_well(i=4, pwf=3899, s=1.5, r=3.5)
        model.set_boundaries({0: {'pressure': 4000}, -1: {'rate': 0}})
        return model
    unittest.main()

