from venv import create
from openresim import grids, fluids, models
import numpy as np
import unittest


class TestApp(unittest.TestCase):
    def test_trans(self):
        trans_desired = np.array([28.4004, 28.4004, 28.4004, 28.4004, 28.4004])
        model = create_model()
        np.testing.assert_array_equal(model.T["x"], trans_desired)

    def test_RHS(self):
        RHS_desired = np.array(
            [
                2.221714417615699,
                2.221714417615699,
                2.221714417615699,
                2.221714417615699,
                2.221714417615699,
                2.221714417615699,
            ]
        )
        model = create_model()
        np.testing.assert_array_equal(model.RHS, RHS_desired)

    def test_step_0(self):
        initial_pressures_desired = np.array(
            [4000.0, 4000.0, 4000.0, 4000.0, 4000.0, np.nan]
        )
        model = create_model()
        np.testing.assert_almost_equal(
            model.pressures[0], initial_pressures_desired, decimal=5
        )

    def test_step_1(self):
        pressures_desired = np.array(
            [
                4000.0,
                3993.7457054727147,
                3980.7478537203187,
                3966.2439396955738,
                3949.0993471641805,
                np.nan,
            ]
        )
        rates_desired = np.array([355.24893259, 0.0, 0.0, 0.0, -600.0, 0.0])
        error_desired = 0.99891
        cumulative_error_desired = 0.28973

        model = create_model()
        model.solve(sparse=True, update=True, check_MB=True)
        np.testing.assert_almost_equal(model.pressures[1], pressures_desired, decimal=5)
        np.testing.assert_almost_equal(model.rates[1], rates_desired, decimal=5)
        np.testing.assert_almost_equal(1 - model.error, error_desired, decimal=3)
        np.testing.assert_almost_equal(
            model.cumulative_error, cumulative_error_desired, decimal=5
        )

        model = create_model()
        model.solve(sparse=False, update=True, check_MB=True)
        np.testing.assert_almost_equal(model.pressures[1], pressures_desired, decimal=5)
        np.testing.assert_almost_equal(model.rates[1], rates_desired, decimal=5)
        np.testing.assert_almost_equal(1 - model.error, error_desired, decimal=3)
        np.testing.assert_almost_equal(
            model.cumulative_error, cumulative_error_desired, decimal=5
        )

    def test_step_2(self):
        pressures_desired = np.array(
            [
                4000.0,
                3990.953458875824,
                3972.6419439813862,
                3953.6963177175703,
                3933.7691125796787,
                np.nan,
            ]
        )
        rates_desired = np.array([513.85077309, 0.0, 0.0, 0.0, -600.0, 0.0])
        error_desired = 0.99891
        cumulative_error_desired = 0.177731

        model = create_model()
        model.run(nsteps=2, sparse=True, check_MB=True)
        np.testing.assert_almost_equal(model.pressures[2], pressures_desired, decimal=5)
        np.testing.assert_almost_equal(model.rates[2], rates_desired, decimal=5)
        np.testing.assert_almost_equal(1 - model.error, error_desired, decimal=3)
        np.testing.assert_almost_equal(
            model.cumulative_error, cumulative_error_desired, decimal=5
        )

    def test_well(self):
        model = create_model()
        model.solve(sparse=True, update=True, check_MB=True)
        np.testing.assert_almost_equal(
            model.wells[4]["r_eq"], 64.53681120105021, decimal=5
        )
        np.testing.assert_almost_equal(model.wells[4]["q"], -600, decimal=5)
        np.testing.assert_almost_equal(
            model.wells[4]["G"], 11.08453575337366, decimal=5
        )
        np.testing.assert_almost_equal(
            model.wells[4]["pwf"], 3922.0346142177686, decimal=5
        )


if __name__ == "__main__":

    def create_model():
        grid = grids.Cartesian(
            nx=4,
            ny=1,
            nz=1,
            dx=300,
            dy=350,
            dz=40,
            phi=0.27,
            kx=270,
            comp=1 * 10**-6,
            dtype="double",
        )
        fluid = fluids.SinglePhase(
            mu=0.5, B=1, rho=50, comp=1 * 10**-5, dtype="double"
        )
        model = models.Model(grid, fluid, pi=4000, dtype="double", verbose=False)
        model.set_well(id=4, q=-600, s=1.5, r=3.5)
        model.set_boundaries({0: ("pressure", 4000), 5: ("rate", 0)})
        return model

    unittest.main()
