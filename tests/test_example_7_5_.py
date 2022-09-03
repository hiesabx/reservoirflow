from openresim import grids, fluids, models
import numpy as np
import unittest
import warnings


class TestApp(unittest.TestCase):
    def test_trans(self):
        model = create_model()
        trans_desired = np.array([28.4004, 28.4004, 28.4004, 28.4004, 28.4004])
        np.testing.assert_array_equal(model.T["x"], trans_desired)

    def test_fluid(self):
        model = create_model()
        gravity_desired = 0.347221808
        np.testing.assert_almost_equal(model.fluid.g, gravity_desired, decimal=5)

    def test_pressures(self):
        p_desired = np.array([3978.88469, 3936.65409, 3894.42348, 3852.19288])
        # sparse:
        model = create_model()
        model.solve(sparse=True, update=True, check_MB=True)
        p_sparse = model.pressures[-1, model.cells_id]
        np.testing.assert_almost_equal(p_sparse, p_desired, decimal=5)
        # dense:
        model = create_model()
        model.solve(sparse=False, update=True, check_MB=True)
        p_not_sparse = model.pressures[-1, model.cells_id]
        np.testing.assert_almost_equal(p_not_sparse, p_desired, decimal=5)

    def test_well(self):
        model = create_model()
        model.solve(sparse=True, update=True, check_MB=True)
        np.testing.assert_almost_equal(model.wells[4]["q"], -600, decimal=5)
        np.testing.assert_almost_equal(
            model.wells[4]["r_eq"], 64.53681120105021, decimal=5
        )
        np.testing.assert_almost_equal(
            model.wells[4]["G"], 11.08453575337366, decimal=5
        )
        np.testing.assert_almost_equal(
            model.wells[4]["pwf"], 3825.1281513112767, decimal=5
        )

    def test_rates(self):
        model = create_model()
        model.solve(sparse=True, update=True, check_MB=True)
        rates_desired = np.array(
            [600.0000000000247, 0.0, 0.0, 0.0, -600.0000000000036, 0.0]
        )
        np.testing.assert_almost_equal(model.rates[-1], rates_desired, decimal=5)

    def test_simulation_run(self):
        model = create_model()
        model.run(30, True)
        model = create_model()
        model.run(30, False)


def create_model():
    z = np.array([3212.73, 3182.34, 3121.56, 3060.78, 3000, 2969.62])
    grid = grids.Cartesian(
        nx=4,
        ny=1,
        nz=1,
        dx=300,
        dy=350,
        dz=40,
        z=z,
        phi=0.27,
        kx=270,
        dtype="double",
    )
    fluid = fluids.SinglePhase(mu=0.5, B=1, rho=50, dtype="double")
    model = models.Model(grid, fluid, dtype="double", verbose=False)
    model.set_well(id=4, q=-600, s=1.5, r=3.5)
    model.set_boundaries({0: ("pressure", 4000), 5: ("rate", 0)})
    return model


if __name__ == "__main__":
    unittest.main()
