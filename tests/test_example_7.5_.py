from openresim import grids, fluids, models
import numpy as np
import unittest


class TestApp(unittest.TestCase):
    def test_trans(self):
        model = create_model()
        trans_desired = np.array([28.4004, 28.4004, 28.4004, 28.4004, 28.4004])
        np.testing.assert_array_equal(model.T, trans_desired)

    def test_fluid(self):
        model = create_model()
        gravity_desired = 0.347221808
        np.testing.assert_almost_equal(model.fluid.gravity, gravity_desired, decimal=5)

    def test_pressures(self):
        model = create_model()
        pressures_sparse = model.solve(sparse=True, update=True, check_MB=True)
        model = create_model()
        pressures_not_sparse = model.solve(sparse=False, update=True, check_MB=True)
        pressures_desired = np.array(
            [
                3978.884697751098,
                3936.654093253294,
                3894.423488755491,
                3852.1928842576876,
            ]
        )
        np.testing.assert_almost_equal(
            pressures_sparse, pressures_not_sparse, decimal=5
        )
        np.testing.assert_almost_equal(pressures_sparse, pressures_desired, decimal=5)

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
        np.testing.assert_almost_equal(model.rates, rates_desired, decimal=5)


if __name__ == "__main__":

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
        model = models.Model(grid, fluid, dtype="double")

        # model.grid.pv_grid.points = model.grid.get_pv_grid(show_boundary=True) += z
        model.set_well(id=4, q=-600, s=1.5, r=3.5)
        model.set_boundaries({0: {"pressure": 4000}, -1: {"rate": 0}})
        return model

    model = create_model()
    model.solve(verbose=False)
    model.show_grid("pressures")
    # unittest.main()
