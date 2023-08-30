from reservoirflow import grids, fluids, models
import numpy as np
import unittest


class TestApp(unittest.TestCase):
    def test_trans(self):
        trans_desired = np.array([28.4004, 28.4004, 28.4004, 28.4004, 28.4004])
        model = create_model()
        trans = model.get_cells_T_diag(True, 1)
        np.testing.assert_array_equal(trans, trans_desired)

    def test_pressures(self):
        p_desired = np.array(
            [
                3959.4367684962185,
                3878.3103054886546,
                3797.1838424810912,
                3716.057379473528,
            ]
        )
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
            model.wells[4]["pwf"], 3688.9926465271164, decimal=5
        )

    def test_rates(self):
        rates_desired = np.array(
            [2304.0239999999912, 0.0, 0.0, 0.0, -600.0000000000036, -1704.0240000000003]
        )
        model = create_model()
        model.solve(sparse=True, update=True, check_MB=True)
        np.testing.assert_almost_equal(model.rates[1], rates_desired, decimal=5)

    def test_simulation_run(self):
        model = create_model()
        model.run(30, True)
        model = create_model()
        model.run(30, False)


def create_model():
    grid = grids.Cartesian(
        nx=4, ny=1, nz=1, dx=300, dy=350, dz=40, phi=0.27, kx=270, dtype="double"
    )
    fluid = fluids.SinglePhase(mu=0.5, B=1, dtype="double")
    model = models.Numerical(grid, fluid, dtype="double", verbose=False)
    model.set_well(id=4, q=-600, s=1.5, r=3.5)
    model.set_boundaries({0: ("pressure", 4000), 5: ("gradient", -0.2)})
    return model


if __name__ == "__main__":
    unittest.main()
