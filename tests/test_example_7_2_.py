import unittest
import numpy as np
from reservoirflow import grids, fluids, models


class TestApp(unittest.TestCase):
    def test_trans(self):
        trans_desired = np.array([28.4004, 28.4004, 28.4004, 28.4004, 28.4004])
        model = create_model()
        trans = model.get_cells_T_diag(True, 1)
        np.testing.assert_array_equal(trans, trans_desired)

    def test_pressures(self):
        p_desired = np.array([4000.0, 3989.44, 3968.31, 3947.19, 3926.06, np.nan])
        # sparse:
        model = create_model()
        model.solve(sparse=True, update=True, check_MB=True)
        np.testing.assert_almost_equal(model.pressures[1], p_desired, decimal=2)
        # dense:
        model = create_model()
        model.solve(sparse=False, update=True, check_MB=True)
        np.testing.assert_almost_equal(model.pressures[1], p_desired, decimal=2)

    def test_well(self):
        model = create_model()
        model.solve(sparse=True, update=True, check_MB=True)
        np.testing.assert_almost_equal(model.wells[4]["q"], -599.9, decimal=1)
        np.testing.assert_almost_equal(model.wells[4]["r_eq"], 64.537, decimal=3)
        np.testing.assert_almost_equal(model.wells[4]["G"], 11.0845, decimal=4)
        np.testing.assert_almost_equal(model.wells[4]["pwf"], 3899, decimal=1)

    def test_rates(self):
        rates_desired = np.array([599.9, 0.0, 0.0, 0.0, -599.9, 0.0])
        model = create_model()
        model.solve(sparse=True, update=True, check_MB=True)
        np.testing.assert_almost_equal(model.rates[1], rates_desired, decimal=1)

    def test_simulation_run(self):
        model = create_model()
        model.run(30, False)
        model = create_model()
        model.run(30, True)


def create_model():
    grid = grids.Cartesian(
        nx=4, ny=1, nz=1, dx=300, dy=350, dz=40, phi=0.27, kx=270, dtype="double"
    )
    fluid = fluids.SinglePhase(mu=0.5, B=1, dtype="double")
    model = models.BlackOil(grid, fluid, dtype="double", verbose=False)
    model.set_well(id=4, q=-600, s=1.5, r=3.5)
    model.set_boundaries({0: ("pressure", 4000)})
    return model


if __name__ == "__main__":
    unittest.main()
