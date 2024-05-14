import unittest

import numpy as np

from reservoirflow import fluids, grids, models


class TestApp(unittest.TestCase):
    def test_trans(self):
        trans_desired = np.array([28.4004, 28.4004, 28.4004, 28.4004, 28.4004])
        model = create_model(sparse=False)
        trans = model.get_cells_trans_diag(True, 1)
        np.testing.assert_array_equal(trans, trans_desired)

    def test_pressures(self):
        p_desired = np.array(
            [4000.0, 3989.43676, 3968.310305, 3947.18384, 3926.05737, np.nan]
        )
        # vectorize, sparse:
        model = create_model(sparse=True)
        model.solve(vectorize=True, update=True, check_MB=True)
        np.testing.assert_almost_equal(model.pressures[1], p_desired, decimal=5)
        # vectorize, dense:
        model = create_model(sparse=False)
        model.solve(vectorize=True, update=True, check_MB=True)
        np.testing.assert_almost_equal(model.pressures[1], p_desired, decimal=5)
        # symbolic, sparse:
        model = create_model(sparse=True)
        model.solve(vectorize=False, update=True, check_MB=True)
        np.testing.assert_almost_equal(model.pressures[1], p_desired, decimal=5)
        # symbolic, dense:
        model = create_model(sparse=False)
        model.solve(vectorize=False, update=True, check_MB=True)
        np.testing.assert_almost_equal(model.pressures[1], p_desired, decimal=5)

    def test_well(self):
        model = create_model(sparse=True)
        model.solve(update=True, check_MB=True)
        np.testing.assert_almost_equal(model.wells[4]["q"], -600, decimal=5)
        np.testing.assert_almost_equal(model.wells[4]["r_eq"], 64.53681, decimal=5)
        np.testing.assert_almost_equal(model.wells[4]["G"], 11.08453, decimal=5)
        np.testing.assert_almost_equal(model.wells[4]["pwf"], 3898.99264, decimal=5)

    def test_rates(self):
        rates_desired = np.array(
            [600.0000000000169, 0.0, 0.0, 0.0, -600.0000000000036, 0.0]
        )
        # vectorize, sparse:
        model = create_model(sparse=True)
        model.solve(vectorize=True, update=True, check_MB=True)
        np.testing.assert_almost_equal(model.rates[1], rates_desired, decimal=5)
        # vectorize, dense:
        model = create_model(sparse=False)
        model.solve(vectorize=True, update=True, check_MB=True)
        np.testing.assert_almost_equal(model.rates[1], rates_desired, decimal=5)
        # symbolic, sparse:
        model = create_model(sparse=True)
        model.solve(vectorize=False, update=True, check_MB=True)
        np.testing.assert_almost_equal(model.rates[1], rates_desired, decimal=5)
        # symbolic, dense:
        model = create_model(sparse=False)
        model.solve(vectorize=False, update=True, check_MB=True)
        np.testing.assert_almost_equal(model.rates[1], rates_desired, decimal=5)

    def test_simulation_run(self):
        model = create_model(sparse=False)
        model.run(nsteps=30)
        model = create_model(sparse=True)
        model.run(nsteps=30)


def create_model(sparse):
    grid = grids.RegularCartesian(
        nx=4,
        ny=1,
        nz=1,
        dx=300,
        dy=350,
        dz=40,
        phi=0.27,
        kx=270,
        dtype="double",
    )
    fluid = fluids.SinglePhase(
        mu=0.5,
        B=1,
        dtype="double",
    )
    model = models.BlackOil(
        grid,
        fluid,
        dtype="double",
        sparse=sparse,
        verbose=False,
    )
    model.set_well(
        cell_id=4,
        q=-600,
        s=1.5,
        r=3.5,
    )
    model.set_boundaries(
        {0: ("pressure", 4000), 5: ("rate", 0)},
    )

    model.compile(stype="numerical", method="fdm")

    return model


if __name__ == "__main__":
    unittest.main()
