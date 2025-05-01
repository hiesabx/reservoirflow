import unittest

import numpy as np

from reservoirflow import fluids, grids, models


class TestApp(unittest.TestCase):
    def test_trans(self):
        trans_desired = np.array([28.4004, 28.4004, 28.4004, 28.4004, 28.4004])
        model = create_model(sparse=True)
        trans = model.get_cells_trans_diag(True, 1)
        np.testing.assert_array_equal(trans, trans_desired)

    def test_pressures(self):
        p_desired = np.array(
            [
                4000.0,
                3989.43753,
                3968.31261,
                3947.18768,
                3926.06276,
                np.nan,
            ]
        )
        # sparse:
        model = create_model(sparse=True)
        model.solve(update=True, check_MB=True)
        np.testing.assert_almost_equal(model.solution.pressures[1], p_desired, decimal=5)
        # dense:
        model = create_model(sparse=False)
        model.solve(update=True, check_MB=True)
        np.testing.assert_almost_equal(model.solution.pressures[1], p_desired, decimal=5)

    def test_well(self):
        model = create_model(sparse=True)
        model.solve(update=True, check_MB=True)
        np.testing.assert_almost_equal(model.wells[4]["q"], -599.95631, decimal=5)
        np.testing.assert_almost_equal(model.wells[4]["r_eq"], 64.53681, decimal=5)
        np.testing.assert_almost_equal(model.wells[4]["G"], 11.08453, decimal=5)
        np.testing.assert_almost_equal(model.wells[4]["pwf"], 3899, decimal=5)

    def test_rates(self):
        rates_desired = np.array(
            [599.9563191829617, 0.0, 0.0, 0.0, -599.9563191829004, 0.0]
        )
        # sparse:
        model = create_model(sparse=True)
        model.solve(update=True, check_MB=True)
        np.testing.assert_almost_equal(model.solution.rates[1], rates_desired, decimal=5)
        # dense:
        model = create_model(sparse=False)
        model.solve(update=True, check_MB=True)
        np.testing.assert_almost_equal(model.solution.rates[1], rates_desired, decimal=5)

    def test_simulation_run(self):
        model = create_model(sparse=False)
        model.run(nsteps=30)
        model = create_model(sparse=True)
        model.run(nsteps=30)


def create_model(sparse):
    grid = grids.RegularCartesian(
        nx=4, ny=1, nz=1, dx=300, dy=350, dz=40, phi=0.27, kx=270, dtype="double"
    )
    fluid = fluids.SinglePhase(mu=0.5, B=1, dtype="double")
    model = models.BlackOil(grid, fluid, dtype="double", verbose=False)
    model.set_well(cell_id=4, pwf=3899, s=1.5, r=3.5)
    model.set_boundaries({0: ("pressure", 4000), 5: ("rate", 0)})
    model.compile(stype="numerical", method="fdm", sparse=sparse)
    return model


def print_matrices():
    sparse = True
    model = create_model(sparse)
    for s in range(20):
        A, d = model.get_matrices_symb(sparse)
        A_, d_ = model.get_matrices_vect(sparse)
        if sparse:
            A, d = model.A.toarray(), model.d.toarray()
            A_, d_ = model.A_.toarray(), model.d_.toarray()
        print("step:", s)
        print(np.concatenate([A, A_, A - A_], axis=0))
        print(np.concatenate([d, d_, d - d_], axis=1))
        model.solve()
        print()


if __name__ == "__main__":
    unittest.main()
