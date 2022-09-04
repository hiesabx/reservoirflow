from openresim import grids, fluids, models
import numpy as np
import unittest


class TestApp(unittest.TestCase):
    def test_trans(self):
        trans_desired = np.array([28.4004, 28.4004, 28.4004, 28.4004, 28.4004])
        model = create_model()
        np.testing.assert_array_equal(model.T["x"], trans_desired)

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
        model = create_model()
        model.solve(sparse=True, update=True, check_MB=True)
        np.testing.assert_almost_equal(model.pressures[1], p_desired, decimal=5)
        # dense:
        model = create_model()
        model.solve(sparse=False, update=True, check_MB=True)
        np.testing.assert_almost_equal(model.pressures[1], p_desired, decimal=5)

    def test_well(self):
        model = create_model()
        model.solve(sparse=True, update=True, check_MB=True)
        np.testing.assert_almost_equal(model.wells[4]["q"], -599.95631, decimal=5)
        np.testing.assert_almost_equal(model.wells[4]["r_eq"], 64.53681, decimal=5)
        np.testing.assert_almost_equal(model.wells[4]["G"], 11.08453, decimal=5)
        np.testing.assert_almost_equal(model.wells[4]["pwf"], 3899, decimal=5)

    def test_rates(self):
        rates_desired = np.array(
            [599.9563191829617, 0.0, 0.0, 0.0, -599.9563191829004, 0.0]
        )
        # sparse:
        model = create_model()
        model.solve(sparse=True, update=True, check_MB=True)
        np.testing.assert_almost_equal(model.rates[1], rates_desired, decimal=5)
        # dense:
        model = create_model()
        model.solve(sparse=False, update=True, check_MB=True)
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
    model = models.Model(grid, fluid, dtype="double", verbose=False)
    model.set_well(id=4, pwf=3899, s=1.5, r=3.5)
    model.set_boundaries({0: ("pressure", 4000), 5: ("rate", 0)})
    return model


def print_matrices():
    sparse = True
    model = create_model()
    for s in range(20):
        A, d = model.init_matrices(sparse)
        A_, d_ = model.get_matrices_vectorized(sparse)
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
