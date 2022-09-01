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
                3989.437537513856,
                3968.3126125415683,
                3947.1876875692815,
                3926.0627625969946,
            ]
        )

        model = create_model()
        p_sparse = model.solve(sparse=True, update=True, check_MB=True)
        np.testing.assert_almost_equal(p_sparse, p_desired, decimal=5)

        model = create_model()
        p_not_sparse = model.solve(sparse=False, update=True, check_MB=True)
        np.testing.assert_almost_equal(p_not_sparse, p_desired, decimal=5)

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
        model = create_model()
        model.solve(sparse=True, update=True, check_MB=True)
        np.testing.assert_almost_equal(model.rates[1], rates_desired, decimal=5)


def create_model():
    grid = grids.Cartesian(
        nx=4, ny=1, nz=1, dx=300, dy=350, dz=40, phi=0.27, kx=270, dtype="double"
    )
    fluid = fluids.SinglePhase(mu=0.5, B=1, dtype="double")
    model = models.Model(grid, fluid, dtype="double", verbose=False)
    model.set_well(id=4, pwf=3899, s=1.5, r=3.5)
    model.set_boundaries({0: ("pressure", 4000), 5: ("rate", 0)})
    return model


if __name__ == "__main__":
    # unittest.main()
    model = create_model()
    sparse = True
    for step in range(20):
        print("step:", step)
        A, d = model.init_matrices(sparse)
        A_, d_ = model.init_matrices_parallel(sparse)
        if sparse:
            A, d = A.toarray(), d.toarray()
            A_, d_ = A_.toarray(), d_.toarray()
        print("before solve:")
        print(
            np.concatenate(
                [
                    A,
                    A_,
                    np.subtract(A, A_),
                ],
                axis=0,
            )
        )
        print(
            np.concatenate(
                [
                    d,
                    d_,
                    np.subtract(d, d_),
                ],
                axis=1,
            )
        )
        # print(A)
        # print(A)
        model.solve(sparse)
        print()
