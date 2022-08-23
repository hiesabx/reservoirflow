from venv import create
from openresim import grids, fluids, models
import numpy as np
import pandas as pd
import unittest


class TestApp(unittest.TestCase):
    def test_data(self):
        df_desired = pd.read_csv(
            "tests/test_example_7.9.csv",
            index_col=0,
            dtype={None: "float64", "Time [Days]": "int32"},
        )
        model = create_model()
        model.run(nsteps=10, sparse=False)
        df = model.data(*6 * [True], False)
        pd.testing.assert_frame_equal(df, df_desired)
        np.testing.assert_almost_equal(model.error, 3.320340669077382e-10)
        self.assertLess(model.ctime, 5)

    def test_trans(self):
        trans_desired = np.array(
            [
                28.4004,
                28.4004,
                28.4004,
                28.4004,
                28.4004,
            ]
        )
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

    def test_error(self):
        error_desired = 0.99891
        cumulative_error_desired = 0.999999

        model = create_model()
        model.solve(sparse=False, update=True, check_MB=True)
        np.testing.assert_almost_equal(1 - model.error, error_desired, decimal=3)
        np.testing.assert_almost_equal(
            model.cumulative_error, cumulative_error_desired, decimal=5
        )

        model_sp = create_model()
        model_sp.solve(sparse=True, update=True, check_MB=True)
        np.testing.assert_almost_equal(1 - model_sp.error, error_desired, decimal=3)
        np.testing.assert_almost_equal(
            model_sp.cumulative_error, cumulative_error_desired, decimal=5
        )

        error_desired = 0.99891
        cumulative_error_desired = 0.499999

        model.solve(sparse=False, update=True, check_MB=True)
        np.testing.assert_almost_equal(1 - model.error, error_desired, decimal=3)
        np.testing.assert_almost_equal(
            model.cumulative_error, cumulative_error_desired, decimal=5
        )

        model_sp.solve(sparse=True, update=True, check_MB=True)
        np.testing.assert_almost_equal(1 - model_sp.error, error_desired, decimal=3)
        np.testing.assert_almost_equal(
            model_sp.cumulative_error, cumulative_error_desired, decimal=5
        )

    def test_well(self):
        model = create_model()
        model.solve(sparse=True, update=True, check_MB=True)
        np.testing.assert_almost_equal(model.wells[4]["r_eq"], 64.536811, decimal=5)
        np.testing.assert_almost_equal(model.wells[4]["q"], -600, decimal=5)
        np.testing.assert_almost_equal(model.wells[4]["G"], 11.084535, decimal=5)
        np.testing.assert_almost_equal(model.wells[4]["pwf"], 3922.034614, decimal=5)


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
