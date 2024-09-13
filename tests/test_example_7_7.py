import unittest

import numpy as np
import pandas as pd

from reservoirflow import fluids, grids, models


class TestApp(unittest.TestCase):
    def test_data(self):
        df_desired = pd.read_csv(
            "tests/test_example_7_7.csv",
            index_col=0,
            dtype={"Step": "int32", "Time [days]": "int32"},
        )
        df_desired.index = df_desired.index.astype(int)

        model = create_model(sparse=True)
        model.solve(threading=True, vectorize=True)
        df = model.get_df(
            columns=["time", "cells_rate", "cells_pressure", "wells"],
            boundary=True,
            units=True,
            melt=False,
            scale=False,
            save=False,
            drop_nan=False,
            drop_zero=False,
        )
        df.index = df.index.astype(int)
        # df.to_csv("tests/test_example_7_7_.csv")
        pd.testing.assert_frame_equal(df, df_desired)
        np.testing.assert_almost_equal(model.solution.error, 3.320340669077382e-10)
        self.assertLess(model.solution.ctime, 5)

    def test_trans(self):
        T_x_desired = np.array([1.035, 0, 1.035])
        T_y_desired = np.array([1.3524, 1.3524])
        model = create_model(sparse=True)
        T_x = model.get_cells_trans_diag(False, 1)
        T_y = model.get_cells_trans_diag(False, 2)
        # T_x = model.get_cells_T_vect("x", False, False)
        # T_y = model.get_cells_T_vect("y", False, False)
        np.testing.assert_almost_equal(T_x, T_x_desired, decimal=5)
        np.testing.assert_almost_equal(T_y, T_y_desired, decimal=5)

    def test_pressures(self):
        p_desired = np.array([3772.36025, 3354.19841, 3267.38946, 3187.2711])
        # sparse:
        model = create_model(sparse=True)
        model.solve(update=True, check_MB=True)
        p_sparse = model.solution.pressures[-1, model.cells_id]
        np.testing.assert_almost_equal(p_sparse, p_desired, decimal=5)
        # dense:
        model = create_model(sparse=False)
        model.solve(update=True, check_MB=True)
        p_not_sparse = model.solution.pressures[-1, model.cells_id]
        np.testing.assert_almost_equal(p_not_sparse, p_desired, decimal=5)

    def test_well(self):
        model = create_model(sparse=True)
        model.solve(update=True, check_MB=True)
        np.testing.assert_almost_equal(model.wells[6]["r_eq"], 58.527, decimal=3)
        np.testing.assert_almost_equal(model.wells[9]["r_eq"], 58.527, decimal=3)
        np.testing.assert_almost_equal(model.wells[6]["G"], 4.7688, decimal=4)
        np.testing.assert_almost_equal(model.wells[9]["G"], 4.7688, decimal=4)
        np.testing.assert_almost_equal(model.wells[6]["pwf"], 2000)
        np.testing.assert_almost_equal(model.wells[9]["q"], -600)

    def test_rates(self):
        rates_desired = np.array(
            [
                0.0,  # 0
                615.72,  # 1
                1746.76413,  # 2
                0.0,  # 3
                500.0,  # 4
                -108.675,  # 7
                0.0,  # 8
                -108.675,  # 11
                0.0,  # 12
                0.0,  # 13
                -200.0,  # 14
                0.0,  # 15
            ]
        )

        model = create_model(sparse=True)
        boundaries = model.grid.get_boundaries("id", "list")

        model.solve(update=True, check_MB=True)
        q_sparse = model.solution.rates[1][boundaries]
        np.testing.assert_almost_equal(q_sparse, rates_desired, decimal=5)

        model = create_model(sparse=False)
        model.solve(update=True, check_MB=True)
        q_not_sparse = model.solution.rates[1][boundaries]
        np.testing.assert_almost_equal(q_not_sparse, rates_desired, decimal=5)

    def test_error(self):
        model = create_model(sparse=True)
        model.solve(update=True, check_MB=True)
        error_desired = -4.547e-12
        np.testing.assert_almost_equal(model.solution.error, error_desired, decimal=5)

    def test_simulation_run(self):
        model = create_model(sparse=False)
        model.run(nsteps=30)
        model = create_model(sparse=True)
        model.run(nsteps=30)


def create_model(sparse):
    grid = grids.RegularCartesian(
        nx=2,
        ny=2,
        nz=1,
        dx=350,
        dy=250,
        dz=30,
        phi=0.27,
        kx=150,
        ky=100,
        dtype="double",
        unify=False,
    )
    fluid = fluids.SinglePhase(
        mu=3.5,
        B=1,
        rho=50,
        dtype="double",
    )
    model = models.BlackOil(grid, fluid, dtype="double", verbose=False)
    model.set_well(cell_id=6, pwf=2000, s=0, r=3)
    model.set_well(cell_id=9, q=-600, s=0, r=3)
    model.set_boundaries(
        {
            1: ("pressure", 4000),
            2: ("pressure", 4000),
            4: ("rate", 500),
            7: ("gradient", -0.3),
            11: ("gradient", -0.3),
            14: ("rate", -200),
        }
    )
    model.compile(stype="numerical", method="fdm", sparse=sparse)
    return model


if __name__ == "__main__":
    unittest.main()
