import unittest

import numpy as np
import pandas as pd

from reservoirflow import fluids, grids, models


class TestApp(unittest.TestCase):
    def test_data(self):
        df_desired = pd.read_csv(
            "tests/test_example_7_11.csv",
            index_col=0,
            dtype={"Step": "int32", "Time [days]": "int32"},
        )
        # vectorize, dense:
        model = create_model(sparse=False)
        model.run(nsteps=27, vectorize=True, threading=True)
        df = model.get_df(
            columns=["time", "cells_pressure", "wells"],
            boundary=False,
            scale=False,
            units=True,
            melt=False,
            save=False,
            drop_nan=True,
            drop_zero=True,
        )
        df = df.astype({"Time [days]": "int32"})
        # df.to_csv("tests/test_example_7_11_.csv")
        pd.testing.assert_frame_equal(df, df_desired)
        np.testing.assert_almost_equal(model.solution.error, 3.450062457943659e-11)
        # vectorize, sparse:
        model = create_model(sparse=True)
        model.run(nsteps=27, vectorize=True, threading=True)
        df = model.get_df(
            columns=["time", "cells_pressure", "wells"],
            boundary=False,
            scale=False,
            units=True,
            melt=False,
            save=False,
            drop_nan=True,
            drop_zero=True,
        )
        df = df.astype({"Time [days]": "int32"})
        pd.testing.assert_frame_equal(df, df_desired)
        np.testing.assert_almost_equal(model.solution.error, 3.450062457943659e-11)
        # symbolic, dense:
        model = create_model(sparse=False)
        model.run(nsteps=27, vectorize=False, threading=True)
        df = model.get_df(
            columns=["time", "cells_pressure", "wells"],
            boundary=False,
            scale=False,
            units=True,
            melt=False,
            save=False,
            drop_nan=True,
            drop_zero=True,
        )
        df = df.astype({"Time [days]": "int32"})
        # df.to_csv("tests/test_example_7_11_.csv")
        pd.testing.assert_frame_equal(df, df_desired)
        np.testing.assert_almost_equal(model.solution.error, 3.450062457943659e-11)
        # symbolic, sparse:
        model = create_model(sparse=True)
        model.run(nsteps=27, vectorize=False, threading=True)
        df = model.get_df(
            columns=["time", "cells_pressure", "wells"],
            boundary=False,
            scale=False,
            units=True,
            melt=False,
            save=False,
            drop_nan=True,
            drop_zero=True,
        )
        df = df.astype({"Time [days]": "int32"})
        pd.testing.assert_frame_equal(df, df_desired)
        np.testing.assert_almost_equal(model.solution.error, 3.450062457943659e-11)

    def test_trans(self):
        trans_desired = np.array(
            [12.81962, 14.04424, 15.71314, 21.08469, 20.16215, 14.8764]
        )
        model = create_model(sparse=True)
        Tx = model.get_cells_trans_diag(True, 1)
        np.testing.assert_almost_equal(Tx, trans_desired, 5)

    def test_RHS(self):
        RHS_desired = np.array(
            [
                1.87013,
                1.87013,
                1.135436,
                0.333952,
                1.113173,
                0.723562,
                0.723562,
            ]
        )
        model = create_model(sparse=True)
        np.testing.assert_almost_equal(model.RHS, RHS_desired, 5)

    def test_well(self):
        model = create_model(sparse=True)
        model.solve(update=True, check_MB=True)
        np.testing.assert_almost_equal(model.wells[4]["r_eq"], 75.392307, 5)
        np.testing.assert_almost_equal(model.wells[4]["q"], -400, 5)
        np.testing.assert_almost_equal(model.wells[4]["G"], 20.651804, 5)

    def test_simulation_run(self):
        # with self.assertWarns(Warning):
        model = create_model(sparse=False)
        model.run(nsteps=30, vectorize=True)
        model = create_model(sparse=True)
        model.run(nsteps=30, vectorize=True)
        model = create_model(sparse=False)
        model.run(nsteps=30, vectorize=False)
        model = create_model(sparse=True)
        model.run(nsteps=30, vectorize=False)


def create_model(sparse):
    dx = np.array([400, 400, 300, 150, 200, 250, 250])
    phi = np.array([0.21, 0.21, 0.17, 0.10, 0.25, 0.13, 0.13])
    kx = np.array([273, 273, 248, 127, 333, 198, 198])
    grid = grids.RegularCartesian(
        nx=5,
        ny=1,
        nz=1,
        dx=dx,
        dy=500,
        dz=50,
        phi=phi,
        kx=kx,
        comp=0,
        dtype="double",
        unify=False,
    )
    fluid = fluids.SinglePhase(mu=1.5, B=1, rho=50, comp=2.5 * 10**-5, dtype="double")
    model = models.BlackOil(grid, fluid, pi=3000, dt=5, dtype="double")
    model.set_well(cell_id=4, q=-400, pwf=1500, s=0, r=3)
    model.set_boundaries({0: ("rate", 0), 6: ("rate", 0)})
    model.compile(stype="numerical", method="fdm", sparse=sparse)
    return model


if __name__ == "__main__":
    unittest.main()
