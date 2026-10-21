# ToDo: minimize MBE at late time steps.
# model = create_model(sparse=False)
# model.run(30, True, True, True, print_arrays=True)
# model.get_df(True, True, save=True)
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from reservoirflow import fluids, grids, models

TESTS_DIR = Path(__file__).parent
SPARSE = pytest.mark.parametrize("sparse", [False, True], ids=["dense", "sparse"])


def create_model(sparse):
    grid = grids.RegularCartesianGrid(
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
    fluid = fluids.SinglePhaseFluid(
        mu=0.5, B=1, rho=50, comp=1 * 10**-5, dtype="double"
    )
    model = models.BlackOilModel(
        grid, fluid, pi=4000, dt=1, dtype="double", verbose=False
    )
    model.set_well(cell_id=4, q=-600, s=1.5, r=3.5)
    model.set_boundaries({0: ("pressure", 4000), 5: ("rate", 0)})
    model.compile(stype="numerical", method="FdmSolution", sparse=sparse)
    return model


@pytest.fixture(scope="module")
def df_desired():
    return pd.read_csv(
        TESTS_DIR / "test_example_7_9.csv",
        index_col=0,
        dtype={"Step": "int32", "Time [days]": "int32"},
    )


def test_data(df_desired):
    model = create_model(sparse=False)
    model.run(nsteps=10)
    df = model.get_df(
        columns=["time", "cells_rate", "cells_pressure", "wells"],
        boundary=True,
        scale=False,
        units=True,
        melt=False,
        save=False,
        drop_nan=False,
        drop_zero=False,
    )
    df = df.astype({"Time [days]": "int32"})
    # df.to_csv(TESTS_DIR / "test_example_7_9_.csv")
    pd.testing.assert_frame_equal(df, df_desired)
    np.testing.assert_almost_equal(model.solution.tstep_error, 3.320340669077382e-10)


def test_trans():
    trans_desired = np.array([28.4004, 28.4004, 28.4004, 28.4004, 28.4004])
    model = create_model(sparse=False)
    Tx = model.get_cells_trans_diag(True, 1)
    np.testing.assert_array_equal(Tx, trans_desired)


def test_RHS():
    RHS_desired = np.array([2.221714417615699] * 6)
    model = create_model(sparse=True)
    np.testing.assert_array_equal(model.RHS, RHS_desired)


@SPARSE
def test_error(sparse):
    model = create_model(sparse=sparse)

    # first time step:
    model.solve(update=True, check_MB=True)
    np.testing.assert_almost_equal(1 - model.solution.tstep_error, 0.99891, decimal=3)
    np.testing.assert_almost_equal(model.solution.ctime_error, 0.999999, decimal=5)

    # second time step:
    model.solve(update=True, check_MB=True)
    np.testing.assert_almost_equal(1 - model.solution.tstep_error, 0.99891, decimal=3)
    np.testing.assert_almost_equal(model.solution.ctime_error, 0.499999, decimal=5)


def test_well():
    model = create_model(sparse=True)
    model.solve(update=True, check_MB=True)
    np.testing.assert_almost_equal(model.wells[4]["r_eq"], 64.536811, decimal=5)
    np.testing.assert_almost_equal(model.wells[4]["q"], -600, decimal=5)
    np.testing.assert_almost_equal(model.wells[4]["G"], 11.084535, decimal=5)
    np.testing.assert_almost_equal(model.wells[4]["pwf"], 3922.034614, decimal=5)


@SPARSE
def test_simulation_run(sparse):
    model = create_model(sparse=sparse)
    model.run(nsteps=30)
