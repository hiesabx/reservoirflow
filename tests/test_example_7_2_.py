import numpy as np
import pytest

from reservoirflow import fluids, grids, models

SPARSE = pytest.mark.parametrize("sparse", [False, True], ids=["dense", "sparse"])


def create_model(sparse):
    grid = grids.RegularCartesianGrid(
        nx=4, ny=1, nz=1, dx=300, dy=350, dz=40, phi=0.27, kx=270, dtype="double"
    )
    fluid = fluids.SinglePhaseFluid(mu=0.5, B=1, dtype="double")
    model = models.BlackOilModel(grid, fluid, dtype="double", verbose=False)
    model.set_well(cell_id=4, q=-600, s=1.5, r=3.5)
    model.set_boundaries({0: ("pressure", 4000)})
    model.compile(stype="numerical", method="FdmSolution", sparse=sparse)
    return model


def test_trans():
    trans_desired = np.array([28.4004, 28.4004, 28.4004, 28.4004, 28.4004])
    model = create_model(sparse=True)
    trans = model.get_cells_trans_diag(True, 1)
    np.testing.assert_array_equal(trans, trans_desired)


@SPARSE
def test_pressures(sparse):
    p_desired = np.array([4000.0, 3989.44, 3968.31, 3947.19, 3926.06, np.nan])
    model = create_model(sparse=sparse)
    model.solve(update=True, check_MB=True)
    np.testing.assert_almost_equal(model.solution.pressures[1], p_desired, decimal=2)


def test_well():
    model = create_model(sparse=True)
    model.solve(update=True, check_MB=True)
    np.testing.assert_almost_equal(model.wells[4]["q"], -599.9, decimal=1)
    np.testing.assert_almost_equal(model.wells[4]["r_eq"], 64.537, decimal=3)
    np.testing.assert_almost_equal(model.wells[4]["G"], 11.0845, decimal=4)
    np.testing.assert_almost_equal(model.wells[4]["pwf"], 3899, decimal=1)


def test_rates():
    rates_desired = np.array([599.9, 0.0, 0.0, 0.0, -599.9, 0.0])
    model = create_model(sparse=True)
    model.solve(update=True, check_MB=True)
    np.testing.assert_almost_equal(model.solution.rates[1], rates_desired, decimal=1)


@SPARSE
def test_simulation_run(sparse):
    model = create_model(sparse=sparse)
    model.run(nsteps=30)
