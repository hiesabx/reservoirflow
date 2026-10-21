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
    model.set_boundaries({0: ("pressure", 4000), 5: ("gradient", -0.2)})
    model.compile(stype="numerical", method="FdmSolution", sparse=sparse)
    return model


def test_trans():
    trans_desired = np.array([28.4004, 28.4004, 28.4004, 28.4004, 28.4004])
    model = create_model(sparse=True)
    trans = model.get_cells_trans_diag(True, 1)
    np.testing.assert_array_equal(trans, trans_desired)


@SPARSE
def test_pressures(sparse):
    p_desired = np.array(
        [
            3959.4367684962185,
            3878.3103054886546,
            3797.1838424810912,
            3716.057379473528,
        ]
    )
    model = create_model(sparse=sparse)
    model.solve(update=True, check_MB=True)
    p = model.solution.pressures[-1, model.cells_id]
    np.testing.assert_almost_equal(p, p_desired, decimal=5)


def test_well():
    model = create_model(sparse=True)
    model.solve(update=True, check_MB=True)
    np.testing.assert_almost_equal(model.wells[4]["q"], -600, decimal=5)
    np.testing.assert_almost_equal(model.wells[4]["r_eq"], 64.53681120105021, decimal=5)
    np.testing.assert_almost_equal(model.wells[4]["G"], 11.08453575337366, decimal=5)
    np.testing.assert_almost_equal(model.wells[4]["pwf"], 3688.9926465271164, decimal=5)


def test_rates():
    rates_desired = np.array(
        [2304.0239999999912, 0.0, 0.0, 0.0, -600.0000000000036, -1704.0240000000003]
    )
    model = create_model(sparse=True)
    model.solve(update=True, check_MB=True)
    np.testing.assert_almost_equal(model.solution.rates[1], rates_desired, decimal=5)


@SPARSE
def test_simulation_run(sparse):
    model = create_model(sparse=sparse)
    model.run(nsteps=30)
