import numpy as np
import pytest

from reservoirflow import fluids, grids, models

SPARSE = pytest.mark.parametrize("sparse", [False, True], ids=["dense", "sparse"])
VECTORIZE = pytest.mark.parametrize(
    "vectorize", [True, False], ids=["vectorized", "symbolic"]
)


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
        dtype="double",
    )
    fluid = fluids.SinglePhaseFluid(
        mu=0.5,
        B=1,
        dtype="double",
    )
    model = models.BlackOilModel(
        grid,
        fluid,
        dtype="double",
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

    model.compile(stype="numerical", method="FdmSolution", sparse=sparse)

    return model


def test_trans():
    trans_desired = np.array([28.4004, 28.4004, 28.4004, 28.4004, 28.4004])
    model = create_model(sparse=False)
    trans = model.get_cells_trans_diag(True, 1)
    np.testing.assert_array_equal(trans, trans_desired)


@SPARSE
@VECTORIZE
def test_pressures(sparse, vectorize):
    p_desired = np.array(
        [4000.0, 3989.43676, 3968.310305, 3947.18384, 3926.05737, np.nan]
    )
    model = create_model(sparse=sparse)
    model.solve(vectorize=vectorize, update=True, check_MB=True)
    np.testing.assert_almost_equal(model.solution.pressures[1], p_desired, decimal=5)


def test_well():
    model = create_model(sparse=True)
    model.solve(update=True, check_MB=True)
    np.testing.assert_almost_equal(model.wells[4]["q"], -600, decimal=5)
    np.testing.assert_almost_equal(model.wells[4]["r_eq"], 64.53681, decimal=5)
    np.testing.assert_almost_equal(model.wells[4]["G"], 11.08453, decimal=5)
    np.testing.assert_almost_equal(model.wells[4]["pwf"], 3898.99264, decimal=5)


@SPARSE
@VECTORIZE
def test_rates(sparse, vectorize):
    rates_desired = np.array(
        [600.0000000000169, 0.0, 0.0, 0.0, -600.0000000000036, 0.0]
    )
    model = create_model(sparse=sparse)
    model.solve(vectorize=vectorize, update=True, check_MB=True)
    np.testing.assert_almost_equal(model.solution.rates[1], rates_desired, decimal=5)


@SPARSE
def test_simulation_run(sparse):
    model = create_model(sparse=sparse)
    model.run(nsteps=30)
