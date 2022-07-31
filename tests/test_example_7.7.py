from venv import create
from openresim import grids, fluids, models
import numpy as np
import unittest


class TestApp(unittest.TestCase):
    def test_trans(self):
        T_x_desired = np.array([1.035, 1.035, 1.035])
        T_y_desired = np.array([1.3524, 1.3524, 1.3524])
        model = create_model()
        T_x = model.get_cells_T("x", False, False)
        T_y = model.get_cells_T("y", False, False)
        np.testing.assert_almost_equal(T_x, T_x_desired, decimal=5)
        np.testing.assert_almost_equal(T_y, T_y_desired, decimal=5)

    def test_pressures(self):
        p_desired = np.array([3772.36025, 3354.19841, 3267.38946, 3187.2711])
        model = create_model()
        p_sparse = model.solve(sparse=True, update=True, check_MB=True)
        boundaries = model.grid.get_boundaries("id", "list")
        np.testing.assert_almost_equal(p_sparse, p_desired, decimal=5)

        model = create_model()
        pressures_not_sparse = model.solve(sparse=False, update=True, check_MB=True)
        q_not_sparse = model.rates[1][boundaries]
        np.testing.assert_almost_equal(pressures_not_sparse, p_desired, decimal=5)

    def test_well(self):
        model = create_model()
        model.solve(sparse=True, update=True, check_MB=True)
        np.testing.assert_almost_equal(
            model.wells[6]["r_eq"], 58.52694766871065, decimal=5
        )
        np.testing.assert_almost_equal(
            model.wells[9]["r_eq"], 58.52694766871065, decimal=5
        )
        np.testing.assert_almost_equal(model.wells[6]["pwf"], 2000, decimal=5)
        np.testing.assert_almost_equal(model.wells[9]["q"], -600, decimal=5)
        np.testing.assert_almost_equal(
            model.wells[6]["G"], 4.768850278197406, decimal=5
        )
        np.testing.assert_almost_equal(
            model.wells[9]["G"], 4.768850278197406, decimal=5
        )

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

        model = create_model()
        model.solve(sparse=True, update=True, check_MB=True)
        boundaries = model.grid.get_boundaries("id", "list")
        q_sparse = model.rates[1][boundaries]
        np.testing.assert_almost_equal(q_sparse, rates_desired, decimal=5)

        model = create_model()
        model.solve(sparse=False, update=True, check_MB=True)
        q_not_sparse = model.rates[1][boundaries]
        np.testing.assert_almost_equal(q_not_sparse, rates_desired, decimal=5)

    def test_error(self):
        model = create_model()
        model.solve(sparse=True, update=True, check_MB=True)
        error_desired = -4.547e-12
        np.testing.assert_almost_equal(model.error, error_desired, decimal=5)


if __name__ == "__main__":

    def create_model():

        grid = grids.Cartesian(
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
        model = models.Model(grid, fluid, dtype="double", verbose=False)
        model.set_well(id=6, pwf=2000, s=0, r=3)
        model.set_well(id=9, q=-600, s=0, r=3)
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
        return model

    model = create_model()
    # for id_b in model.grid.get_boundaries("id", "tuple"):
    #     print(id_b, ":", end=" ")
    #     print(model.grid.get_cell_neighbors(id_b, None, False, "list"))

    # print(model.solve())
    # print(model.get_cells_T("x", False, False))
    # print(model.get_cells_T("y", False, False))
    # model.grid.show("id")
    # l, r = model.get_cell_eq(5)
    # model.init_matrices(False)
    # model.get_d(False)
    # model.solve(True)
    # model.get_cells_eq()
    # print(model.get_T("x", False))
    # print(model.get_T("y", False))
    # model.get_T("z", False)
    # print(model.RHS)
    # print(model.T)
    # print(model.get_d(False))
    # model.get_A()
    unittest.main()
