from openresim import grids
import numpy as np
import unittest


class Test_2D_XY(unittest.TestCase):
    """_summary_

    ToDo:
    ----------
    test 2d
    """

    def test_cell_id(self):
        id = 6
        coords = grid_2d_xy.get_cell_coords(id)
        self.assertEqual(coords, (1, 1, 0))
        neighbors = grid_2d_xy.get_cell_neighbors(id=id, boundary=False, fmt="list")
        self.assertEqual(neighbors, [7, 11])
        neighbors = grid_2d_xy.get_cell_neighbors(id=id, boundary=False, fmt="dict")
        self.assertEqual(neighbors, {"x": [7], "y": [11], "z": []})
        boundaries = grid_2d_xy.get_cell_boundaries(id=id, fmt="list")
        self.assertEqual(boundaries, [1, 5])
        boundaries = grid_2d_xy.get_cell_boundaries(id=id, fmt="dict")
        self.assertEqual(boundaries, {"x": [5], "y": [1], "z": []})

    def test_cell_coords(self):
        coords = (1, 1, 0)
        id = grid_2d_xy.get_cell_id(coords)
        self.assertEqual(id, 6)
        neighbors = grid_2d_xy.get_cell_neighbors(
            coords=coords, boundary=False, fmt="list"
        )
        self.assertEqual(neighbors, [(2, 1, 0), (1, 2, 0)])
        neighbors = grid_2d_xy.get_cell_neighbors(
            coords=coords, boundary=False, fmt="dict"
        )
        self.assertEqual(neighbors, {"x": [(2, 1, 0)], "y": [(1, 2, 0)], "z": []})

        boundaries = grid_2d_xy.get_cell_boundaries(coords=coords, fmt="list")
        self.assertEqual(boundaries, [(1, 0, 0), (0, 1, 0)])
        boundaries = grid_2d_xy.get_cell_boundaries(coords=coords, fmt="dict")
        self.assertEqual(boundaries, {"x": [(0, 1, 0)], "y": [(1, 0, 0)], "z": []})


class Test_3D(unittest.TestCase):
    """_summary_

    ToDo:
    ----------
    test 2d
    """

    def test_cell_id_3d_xyz(self):
        id = 31
        coords = grid_3d.get_cell_coords(id)
        self.assertEqual(coords, (1, 1, 1))
        neighbors = grid_3d.get_cell_neighbors(id=id, boundary=False, fmt="list")
        self.assertEqual(neighbors, [32, 36, 56])
        neighbors = grid_3d.get_cell_neighbors(id=id, boundary=False, fmt="dict")
        self.assertEqual(neighbors, {"x": [32], "y": [36], "z": [56]})
        boundaries = grid_3d.get_cell_boundaries(id=id, fmt="list")
        self.assertEqual(boundaries, [26, 6, 30])
        boundaries = grid_3d.get_cell_boundaries(id=id, fmt="dict")
        self.assertEqual(boundaries, {"x": [30], "y": [26], "z": [6]})

    def test_cell_coords_3d_xyz(self):
        coords = (1, 1, 1)
        id = grid_3d.get_cell_id(coords)
        self.assertEqual(id, 31)
        neighbors = grid_3d.get_cell_neighbors(
            coords=coords, boundary=False, fmt="list"
        )
        self.assertEqual(neighbors, [(2, 1, 1), (1, 2, 1), (1, 1, 2)])
        neighbors = grid_3d.get_cell_neighbors(
            coords=coords, boundary=False, fmt="dict"
        )
        self.assertEqual(
            neighbors, {"x": [(2, 1, 1)], "y": [(1, 2, 1)], "z": [(1, 1, 2)]}
        )
        boundaries = grid_3d.get_cell_boundaries(coords=coords, fmt="list")
        self.assertEqual(boundaries, [(0, 1, 1), (1, 1, 0), (1, 0, 1)])
        boundaries = grid_3d.get_cell_boundaries(coords=coords, fmt="dict")
        self.assertEqual(
            boundaries, {"x": [(0, 1, 1)], "y": [(1, 0, 1)], "z": [(1, 1, 0)]}
        )


if __name__ == "__main__":
    grid_2d_xy = grids.CartGrid(
        nx=3,
        ny=3,
        nz=1,
        dx=10,
        dy=10,
        dz=10,
        phi=0.27,
        kx=270,
        ky=270,
        kz=270,
        verbose=False,
    )
    grid_2d_xz = grids.CartGrid(
        nx=3,
        ny=1,
        nz=3,
        dx=10,
        dy=10,
        dz=10,
        phi=0.27,
        kx=270,
        ky=270,
        kz=270,
        verbose=False,
    )
    grid_2d_yz = grids.CartGrid(
        nx=1,
        ny=3,
        nz=3,
        dx=10,
        dy=10,
        dz=10,
        phi=0.27,
        kx=270,
        ky=270,
        kz=270,
        verbose=False,
    )
    grid_3d = grids.CartGrid(
        nx=3,
        ny=3,
        nz=3,
        dx=10,
        dy=10,
        dz=10,
        phi=0.27,
        kx=270,
        ky=270,
        kz=270,
        verbose=False,
    )
    unittest.main()
