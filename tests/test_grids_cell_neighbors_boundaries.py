from reservoirflow import grids
import numpy as np
import unittest


class Test_1D_X(unittest.TestCase):
    def test_cell_id(self):
        id = 1

        coords = grid_1d_x.get_cell_coords(id)
        self.assertEqual(coords, (1, 0, 0))

        neighbors = grid_1d_x.get_cell_neighbors(id=id, boundary=False, fmt="list")
        self.assertEqual(neighbors, [2])

        neighbors = grid_1d_x.get_cell_neighbors(id=id, boundary=False, fmt="dict")
        self.assertEqual(neighbors, {"x": [2], "y": [], "z": []})

        boundaries = grid_1d_x.get_cell_boundaries(id=id, fmt="list")
        self.assertEqual(boundaries, [0])

        boundaries = grid_1d_x.get_cell_boundaries(id=id, fmt="dict")
        self.assertEqual(boundaries, {"x": [0], "y": [], "z": []})

    def test_cell_coords(self):
        coords = (1, 0, 0)

        id = grid_1d_x.get_cell_id(coords)
        self.assertEqual(id, 1)

        neighbors = grid_1d_x.get_cell_neighbors(
            coords=coords, boundary=False, fmt="list"
        )
        self.assertEqual(neighbors, [(2, 0, 0)])

        neighbors = grid_1d_x.get_cell_neighbors(
            coords=coords, boundary=False, fmt="dict"
        )
        self.assertEqual(neighbors, {"x": [(2, 0, 0)], "y": [], "z": []})

        boundaries = grid_1d_x.get_cell_boundaries(coords=coords, fmt="list")
        self.assertEqual(boundaries, [(0, 0, 0)])

        boundaries = grid_1d_x.get_cell_boundaries(coords=coords, fmt="dict")
        self.assertEqual(boundaries, {"x": [(0, 0, 0)], "y": [], "z": []})


class Test_1D_Y(unittest.TestCase):
    def test_cell_id(self):
        id = 1

        coords = grid_1d_y.get_cell_coords(id)
        self.assertEqual(coords, (0, 1, 0))

        neighbors = grid_1d_y.get_cell_neighbors(id=id, boundary=False, fmt="list")
        self.assertEqual(neighbors, [2])

        neighbors = grid_1d_y.get_cell_neighbors(id=id, boundary=False, fmt="dict")
        self.assertEqual(neighbors, {"x": [], "y": [2], "z": []})

        boundaries = grid_1d_y.get_cell_boundaries(id=id, fmt="list")
        self.assertEqual(boundaries, [0])

        boundaries = grid_1d_y.get_cell_boundaries(id=id, fmt="dict")
        self.assertEqual(boundaries, {"x": [], "y": [0], "z": []})

    def test_cell_coords(self):
        coords = (0, 1, 0)

        id = grid_1d_y.get_cell_id(coords)
        self.assertEqual(id, 1)

        neighbors = grid_1d_y.get_cell_neighbors(
            coords=coords, boundary=False, fmt="list"
        )
        self.assertEqual(neighbors, [(0, 2, 0)])

        neighbors = grid_1d_y.get_cell_neighbors(
            coords=coords, boundary=False, fmt="dict"
        )
        self.assertEqual(neighbors, {"x": [], "y": [(0, 2, 0)], "z": []})

        boundaries = grid_1d_y.get_cell_boundaries(coords=coords, fmt="list")
        self.assertEqual(boundaries, [(0, 0, 0)])

        boundaries = grid_1d_y.get_cell_boundaries(coords=coords, fmt="dict")
        self.assertEqual(boundaries, {"x": [], "y": [(0, 0, 0)], "z": []})


class Test_1D_Z(unittest.TestCase):
    def test_cell_id(self):
        id = 1

        coords = grid_1d_z.get_cell_coords(id)
        self.assertEqual(coords, (0, 0, 1))

        neighbors = grid_1d_z.get_cell_neighbors(id=id, boundary=False, fmt="list")
        self.assertEqual(neighbors, [2])

        neighbors = grid_1d_z.get_cell_neighbors(id=id, boundary=False, fmt="dict")
        self.assertEqual(neighbors, {"x": [], "y": [], "z": [2]})

        boundaries = grid_1d_z.get_cell_boundaries(id=id, fmt="list")
        self.assertEqual(boundaries, [0])

        boundaries = grid_1d_z.get_cell_boundaries(id=id, fmt="dict")
        self.assertEqual(boundaries, {"x": [], "y": [], "z": [0]})

    def test_cell_coords(self):
        coords = (0, 0, 1)

        id = grid_1d_z.get_cell_id(coords)
        self.assertEqual(id, 1)

        neighbors = grid_1d_z.get_cell_neighbors(
            coords=coords, boundary=False, fmt="list"
        )
        self.assertEqual(neighbors, [(0, 0, 2)])

        neighbors = grid_1d_z.get_cell_neighbors(
            coords=coords, boundary=False, fmt="dict"
        )
        self.assertEqual(neighbors, {"x": [], "y": [], "z": [(0, 0, 2)]})

        boundaries = grid_1d_z.get_cell_boundaries(coords=coords, fmt="list")
        self.assertEqual(boundaries, [(0, 0, 0)])

        boundaries = grid_1d_z.get_cell_boundaries(coords=coords, fmt="dict")
        self.assertEqual(boundaries, {"x": [], "y": [], "z": [(0, 0, 0)]})


class Test_2D_XY(unittest.TestCase):
    def test_cell_id(self):
        id = 6

        coords = grid_2d_xy.get_cell_coords(id)
        self.assertEqual(coords, (1, 1, 0))

        neighbors = grid_2d_xy.get_cell_neighbors(id=id, boundary=False, fmt="list")
        self.assertEqual(neighbors, [7, 11])

        neighbors = grid_2d_xy.get_cell_neighbors(id=id, boundary=False, fmt="dict")
        self.assertEqual(neighbors, {"x": [7], "y": [11], "z": []})

        boundaries = grid_2d_xy.get_cell_boundaries(id=id, fmt="list")
        self.assertEqual(boundaries, [5, 1])

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
        self.assertEqual(boundaries, [(0, 1, 0), (1, 0, 0)])

        boundaries = grid_2d_xy.get_cell_boundaries(coords=coords, fmt="dict")
        self.assertEqual(boundaries, {"x": [(0, 1, 0)], "y": [(1, 0, 0)], "z": []})


class Test_2D_XZ(unittest.TestCase):
    def test_cell_id(self):
        id = 6

        coords = grid_2d_xz.get_cell_coords(id)
        self.assertEqual(coords, (1, 0, 1))

        neighbors = grid_2d_xz.get_cell_neighbors(id=id, boundary=False, fmt="list")
        self.assertEqual(neighbors, [7, 11])

        neighbors = grid_2d_xz.get_cell_neighbors(id=id, boundary=False, fmt="dict")
        self.assertEqual(neighbors, {"x": [7], "y": [], "z": [11]})

        boundaries = grid_2d_xz.get_cell_boundaries(id=id, fmt="list")
        self.assertEqual(boundaries, [5, 1])

        boundaries = grid_2d_xz.get_cell_boundaries(id=id, fmt="dict")
        self.assertEqual(boundaries, {"x": [5], "y": [], "z": [1]})

    def test_cell_coords(self):
        coords = (1, 0, 1)

        id = grid_2d_xz.get_cell_id(coords)
        self.assertEqual(id, 6)

        neighbors = grid_2d_xz.get_cell_neighbors(
            coords=coords, boundary=False, fmt="list"
        )
        self.assertEqual(neighbors, [(2, 0, 1), (1, 0, 2)])

        neighbors = grid_2d_xz.get_cell_neighbors(
            coords=coords, boundary=False, fmt="dict"
        )
        self.assertEqual(neighbors, {"x": [(2, 0, 1)], "y": [], "z": [(1, 0, 2)]})

        boundaries = grid_2d_xz.get_cell_boundaries(coords=coords, fmt="list")
        self.assertEqual(boundaries, [(0, 0, 1), (1, 0, 0)])

        boundaries = grid_2d_xz.get_cell_boundaries(coords=coords, fmt="dict")
        self.assertEqual(boundaries, {"x": [(0, 0, 1)], "y": [], "z": [(1, 0, 0)]})


class Test_2D_YZ(unittest.TestCase):
    def test_cell_id(self):
        id = 6

        coords = grid_2d_yz.get_cell_coords(id)
        self.assertEqual(coords, (0, 1, 1))

        neighbors = grid_2d_yz.get_cell_neighbors(id=id, boundary=False, fmt="list")
        self.assertEqual(neighbors, [7, 11])

        neighbors = grid_2d_yz.get_cell_neighbors(id=id, boundary=False, fmt="dict")
        self.assertEqual(neighbors, {"x": [], "y": [7], "z": [11]})

        boundaries = grid_2d_yz.get_cell_boundaries(id=id, fmt="list")
        self.assertEqual(boundaries, [5, 1])

        boundaries = grid_2d_yz.get_cell_boundaries(id=id, fmt="dict")
        self.assertEqual(boundaries, {"x": [], "y": [5], "z": [1]})

    def test_cell_coords(self):
        coords = (0, 1, 1)

        id = grid_2d_yz.get_cell_id(coords)
        self.assertEqual(id, 6)

        neighbors = grid_2d_yz.get_cell_neighbors(
            coords=coords, boundary=False, fmt="list"
        )
        self.assertEqual(neighbors, [(0, 2, 1), (0, 1, 2)])

        neighbors = grid_2d_yz.get_cell_neighbors(
            coords=coords, boundary=False, fmt="dict"
        )
        self.assertEqual(neighbors, {"x": [], "y": [(0, 2, 1)], "z": [(0, 1, 2)]})

        boundaries = grid_2d_yz.get_cell_boundaries(coords=coords, fmt="list")
        self.assertEqual(boundaries, [(0, 0, 1), (0, 1, 0)])

        boundaries = grid_2d_yz.get_cell_boundaries(coords=coords, fmt="dict")
        self.assertEqual(boundaries, {"x": [], "y": [(0, 0, 1)], "z": [(0, 1, 0)]})


class Test_3D(unittest.TestCase):
    def test_cell_id(self):
        id = 31
        coords = grid_3d.get_cell_coords(id)
        self.assertEqual(coords, (1, 1, 1))
        neighbors = grid_3d.get_cell_neighbors(id=id, boundary=False, fmt="list")
        self.assertEqual(neighbors, [32, 36, 56])
        neighbors = grid_3d.get_cell_neighbors(id=id, boundary=False, fmt="dict")
        self.assertEqual(neighbors, {"x": [32], "y": [36], "z": [56]})
        boundaries = grid_3d.get_cell_boundaries(id=id, fmt="list")
        self.assertEqual(boundaries, [30, 26, 6])
        boundaries = grid_3d.get_cell_boundaries(id=id, fmt="dict")
        self.assertEqual(boundaries, {"x": [30], "y": [26], "z": [6]})

    def test_cell_coords(self):
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
        self.assertEqual(boundaries, [(0, 1, 1), (1, 0, 1), (1, 1, 0)])
        boundaries = grid_3d.get_cell_boundaries(coords=coords, fmt="dict")
        self.assertEqual(
            boundaries, {"x": [(0, 1, 1)], "y": [(1, 0, 1)], "z": [(1, 1, 0)]}
        )


unify = False
grid_1d_x = grids.Cartesian(
    nx=3,
    ny=1,
    nz=1,
    dx=10,
    dy=10,
    dz=10,
    phi=0.27,
    kx=270,
    ky=270,
    kz=270,
    verbose=False,
    unify=unify,
)

grid_1d_y = grids.Cartesian(
    nx=1,
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
    unify=unify,
)

grid_1d_z = grids.Cartesian(
    nx=1,
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
    unify=unify,
)

grid_2d_xy = grids.Cartesian(
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
    unify=unify,
)
grid_2d_xz = grids.Cartesian(
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
    unify=unify,
)
grid_2d_yz = grids.Cartesian(
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
    unify=unify,
)
grid_3d = grids.Cartesian(
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
    unify=unify,
)


if __name__ == "__main__":
    unittest.main()
