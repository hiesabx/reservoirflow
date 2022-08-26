from openresim import grids
import numpy as np
import unittest


class TestApp(unittest.TestCase):
    def test_cell_d(self):
        # cell: coords=(1,1,1), id=21
        self.assertEqual(grid.get_cell_d(dir="x", coords=(1, 1, 1)), 21)
        self.assertEqual(grid.get_cell_d(dir="y", coords=(1, 1, 1)), 22)
        self.assertEqual(grid.get_cell_d(dir="z", coords=(1, 1, 1)), 23)
        self.assertEqual(grid.get_cell_d(dir="x", id=21), 21)
        self.assertEqual(grid.get_cell_d(dir="y", id=21), 22)
        self.assertEqual(grid.get_cell_d(dir="z", id=21), 23)
        self.assertEqual(grid.get_cell_dx(coords=(1, 1, 1)), 21)
        self.assertEqual(grid.get_cell_dy(coords=(1, 1, 1)), 22)
        self.assertEqual(grid.get_cell_dz(coords=(1, 1, 1)), 23)
        self.assertEqual(grid.get_cell_dx(id=21), 21)
        self.assertEqual(grid.get_cell_dy(id=21), 22)
        self.assertEqual(grid.get_cell_dz(id=21), 23)
        # cell: coords=(2,2,2), id=42
        self.assertEqual(grid.get_cell_d(dir="x", coords=(2, 2, 2)), 31)
        self.assertEqual(grid.get_cell_d(dir="y", coords=(2, 2, 2)), 32)
        self.assertEqual(grid.get_cell_d(dir="z", coords=(2, 2, 2)), 33)
        self.assertEqual(grid.get_cell_d(dir="x", id=42), 31)
        self.assertEqual(grid.get_cell_d(dir="y", id=42), 32)
        self.assertEqual(grid.get_cell_d(dir="z", id=42), 33)
        self.assertEqual(grid.get_cell_dx(coords=(2, 2, 2)), 31)
        self.assertEqual(grid.get_cell_dy(coords=(2, 2, 2)), 32)
        self.assertEqual(grid.get_cell_dz(coords=(2, 2, 2)), 33)
        self.assertEqual(grid.get_cell_dx(id=42), 31)
        self.assertEqual(grid.get_cell_dy(id=42), 32)
        self.assertEqual(grid.get_cell_dz(id=42), 33)

    def test_cell_area(self):
        # cell: coords=(1,1,1), id=21
        self.assertEqual(grid.get_cell_A(dir="x", coords=(1, 1, 1)), 506)
        self.assertEqual(grid.get_cell_A(dir="y", coords=(1, 1, 1)), 483)
        self.assertEqual(grid.get_cell_A(dir="z", coords=(1, 1, 1)), 462)
        self.assertEqual(grid.get_cell_A(dir="x", id=21), 506)
        self.assertEqual(grid.get_cell_A(dir="y", id=21), 483)
        self.assertEqual(grid.get_cell_A(dir="z", id=21), 462)
        self.assertEqual(grid.get_cell_Ax(coords=(1, 1, 1)), 506)
        self.assertEqual(grid.get_cell_Ay(coords=(1, 1, 1)), 483)
        self.assertEqual(grid.get_cell_Az(coords=(1, 1, 1)), 462)
        self.assertEqual(grid.get_cell_Ax(id=21), 506)
        self.assertEqual(grid.get_cell_Ay(id=21), 483)
        self.assertEqual(grid.get_cell_Az(id=21), 462)
        # cell: coords=(2,2,2), id=42
        self.assertEqual(grid.get_cell_A(dir="x", coords=(2, 2, 2)), 1056)
        self.assertEqual(grid.get_cell_A(dir="y", coords=(2, 2, 2)), 1023)
        self.assertEqual(grid.get_cell_A(dir="z", coords=(2, 2, 2)), 992)
        self.assertEqual(grid.get_cell_A(dir="x", id=42), 1056)
        self.assertEqual(grid.get_cell_A(dir="y", id=42), 1023)
        self.assertEqual(grid.get_cell_A(dir="z", id=42), 992)
        self.assertEqual(grid.get_cell_Ax(coords=(2, 2, 2)), 1056)
        self.assertEqual(grid.get_cell_Ay(coords=(2, 2, 2)), 1023)
        self.assertEqual(grid.get_cell_Az(coords=(2, 2, 2)), 992)
        self.assertEqual(grid.get_cell_Ax(id=42), 1056)
        self.assertEqual(grid.get_cell_Ay(id=42), 1023)
        self.assertEqual(grid.get_cell_Az(id=42), 992)

    def test_cell_coords(self):
        # cell: coords=(1,1,1), id=21
        self.assertEqual(grid.get_cell_V(coords=(1, 1, 1)), 10626)
        self.assertEqual(grid.get_cell_V(id=21), 10626)
        # cell: coords=(2,2,2), id=42
        self.assertEqual(grid.get_cell_V(coords=(2, 2, 2)), 32736)
        self.assertEqual(grid.get_cell_V(id=42), 32736)


if __name__ == "__main__":
    dx = [11, 21, 31, 41]
    dy = [12, 22, 32, 42]
    dz = [13, 23, 33, 43]
    grid = grids.Cartesian(
        nx=2,
        ny=2,
        nz=2,
        dx=dx,
        dy=dy,
        dz=dz,
        kx=270,
        ky=270,
        kz=270,
        phi=0.27,
        unify=True,
    )

    unittest.main()
