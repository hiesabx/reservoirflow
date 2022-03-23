from openresim import grids
import numpy as np
import unittest


class TestApp(unittest.TestCase):

    def test_cell_id(self):
        self.assertEqual(grid.cell_id(0,0,0, boundary=False),  0)
        self.assertEqual(grid.cell_id(2,0,0, boundary=False),  2)
        self.assertEqual(grid.cell_id(0,2,0, boundary=False),  6)
        self.assertEqual(grid.cell_id(0,0,2, boundary=False), 18)
        self.assertEqual(grid.cell_id(2,2,0, boundary=False),  8)
        self.assertEqual(grid.cell_id(2,2,2, boundary=False), 26)
        self.assertEqual(grid.cell_id(1,1,1, boundary=True),  31)
        self.assertEqual(grid.cell_id(3,1,1, boundary=True),  33)
        self.assertEqual(grid.cell_id(1,3,1, boundary=True),  41)
        self.assertEqual(grid.cell_id(1,1,3, boundary=True),  81)
        self.assertEqual(grid.cell_id(3,3,1, boundary=True),  43)
        self.assertEqual(grid.cell_id(3,3,3, boundary=True),  93)
        
    def test_cell_coords(self):
        self.assertEqual(grid.cell_coords( 0, boundary=False), (0,0,0))
        self.assertEqual(grid.cell_coords( 2, boundary=False), (2,0,0))
        self.assertEqual(grid.cell_coords( 6, boundary=False), (0,2,0))
        self.assertEqual(grid.cell_coords(18, boundary=False), (0,0,2))
        self.assertEqual(grid.cell_coords( 8, boundary=False), (2,2,0))
        self.assertEqual(grid.cell_coords(26, boundary=False), (2,2,2))
        self.assertEqual(grid.cell_coords(31, boundary=True),  (1,1,1))
        self.assertEqual(grid.cell_coords(33, boundary=True),  (3,1,1))
        self.assertEqual(grid.cell_coords(41, boundary=True),  (1,3,1))
        self.assertEqual(grid.cell_coords(81, boundary=True),  (1,1,3))
        self.assertEqual(grid.cell_coords(43, boundary=True),  (3,3,1))
        self.assertEqual(grid.cell_coords(93, boundary=True),  (3,3,3))
        
if __name__ == '__main__':        
    grid = grids.CartGrid(nx=3, ny=3, nz=3, dx=10, dy=10, dz=10, phi=0.27, k=270)
    unittest.main()
