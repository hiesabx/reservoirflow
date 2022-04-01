from openresim import grids
import numpy as np
import unittest


class TestApp(unittest.TestCase):

    def test_cell_id(self):
        id = 31
        coords = grid.get_cell_coords(id)
        self.assertEqual(coords, (1, 1, 1))
        neighbors = grid.get_cell_neighbors(id=id, boundary=False)
        self.assertEqual(neighbors, [32, 36, 56])
        boundaries = grid.get_cell_boundaries(id=id)
        self.assertEqual(boundaries, [26, 6, 30])
        
    def test_cell_coords(self):
        coords = (1,1,1)
        id = grid.get_cell_id(coords)
        self.assertEqual(id, 31)
        neighbors = grid.get_cell_neighbors(coords=coords, boundary=False)
        self.assertEqual(neighbors, [(2, 1, 1), (1, 2, 1), (1, 1, 2)])
        boundaries = grid.get_cell_boundaries(coords=coords)
        self.assertEqual(boundaries, [(0, 1, 1), (1, 1, 0), (1, 0, 1)])
        
if __name__ == '__main__':        
    grid = grids.CartGrid(nx=3, ny=3, nz=3, 
                          dx=10, dy=10, dz=10, 
                          phi=0.27, k=270)
    unittest.main()
