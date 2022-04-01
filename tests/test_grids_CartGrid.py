from openresim import grids
import numpy as np
import unittest


class TestApp(unittest.TestCase):

    def test_CartGrid(self):
        test_grid(nx=1, ny=1, nz=1)
        test_grid(nx=2, ny=1, nz=1)
        test_grid(nx=1, ny=2, nz=1)
        test_grid(nx=1, ny=1, nz=2)
        test_grid(nx=2, ny=2, nz=1)
        test_grid(nx=2, ny=1, nz=2)
        test_grid(nx=1, ny=2, nz=2)
        test_grid(nx=2, ny=2, nz=2)
        test_grid(nx=3, ny=1, nz=1)
        test_grid(nx=1, ny=3, nz=1)
        test_grid(nx=1, ny=1, nz=3)
        test_grid(nx=3, ny=3, nz=1)
        test_grid(nx=3, ny=1, nz=3)
        test_grid(nx=1, ny=3, nz=3)
        test_grid(nx=3, ny=3, nz=3)
        
if __name__ == '__main__':
    def test_grid(nx,ny,nz):
        grid = grids.CartGrid(nx=nx, ny=ny, nz=nz, dx=10, dy=10, dz=10, phi=0.27, k=270)
        grid.show(boundary=True, centers='id')
        grid.show(boundary=False, centers='id')
    unittest.main()
