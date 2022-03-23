from openresim import grids
import numpy as np
import unittest


class TestApp(unittest.TestCase):

    def test_CartGrid(self):
        grids.CartGrid(nx=1, ny=1, nz=1, dx=10, dy=10, dz=10, phi=0.27, k=270)
        grids.CartGrid(nx=2, ny=1, nz=1, dx=10, dy=10, dz=10, phi=0.27, k=270)
        grids.CartGrid(nx=1, ny=2, nz=1, dx=10, dy=10, dz=10, phi=0.27, k=270)
        grids.CartGrid(nx=1, ny=1, nz=2, dx=10, dy=10, dz=10, phi=0.27, k=270)
        grids.CartGrid(nx=2, ny=2, nz=1, dx=10, dy=10, dz=10, phi=0.27, k=270)
        grids.CartGrid(nx=2, ny=1, nz=2, dx=10, dy=10, dz=10, phi=0.27, k=270)
        grids.CartGrid(nx=1, ny=2, nz=2, dx=10, dy=10, dz=10, phi=0.27, k=270)
        grids.CartGrid(nx=2, ny=2, nz=2, dx=10, dy=10, dz=10, phi=0.27, k=270)
        grids.CartGrid(nx=3, ny=1, nz=1, dx=10, dy=10, dz=10, phi=0.27, k=270)
        grids.CartGrid(nx=1, ny=3, nz=1, dx=10, dy=10, dz=10, phi=0.27, k=270)
        grids.CartGrid(nx=1, ny=1, nz=3, dx=10, dy=10, dz=10, phi=0.27, k=270)
        grids.CartGrid(nx=3, ny=3, nz=1, dx=10, dy=10, dz=10, phi=0.27, k=270)
        grids.CartGrid(nx=3, ny=1, nz=3, dx=10, dy=10, dz=10, phi=0.27, k=270)
        grids.CartGrid(nx=1, ny=3, nz=3, dx=10, dy=10, dz=10, phi=0.27, k=270)
        grids.CartGrid(nx=3, ny=3, nz=3, dx=10, dy=10, dz=10, phi=0.27, k=270)
        
if __name__ == '__main__':        
    unittest.main()
