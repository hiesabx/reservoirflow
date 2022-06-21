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


if __name__ == "__main__":

    def get_d(d_0, n):
        if n > 1:
            return [d_0] + [d_0 + (i * d_0) for i in range(1, n + 1)] + [d_0]
        else:
            return d_0

    def test_grid(nx, ny, nz):
        grid = grids.CartGrid(
            nx=nx,
            ny=ny,
            nz=nz,
            dx=get_d(10, nx),
            dy=get_d(10, ny),
            dz=get_d(10, nz),
            kx=270,
            ky=270,
            kz=270,
            phi=0.27,
        )
        grid.show(boundary=True, label="id")
        grid.show(boundary=False, label="id")

    unittest.main()
