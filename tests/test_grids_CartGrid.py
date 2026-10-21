import pytest
import pyvista as pv

from reservoirflow import grids

pv.OFF_SCREEN = True


def get_d(d_0, n):
    if n > 1:
        return [d_0] + [d_0 + (i * d_0) for i in range(1, n + 1)] + [d_0]
    else:
        return d_0


@pytest.mark.parametrize(
    "nx, ny, nz",
    [
        (1, 1, 1),
        (2, 1, 1),
        (1, 2, 1),
        (1, 1, 2),
        (2, 2, 1),
        (2, 1, 2),
        (1, 2, 2),
        (2, 2, 2),
        (3, 1, 1),
        (1, 3, 1),
        (1, 1, 3),
        (3, 3, 1),
        (3, 1, 3),
        (1, 3, 3),
        (3, 3, 3),
    ],
    ids=lambda v: str(v),
)
def test_grid(nx, ny, nz):
    grid = grids.RegularCartesianGrid(
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
