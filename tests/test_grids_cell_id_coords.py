import pytest

from reservoirflow import grids

# (coords, boundary, id)
CELLS = [
    ((0, 0, 0), False, 0),
    ((2, 0, 0), False, 2),
    ((0, 2, 0), False, 6),
    ((0, 0, 2), False, 18),
    ((2, 2, 0), False, 8),
    ((2, 2, 2), False, 26),
    ((1, 1, 1), True, 31),
    ((3, 1, 1), True, 33),
    ((1, 3, 1), True, 41),
    ((1, 1, 3), True, 81),
    ((3, 3, 1), True, 43),
    ((3, 3, 3), True, 93),
]


@pytest.fixture(scope="module")
def grid():
    return grids.RegularCartesianGrid(
        nx=3, ny=3, nz=3, dx=10, dy=10, dz=10, phi=0.27, kx=270, ky=270, kz=270
    )


@pytest.mark.parametrize("coords, boundary, cell_id", CELLS)
def test_cell_id(grid, coords, boundary, cell_id):
    assert grid.get_cell_id(coords, boundary=boundary) == cell_id


@pytest.mark.parametrize("coords, boundary, cell_id", CELLS)
def test_cell_coords(grid, coords, boundary, cell_id):
    assert grid.get_cell_coords(cell_id, boundary=boundary) == coords
