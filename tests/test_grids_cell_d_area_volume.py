import pytest

from reservoirflow import grids

# (coords, id, (dx, dy, dz), (Ax, Ay, Az), V)
CELLS = [
    pytest.param((1, 1, 1), 21, (21, 22, 23), (506, 483, 462), 10626, id="cell21"),
    pytest.param((2, 2, 2), 42, (31, 32, 33), (1056, 1023, 992), 32736, id="cell42"),
]


@pytest.fixture(scope="module")
def grid():
    return grids.RegularCartesianGrid(
        nx=2,
        ny=2,
        nz=2,
        dx=[11, 21, 31, 41],
        dy=[12, 22, 32, 42],
        dz=[13, 23, 33, 43],
        kx=270,
        ky=270,
        kz=270,
        phi=0.27,
        unify=True,
    )


@pytest.mark.parametrize("coords, cell_id, d, A, V", CELLS)
def test_cell_d(grid, coords, cell_id, d, A, V):
    for direction, expected in zip("xyz", d):
        get_d = getattr(grid, f"get_cell_d{direction}")
        assert grid.get_cell_d(dir=direction, coords=coords) == expected
        assert grid.get_cell_d(dir=direction, id=cell_id) == expected
        assert get_d(coords=coords) == expected
        assert get_d(id=cell_id) == expected


@pytest.mark.parametrize("coords, cell_id, d, A, V", CELLS)
def test_cell_area(grid, coords, cell_id, d, A, V):
    for direction, expected in zip("xyz", A):
        get_A = getattr(grid, f"get_cell_A{direction}")
        assert grid.get_cell_A(dir=direction, coords=coords) == expected
        assert grid.get_cell_A(dir=direction, id=cell_id) == expected
        assert get_A(coords=coords) == expected
        assert get_A(id=cell_id) == expected


@pytest.mark.parametrize("coords, cell_id, d, A, V", CELLS)
def test_cell_volume(grid, coords, cell_id, d, A, V):
    assert grid.get_cell_V(coords=coords) == V
    assert grid.get_cell_V(id=cell_id) == V
