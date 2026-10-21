from functools import cache

import pytest

from reservoirflow import grids

AXES = ("x", "y", "z")


def full(d):
    """Dict with all three axes, missing ones empty."""
    return {a: d.get(a, []) for a in AXES}


def flat(d):
    """List form of a per-axis dict, in x, y, z order."""
    return [v for a in AXES for v in d.get(a, [])]


# id, coords, neighbors and boundaries (by id and by coords), per axis
CASES = {
    "1d_x": {
        "n": (3, 1, 1),
        "id": 1,
        "coords": (1, 0, 0),
        "nbr": {"x": [2]},
        "bnd": {"x": [0]},
        "nbr_c": {"x": [(2, 0, 0)]},
        "bnd_c": {"x": [(0, 0, 0)]},
    },
    "1d_y": {
        "n": (1, 3, 1),
        "id": 1,
        "coords": (0, 1, 0),
        "nbr": {"y": [2]},
        "bnd": {"y": [0]},
        "nbr_c": {"y": [(0, 2, 0)]},
        "bnd_c": {"y": [(0, 0, 0)]},
    },
    "1d_z": {
        "n": (1, 1, 3),
        "id": 1,
        "coords": (0, 0, 1),
        "nbr": {"z": [2]},
        "bnd": {"z": [0]},
        "nbr_c": {"z": [(0, 0, 2)]},
        "bnd_c": {"z": [(0, 0, 0)]},
    },
    "2d_xy": {
        "n": (3, 3, 1),
        "id": 6,
        "coords": (1, 1, 0),
        "nbr": {"x": [7], "y": [11]},
        "bnd": {"x": [5], "y": [1]},
        "nbr_c": {"x": [(2, 1, 0)], "y": [(1, 2, 0)]},
        "bnd_c": {"x": [(0, 1, 0)], "y": [(1, 0, 0)]},
    },
    "2d_xz": {
        "n": (3, 1, 3),
        "id": 6,
        "coords": (1, 0, 1),
        "nbr": {"x": [7], "z": [11]},
        "bnd": {"x": [5], "z": [1]},
        "nbr_c": {"x": [(2, 0, 1)], "z": [(1, 0, 2)]},
        "bnd_c": {"x": [(0, 0, 1)], "z": [(1, 0, 0)]},
    },
    "2d_yz": {
        "n": (1, 3, 3),
        "id": 6,
        "coords": (0, 1, 1),
        "nbr": {"y": [7], "z": [11]},
        "bnd": {"y": [5], "z": [1]},
        "nbr_c": {"y": [(0, 2, 1)], "z": [(0, 1, 2)]},
        "bnd_c": {"y": [(0, 0, 1)], "z": [(0, 1, 0)]},
    },
    "3d": {
        "n": (3, 3, 3),
        "id": 31,
        "coords": (1, 1, 1),
        "nbr": {"x": [32], "y": [36], "z": [56]},
        "bnd": {"x": [30], "y": [26], "z": [6]},
        "nbr_c": {"x": [(2, 1, 1)], "y": [(1, 2, 1)], "z": [(1, 1, 2)]},
        "bnd_c": {"x": [(0, 1, 1)], "y": [(1, 0, 1)], "z": [(1, 1, 0)]},
    },
}

ALL_CASES = pytest.mark.parametrize("case", CASES.values(), ids=CASES.keys())


@cache
def make_grid(nx, ny, nz):
    return grids.RegularCartesianGrid(
        nx=nx,
        ny=ny,
        nz=nz,
        dx=10,
        dy=10,
        dz=10,
        phi=0.27,
        kx=270,
        ky=270,
        kz=270,
        verbose=False,
        unify=False,
    )


@ALL_CASES
def test_cell_id(case):
    grid = make_grid(*case["n"])
    id = case["id"]

    assert grid.get_cell_coords(id) == case["coords"]

    neighbors = grid.get_cell_neighbors(id=id, boundary=False, fmt="list")
    assert neighbors == flat(case["nbr"])
    neighbors = grid.get_cell_neighbors(id=id, boundary=False, fmt="dict")
    assert neighbors == full(case["nbr"])

    boundaries = grid.get_cell_boundaries(id=id, fmt="list")
    assert boundaries == flat(case["bnd"])
    boundaries = grid.get_cell_boundaries(id=id, fmt="dict")
    assert boundaries == full(case["bnd"])


@ALL_CASES
def test_cell_coords(case):
    grid = make_grid(*case["n"])
    coords = case["coords"]

    assert grid.get_cell_id(coords) == case["id"]

    neighbors = grid.get_cell_neighbors(coords=coords, boundary=False, fmt="list")
    assert neighbors == flat(case["nbr_c"])
    neighbors = grid.get_cell_neighbors(coords=coords, boundary=False, fmt="dict")
    assert neighbors == full(case["nbr_c"])

    boundaries = grid.get_cell_boundaries(coords=coords, fmt="list")
    assert boundaries == flat(case["bnd_c"])
    boundaries = grid.get_cell_boundaries(coords=coords, fmt="dict")
    assert boundaries == full(case["bnd_c"])
