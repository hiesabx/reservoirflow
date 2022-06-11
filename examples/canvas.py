#%% # Import
from openresim import grids

#%% # Create Grid:
g = grids.CartGrid(
    nx=3,
    ny=3,
    nz=3,
    dx=10,
    dy=10,
    dz=10,
)

#%% Show:
def show(label="id"):
    g.show(label=label, boundary=True)
    g.show(label=label, boundary=False)


#%% Check Neighbors:
def check_neighbors():
    fmt1 = "set"
    fmt2 = "set"
    id = g.get_cells_id(False, False, fmt1)
    coords = g.get_cells_coords(False, False, fmt1)

    for i, c in zip(id, coords):
        print(f"id: {i} - coords: {c}")
        neighbors_id = g.get_cell_neighbors(id=i, boundary=False, fmt=fmt2)
        neighbors_coords = g.get_cell_neighbors(coords=c, boundary=False, fmt=fmt2)
        print(" - neighbors id:", neighbors_id)
        print(" - neighbors coords:", neighbors_coords)


def boundaries():
    fmt1 = "array"
    # fmt2 = "set"
    id = g.get_cells_id(True, True, fmt1)
    coords = g.get_cells_coords(True, True, fmt1)
    print(g.remove_boundaries(id))
    print(g.remove_boundaries(coords))
    print(g.keep_boundaries(id, fmt="list"))
    print(g.keep_boundaries(coords, fmt="list"))


# %%
if __name__ == "__main__":
    # check_neighbors()
    boundaries()
    # show("id")
