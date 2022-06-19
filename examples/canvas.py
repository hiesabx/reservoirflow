#%% # Import
from openresim import grids

#%% # Create Grid:
def model_1():

    g = grids.CartGrid(nx=3, ny=3, nz=1, dx=10, dy=10, dz=10, verbose=False, unify=True)

    return g


#%% Show:
def show(label="id"):
    g.show(label=label, boundary=True)
    # g.show(label=label, boundary=False)


#%% Check Neighbors:
def neighbors():
    fmt = "set"
    id = g.get_cells_id(False, False)
    coords = g.get_cells_coords(False, False)
    print(id, coords)
    for i, c in zip(id, coords):
        print(f"id: {i} - coords: {c}")
        neighbors_id = g.get_cell_neighbors(id=i, boundary=True, fmt="array")
        print(" - neighbors id:", neighbors_id)
        neighbors_coords = g.get_cell_neighbors(coords=c, boundary=True, fmt="array")
        print(" - neighbors coords:", neighbors_coords)


def boundaries():
    fmt = "tuple"
    id = g.get_cells_id(False, False, fmt)
    coords = g.get_cells_coords(False, False, fmt)
    # print(g.remove_boundaries(coords))
    # print(g.extract_boundaries(id, fmt="tuple"))
    # print(g.extract_boundaries(coords, fmt="set"))
    # print(g.get_boundaries("coords", "set"))
    for i, c in zip(id, coords):
        # print(f"- id: {i}")
        print(f"- coords: {c}")
        # print(g.get_cell_boundaries(id=i, fmt="dict"))
        print(g.get_cell_boundaries(coords=c, fmt="array"))
        # break€€


def props():
    g.set_props(phi=0.30)
    g.set_props(phi=0.20, id=(6))
    # print(g.phi)
    g.set_props(kx=1)
    g.set_props(kx=100000, id=6)
    print(g.kx)


def prop():
    g.verbose = True
    g.set_prop(name="phi", value=0.30)
    g.set_prop(name="phi", value=0, id=11)
    # print(g.props["phi"])
    print(g.get_prop("kx"))


# %%
if __name__ == "__main__":
    g = model_1()
    prop()
    # props()
    # neighbors()
    # boundaries()
    # show("id")
    # show("coords")
