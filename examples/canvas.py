#%% # Import
from openresim import grids

#%% # Create Grid:
def get_d(d_0, n):
    if n > 1:
        return [d_0] + [d_0 + (i * d_0) for i in range(1, n + 1)] + [d_0]
    else:
        return d_0


def model_1(n, unify):

    nx, ny, nz = n
    g = grids.Cartesian(
        nx=nx,
        ny=ny,
        nz=nz,
        dx=get_d(10, nx),
        dy=get_d(10, ny),
        dz=get_d(10, nz),
        verbose=False,
        unify=unify,
    )

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
    # g.set_props(phi=0.20, id=(6))
    g.set_prop("phi", value=0.20, id=6)
    print(id(g.phi))
    print(id(g.props["phi"]))
    g.set_props(kx=1)
    g.set_prop("kx", value=100000, id=6)
    print(id(g.kx))
    print(id(g.props["kx"]))


def prop():
    g.verbose = True
    g.set_prop(name="phi", value=0.30)
    g.set_prop(name="phi", value=0, id=11)
    g.set_prop(name="kx", value=300)
    g.set_prop(name="kx", value=30, id=11)
    # print(g.props["phi"])
    print(g.is_homogeneous)
    print(g.is_heterogeneous)


def G():
    g.get_G(dir="x")


def cells_D():
    print(g.get_cells_d("x", False, True))
    print(g.get_cells_d("y", False, False))
    print(g.get_cells_d("y", True, False))
    print(g.get_cells_d("z", True, True))


def cell_d():
    print(g.get_cells_id(True, True))
    print(g.get_cell_d(dir="x", id=13))
    print(g.dx)
    print(g.get_cell_d(dir="x", coords=(3, 2, 0)))


def cells_A():
    print(g.get_cells_id(True, True))
    print(g.get_cells_Ax(True, True))
    print(g.get_cell_Ax(id=13))
    coords = g.get_cell_coords(13, True)
    print(g.get_cell_Ax(coords=coords))


def volume():
    # print(g.dx)
    # print(g.dy)
    # print(g.dz)
    # print(g.Ax)
    # print(g.Ay)
    # print(g.Az)
    # print(g.get_Vt(boundary=True, pyvista=True))
    # print(g.get_Vt(boundary=False, pyvista=True))
    # print(g.get_Vt(boundary=True, pyvista=False))
    # print(g.get_Vt(boundary=False, pyvista=False))
    print(g.get_cells_V(boundary=True, fshape=True, pyvista=True))
    # print(g.get_cells_V(boundary=True, fshape=False, pyvista=True))
    # print(g.get_cells_V(boundary=False, fshape=True, pyvista=True))
    # print(g.get_cells_V(boundary=False, fshape=False, pyvista=True))

    # print(g.get_cells_V(boundary=True, fshape=True, pyvista=False))
    # print(g.get_cells_V(boundary=True, fshape=False, pyvista=False))
    # print(g.get_cells_V(boundary=False, fshape=True, pyvista=False))
    # print(g.get_cells_V(boundary=False, fshape=False, pyvista=False))
    print(g.get_cell_V(id=6))
    print(g.get_cell_V(id=13))
    print(g.get_cell_V(coords=(1, 1, 0)))
    print(g.get_cell_V(coords=(3, 2, 0)))


def center():
    print(g.get_cells_center(True, False, True))
    print(g.get_cells_center(True, False, False))
    print(g.get_cell_center(id=6), g.get_cell_center(coords=(2, 1, 0)))
    print(g.get_cell_center(id=10), g.get_cell_center(coords=(2, 2, 0)))


def G():
    g.set_prop("kx", 100)
    g.set_prop("ky", 100)
    g.set_prop("kz", 100)
    print(g.get_Gx())
    print(g.get_Gy())
    print(g.get_Gz())


# %%
if __name__ == "__main__":
    g = model_1((2, 2, 1), unify=False)
    # G()
    # center()
    # volume()
    # cell_d()
    # cells_A()
    prop()
    props()
    # neighbors()
    # boundaries()
    show("id")
    show("coords")
    # show("dx")
    # show("center")
