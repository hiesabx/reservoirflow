from re import S
from numpy import size
from openresim import grids
import pyvista as pv
import numpy as np

z = [3212.73, 3182.34, 3121.56, 3060.78, 3000, 2969.62]
grid = grids.Cartesian(nx=4, ny=1, nz=1, dx=300, dy=350, dz=40, z=z, phi=0.27, kx=270)

pv_grid = grid.get_pyvista_grid()
pv_grid.cell_data["dx"] = grid.Ax

# https://docs.pyvista.org/api/core/_autosummary/pyvista.DataSet.html
print("cell_data:", pv_grid.cell_data)
print("dx:", pv_grid.cell_data["dx"])  # or pv_grid.cell_arrays['dx']
print("n_arrays:", pv_grid.n_arrays)

print("n_cells:", pv_grid.n_cells)
print("cells_id", np.arange(pv_grid.n_cells))
print("n_points:", pv_grid.n_points)
print("cell indecies:", pv_grid)
print("number of points for cell 0:", pv_grid.cell_n_points(0))
print("center:", pv_grid.center)
print("cell_centers points:", pv_grid.cell_centers().points)
print("points:", pv_grid.points)
print("cell points of cell 0:", pv_grid.cell_points(0))
print("neighbors cell 0:", pv_grid.neighbors(0))
print("type of cell 0:", pv_grid.cell_type(0))
print("bounds:", pv_grid.bounds)
print("bounds for cell 0:", pv_grid.cell_bounds(0))
print("area", pv_grid.compute_cell_quality("area").cell_data["CellQuality"])

show = False
if show:
    # 1. Create Plotter:
    p = pv.Plotter()
    # 2. Plot Grid:
    p.add_mesh(pv_grid, show_edges=True, color="white", ambient=0.2, opacity=0.7)
    # 3. Plot centers:
    p.add_mesh(
        pv_grid.cell_centers(),
        render_points_as_spheres=True,
        show_edges=True,
        color="black",
    )
    # 4. Plot main center:
    center = pv_grid.center
    sphere = pv.Sphere(center=center, radius=100)
    p.add_mesh(sphere, render_points_as_spheres=True, color="green")
    # 5. Plot specific grid center:
    center_1 = pv_grid.extract_cells(1).GetCenter()
    well = pv.Cylinder(center=center_1, height=100, radius=10, direction=(0, 0, 1))
    p.add_mesh(well)
    # 6. Show:
    p.show(full_screen=True)
