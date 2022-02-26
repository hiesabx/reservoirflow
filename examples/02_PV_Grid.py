from numpy import size
from openresim import grids
import pyvista as pv


grid = grids.Grid1D(nx=3, ny=1, nz=1, dx=300, dy=350, dz=40, phi=0.27, k=270)

pv_grid = grid.get_pv_grid()
pv_grid.cell_data['dx'] = grid.dx

print('n_cells:', pv_grid.n_cells)
print('n_points:', pv_grid.n_points)
print('n_arrays:', pv_grid.n_arrays)
print('dx:', pv_grid.cell_data['dx']) # or pv_grid.cell_arrays['dx']
print('points:', pv_grid.points)
print('cell_centers points:', pv_grid.cell_centers().points)
print('neighbors:', pv_grid.neighbors(2))

# 1. Create Plotter:
p = pv.Plotter()
# 2. Plot Grid:
p.add_mesh(pv_grid, show_edges=True, color='white', ambient=0.2, opacity=0.7)
# 3. Plot centers:
p.add_mesh(pv_grid.cell_centers(), render_points_as_spheres=True, show_edges=True, color='black')
# 4. Plot main center:
center = pv_grid.center
sphere = pv.Sphere(center=center, radius=100)
p.add_mesh(sphere, render_points_as_spheres=True, color='green')
# 5. Plot specific grid center:
center_1 = pv_grid.extract_cells(1).GetCenter()
well = pv.Cylinder(center=center_1, height=100, radius=10, direction=(0,0,1))
p.add_mesh(well)
# 6. Show:
p.show(full_screen=True)