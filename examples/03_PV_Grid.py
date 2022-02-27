import numpy as np
import pyvista as pv
from pyvista import examples


dem = examples.download_crater_topo()
subset = dem.extract_subset((500, 900, 400, 800, 0, 0), (5,5,1))
# subset.plot(cpos="xy")
terrain = subset.warp_by_scalar()
print(terrain.points)
z_cells = np.array([10, 100, 10])
print(np.cumsum(z_cells).reshape((1, 1, -1)))

xx = np.repeat(terrain.x, len(z_cells), axis=-1)
yy = np.repeat(terrain.y, len(z_cells), axis=-1)
zz = np.repeat(terrain.z, len(z_cells), axis=-1) - np.cumsum(z_cells).reshape((1, 1, -1))

mesh = pv.StructuredGrid(xx, yy, zz)
mesh["Elevation"] = zz.ravel(order="F")
cpos = [(1826736.796308761, 5655837.275274233, 4676.8405505181745),
 (1821066.1790519988, 5649248.765538796, 943.0995128226014),
 (-0.2797856225380979, -0.27966946337594883, 0.9184252809434081)]

mesh.plot(show_edges=True, 
    lighting=False, 
    cpos=cpos,
    full_screen=True
)