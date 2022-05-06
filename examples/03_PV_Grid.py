import pyvista
sphere = pyvista.Sphere(theta_resolution=20, phi_resolution=20)
cqual = sphere.compute_cell_quality('area')
# cqual.plot(show_edges=True)
print(cqual['CellQuality'])