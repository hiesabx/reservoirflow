import numpy as np
import pyvista as pv
from reservoirflow import models
import math
import os


# Default settings color bar:
# cbar_opt = dict(height=0.40,
#                 vertical=True,
#                 position_x=0.85,
#                 position_y=0.60,
#                 fmt="%.f",
#                 interactive=True,
# )


# Default settings annotations:
# annotations = {
#     model.pressures[0]: 'Boundary',
#     model.wells[4]['pwf']: "Pwf",
# }


def show(model, property: str, centers=False, boundary=False, bounds=False):
    # Extract property values:
    props = list(vars(model))
    assert (
        property in props
    ), f"""Unknown property! Available properties are : {props}"""
    local_dict = locals()
    exec(f"values = model.{property}", local_dict)
    values = local_dict["values"]

    # Show boundary:
    if not boundary:
        cells_id = model.grid.get_cells_id(False, False, "list")
        values = values[:, cells_id]

    # Define limits:
    max_v = np.nanmax(values)
    min_v = math.floor(np.nanmin(values) / 1000) * 1000
    limits = [min_v, max_v]  # [min_v - (min_v*0.2), max_v + (max_v*0.2)]

    # Number Format:
    if min_v < 1:
        fmt = "%.2f"
    else:
        fmt = "%.f"

    # Color bar options:
    cbar_opt = dict(
        height=0.07,
        width=0.7,
        # horizontal=True,
        position_x=0.25,
        position_y=0.10,
        fmt=fmt,
        title=property,
        # interactive=True,
        # full_screen = True,
    )

    grid = model.grid.get_pyvista_grid(boundary)
    grid.cell_data[property] = values[1]

    # Setup plotter:
    pl = pv.Plotter()

    # Show grid:
    pl.add_mesh(
        grid,
        clim=limits,
        # style='wireframe',
        show_edges=True,
        opacity=0.7,
        lighting=True,
        ambient=0.2,
        n_colors=5,
        colormap="Blues",
        label=property,
        categories=True,
        # nan_color='gray',
        nan_opacity=0.7,
        # use_transparency=True,
        scalars=values[-1],  # or str 'pressures'
        scalar_bar_args=cbar_opt,
        show_scalar_bar=True,
        # annotations=annotations,
    )

    # Show centers:
    if centers:
        pl.add_points(
            grid.cell_centers(),
            point_size=10,
            render_points_as_spheres=True,
            show_edges=True,
            color="black",
        )

    # Show wells:
    def show_wells(pl):
        for w in model.wells:
            # x = model.grid.dx[1:w+1].sum() + model.grid.dx[w]//2
            # y = model.grid.dy[w]//2
            # z = 100
            height = model.grid.dz[w] * 10
            # well_cell_i = w if boundary else w - 1
            well_cell_center = list(
                model.grid.get_pyvista_grid(True).extract_cells(w).GetCenter()
            )
            well_cell_center[2] = height // 2
            well = pv.Cylinder(
                center=well_cell_center,
                height=height,
                radius=model.wells[w]["r"],
                direction=(0, 0, 1),
            )
            pl.add_mesh(well)

    show_wells(pl)
    # Plot configuration:
    pl.camera_position = "xz"
    # p.show_bounds(grid='front', location='outer', all_edges=True)
    # _ = p.update_scalar_bar_range(0, model.pressures.max())
    pl.add_camera_orientation_widget()
    pl.enable_fly_to_right_click()
    # pl.enable_zoom_style()
    pl.show_axes()
    pl.set_background("black", top="gray")
    if bounds:
        pl.show_bounds(
            grid="front", location="outer", all_edges=True
        )  # or pl.show_axes() # or pl.show_grid()

    pl.show(full_screen=True)

    pl = pv.Plotter(notebook=False, off_screen=True)
    show_wells(pl)

    pl.add_mesh(
        grid,
        clim=limits,
        # style='wireframe',
        show_edges=True,
        opacity=0.7,
        lighting=True,
        ambient=0.2,
        n_colors=5,
        colormap="Blues",
        label=property,
        categories=True,
        # nan_color='gray',
        nan_opacity=0.7,
        # use_transparency=True,
        scalars=values[-1],  # or str 'pressures'
        scalar_bar_args=cbar_opt,
        show_scalar_bar=True,
        # annotations=annotations,
    )

    if not os.path.exists("images/"):
        os.makedirs("images/")

    pl.open_gif("images/grid.gif")

    pts = grid.points.copy()
    for step in range(model.nsteps):
        pl.update_coordinates(pts, render=False)
        pl.update_scalars(values[step], render=False)
        pl.render()
        pl.write_frame()
        # time.sleep(2)

    pl.close()

    # To show time with 3D
    # https://docs.pyvista.org/api/plotting/plotting.html
