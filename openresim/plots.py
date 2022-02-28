import numpy as np
import pyvista as pv
from openresim import models



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

def show_grid(model: models, property:str, show_centers=True, show_boundary=False, show_bounds=False):
    # Extract property values:
    try:
        local_dict = locals()
        exec(f"values = model.{property}", local_dict)
        values = local_dict['values'].reshape(model.grid.nx+2, model.grid.ny, model.grid.nz)
    except:
        props = list(vars(model))
        raise ValueError(f"""Unknown property. 
                            Available properties are : {props}""")
    if not show_boundary:
        values = values[1:-1]
    # Define limits: 
    max_v = max(values)
    min_v = min(values)
    limits = [min_v, max_v] # [min_v - (min_v*0.2), max_v + (max_v*0.2)]
    # Number Format: 
    if min_v < 1:
        fmt="%.2f"
    else:
        fmt="%.f"
    # Color bar options:
    cbar_opt = dict(height=0.07,
                    width=0.7,
                    #horizontal=True,
                    position_x=0.25, 
                    position_y=0.10,
                    fmt=fmt,
                    title=property,
                    #interactive=True,
                    # full_screen = True,
    )

    grid = model.grid.get_pv_grid(show_boundary)
    grid.cell_data[property] = values
    # Setup plotter:
    pl = pv.Plotter()
    pl.add_camera_orientation_widget()
    pl.enable_fly_to_right_click()
    # pl.enable_zoom_style()
    # pl.add_mesh(grid0)
    pl.add_mesh(grid,
        clim=limits,
        # style='wireframe',
        show_edges=True,
        opacity=0.7,
        lighting=True,
        ambient=0.2,
        n_colors=5,
        colormap='Blues',
        label=property,
        categories=True,
        #nan_color='gray',
        nan_opacity=0.7,
        # use_transparency=True,
        # scalars='Pressures', 
        scalar_bar_args=cbar_opt, 
        show_scalar_bar=True,
        # annotations=annotations,
    )
    if show_centers:
        pl.add_points(grid.cell_centers(),
            point_size=10, 
            render_points_as_spheres=True, 
            show_edges=True, 
            color='black'
        )
    for w in model.wells:
        # x = model.grid.dx[1:w+1].sum() + model.grid.dx[w]//2
        # y = model.grid.dy[w]//2
        # z = 100
        # well_center = (x, y, z)
        height = model.grid.dz[w] * 2
        well_cell_i = w if show_boundary else w - 1
        well_cell_center = list(model.grid.pv_grid.extract_cells(well_cell_i).GetCenter())
        well_cell_center[2] = height//2
        well = pv.Cylinder(center=well_cell_center, height=height, radius=model.wells[w]['r'], direction=(0,0,1))
        pl.add_mesh(well)
    pl.camera_position = 'xz'
    #p.show_bounds(grid='front', location='outer', all_edges=True)
    #_ = p.update_scalar_bar_range(0, model.pressures.max())
    pl.show_axes()
    if show_bounds:
        pl.show_bounds(grid='front', location='outer', all_edges=True) # or pl.show_axes() # or pl.show_grid()
    pl.set_background('black', top='gray')
    pl.show(full_screen=True)