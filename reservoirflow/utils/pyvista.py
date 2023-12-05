"""
pyvista
=======

This module is used to provide 3D shows for grids and model modules
using Pyvista. These functions are used as class methods to provide
direct accessibility to 3D visualizations using ``show`` attribute.
"""

import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv

from . import helpers

# -----------------------------------------------------------------------------
# BACKGROUND:
# -----------------------------------------------------------------------------

WINDOW_TITLE = "ReservoirFlow 3D Show"
BACKGROUND_COLOR = "black"
TEXT_COLOR = "white"


def set_background_color(color):
    global BACKGROUND_COLOR, TEXT_COLOR


def set_text_color(color):
    global BACKGROUND_COLOR, TEXT_COLOR


def set_mode(mode="dark"):
    """Set dark/light mode.

    Parameters
    ----------
    mode : str, optional
        mode as str in ["dark", "light"].

    Raises
    ------
    ValueError
        Unknown mode is used.
    """
    global BACKGROUND_COLOR, TEXT_COLOR
    if mode in ["dark", "black"]:
        BACKGROUND_COLOR = "black"
        TEXT_COLOR = "white"
    elif mode in ["light", "white"]:
        BACKGROUND_COLOR = "white"
        TEXT_COLOR = "black"
    else:
        raise ValueError("Unknown mode is used.")


def set_plotter_backend(static: bool) -> None:
    """Set backend for a pyvista plotter.

    Parameters
    ----------
    static : bool, optional
        show as a static image in a jupyter notebook. True
        value is used to render images for the documentation.
    """
    if static:
        pv.set_jupyter_backend("static")
    else:
        pv.set_jupyter_backend("trame")


def decide_widget(static: bool, notebook: bool) -> bool:
    """Decide to add widget or not.

    Parameters
    ----------
    static : bool
        show as a static image in a jupyter notebook. This argument
        is ignored when notebook argument is set to False. True
        value is used to render images for the documentation.
    notebook : bool
        show plot is placed inline a jupyter notebook. If False,
        then an interactive window will be opened outside of jupyter
        notebook.

    Returns
    -------
    bool
        add widget bool.
    """
    if static and notebook:
        return False
    else:
        return True


def set_plotter_config(
    pl: pv.Plotter,
    static: bool = False,
    notebook: bool = False,
) -> None:
    """Set the configuration for a pyvista plotter.

    Parameters
    ----------
    pl : pyvista.Plotter
        plotter object from pyvista.
    static : bool, optional
        show as a static image in a jupyter notebook. This argument
        is ignored when notebook argument is set to False. True
        value is used to render images for the documentation.
    notebook : bool, optional
        show plot is placed inline a jupyter notebook. If False,
        then an interactive window will be opened outside of jupyter
        notebook.
    """
    pl.set_background(BACKGROUND_COLOR)
    pl.enable_fly_to_right_click()
    pl.add_axes(color=TEXT_COLOR)
    if decide_widget(static, notebook):
        pl.add_camera_orientation_widget()


def get_text_locs(
    n: int,
    hloc: float,
    vloc: float,
    space: float,
) -> list:
    """Returns text locations for a pyvista plotter.

    Parameters
    ----------
    n : int
        number of text lines.
    hloc : float
        horizontal location of text lines.
    vloc : float
        vertical location of text lines.
    space : float
        spacing between text lines.

    Returns
    -------
    list
        list of tuples as [(hloc, vloc), ..] based on n.
    """
    vlocs = vloc - np.arange(0, space * n, space)
    return [(hloc, vloc) for vloc in vlocs]


def get_limits_fmt(values: np.ndarray) -> tuple:
    """Get values limits and format.

    Parameters
    ----------
    values : ndarray
        numpy array where values are used to define min and max.

    Returns
    -------
    tuple
        tuple(list, str) > ([min, max], "%.#f")
    """
    min_v = np.nanmin(values)
    max_v = np.nanmax(values)
    if min_v < 1:
        fmt = "%.2f"
    else:
        fmt = "%.f"
    return [min_v, max_v], fmt


def get_cdir(grid) -> str:
    """Get camera direction.

    Parameters
    ----------
    grid : Grid
        grid object from reservoirflow grids module.

    Returns
    -------
    str
        a string to indicate the flow direction based on xyz axis.
    """

    if grid.D == 1:
        if grid.fdir == "y":
            cdir = "yz"
        else:
            cdir = "xz"
    elif grid.D == 2:
        cdir = grid.fdir
    else:
        cdir = "xz"

    return cdir


def align_camera(
    pl: pv.Plotter,
    cdir: str = "xz",
    azimuth: float = 30,
    elevation: float = 30,
    zoom: float = 1.3,
) -> None:
    """Align the camera for a pyvista plotter.

    Parameters
    ----------
    pl : pv.Plotter
        plotter object from pyvista.
    cdir : str, optional
        plotter camera direction.
    azimuth : float, optional
        adjust camera azimuth which is a horizontal rotation around the
        central focal point, see `pyvista.Camera
        <https://docs.pyvista.org/version/stable/api/core/camera.html>`_.
    elevation : float, optional
        adjust camera elevation which is a vertical rotation around the
        central focal point, see `pyvista.Camera
        <https://docs.pyvista.org/version/stable/api/core/camera.html>`_.
    zoom : float, optional
        adjust camera azimuth which is a horizontal rotation around the
        central focal point, see `pyvista.Camera
        <https://docs.pyvista.org/version/stable/api/core/camera.html>`_.
    """
    pl.camera_position = cdir
    pl.camera.azimuth = azimuth
    pl.camera.elevation = elevation
    pl.camera.zoom(zoom)


def add_wells(
    pl: pv.Plotter,
    model,
) -> None:
    """Add wells to a pyvista plotter.

    Parameters
    ----------
    pl : pv.Plotter
        plotter object from pyvista.
    model : rf.models.Model
        a model object from models module.

    Notes
    -----
    Getting cells center using pyvista (you may need to convert to list):

    >>> model.grid.get_pyvista_grid(True).extract_cells(w).GetCenter()
    """
    cells_center = model.grid.get_cells_center(True, False, False).copy()
    for w in model.wells:
        height = model.grid.dz[w] * 10
        center = cells_center[w]
        center[2] += height // 2
        well = pv.Cylinder(
            center=center,
            direction=(0, 0, 1),
            radius=model.wells[w]["r"],
            height=height,
        )
        pl.add_mesh(well, render=False)


def get_cbar_dict(
    prop: str = "pressures",
    n_colors: int = 10,
    fmt: str = "%.f",
) -> dict:
    """Get color bar dictionary for a pyvista plotter.

    Parameters
    ----------
    prop : str, optional
        name of the property as a string.
    n_colors : int, optional
        number of colors of the color bar.
    fmt : str, optional
        text format of the color bar.

    Returns
    -------
    dict
        color bar configuration as a dictionary.
    """
    n_bins = n_colors + 1
    cbar_dict = dict(
        title=prop,
        n_labels=n_bins,
        title_font_size=24,
        label_font_size=18,
        color=TEXT_COLOR,
        font_family="arial",
        width=0.07,
        height=0.7,
        position_x=0.90,
        position_y=0.03,
        vertical=True,
        fmt=fmt,
        use_opacity=False,
        outline=False,
    )

    return cbar_dict


def get_colormap(
    cmap: str,
    gamma: float,
    n_colors: int,
):  # -> plt.cm:
    """Get colormap for a pyvista plotter from plt.cm

    Parameters
    ----------
    cmap : str, optional
        color map name based on Matplotlib, see
        `Choosing Colormaps in Matplotlib <https://matplotlib.org/stable
        /users/explain/colors/colormaps.html>`_.
    gamma : float, optional
        shift color map distribution to left when values less than 1 and
        to right when values larger than 1. In case of qualitative
        colormaps, this argument is ignored.
    n_colors : int, optional
        number of colors. In case of qualitative colormaps, n_colors
        should not exceed the total number of available colors (e.g. 10
        for cmap="tab10" and 20 for cmap="tab20").

    Returns
    -------
    plt.cm
        color map as a cm object from plt (matplotlib.pyplot).
    """
    qualitative_cmaps = [
        "Pastel1",
        "Pastel2",
        "Paired",
        "Accent",
        "Dark2",
        "Set1",
        "Set2",
        "Set3",
        "tab10",
        "tab20",
        "tab20b",
        "tab20c",
    ]

    colormap = plt.cm.get_cmap(cmap, n_colors)
    if cmap in qualitative_cmaps:
        cshape = np.unique(colormap.colors, axis=0).shape
        if cshape[0] < n_colors:
            print(f"Some colors are repetitive, consider n_colors<={cshape[0]}.")
    else:
        colormap.set_gamma(gamma)

    return colormap


def get_annotations(values: np.ndarray) -> dict:
    """Get

    Parameters
    ----------
    values : np.ndarray
        array of values.

    Returns
    -------
    dict
        annotations with values as keys and str as values.
    """
    return {
        np.nanmin(values): "min",
        np.nanmax(values): "max",
        np.nanmean(values): "avg",
    }


def get_window_size(resolution: str = "HD") -> tuple:
    """Get window size in pixels

    Parameters
    ----------
    resolution : str, optional
        resolution str in ["1K", "HD", "FHD", "2K", "4K", "8K", "10K"].

    Returns
    -------
    tuple
        window size in pixels.

    Raises
    ------
    ValueError
        Unknown resolution value.
    """
    resolution = resolution.lower()
    if resolution in ["1k"]:
        return (1024, 768)
    elif resolution in ["hd"]:
        return (1280, 720)
    elif resolution in ["full hd", "fhd"]:
        return (1920, 1080)
    elif resolution in ["2k"]:
        return (2048, 1080)
    elif resolution in ["4k"]:
        return (3840, 2160)
    elif resolution in ["8k"]:
        return (7680, 3420)
    elif resolution in ["10k"]:
        return (10240, 4320)
    else:
        raise ValueError("Unknown resolution value.")


def add_ruler(
    pl: pv.Plotter,
    mesh,
) -> None:
    """Add a ruler to a pyvista plotter

    Parameters
    ----------
    pl : pv.Plotter
        plotter object from pyvista.
    mesh :
        mesh object from pyvista.
    """
    pl.show_bounds(
        mesh=mesh,
        font_family="arial",
        grid="front",
        location="outer",
        ticks="outside",
        all_edges=False,
        padding=0.001,
    )


def add_title(pl: pv.Plotter) -> None:
    """Add a title to a pyvista plotter

    Parameters
    ----------
    pl : pv.Plotter
        plotter object from pyvista.
    """
    pl.add_text(
        text="ReservoirFlow",
        position=(0.45, 0.95),
        font_size=14,
        font="arial",
        color=TEXT_COLOR,
        viewport=True,
    )


def add_desc(
    pl: pv.Plotter,
    label: str,
    D: int,
    fdir: str,
    boundary: bool,
) -> None:
    """Add a description to a pyvista plotter

    Parameters
    ----------
    pl : pv.Plotter
        plotter object from pyvista.
    label : str
        label as str.
    D : int
        number of flow dimensions.
    fdir : str
        flow direction.
    boundary : bool
        add boundary string.
    """
    s = helpers.get_boundary_str(boundary)
    desc = "{}D grid model by {} (flow at {}-direction {})".format(
        D,
        label,
        fdir,
        s,
    )
    pl.add_text(
        text=desc,
        position=(0.01, 0.01),
        font_size=10,
        font="arial",
        color=TEXT_COLOR,
        viewport=True,
    )


def add_grid_labels(
    pl: pv.Plotter,
    grid,
    label: str,
    opacity: float,
    boundary: bool,
) -> None:
    """Add labels to a grid

    Parameters
    ----------
    pl : pv.Plotter
        plotter object from pyvista.
    grid :
        grid object from reservoirflow grids module.
    label : str
        label of grid centers as str in ['id', 'coords', 'icoords',
        'dx', 'dy', 'dz', 'Ax', 'Ay', 'Az', 'V', 'center', 'sphere']. If
        None, this nothing will be added.
    opacity : float
        adjust transparency between 0 and 1 where 0 is fully transparent
        and 1 is fully nontransparent.
    boundary : bool
        include boundary cells.

    Raises
    ------
    ValueError
        label is not recognized.
    """
    points = grid.get_cells_center(boundary, False, False)
    if opacity < 1 and not label in [None, False]:
        if label == "coords":
            labels = grid.get_cells_coords(boundary, False, "tuple")
        elif label == "icoords":
            labels = grid.get_cells_icoords(boundary, False, "tuple")
        elif label == "id":
            labels = grid.get_cells_id(boundary, False, "tuple")
        elif label == "i":
            labels = grid.get_cells_i(boundary, False, "tuple")
        elif label == "dx":
            labels = grid.get_cells_dx(boundary, False)
        elif label == "dy":
            labels = grid.get_cells_dy(boundary, False)
        elif label == "dz":
            labels = grid.get_cells_dz(boundary, False)
        elif label in ["area_x", "Ax"]:
            labels = grid.get_cells_Ax(boundary, False)
        elif label in ["area_y", "Ay"]:
            labels = grid.get_cells_Ay(boundary, False)
        elif label in ["area_z", "Az"]:
            labels = grid.get_cells_Az(boundary, False)
        elif label in ["volume", "V"]:
            labels = grid.get_cells_V(boundary, False, False)
        elif label in ["center", "centers", "sphere"]:
            labels = points
        else:
            raise ValueError(f"label='{label}' is not recognized.")

        if label == "sphere":
            pl.add_points(
                points,
                point_size=10,
                render_points_as_spheres=True,
                show_edges=True,
                color=TEXT_COLOR,
            )
        else:
            pl.add_point_labels(
                points=points,
                labels=labels,
                font_size=10,
                text_color=TEXT_COLOR,
                point_size=10,
            )


def get_grid_plotter(
    grid,
    prop,
    label,
    boundary,
    corners,
    desc,
    # cbar,
    cmap,
    gamma,
    n_colors,
    opacity,
    azimuth,
    elevation,
    zoom,
    static,
    notebook,
    window_size,
    # prop: str = "pressures",
    # label: str = None,
    # boundary: bool = False,
    # corners: bool = False,
    # # cbar: bool = True,
    # cmap: str = "Blues",
    # gamma: float = 0.7,
    # n_colors: int = 10,
    # opacity: float = 0.9,
    # azimuth: float = 45,
    # elevation: float = 40,
    # zoom: float = 1,
    # static: bool = False,
    # notebook: bool = False,
    # window_size: tuple = None,
    **kwargs,
):
    # values = get_values(model, prop, boundary)
    # limits, fmt = get_limits_fmt(values)
    colormap = get_colormap(cmap, gamma, n_colors)
    # annotations = get_annotations(values)

    if static and window_size is None:
        window_size = get_window_size("fhd")

    pl = pv.Plotter(
        notebook=notebook,
        window_size=window_size,
        **kwargs,
    )

    grid_pv = grid.get_pyvista_grid(boundary)

    grid_mesh = dict(
        mesh=grid_pv,
        # clim=limits,
        show_edges=True,
        opacity=opacity,
        nan_color="gray",
        nan_opacity=0.05,
        lighting=True,
        colormap=colormap,
        color="white",
        # scalars=values[0].copy(),
        show_scalar_bar=False,
        # annotations=annotations,
    )

    if corners:
        points = grid_pv.points
        pl.add_point_labels(
            points=points,
            labels=points.tolist(),
            font_size=10,
            text_color=TEXT_COLOR,
            point_size=10,
        )
    pl.add_mesh(**grid_mesh)

    # if cbar:
    #     cbar_dict = get_cbar_dict(property.capitalize(), n_colors, fmt)
    #     pl.add_scalar_bar(**cbar_dict)

    cdir = get_cdir(grid)
    align_camera(pl, cdir, azimuth=azimuth, elevation=elevation, zoom=zoom)
    set_plotter_config(pl, static=static, notebook=notebook)
    add_grid_labels(pl, grid, label, opacity, boundary)
    add_title(pl)
    if desc:
        add_desc(pl, label, grid.D, grid.fdir, boundary)

    return pl, grid_pv


def show_grid(
    grid,
    prop: str = "pressures",
    label: str = None,
    boundary: bool = False,
    # wells: bool = True,
    corners: bool = False,
    ruler: bool = False,
    desc: bool = False,
    # info: bool = True,
    # cbar: bool = True,
    cmap: str = "Blues",
    gamma: float = 0.7,
    n_colors: int = 10,
    opacity: float = 0.9,
    azimuth: float = 45,
    elevation: float = 40,
    zoom: float = 1,
    static: bool = False,
    notebook: bool = False,
    window_size: tuple = None,
    **kwargs,
):
    """Shows pyvista grid.

    This function shows the grid using pyvista object in 3D. Only if
    the total number of cells is lower than 20, then the grid will
    be transparent. Therefore, to be able debug your model, try to
    first test a small model.

    Parameters
    ----------
    label : str, optional
        label of grid centers as str in ['id', 'coords', 'icoords',
        'dx', 'dy', 'dz', 'Ax', 'Ay', 'Az', 'V', 'center', 'sphere']. If
        None, this nothing will be added.
    boundary : bool, optional
        include boundary cells.
    corners : bool, optional
        show corners ijk values.
    static : bool, optional
        show as a static image in a jupyter notebook. This argument
        is ignored when notebook argument is set to False. True
        value is used to render images for the documentation.
    notebook : bool, optional
        show plot is placed inline a jupyter notebook. If False,
        then an interactive window will be opened outside of jupyter
        notebook.
    azimuth : float, optional
        adjust camera azimuth which is a horizontal rotation around the
        central focal point, see
        `pyvista.Camera <https://docs.pyvista.org/version/stable/api/core/camera.html>`_.
    elevation : float, optional
        adjust camera elevation which is a vertical rotation around the
        central focal point, see
        `pyvista.Camera <https://docs.pyvista.org/version/stable/api/core/camera.html>`_.
    zoom : float, optional
        adjust camera zoom which is direct zooming into the central
        focal point. see
        `pyvista.Camera.zoom <https://docs.pyvista.org/version/stable/api/core/_autosummary/pyvista.Camera.zoom.html>`_.
    **kwargs :
        you can pass any argument for pyvista Plotter except if it is
        defined as argument in this function (e.g. window_size), see
        `pyvista.Plotter <https://docs.pyvista.org/version/stable/api/plotting/_autosummary/pyvista.Plotter.html>`_.

    Raises
    ------
    ValueError
        label is not recognized.
    """
    set_plotter_backend(static)
    pl, grid_pv = get_grid_plotter(
        grid,
        prop=prop,
        label=label,
        boundary=boundary,
        # wells=wells,
        corners=corners,
        desc=desc,
        # cbar=cbar,
        cmap=cmap,
        gamma=gamma,
        n_colors=n_colors,
        opacity=opacity,
        azimuth=azimuth,
        elevation=elevation,
        zoom=zoom,
        static=static,
        notebook=notebook,
        window_size=window_size,
        **kwargs,
    )

    if ruler:
        add_ruler(pl, grid_pv)

    pl.show(title=WINDOW_TITLE)


def get_model_values(model, prop, boundary):
    prop = prop.lower()
    if prop in ["p", "pressure", "pressures"]:
        values = model.pressures
    elif prop in ["q", "rate", "rates"]:
        values = model.rates
    else:
        raise ValueError("Unknown property.")

    if not boundary:
        cells_id = model.grid.get_cells_id(False, False, "list")
        return values[:, cells_id]
    return values


def get_model_plotter(
    model,
    prop,
    label,
    boundary,
    wells,
    desc,
    cbar,
    cmap,
    gamma,
    n_colors,
    opacity,
    azimuth,
    elevation,
    zoom,
    static,
    notebook,
    window_size,
    # prop: str = "pressures",
    # label: str = None,
    # boundary: bool = False,
    # wells: bool = True,
    # cbar: bool = True,
    # cmap: str = "Blues",
    # gamma: float = 0.7,
    # n_colors: int = 10,
    # opacity: float = 0.9,
    # azimuth: float = 45,
    # elevation: float = 40,
    # zoom: float = 1,
    # static: bool = False,
    # notebook: bool = False,
    # window_size: tuple = None,
    **kwargs,
):
    values = get_model_values(model, prop, boundary)
    limits, fmt = get_limits_fmt(values)
    colormap = get_colormap(cmap, gamma, n_colors)
    annotations = get_annotations(values)

    if static and window_size is None:
        window_size = get_window_size("fhd")

    pl = pv.Plotter(
        notebook=notebook,
        window_size=window_size,
        **kwargs,
    )

    grid_pv = model.grid.get_pyvista_grid(boundary)

    grid_mesh = dict(
        mesh=grid_pv,
        clim=limits,
        show_edges=True,
        opacity=opacity,
        nan_color="gray",
        nan_opacity=0.05,
        lighting=True,
        colormap=colormap,
        color="white",
        scalars=values[0].copy(),
        show_scalar_bar=False,
        annotations=annotations,
    )

    if wells:
        add_wells(pl, model)
    pl.add_mesh(**grid_mesh)

    if cbar:
        cbar_dict = get_cbar_dict(prop.capitalize(), n_colors, fmt)
        pl.add_scalar_bar(**cbar_dict)

    cdir = get_cdir(model.grid)
    align_camera(pl, cdir, azimuth=azimuth, elevation=elevation, zoom=zoom)
    set_plotter_config(pl, static=static, notebook=notebook)
    add_grid_labels(pl, model.grid, label, opacity, boundary)
    add_title(pl)
    if desc:
        add_desc(pl, label, model.grid.D, model.grid.fdir, boundary)

    return pl, grid_pv


def save_gif(
    model,
    prop: str = "pressures",
    label: str = None,
    boundary: bool = False,
    wells: bool = True,
    ruler: bool = False,
    desc: bool = False,
    info: bool = True,
    cbar: bool = True,
    cmap: str = "Blues",
    gamma: float = 0.7,
    n_colors: int = 10,
    opacity: float = 0.9,
    azimuth: float = 45,
    elevation: float = 40,
    zoom: float = 1,
    fps: int = 10,
    file_name: str = "grid_animated.gif",
    window_size: tuple = None,  # (2048, 1080),
    **kwargs,
):
    """Saves pyvista show as gif

    Parameters
    ----------
    prop : str, optional
        property to be visualized in ["pressures", "rates"].
    boundary : bool, optional
        include boundary cells.
    wells : bool, optional
        show wells.
    ruler : bool, optional
        show ruler grid lines.
    info : bool, optional
        show simulation information.
    cbar : bool, optional
        show color bar.
    cmap : str, optional
        color map name based on Matplotlib, see
        `Choosing Colormaps in Matplotlib <https://matplotlib.org/stable/users/explain/colors/colormaps.html>`_.
    gamma : float, optional
        shift color map distribution to left when values less than 1 and
        to right when values larger than 1. In case of qualitative
        colormaps, this argument is ignored.
    n_colors : int, optional
        number of colors. In case of qualitative colormaps, n_colors
        should not exceed the total number of available colors (e.g. 10
        for cmap="tab10" and 20 for cmap="tab20").
    opacity : float, optional
        adjust transparency between 0 and 1 where 0 is fully transparent
        and 1 is fully nontransparent.
    azimuth : float, optional
        adjust camera azimuth which is a horizontal rotation around the
        central focal point, see
        `pyvista.Camera <https://docs.pyvista.org/version/stable/api/core/camera.html>`_.
    elevation : float, optional
        adjust camera elevation which is a vertical rotation around the
        central focal point, see
        `pyvista.Camera <https://docs.pyvista.org/version/stable/api/core/camera.html>`_.
    zoom : float, optional
        adjust camera zoom which is direct zooming into the central
        focal point. see
        `pyvista.Camera.zoom <https://docs.pyvista.org/version/stable/api/core/_autosummary/pyvista.Camera.zoom.html>`_.
    fps : int, optional
        the number of frames per second starting at 1 and higher. Use
        this to control the speed of the gif image. Note that at some
        point, higher fps will not have affect on speed.
    fine_name : str, optional
        file name of the gif file (including the extension .gif).
    window_size : tuple, optional
        pyvista plotter window size in pixels.
    **kwargs :
        you can pass any argument for pyvista Plotter except if it is
        defined as argument in this function (e.g. window_size), see
        `pyvista.Plotter <https://docs.pyvista.org/version/stable/api/plotting/_autosummary/pyvista.Plotter.html>`_.

    Raises
    ------
    ValueError
        label is not recognized.

    Warnings
    --------
    Jittering :

    There is still jittering affect especially when a
    continues colors are used. In addition, the same affect appears
    on side 0 of cell 0 even with discretized colors.
    This is a common issue in pyvista also causing issus with color
    bar. So far there is no solution to this issue. The affect might
    be mitigated by using discretized colors or by making the first
    cell or layer nontransparent.

    .. highlight:: python
    .. code-block:: python

        # making cell 0 nontransparent:
        base_cells = 0
        # or layer 0:
        base_cells = model.grid.get_cells_i(boundary, True)[0, :, :]
        # add mesh:
        pl.add_mesh(
            grid_pv.extract_cells(base_cells),
            clim=limits,
            show_edges=True,
            opacity=1,
            lighting=True,
            colormap=colormap,
            show_scalar_bar=False,
        )
    """

    values = get_model_values(model, prop, boundary)

    pl, grid_pv = get_model_plotter(
        model,
        prop=prop,
        label=label,
        boundary=boundary,
        wells=wells,
        desc=desc,
        cbar=cbar,
        cmap=cmap,
        gamma=gamma,
        n_colors=n_colors,
        opacity=opacity,
        azimuth=azimuth,
        elevation=elevation,
        zoom=zoom,
        static=True,
        notebook=True,
        window_size=window_size,
        **kwargs,
    )

    if info:
        font = dict(font="courier", color=TEXT_COLOR)
        conf = dict(font_size=11, viewport=True, render=False)
        locs = get_text_locs(6, 0.01, 0.95, 0.03)
        dates = model.get_df(columns=["Date"])["Date"]
        date = pl.add_text("", position=locs[0], **font, **conf)
        time_step = pl.add_text("", position=locs[1], **font, **conf)
        value_init = f"initial pressure: {model.pi}"
        pl.add_text(value_init, position=locs[2], **font, **conf)
        value_avg = pl.add_text("", position=locs[3], **font, **conf)
        value_min = pl.add_text("", position=locs[4], **font, **conf)
        value_max = pl.add_text("", position=locs[5], **font, **conf)
    if ruler:
        add_ruler(pl, grid_pv)

    pl.open_gif(
        filename=file_name,
        loop=0,
        fps=fps,
        palettesize=256,
        subrectangles=False,
    )

    for tstep in range(model.nsteps):
        if info:
            date.SetInput(f"    current date: {dates[tstep]}")
            time_step.SetInput(f"current timestep: {tstep:02d}")
            value_avg.SetInput(f"    avg pressure: {np.nanmean(values[tstep]):.0f}")
            value_min.SetInput(f"    min pressure: {np.nanmin(values[tstep]):.0f}")
            value_max.SetInput(f"    max pressure: {np.nanmax(values[tstep]):.0f}")
        pl.update_scalars(values[tstep], grid_pv, render=False)
        pl.write_frame()

    pl.close()


def show_model(
    model,
    prop: str = "pressures",
    label: str = None,
    boundary: bool = False,
    wells: bool = True,
    ruler: bool = False,
    desc: bool = False,
    info: bool = True,
    cbar: bool = True,
    cmap: str = "Blues",
    gamma: float = 0.7,
    n_colors: int = 10,
    opacity: float = 0.9,
    azimuth: float = 45,
    elevation: float = 40,
    zoom: float = 1,
    static: bool = False,
    notebook: bool = False,
    window_size: tuple = None,
    **kwargs,
):
    """Show pyvista plotter

    Parameters
    ----------
    prop : str, optional
        property to be visualized in ["pressures", "rates"].
    label : str, optional
        label of grid centers as str in ['id', 'coords', 'icoords',
        'dx', 'dy', 'dz', 'Ax', 'Ay', 'Az', 'V', 'center', 'sphere']. If
        None or False, nothing will be added as a label.
    boundary : bool, optional
        include boundary cells.
    wells : bool, optional
        show wells.
    ruler : bool, optional
        show ruler grid lines.
    info : bool, optional
        show simulation information.
    cbar : bool, optional
        show color bar.
    cmap : str, optional
        color map name based on Matplotlib, see
        `Choosing Colormaps in Matplotlib <https://matplotlib.org/stable/users/explain/colors/colormaps.html>`_.
    gamma : float, optional
        shift color map distribution to left when values less than 1 and
        to right when values larger than 1. In case of qualitative
        colormaps, this argument is ignored.
    n_colors : int, optional
        number of colors. In case of qualitative colormaps, n_colors
        should not exceed the total number of available colors (e.g. 10
        for cmap="tab10" and 20 for cmap="tab20").
    opacity : float, optional
        adjust transparency between 0 and 1 where 0 is fully transparent
        and 1 is fully nontransparent.
    azimuth : float, optional
        adjust camera azimuth which is a horizontal rotation around the
        central focal point, see
        `pyvista.Camera <https://docs.pyvista.org/version/stable/api/core/camera.html>`_.
    elevation : float, optional
        adjust camera elevation which is a vertical rotation around the
        central focal point, see
        `pyvista.Camera <https://docs.pyvista.org/version/stable/api/core/camera.html>`_.
    zoom : float, optional
        adjust camera zoom which is direct zooming into the central
        focal point. see
        `pyvista.Camera.zoom <https://docs.pyvista.org/version/stable/api/core/_autosummary/pyvista.Camera.zoom.html>`_.
    static : bool, optional
        show as a static image in a jupyter notebook. This argument
        is ignored when notebook argument is set to False. True
        value is used to render images for the documentation.
    notebook : bool, optional
        show plot is placed inline a jupyter notebook. If False,
        then an interactive window will be opened outside of jupyter
        notebook.
    window_size : tuple, optional
        pyvista plotter window size in pixels.
    **kwargs :
        you can pass any argument for pyvista Plotter except if it is
        defined as argument in this function (e.g. window_size), see
        `pyvista.Plotter <https://docs.pyvista.org/version/stable/api/plotting/_autosummary/pyvista.Plotter.html>`_

    Raises
    ------
    ValueError
        label is not recognized.

    Notes
    -----
    Key callbacks:

    .. highlight:: python
    .. code-block:: python

        def my_cpos_callback():
            pl.add_text(text=str(pl.camera.position), name="cpos")
        pl.add_key_event("p", my_cpos_callback)
    """
    set_plotter_backend(static)
    pl, grid_pv = get_model_plotter(
        model,
        prop=prop,
        label=label,
        boundary=boundary,
        wells=wells,
        desc=desc,
        cbar=cbar,
        cmap=cmap,
        gamma=gamma,
        n_colors=n_colors,
        opacity=opacity,
        azimuth=azimuth,
        elevation=elevation,
        zoom=zoom,
        static=static,
        notebook=notebook,
        window_size=window_size,
        **kwargs,
    )

    if info:
        font = dict(font="courier", color=TEXT_COLOR)
        conf = dict(font_size=11, viewport=True, render=False)
        locs = get_text_locs(6, 0.01, 0.95, 0.03)
        dates = model.get_df(columns=["Date"])["Date"]
        date = pl.add_text("", position=locs[0], **font, **conf)
        time_step = pl.add_text("", position=locs[1], **font, **conf)
        value_init = f"initial pressure: {model.pi}"
        pl.add_text(value_init, position=locs[2], **font, **conf)
        value_avg = pl.add_text("", position=locs[3], **font, **conf)
        value_min = pl.add_text("", position=locs[4], **font, **conf)
        value_max = pl.add_text("", position=locs[5], **font, **conf)

    if ruler:
        add_ruler(pl, grid_pv)

    def update_pl(tstep):
        tstep = int(tstep)
        if info:
            date.SetInput(f"    current date: {dates[tstep]}")
            time_step.SetInput(f"current timestep: {tstep:02d}")
            value_avg.SetInput(f"    avg pressure: {np.nanmean(values[tstep]):.0f}")
            value_min.SetInput(f"    min pressure: {np.nanmin(values[tstep]):.0f}")
            value_max.SetInput(f"    max pressure: {np.nanmax(values[tstep]):.0f}")
        pl.update_scalars(values[tstep], grid_pv, render=False)
        return pl

    values = get_model_values(model, prop, boundary)
    last_tstep = model.nsteps - 1
    pl = update_pl(last_tstep)

    if not static or not notebook:
        pl.add_slider_widget(
            callback=update_pl,
            rng=(0, last_tstep),
            value=last_tstep,
            title="Timestep",
            pointa=(0.40, 0.90),
            pointb=(0.60, 0.90),
            color=TEXT_COLOR,
            interaction_event="always",
            style="modern",
            title_height=0.025,
            title_opacity=0.7,
            title_color=TEXT_COLOR,
            fmt="%.0f",
            slider_width=0.02,
            tube_width=0.01,
        )

    pl.show(title=WINDOW_TITLE)
