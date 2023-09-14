"""conf file"""
from datetime import datetime

from tabulate import tabulate

import reservoirflow as rf

project = "ReservoirFlow"
author = "Zakariya Abugrin"
copyright = f"{datetime.now().year}, {author}"
version = rf.__version__
release = f"v{version}"
master_doc = "index"

extensions = [
    "sphinx.ext.autodoc",  # extract docs
    "sphinx.ext.napoleon",  # enhance parameters section
    "myst_nb",  # read md and ipynb files or "myst_parser",  # read md files
    "sphinx_gallery.gen_gallery",  # read py files as sphinx gallery
    # "sphinx.ext.autosummary",
    # "sphinx.ext.mathjax",  # read latex
]


# source_suffix = {
#     ".rst": "restructuredtext",
#     ".ipynb": "myst-nb",
#     ".myst": "myst-nb",
# }

templates_path = ["_templates"]
exclude_patterns = []

# myst:
myst_all_links_external = True
myst_admonition_enable = True
myst_amsmath_enable = True
myst_html_img_enable = True
myst_url_schemes = ("http", "https", "mailto")
# myst_enable_extensions = [
#     "amsmath",
#     "colon_fence",
#     "deflist",
#     "dollarmath",
#     "html_image",
# ]
# myst_url_schemes = ("http", "https", "mailto")

# myst-nb: https://myst-nb.readthedocs.io/en/latest/render/format_code_cells.html
# nb_number_source_lines = True
# nb_output_stderr = (
#     "show"  # ['show', 'remove', 'remove-warn', 'warn', 'error', 'severe']
# )
# nb_merge_streams = True
# nb_mime_priority_overrides = [
#     ("html", "text/plain", 0),
#     ("latex", "image/jpeg", None),
#     ("*", "customtype", 20),
# ]

# switcher: (rc: release candidate)
switcher_version = version
if ".dev" in version:
    switcher_version = "dev"
elif "rc" in version:
    switcher_version = version.split("rc", maxsplit=1)[0] + " (rc)"
html_static_path = ["_static", "_static/versions.json"]

# html:
html_theme = "pydata_sphinx_theme"
html_title = project  # + release
html_logo = "_static/logo.png"
html_favicon = "_static/logo_grid.png"
html_static_path = ["_static"]
html_show_sourcelink = False
html_theme_options = {
    # icons:
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/zakgrin/reservoirflow",
            "icon": "fab fa-github",
            "type": "fontawesome",
        },
        {
            "name": "LinkedIn",
            "url": "https://www.linkedin.com/in/zakariya-abugrin/",
            "icon": "fab fa-linkedin-in",
        },
    ],
    "footer_start": ["copyright"],  # remove PyData and Sphinx notes.
    "navbar_start": ["navbar-logo"],
    # switcher:
    "navbar_end": ["version-switcher", "theme-switcher", "navbar-icon-links"],
    "switcher": {
        "json_url": "_static/versions.json",
        "version_match": switcher_version,
    },
    "check_switcher": False,
}

# autodoc:
autoclass_content = "init"
autodoc_class_signature = "mixed"
autodoc_member_order = "bysource"
autodoc_docstring_signature = True
autodoc_typehints = "signature"
autodoc_typehints_description_target = "all"
autodoc_typehints_format = "short"
autodoc_inherit_docstrings = False
autodoc_preserve_defaults = False

# sphinx gallery
examples_dirs = ["user_guide/tutorials/example_sphinx_gallery"]
gallery_dirs = [d + "/build" for d in examples_dirs]
exclude_patterns.extend([d + "/*.ipynb" for d in gallery_dirs])
sphinx_gallery_conf = {
    "examples_dirs": examples_dirs,
    "gallery_dirs": gallery_dirs,
    "line_numbers": True,
    # "show_memory": True,
    # "notebook_images": True,
    "remove_config_comments": True,
    "first_notebook_cell": f"# {project} ({release}),  {copyright}",
    "last_notebook_cell": "# The End.",
    "capture_repr": ("_repr_html_", "__repr__", "__str__"),
    "default_thumb_file": "source/_static/logo_grid.png",
}


# Units and Factors:
def store_dict(in_dict, name="FACTORS"):
    name = name.upper()
    if name == "UNITS":
        label = "property"
    elif name == "FACTORS":
        label = "factor"
    else:
        raise ValueError("name is unknown")
    units = list(in_dict.keys())
    headers = [label] + units
    # headers = [h.capitalize() for h in headers]
    props = in_dict[units[0]]
    report = []
    for prop in props:
        if label == "property":
            prop_field = ":math:`" + in_dict[units[0]][prop] + "`"
            prop_metric = ":math:`" + in_dict[units[1]][prop] + "`"
            prop_lab = ":math:`" + in_dict[units[2]][prop] + "`"
        elif label == "factor":
            prop_field = in_dict[units[0]][prop]
            prop_metric = in_dict[units[1]][prop]
            prop_lab = in_dict[units[2]][prop]
        report.append(
            [
                prop,
                prop_field,
                prop_metric,
                prop_lab,
            ]
        )
    table = tabulate(report, headers=headers, showindex=False, tablefmt="rst")
    path = f"user_guide/units_factors/{name}.rst"
    with open(path, "w") as f:
        # name = name.capitalize()
        # f.write(f"{name}\n{'#'*len(name)}\n")
        f.writelines(table)


store_dict(rf.UNITS, "UNITS")
store_dict(rf.FACTORS, "FACTORS")
