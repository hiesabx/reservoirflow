from datetime import datetime
from io import StringIO

from tabulate import SEPARATING_LINE, tabulate

import reservoirflow as rf

project = "ReservoirFlow"
copyright = f"{datetime.now().year}, Zakariya Abugrin"
author = "Zakariya Abugrin"
release = f" v{rf.__version__}"
master_doc = "index"

extensions = [
    "sphinx.ext.autodoc",  # extract docs
    "sphinx.ext.napoleon",  # enhance parameters section
    "myst_parser",  # add md files
]

templates_path = ["_templates"]
exclude_patterns = []

# myst:
myst_all_links_external = True

# html:
html_theme = "pydata_sphinx_theme"
html_title = project  # + release
html_logo = "_static/logo.png"
html_favicon = "_static/logo_grid.png"
html_static_path = ["_static"]
html_show_sourcelink = False
html_theme_options = {
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
        report.append(
            [
                prop,
                in_dict[units[0]][prop],
                in_dict[units[1]][prop],
                in_dict[units[2]][prop],
            ]
        )
    table = tabulate(report, headers=headers, showindex=False, tablefmt="rst")
    path = f"user_guide/{name}.rst"
    with open(path, "w") as f:
        f.writelines(table)


store_dict(rf.UNITS, "UNITS")
store_dict(rf.FACTORS, "FACTORS")
