"""conf file"""
# import os
# import sys
from datetime import datetime

from tabulate import tabulate

import reservoirflow as rf

# sys.path.insert(0, os.path.abspath("../"))

project = "ReservoirFlow"
author = "Bayanatics"
copyright = f"{datetime.now().year}, {author}"
version = rf.__version__
release = f"v{version}"
master_doc = "index"

# remove warnings:
suppress_warnings = [
    "myst.header",
    "myst.footnote",
    "toc.excluded",
]

# source_suffix = {
#     ".rst": "restructuredtext",
#     # ".ipynb": "myst-nb",
#     # ".myst": "myst-nb",
#     ".py": "myst-nb",
# }

templates_path = ["_templates"]
exclude_patterns = [
    "build",
    "**/example_*",
    "_api_backup/*",
    "_api_backup1/*",
    "api/api.rst",
]

# switcher: (rc: release candidate)
switcher_version = version
json_url = "_static/versions.json"
if ".dev" in version:
    switcher_version = "dev"
elif "rc" in version:
    switcher_version = version.split("rc", maxsplit=1)[0] + " (rc)"

extensions = [
    "sphinx.ext.autodoc",  # build docstring contents.
    "sphinx.ext.autosummary",  # build docstring pages.
    "sphinx.ext.napoleon",  # using sections in docstring.
    "sphinx.ext.doctest",
    # "sphinx.ext.autosectionlabel",
    "sphinx.ext.todo",  # allow .. todo:: directive.
    "myst_nb",  # read md and ipynb files or "myst_parser",  # read md files
    # "sphinx_gallery.gen_gallery",  # read py files as sphinx gallery
    "sphinx_comments",  # allow comments
    # "numpydoc",  # numpydoc: https://numpydoc.readthedocs.io/en/latest/format.html
    # "autodoc2",  # markdown in docstring: https://sphinx-autodoc2.readthedocs.io/en/latest/quickstart.html
    "sphinx_design",  # https://sphinx-design.readthedocs.io/en/latest/get_started.html
    "sphinxcontrib.bibtex",  # https://sphinxcontrib-bibtex.readthedocs.io/en/latest/quickstart.html
]

# config: https://www.sphinx-doc.org/en/master/usage/configuration.html
add_module_names = True
toc_object_entries_show_parents = "hide"  # domain, hide, all "On this page"

# todo:
todo_include_todos = True

# bibtex:
bibtex_bibfiles = [
    "research_development/references/books.bib",
    "research_development/references/papers.bib",
]
bibtex_default_style = "unsrt"
bibtex_reference_style = "super"
bibtex_foot_reference_style = "foot"

# autosummary:
autosummary_generate = True
autosummary_generate_overwrite = True
autosummary_imported_members = True
autosummary_ignore_module_all = False
autosummary_filename_map = {"reservoirflow": "API"}
autoclass_content = "both"

# autodoc: managed by autosummary templates.
# autodoc_class_signature = "mixed"
# autodoc_member_order = "bysource"
# autodoc_docstring_signature = False  # change to true
# autodoc_typehints = "both"
# autodoc_typehints_description_target = "all"
# autodoc_typehints_format = "short"  # short, fully-qualified
# autodoc_preserve_defaults = True
# autodoc_inherit_docstrings = True
# autodoc_default_options = {
#     "members": True,
#     "member-order": "bysource",
#     "undoc-members": True,
#     "inherited-members": True,
#     "show-inheritance": False,
#     "imported-members": False,
#     # "special-members": "__init__",
#     # "exclude-members": "__weakref__",
# }


# myst and myst-nb: https://myst-nb.readthedocs.io/en/latest/render/format_code_cells.html
myst_all_links_external = True
myst_admonition_enable = True
myst_amsmath_enable = True
myst_html_img_enable = True
myst_url_schemes = ("http", "https", "mailto")
myst_enable_extensions = [
    "amsmath",
    "colon_fence",
    "deflist",
    "dollarmath",
    "html_image",
    "fieldlist",
]
myst_footnote_transition = False

# myst-nb: https://myst-nb.readthedocs.io/en/latest/configuration.html
# nb_output_stderr = "remove"  # remove progress bar
nb_merge_streams = True  # combine print output in one cell

# numpydoc: https://numpydoc.readthedocs.io/en/latest/install.html
# numpydoc_use_plots = True
# numpydoc_show_class_members = True  # table for attributes and methods.
# numpydoc_class_members_toctree = False
# numpydoc_xref_param_type = False

# autodoc2
# autodoc2_packages = [
#     r"..\..\\reservoirflow",
# ]
# autodoc2_docstring_parser_regexes = [
#     # this will render all docstrings as Markdown
#     (r".*", "myst"),
# ]

# sphinx gallery:
# examples_dirs = ["user_guide/tutorials/example_sphinx_gallery"]
# gallery_dirs = [d + "/build" for d in examples_dirs]
# exclude_patterns.extend([d + "/*.ipynb" for d in gallery_dirs])
# sphinx_gallery_conf = {
#     "examples_dirs": examples_dirs,
#     "gallery_dirs": gallery_dirs,
#     "notebook_images": True,
#     "remove_config_comments": True,
#     "first_notebook_cell": f"# {project} ({release}),  {copyright}",
#     "last_notebook_cell": "# The End.",
#     "capture_repr": ("_repr_html_", "__repr__"),
#     "default_thumb_file": "source/_static/logo_grid.png",
# }

# html:
html_theme = "pydata_sphinx_theme"
html_title = project  # + release
html_logo = "_static/RF_logo.png"
html_favicon = "_static/RF_logo.png"
html_show_sourcelink = False
html_static_path = [
    "_static",
]
html_css_files = [
    "css/custom.css",
]
html_theme_options = {
    # "header_links_before_dropdown": 7,
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
    "footer_start": ["copyright"],
    "footer_center": "",
    "footer_end": "",
    "navbar_start": [
        "navbar-logo",
        "version-switcher",
    ],
    "navbar_center": [
        "navbar-nav",
    ],
    "navbar_end": [
        "theme-switcher",
        "navbar-icon-links",
    ],
    "switcher": {
        "json_url": json_url,
        "version_match": switcher_version,
    },
    "check_switcher": True,
    "show_version_warning_banner": True,
    "navigation_with_keys": False,
    "show_toc_level": 2,
    "secondary_sidebar_items": ["page-toc"],
}

# utteranc.es:
comments_config = {
    "utterances": {
        "repo": "zakgrin/reservoirflow_utterances",
    }
}


# Units and Factors:
def store_dict(in_dict, name="FACTORS", folder=""):
    name = name.upper()
    if name == "UNITS":
        label = "property"
    elif name == "FACTORS":
        label = "factor"
    elif name == "NOMENCLATURE":
        label = "nomenclature"
    else:
        raise ValueError("name is unknown")
    columns = list(in_dict.keys())
    headers = [label] + columns
    # headers = [h.capitalize() for h in headers]
    props = in_dict[columns[0]]
    report = []

    for prop in props:
        row = []
        if label == "property":
            for i in range(len(columns)):
                row.append(":math:`" + in_dict[columns[i]][prop] + "`")
            # row.append(":math:`" + in_dict[columns[0]][prop] + "`")
            # row.append(":math:`" + in_dict[columns[1]][prop] + "`")
            # row.append(":math:`" + in_dict[columns[2]][prop] + "`")
            # prop_field = ":math:`" + in_dict[columns[0]][prop] + "`"
            # prop_metric = ":math:`" + in_dict[columns[1]][prop] + "`"
            # prop_lab = ":math:`" + in_dict[columns[2]][prop] + "`"
        elif label == "factor":
            for i in range(len(columns)):
                row.append(in_dict[columns[i]][prop])
            # prop_field = in_dict[columns[0]][prop]
            # prop_metric = in_dict[columns[1]][prop]
            # prop_lab = in_dict[columns[2]][prop]
        elif label == "nomenclature":
            row.append(in_dict[columns[0]][prop])
            row.append(":math:`" + in_dict[columns[1]][prop] + "`")

            # prop_field = ":math:`" + in_dict[columns[0]][prop] + "`"
            # prop_metric = ":math:`" + in_dict[columns[1]][prop] + "`"

        report.append(
            [
                prop,
                *row,
                # prop_field,
                # prop_metric,
                # prop_lab,
            ]
        )
    table = tabulate(report, headers=headers, showindex=False, tablefmt="rst")
    path = f"{folder+name.lower()}_table.rst"
    with open(path, "w") as f:
        f.writelines(table)


store_dict(rf.UNITS, "UNITS", "user_guide/units_factors/")
store_dict(rf.FACTORS, "FACTORS", "user_guide/units_factors/")
store_dict(rf.NOMENCLATURE, "NOMENCLATURE", "user_guide/nomenclature/")
