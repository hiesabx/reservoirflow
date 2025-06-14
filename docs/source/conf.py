"""conf file"""

# import os
# import sys
from datetime import datetime

from tabulate import tabulate

import reservoirflow as rf

# sys.path.insert(0, os.path.abspath("../"))

project = "ReservoirFlow"
author = "Developed by Hiesab"
disclaimer = "Third-party components are copyrighted by their respective authors"
copyright = f"{datetime.now().year}, {author}. {disclaimer}"
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
]

# switcher:
json_url = "_static/versions.json"
if any(v in version for v in ["dev", "a", "b", "rc"]):
    # all pre-release stages are called dev
    switcher_version = "dev"
else:
    # only keep major.minor version number to match versions.json
    switcher_version = ".".join(version.split(".")[:2])

extensions = [
    "sphinx.ext.autodoc",  # build docstring contents.
    "sphinx.ext.autosummary",  # build docstring pages.
    "sphinx.ext.napoleon",  # using sections in docstring.
    "sphinx.ext.doctest",
    # "sphinx.ext.autosectionlabel",
    "sphinx.ext.todo",  # allow .. todo:: directive.
    "myst_nb",  # read md and ipynb files or "myst_parser",  # read md files
    # "sphinx_gallery.gen_gallery",  # read py files as sphinx gallery
    # Admonition: https://myst-parser.readthedocs.io/en/latest/syntax/admonitions.html
    "sphinx_comments",  # allow comments
    # "numpydoc",  # numpydoc: https://numpydoc.readthedocs.io/en/latest/format.html
    # "autodoc2",  # markdown in docstring: https://sphinx-autodoc2.readthedocs.io/en/latest/quickstart.html
    "sphinx_design",  # https://sphinx-design.readthedocs.io/en/latest/get_started.html
    "sphinxcontrib.bibtex",  # https://sphinxcontrib-bibtex.readthedocs.io/en/latest/quickstart.html
    "sphinx_copybutton",  # https://sphinx-copybutton.readthedocs.io/en/latest/
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
myst_footnote_transition = False  # solves the issue with references in nb.
myst_highlight_code_blocks = False

# myst-nb: https://myst-nb.readthedocs.io/en/latest/configuration.html
nb_output_stderr = "remove"  # 'remove' progress bar
nb_merge_streams = True  # combine print output in one cell
nb_execution_mode = "auto"
nb_number_source_lines = True

# copybutton
copybutton_selector = ":not(.prompt) > div.highlight pre"
copybutton_exclude = ".linenos, .gp, .go"

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
html_title = project + f" {version}"
html_logo = "_static/RF_logo.png"
html_favicon = "_static/RF_logo.png"
html_show_sourcelink = False
html_static_path = [
    "_static",
]
html_css_files = [
    "custom.css",
]

# bug in pydata theme: remove in future version starting at 0.15.4
html_sidebars = {
    "support_us": [],
    "capabilities": [],
    "about_us": [],
    "release_notes/*": [],
}
pygments_style = "none"

# announcement = """
# Are you new to ReservoirFlow? Please find some time to read the first release note
# (see <a href='/release_notes/release_note_v0.1.0.html'>Release Note v0.1.0</a>).
# <br>
# In addition, please consider introducing yourself to our community
# (see <a href='/community/forum/introduce_yourself.html'>Introduce Yourself 👋</a>).
# """
html_theme_options = {
    # "header_links_before_dropdown": 7,
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/hiesabx/reservoirflow",
            "icon": "fab fa-github",
            "type": "fontawesome",
        },
        {
            "name": "LinkedIn",
            "url": "https://www.linkedin.com/company/hiesab/",
            "icon": "fab fa-linkedin-in",
        },
    ],
    "footer_start": ["copyright"],
    "footer_center": "",
    "footer_end": "",
    "navbar_start": [
        "navbar-logo",
    ],
    "navbar_center": [
        "navbar-nav",
    ],
    "navbar_end": [
        "version-switcher",
        "theme-switcher",
        "navbar-icon-links",
    ],
    "switcher": {
        "json_url": json_url,
        "version_match": switcher_version,
    },
    "check_switcher": True,
    "show_version_warning_banner": False,  # for multiple versions
    "navigation_with_keys": False,
    "show_toc_level": 2,
    "secondary_sidebar_items": ["page-toc"],
    # "announcement": announcement,
    "pygments_light_style": "default",
    "pygments_dark_style": "monokai",  # 'native'
    "back_to_top_button": True,
    # "primary_sidebar_end": ["sidebar-ads.html"],
}

# comments_config = {
#     # "utterances": {
#     #     "repo": "hiesabx/reservoirflow_comments",
#     # },
#     "giscus": {
#         "repo": "hiesabx/reservoirflow_comments",
#     }
# }


# Units and Factors:
def store_dict(in_dict, name="FACTORS", folder=""):
    name = name.upper()
    if name == "UNITS":
        label = "property"
    elif name == "FACTORS":
        label = "factor"
    elif name == "NOMENCLATURE":
        label = "property"
    else:
        raise ValueError("name is unknown")
    columns = list(in_dict.keys())
    headers = [label] + columns
    # headers = [h.capitalize() for h in headers]
    props = in_dict[columns[0]]
    report = []

    for prop in props:
        row = []
        if name == "UNITS":
            for column in columns:
                row.append(":math:`" + in_dict[column][prop] + "`")
        elif name == "FACTORS":
            for column in columns:
                row.append(in_dict[column][prop])
        elif name == "NOMENCLATURE":
            row.append(in_dict[columns[0]][prop])
            row.append(":math:`" + in_dict[columns[1]][prop] + "`")

        report.append(
            [
                prop,
                *row,
            ]
        )
    table = tabulate(report, headers=headers, showindex=False, tablefmt="rst")
    path = f"{folder+name.lower()}_table.rst"
    with open(path, "w") as f:
        f.writelines(table)


store_dict(rf.UNITS, "UNITS", "user_guide/units_factors/")
store_dict(rf.FACTORS, "FACTORS", "user_guide/units_factors/")
store_dict(rf.NOMENCLATURE, "NOMENCLATURE", "user_guide/nomenclature/")
