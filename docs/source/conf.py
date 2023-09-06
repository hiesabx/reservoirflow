from datetime import datetime

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
    # "sphinx_mdinclude",
]

templates_path = ["_templates"]
exclude_patterns = []

# myst:
myst_all_links_external = True

# html:
html_theme = "pydata_sphinx_theme"
html_title = project  # + release
html_static_path = ["_static"]

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
