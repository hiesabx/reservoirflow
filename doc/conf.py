# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import sys

import reservoirflow as rf

# sys.path.insert(0, os.path.abspath(".."))
sys.path.insert(0, os.path.abspath(os.path.join("..")))
# os.path.abspath(os.path.join(".."))

project = "ReservoirFlow"
copyright = "2023, Zakariya Abugrin"
author = "Zakariya Abugrin"
release = "v" + rf.__version__

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx_copybutton",
    # "sphinx_automodapi.automodapi",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", ".denv"]
master_doc = "index"

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "pydata_sphinx_theme"  # "sphinx_rtd_theme", "'pydata_sphinx_theme'"
html_title = project  # + f" ({release})"
html_static_path = ["_static"]
myst_enable_extension = []
myst_all_links_external = True

# automodeapi:
# https://sphinx-automodapi.readthedocs.io/en/latest/automodapi.html
# automodsumm_inherited_members = False
# numpydoc_show_class_members = True

# https://www.sphinx-doc.org/en/master/usage/extensions/autodoc.html#confval-autodoc_inherit_docstrings

autoclass_content = "init"  # "class", "init", "both"
autodoc_class_signature = "mixed"  # "mixed", "separated"
autodoc_member_order = "bysource"  # 'alphabetical', 'groupwise', 'bysource'
autodoc_default_options = {
    "members": True,
    "member-order": "bysource",
    "undoc-members": True,
    "private-members": False,
    # "special-members": "__init__",
    # "inherited-members": False,
    # "show-inheritance": False,
    # "ignore-module-all": False,
    # "imported-members": False,
    # # "exclude-members": ,
    # "class-doc-from": "__init__",
    # "no-value": False,
}

autodoc_docstring_signature = True
autodoc_typehints = "signature"  # 'signature', 'description', 'none', 'both'
autodoc_typehints_description_target = "all"  # "all", "documented", "documented_params"
autodoc_typehints_format = "short"  # 'fully-qualified', 'short'
autodoc_inherit_docstrings = False
autodoc_preserve_defaults = False
