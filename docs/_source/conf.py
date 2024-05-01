# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys

from datetime import date

sys.path.insert(0, os.path.abspath("../../kosmatau3d/"))


# -- Project information -----------------------------------------------------

project = "kosmatau3d"
copyright = f"2022-{date.today().year}, Craig Yanitski"
author = "Craig Yanitski"

# The full version, including alpha/beta/rc tags
release = "1.0.10"


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    # "sphinxcontrib.embedpdf",  # enable rendering of pdf images in documentation
    # 'sphinx.ext.duration',
    # 'sphinx.ext.doctest',
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.coverage",
    "sphinx.ext.mathjax",
    "sphinx.ext.intersphinx",  # Link to other project's documentation (see mapping below)
    "sphinx.ext.viewcode",  # Add a link to the Python source code for classes, functions etc.
    "sphinx.ext.todo",  # Allow ToDo admonition
    # 'sphinx_autodoc_typehints',   # Automatically document param types (less noise in class signature)
    "nbsphinx",  # Integrate Jupyter Notebooks and Sphinx
    "IPython.sphinxext.ipython_console_highlighting",
]

# Mappings for sphinx.ext.intersphinx. Projects have to have Sphinx-generated doc! (.inv file)
# intersphinx_mapping = {
#     "python": ("https://docs.python.org/3/", None),
# }

autosummary_generate = True  # Turn on sphinx.ext.autosummary
# autosummary_imported_members = True # Used for function template
autoclass_content = "both"  # Add __init__ doc (params) to class summaries
html_show_sourcelink = (
    True  # Remove 'view source code' from top of page (for html, not python)
)
autodoc_inherit_docstrings = True  # Inherit docstring from base class if no docstring
set_type_checking_flag = True  # Enable 'extensive' imports for sphinx_autodoc_typehints
nbsphinx_allow_errors = True  # Continue through Jupyter errors
add_module_names = False  # Remove namespaces from class/method signatures

# Add any paths that contain templates here, relative to this directory.
templates_path = ["../_templates"]
# autosummary_template_path = '../_templates/'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
# exclude_patterns = ['kosmatau3d/models/cyplot*']


def skip_member(app, what, name, obj, skip, options):
    # List of external module names to skip
    external_modules_to_skip = [
        "scipy",
        "scipy.interpolate",
        "kosmatau3d.models.interpolations.interpolate.interpolate",
    ]

    if what == "module" and obj.__name__ in external_modules_to_skip:
        return True
    return None


def setup(app):
    app.connect("autodoc-skip-member", skip_member)


# -- Options for HTML output -------------------------------------------------

html_theme = "pydata_sphinx_theme"
html_static_path = ["../_static"]
html_title = "kosmatau3d"
html_last_updated_fmt = "%d %b %Y"
html_theme_options = {
    "switcher": {
        "json_url": "https://kosmatau3d.readthedocs.io/en/latest/_static/version_switcher.json",
        "version_match": release,
    },
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/CraigYanitski/kosmatau3d/",
            "icon": "fa-brands fa-github",  # icons from `fontawesome`
        },
        {
            "name": "PyPI",
            "url": "https://pypi.org/project/kosmatau3d/",
            "icon": "fa-custom fa-pypi",
        },
    ],
    "navbar_start": ["navbar-logo", "version-switcher"],
    # "navbar_persistent": [],
    "show_version_warning_banner": True,
    # "secondary_sidebar_items": ["page-toc"],
    "show_prev_next": False,
}
html_sidebars = {
    "installation": [],
    "getting-started": [],
    "theory": [],
    "development": [],
    "acknowledgements": [],
}
html_css_files = ["pydata-custom.css"]  # Override some CSS settings


# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
# exclude_patterns = ["_build", "_templates", "Thumbs.db", "jupyter_execute", ".*"]

# Additional stuff
master_doc = "index"

# Cross-reference existing Python objects
# intersphinx_mapping = {
#     "python": ("https://docs.python.org/3/", None),
#     "pandas": ("https://pandas.pydata.org/pandas-docs/stable", None),
#     "numpy": ("https://numpy.org/doc/stable", None),
#     "scipy": ("https://docs.scipy.org/doc/scipy", None),
#     "numba": ("https://numba.pydata.org/numba-doc/latest", None),
# }
