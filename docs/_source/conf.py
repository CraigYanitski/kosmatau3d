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
sys.path.insert(0, os.path.abspath('..'))


# -- Project information -----------------------------------------------------

project = 'kosmatau3d'
copyright = '2023, Craig Yanitski'
author = 'Craig Yanitski'

# The full version, including alpha/beta/rc tags
release = '1.0.3'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    # 'sphinx.ext.duration',
    # 'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',  # Link to other project's documentation (see mapping below)
    'sphinx.ext.viewcode',  # Add a link to the Python source code for classes, functions etc.
    'sphinx_autodoc_typehints', # Automatically document param types (less noise in class signature)
    'nbsphinx',  # Integrate Jupyter Notebooks and Sphinx
    'IPython.sphinxext.ipython_console_highlighting',
]

# Mappings for sphinx.ext.intersphinx. Projects have to have Sphinx-generated doc! (.inv file)
intersphinx_mapping = {
    "python": ("https://docs.python.org/3/", None),
}

autosummary_generate = True   # Turn on sphinx.ext.autosummary
autoclass_content = 'both'    # Add __init__ doc (params) to class summaries
html_show_sourcelink = False  # Remove 'view source code' from top of page (for html, not python)
autodoc_inherit_docstrings = True   # Inherit docstring from bas class if no docstring
set_type_checking_flag = True   # Enable 'extensive' imports for sphinx_autodoc_typehints
nbsphinx_allow_errors = True  # Continue through Jupyter errors
add_module_names = False   # Reemove namespaces from class/method signatures

# Add any paths that contain templates here, relative to this directory.
templates_path = ['../_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
# exclude_patterns = ['kosmatau3d.models.cyplot']

def skip_member(app, what, name, obj, skip, options):
    # List of external module names to skip
    external_modules_to_skip = ['scipy', 'scipy.interpolate']

    if what == 'module' and obj.__name__ in external_modules_to_skip:
        return True
    return None

def setup(app):
    app.connect('autodoc-skip-member', skip_member)


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
# Default: alabaster
html_theme = 'pydata_sphinx_theme'
html_css_files = ["readthedocs-custom.css"] # Override some CSS settings

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['../_static']


# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
# exclude_patterns = ["_build", "_templates", "Thumbs.db", "jupyter_execute", ".*"]

# Additional stuff
master_doc = "index"

# Cross-reference existing Python objects
intersphinx_mapping = {
    "python": ("https://docs.python.org/3/", None),
    "pandas": ("https://pandas.pydata.org/pandas-docs/stable", None),
    "numpy": ("https://numpy.org/doc/stable", None),
    "scipy": ("https://docs.scipy.org/doc/scipy", None),
    "numba": ("https://numba.pydata.org/numba-doc/latest", None),
}