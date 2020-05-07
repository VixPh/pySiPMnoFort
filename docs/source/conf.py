# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.


import sphinx_rtd_theme
import os
import sys
sys.path.insert(0, os.path.abspath('../../'))
sys.path.insert(0, os.path.abspath('../../files'))

# -- Project information -----------------------------------------------------

project = 'pySiPM'
copyright = '2020, Edoardo Proserpio'
author = 'Edoardo Proserpio \\and Massimiliano Antonello \\and Romualdo Santoro'

# The full version, including alpha/beta/rc tags
release = '0.1'




# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx_rtd_theme','sphinx.ext.autodoc', 'sphinx.ext.doctest',
              'sphinx.ext.intersphinx','sphinx.ext.viewcode',
              'sphinx.ext.autosummary','sphinx.ext.napoleon','sphinx.ext.imgmath']

# Settings
numfig = True
html_theme = "sphinx_rtd_theme"
html_logo = "logo.png"
html_theme_options = {
    'sticky_navigation':True
}
imgmath_image_format = 'svg'
imgmath_font_size = 16
imgmath_latex_preamble = r'''
\usepackage{amsmath}
\usepackage{amssymb}
'''

# Latex settings
latex_documents = [
  ('index', 'doc.tex', 'pySiPM Documentation',
   r'Edoardo Proserpio\and Massimiliano Antonello\and Romualdo Santoro', 'manual'),
]
latex_elements = {'papersize':'letterpaper','pointsize':'10pt','preamble': r'''
\usepackage{amsmath}
\usepackage{amssymb}
\DeclareRobustCommand{\and}{\end{tabular}\kern-\tabcolsep\\\begin{tabular}[t]{c}}
'''}
latex_logo = 'logo.png'

# Napoleon settings
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = False
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.

# html_theme = "alabaster"
# html_theme_options = {
#     'github_user': 'EdoPro98',
#     'fixed_sidebar': True,
#     'sidebar_collapse': True
# }

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']
