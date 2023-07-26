# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import sys
import os
from os.path import join
#sys.path.append('/usr/lib/python3.11/site-packages/breathe/')

sys.path.insert(0, join(os.path.abspath('.'), '..', 'build', 'wrapper',
                        'python_sphinx'))

project = 'pyfaust'
copyright = '2015-2023'
author = 'Rémi Gribonval, Luc Le Magoarou,  Adrien Leman (2016), Nicolas Bellot(2015-2016), Thomas Gautrais (2015), Hakim Hadj-Djilani (2018-), Pascal Carrivain (2023-).'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

#extensions = []

#extensions = [ 'breathe' ]

extensions = [
    'myst_parser',
    'sphinx.ext.autodoc',
    'sphinx.ext.doctest',
    'sphinx.ext.mathjax',
    'sphinx_autodoc_typehints',
    'sphinx.ext.autosummary',
    'nbsphinx',
    'sphinx_math_dollar'
]
#breathe_projects = {"pyfaust": "../build/doc/xml"}
#
#breathe_default_project = "pyfaust"

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

#html_theme = 'alabaster'
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_show_sourcelink = False

toc_object_entries_show_parents = 'all'

