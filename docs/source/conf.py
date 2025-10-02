import sys
import os
sys.path.insert(0, os.path.abspath('../..'))

import furo  # ensure theme is available during builds

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'Pseudomancer'
copyright = '2025, Adam Dinan'
author = 'Adam Dinan'
release = '0.1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx_autodoc_typehints",
    "myst_parser",
]
try:  # optional extension; skip if missing
    import sphinx_autodoc_typehints  # noqa: F401
except ImportError:  # pragma: no cover - docs build fallback
    extensions.remove("sphinx_autodoc_typehints")


templates_path = ['_templates']
exclude_patterns = []

autodoc_mock_imports = [
    "Bio",
    "numpy",
    "pandas",
    "requests",
]


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "furo"
html_static_path = ['_static']

# -- Warnings ------------------------------------------------------------
suppress_warnings = ["toc.not_included"]
