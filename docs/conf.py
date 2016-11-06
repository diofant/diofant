#
# Diofant documentation build configuration file.
#
# This file is execfile()d with the current directory set to its containing dir.
#
# The contents of this file are pickled, so don't put values in the namespace
# that aren't pickleable (module imports are okay, they're removed automatically).
#
# All configuration values have a default value; values that are commented out
# serve to show the default value.

import os
import warnings
import diofant

# This turns numpydoc's section warnings to exceptions,
# see numpy/numpydoc#58.
warnings.simplefilter('error', UserWarning)

# General configuration
# ---------------------

# Add any Sphinx extension module names here, as strings. They can be extensions
# coming with Sphinx (named 'sphinx.addons.*') or your custom ones.
extensions = ['sphinx.ext.autodoc', 'sphinx.ext.viewcode', 'sphinx.ext.mathjax',
              'numpydoc', 'sphinx.ext.graphviz', 'sphinx.ext.autosummary',
              'sphinx.ext.intersphinx']

numpydoc_show_class_members = False
numpydoc_class_members_toctree = False

# If true, Sphinx will warn about all references where the target cannot be found.
nitpicky = True

# MathJax file, which is free to use.  See http://www.mathjax.org/docs/2.0/start.html
mathjax_path = 'http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS_HTML-full'

# Add any paths that contain templates here, relative to this directory.
templates_path = ['.templates']

# A list of glob-style patterns that should be excluded when looking for source files.
exclude_patterns = ['README.rst']

# The suffix of source filenames.
source_suffix = '.rst'

# The master toctree document.
master_doc = 'index'

# General substitutions.
project = 'Diofant'
copyright = '2016 SymPy Development Team, Sergey B Kirpichev'

# The default replacements for |version| and |release|, also used in various
# other places throughout the built documents.
#
# The short X.Y version.
version = diofant.__version__
# The full version, including alpha/beta/rc tags.
release = version

# Today_fmt is used as the format for a strftime call.
today_fmt = '%B %d, %Y'

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# The name of a reST role (builtin or Sphinx extension) to use as the
# default role, that is, for text marked up `like this`.
default_role = 'math'

# Options for HTML output
# -----------------------

# If not '', a 'Last updated on:' timestamp is inserted at every page bottom,
# using the given strftime format.
html_last_updated_fmt = '%b %d, %Y'

# If true, generate domain-specific indices in addition to the general
# index. For e.g. the Python domain, this is the global module index. Default
# is True.  This value can be a bool or a list of index names that
# should be generated.
html_domain_indices = ['py-modindex']

# Output file base name for HTML help builder.
htmlhelp_basename = 'Diofantdoc'

# Options for LaTeX output
# ------------------------

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title, author, document class [howto/manual], toctree_only).
# toctree_only is set to True so that the start file document itself is not included in the
# output, only the documents referenced by it via TOC trees.  The extra stuff in the master
# document is intended to show up in the HTML, but doesn't really belong in the LaTeX output.
latex_documents = [('index', 'diofant.tex', 'Diofant Documentation',
                    'Diofant Development Team', 'manual', True)]

# Additional stuff for the LaTeX preamble.
latex_elements = {
    'babel':     r'\usepackage[english]{babel}',
    'fontenc': r'''
\usepackage{bm}
\usepackage{amssymb}
''',
    'fontpkg':   '',
    'inputenc':  '',
    'utf8extra': '',
    'preamble':  r'''
% redefine \LaTeX to be usable in math mode
\expandafter\def\expandafter\LaTeX\expandafter{\expandafter\text\expandafter{\LaTeX}}
'''
}

# Show page numbers next to internal references
latex_show_pagerefs = True

# We use False otherwise the module index gets generated twice.
latex_use_modindex = False

# Options for extensions
# ----------------------

texinfo_documents = [
    (master_doc, 'diofant', 'Diofant Documentation', 'Diofant Development Team',
   'Diofant', 'Computer algebra system (CAS) in Python', 'Programming', 1),
]

# Use svg for graphviz
graphviz_output_format = 'svg'

intersphinx_mapping = {
    'python3': ('https://docs.python.org/3/', None),
    'numpy': ('https://docs.scipy.org/doc/numpy', None),
}

html_static_path = ['_static']
html_context = {
    'css_files': [
        '_static/theme_overrides.css',
    ],
}

# RTD theme support:
on_rtd = os.environ.get('READTHEDOCS', None) == 'True'
if not on_rtd:  # only import and set the theme if we're building docs locally
    import sphinx_rtd_theme
    html_theme = 'sphinx_rtd_theme'
    html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]
else:
    html_context = {
        'css_files': [
            'https://media.readthedocs.org/css/sphinx_rtd_theme.css',
            'https://media.readthedocs.org/css/readthedocs-doc-embed.css',
            '_static/theme_overrides.css',
        ],
    }
