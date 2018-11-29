#
# Diofant documentation build configuration file.
#
# This file is execfile()d with the current directory set to its
# containing dir.
#
# The contents of this file are pickled, so don't put values in the
# namespace that aren't pickleable (module imports are okay, they're
# removed automatically).
#

import warnings

import diofant


# Turns numpydoc's section warnings to exceptions, see numpy/numpydoc#58.
warnings.simplefilter('error', UserWarning)

# Add any Sphinx extension module names here, as strings.
extensions = ['sphinx.ext.autodoc', 'sphinx.ext.viewcode', 'sphinx.ext.mathjax',
              'sphinx.ext.graphviz', 'sphinx.ext.intersphinx',
              'sphinx.ext.extlinks', 'sphinx.ext.napoleon']

napoleon_google_docstring = False
napoleon_use_param = False
napoleon_use_rtype = False

# Sphinx will warn about all references where the target cannot be found.
nitpicky = True

# The document name of the "master" document, that is, the document
# that contains the root toctree directive.
master_doc = 'index'

# Project information.
project = 'Diofant'
copyright = '2006-2018 SymPy Development Team, 2013-2018 Sergey B Kirpichev'
version = diofant.__version__
release = version

# The name of default reST role, that is, for text marked up `like this`.
default_role = 'math'

# The theme to use for HTML and HTML Help pages.
html_theme = 'sphinx_rtd_theme'

# The LaTeX engine to build the docs.
latex_engine = 'xelatex'

# This value determines how to group the document tree into LaTeX source
# files. It must be a list of tuples (startdocname, targetname, title,
# author, documentclass, toctree_only),
latex_documents = [('index', 'diofant.tex', 'Diofant Documentation',
                    'Diofant Development Team', 'manual', True)]

# A dictionary that contains LaTeX snippets that override predefined.
latex_elements = {
    'preamble':  r'''
\setmainfont{DejaVu Serif}
\setsansfont{DejaVu Sans}
\setmonofont{DejaVu Sans Mono}
% redefine \LaTeX to be usable in math mode
\expandafter\def\expandafter\LaTeX\expandafter{\expandafter\text\expandafter{\LaTeX}}
'''
}

# Add page references after internal references.
latex_show_pagerefs = True

# The output format for Graphviz when building HTML files.
graphviz_output_format = 'svg'

# Contains mapping the locations and names of other projects that
# should be linked to in this documentation.
intersphinx_mapping = {
    'python3': ('https://docs.python.org/3/', None),
    'numpy': ('https://docs.scipy.org/doc/numpy', None),
    'scipy': ('https://docs.scipy.org/doc/scipy/reference', None),
}

# Dictionary of external sites, mapping unique short alias names to a
# base URL and a prefix.
extlinks = {
    'issue': ('https://github.com/diofant/diofant/issues/%s', '#'),
    'pull': ('https://github.com/diofant/diofant/pull/%s', '#'),
    'commit': ('https://github.com/diofant/diofant/commit/%s', ''),
    'sympyissue': ('https://github.com/sympy/sympy/issues/%s', 'sympy/sympy#'),
    'sympypull': ('https://github.com/sympy/sympy/pull/%s', 'sympy/sympy#'),
}

# The number of times the linkcheck builder will attempt to check a URL
# before declaring it broken.
linkcheck_retries = 3

# A list of regular expressions that match URIs that should not be checked.
linkcheck_ignore = [r'https://primes.utm.edu/notes/gaps.html',
                    r'https://primes.utm.edu/glossary/xpage/BertrandsPostulate.html',
                    r'https://primes.utm.edu/prove/prove2_3.html',
                    r'https://primes.utm.edu/glossary/xpage/Pseudoprime.html']

# This value controls if docstring for classes or methods, if not explicitly
# set, is inherited form parents.
autodoc_inherit_docstrings = False

# A list of paths that contain custom static files. Relative paths are taken as
# relative to the configuration directory. They are copied to the output’s
# _static directory.
html_static_path = ['_static']


# https://docs.readthedocs.io/en/latest/guides/adding-custom-css.html
def setup(app):
    app.add_stylesheet('custom.css')
