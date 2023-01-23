"""
Diofant documentation build configuration file.

This file is execfile()d with the current directory set to its
containing dir.

The contents of this file are pickled, so don't put values in the
namespace that aren't pickleable (module imports are okay, they're
removed automatically).
"""

import inspect
import os
import sys
import warnings

import sphinx.util.texescape

import diofant


# Turns numpydoc's section warnings to exceptions, see numpy/numpydoc#58.
warnings.simplefilter('error', UserWarning)

# Add any Sphinx extension module names here, as strings.
extensions = ['sphinx.ext.autodoc', 'sphinx.ext.linkcode', 'sphinx.ext.mathjax',
              'sphinx.ext.graphviz', 'sphinx.ext.intersphinx',
              'sphinx.ext.extlinks', 'sphinx.ext.napoleon',
              'sphinxcontrib.bibtex', 'sphinxcontrib.autoprogram']

napoleon_google_docstring = False
napoleon_use_param = False
napoleon_use_rtype = False

# List our *.bib files.
bibtex_bibfiles = ['sources.bib']

# Sphinx will warn about all references where the target cannot be found.
nitpicky = True

# Glob-style patterns that should be excluded when looking for sources.
exclude_patterns = ['README.rst']

# Project information.
project = 'Diofant'
copyright = '2006-2022 SymPy Development Team, 2013-2023 Sergey B Kirpichev'
version = diofant.__version__
release = version

# The name of default reST role, that is, for text marked up `like this`.
default_role = 'math'

# The theme to use for HTML and HTML Help pages.
html_theme = 'sphinx_rtd_theme'

# The LaTeX engine to build the docs.
latex_engine = 'xelatex'

# If True, the PDF build from the LaTeX files created by Sphinx will use xindy
# rather than makeindex.
latex_use_xindy = False

# This value determines how to group the document tree into LaTeX source
# files.  It must be a list of tuples (startdocname, targetname, title,
# author, documentclass, toctree_only),
latex_documents = [('index', 'Diofant-' + version + '.tex', 'Diofant Documentation',
                    'Diofant Development Team', 'manual', True)]

# A dictionary that contains LaTeX snippets that override predefined.
latex_elements = {
    'fontpkg': r"""
\setmainfont{DejaVu Serif}
\setsansfont{DejaVu Sans}
\setmonofont{DejaVu Sans Mono}
""",
    'preamble': r"""
% redefine \LaTeX to be usable in math mode
\expandafter\def\expandafter\LaTeX\expandafter{\expandafter\text\expandafter{\LaTeX}}

\fvset{formatcom=\baselineskip10pt\relax\let\strut\empty}
"""
}

# Add page references after internal references.
latex_show_pagerefs = True

# The output format for Graphviz when building HTML files.
graphviz_output_format = 'svg'

# Contains mapping the locations and names of other projects that
# should be linked to in this documentation.
intersphinx_mapping = {'python': ('https://docs.python.org/3/', None),
                       'numpy': ('https://numpy.org/doc/stable/', None),
                       'scipy': ('https://docs.scipy.org/doc/scipy/', None)}

# Dictionary of external sites, mapping unique short alias names to a
# base URL and a prefix.
extlinks = {'issue': ('https://github.com/diofant/diofant/issues/%s', '#%s'),
            'pull': ('https://github.com/diofant/diofant/pull/%s', '#%s'),
            'commit': ('https://github.com/diofant/diofant/commit/%s', '%s'),
            'sympyissue': ('https://github.com/sympy/sympy/issues/%s', 'sympy/sympy#%s'),
            'sympypull': ('https://github.com/sympy/sympy/pull/%s', 'sympy/sympy#%s')}

# The number of times the linkcheck builder will attempt to check a URL
# before declaring it broken.
linkcheck_retries = 3

# A list of regular expressions that match URIs that should not be checked.
linkcheck_ignore = [r'https://primes.utm.edu/notes/gaps.html',
                    r'https://primes.utm.edu/glossary/xpage/BertrandsPostulate.html',
                    r'https://primes.utm.edu/prove/prove2_3.html',
                    r'https://primes.utm.edu/glossary/xpage/Pseudoprime.html',
                    r'https://github.com/.*']

# This value controls if docstring for classes or methods, if not explicitly
# set, is inherited form parents.
autodoc_inherit_docstrings = False

# Create table of contents entries for domain objects (e.g. functions,
# classes, attributes, etc).
toc_object_entries = False

# A list of paths that contain custom static files.  Relative paths are taken as
# relative to the configuration directory.  They are copied to the output’s
# _static directory.
html_static_path = ['_static']

# These paths are either relative to html_static_path
# or fully qualified paths (eg. https://...)
html_css_files = ['custom.css']

# Should we show "Created using Sphinx" in the HTML footer?
html_show_sphinx = False

# Paths to the logo and favicon.ico, relative to the conf.py's directory.
html_logo = '_static/logo.svg'
html_favicon = '_static/favicon.ico'
latex_logo = '_static/logo.png'

# Theme-specific options.
html_theme_options = {'logo_only': True,
                      'display_version': False,
                      'analytics_id': 'UA-147167689-1'}

# The inline configuration options for mathjax.  The value is used as
# a parameter of MathJax.Hub.Config().
mathjax3_config = {'CommonHTML': {'linebreaks': {'automatic': True}},
                   'HTML-CSS': {'linebreaks': {'automatic': True}},
                   'SVG': {'linebreaks': {'automatic': True}}}

# The maximum number of days to cache remote inventories.
intersphinx_cache_limit = 30


def linkcode_resolve(domain, info):
    """Determine the URL corresponding to Python object."""
    if domain != 'py':
        return

    modname = info['module']
    fullname = info['fullname']

    submod = sys.modules.get(modname)
    if submod is None:
        return

    obj = submod
    for part in fullname.split('.'):
        try:
            obj = getattr(obj, part)
        except AttributeError:
            return

    # strip decorators, which would resolve to the source of the decorator
    # possibly an upstream bug in getsourcefile, bpo-1764286
    try:
        unwrap = inspect.unwrap
    except AttributeError:
        pass
    else:
        obj = unwrap(obj)

    try:
        fn = inspect.getsourcefile(obj)
    except TypeError:
        fn = None
    if not fn:
        return

    try:
        source, lineno = inspect.getsourcelines(obj)
    except OSError:
        lineno = None

    if lineno:
        linespec = f'#L{lineno:d}-L{lineno + len(source) - 1:d}'
    else:
        linespec = ''

    fn = os.path.relpath(fn, start=os.path.dirname(diofant.__file__))

    blobpath = 'https://github.com/diofant/diofant/blob/'
    if 'dev' in version:
        return blobpath + f'master/diofant/{fn}{linespec}'
    else:
        return blobpath + f'v{version}/diofant/{fn}{linespec}'


# monkey-patch sphinx
del sphinx.util.texescape.tex_replacements[29:]
