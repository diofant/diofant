How to Build Documentation
==========================

To make the html documentation, install the prerequisites, e.g. on
Debian/Ubuntu (similarly for other distributions)::

    apt-get install python-sphinx texlive-latex-recommended dvipng
    apt-get install python-pip
    pip install numpydoc==0.5 sphinx_rtd_theme

and do::

    cd docs
    make html

and to view it, do::

    epiphany _build/html/index.html
