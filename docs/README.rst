How to Build Documentation
==========================

To make the html documentation, install the prerequisites (optionally,
you can also install the sphinx_rtd_theme package with pip)::

    pip install -r docs/requirements.txt

then enter the project directory and do::

    python setup.py build_sphinx

and to view it, do::

    epiphany build/sphinx/html/index.html
