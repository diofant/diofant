.. _installation:

Installation
------------

The SymPy CAS can be installed on virtually any computer with Python
3.2 or above.  SymPy requires `setuptools`_.  The current recommended
method of installation is directly from the source files.

Source
======

If you are a developer or like to get the latest updates as they come, be
sure to install from git::

    $ git clone git://github.com/skirpichev/omg.git
    $ cd omg
    $ python setup.py develop

.. note::

    You could use `pyvenv`_ (or `virtualenv`_) to create isolated Python
    environment first, instead of installaing everything system-wide.

To update to the latest version, go into your repository and execute::

    $ git pull origin master

You can see old SymPy's history (from Hg and SVN repos) in the
branch sympy-svn-history.  To see this history as part of
master's, simply do::

    $ git fetch origin 'refs/replace/*:refs/replace/*'

Run SymPy
=========

After installation, it is best to verify that your freshly-installed SymPy
works.  To do this, start up the Python interpreter and import the
SymPy libraries::

    >>> from sympy import *

From here, execute some simple SymPy statements like the ones below::

    >>> x = Symbol('x')
    >>> limit(sin(x)/x, x, 0)
    1
    >>> integrate(1/x, x)
    log(x)

For a starter guide on using SymPy effectively, refer to the :ref:`tutorial`.

Questions
=========

If you think there's a bug or you would like to request a feature, please
:ref:`open an issue ticket <reporting-issues>`.

.. _setuptools: https://packaging.python.org/en/latest/projects.html#setuptools
.. _pyvenv: https://docs.python.org/3/library/venv.html
.. _virtualenv: https://virtualenv.pypa.io/
