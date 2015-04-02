.. _installation:

Installation
------------

The SymPy CAS can be installed on virtually any computer with Python
2.7 or above. SymPy use `setuptools`_ and does require `mpmath`_ Python library
to be installed first (`setuptools`_ should handle this dependency
automatically).  The current recommended method of installation
is directly from the source files.

SymPy officially supports Python 2.7, 3.4, and PyPy3.

Source
======

If you are a developer or like to get the latest updates as they come, be sure
to install from git. To download the repository, execute the following from the
command line::

    $ git clone git://github.com/skirpichev/omg.git

From your favorite command line terminal, change directory into that folder and
execute the following::

    $ python setup.py install

Alternatively, if you don't want to install the package onto your computer, you
may run SymPy with the "isympy" console (which automatically imports SymPy
packages and defines common symbols) by executing within that folder::

    $ ./bin/isympy

You may now run SymPy statements directly within the Python shell::

    >>> from __future__ import division
    >>> from sympy import *
    >>> x, y, z, t = symbols('x y z t')
    >>> k, m, n = symbols('k m n', integer=True)
    >>> f, g, h = symbols('f g h', cls=Function)
    >>> diff(x**2/2, x)
    x

To update to the latest version, go into your repository and execute::

    $ git pull origin master

You can see old SymPy's history (from Hg and SVN repos) in the
branch sympy-svn-history.  To see this history as part of master's, simply do::

    $ git fetch origin 'refs/replace/*:refs/replace/*'

Run SymPy
=========

After installation, it is best to verify that your freshly-installed SymPy
works. To do this, start up Python and import the SymPy libraries::

    $ python
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

If you think there's a bug or you would like to request a feature, please open
an `issue ticket`_.

.. _issue ticket: https://github.com/skirpichev/omg/issues
.. _setuptools: https://packaging.python.org/en/latest/projects.html#setuptools
.. _mpmath: http://mpmath.org/
