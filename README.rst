SymPy
=====

A Python library for symbolic mathematics.

See the AUTHORS file for the list of authors.

Documentation and usage
-----------------------

The most recent documentation is available in HTML format online:

http://omg.readthedocs.org/

The online documentation is automatically generated for commits
in the master branch of the repository, build status is available here:

.. image:: https://readthedocs.org/projects/omg/badge/?version=latest
    :target: https://readthedocs.org/projects/omg/?badge=latest
    :alt: Documentation Status

You can generate everything at the above site in your local copy of SymPy by::

    $ cd doc
    $ make html

Then the docs will be in `_build/html`. If you don't want to read that, here
is a short usage:

From this directory, start python and::

    >>> from sympy import Symbol, cos
    >>> x = Symbol('x')
    >>> e = 1/cos(x)
    >>> print e.series(x, 0, 10)
    1 + (1/2)*x**2 + (5/24)*x**4 + (61/720)*x**6 + (277/8064)*x**8 + O(x**10)

SymPy also comes with a console that is a simple wrapper around the
classic python console (or IPython when available) that loads the
sympy namespace and executes some common commands for you.

To start it, issue::

    $ bin/isympy

from this directory if SymPy is not installed or simply::

    $ isympy

if SymPy is installed.

Installation
------------

SymPy has a hard dependency on the `mpmath <http://mpmath.org/>`
library (version >= 0.19).  You should install it first, please refer to
the mpmath installation guide:

https://github.com/fredrik-johansson/mpmath#1-download--installation

To install SymPy itself, then simply run::

    $ python setup.py install

If you install it system-wide, you may need to prefix the previous command with ``sudo``::

    $ sudo python setup.py install

To install the SymPy for Python 3, simply run the above commands with a Python
3 interpreter.

Tests
-----

To execute all tests, run::

    $./setup.py test

in the current directory.

For more fine-grained running of tests or doctest, use ``bin/test`` or
respectively ``bin/doctest``. The master branch is automatically tested by
Travis CI, the results can be seen here:

.. image:: https://secure.travis-ci.org/skirpichev/omg.png?branch=master
    :target: http://travis-ci.org/skirpichev/omg

Bugs
----

Our issue tracker is at https://github.com/skirpichev/omg/issues.  Please report
any bugs that you find.  Or, even better, fork the repository on GitHub and
create a pull request.
