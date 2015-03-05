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

Brief History
-------------

SymPy was started by Ondřej Čertík in 2005, he wrote some code during the
summer, then he wrote some more code during the summer 2006. In February 2007,
Fabian Pedregosa joined the project and helped fixed many things, contributed
documentation and made it alive again. 5 students (Mateusz Paprocki, Brian
Jorgensen, Jason Gedge, Robert Schwarz and Chris Wu) improved SymPy incredibly
during the summer 2007 as part of the Google Summer of Code (GSoC). Pearu Peterson
joined the development during the summer 2007 and he has made SymPy much more
competitive by rewriting the core from scratch, that has made it from 10x to
100x faster. Jurjen N.E. Bos has contributed pretty printing and other patches.
Fredrik Johansson has wrote mpmath and contributed a lot of patches.

SymPy has participated in every GSoC since 2007.  Moderate amount
of SymPy's development has come from GSoC students.

In 2011, Ondřej Čertík stepped down as lead developer, with Aaron Meurer, who
also started as a GSoC student, taking his place.

Ondřej Čertík is still active in the community, but is too busy with work
and family to play a lead development role.  Unfortunately, his remaining
activity neither constructive nor productive anymore and SymPy just
slowly dying now.

This project is a fork of the SymPy, last SymPy's commit is cbdd072 (22
Feb 2015).  The git history goes back to 2007, when development was in svn and
then in hg.   You can see this old history in the branch sympy-svn-history.
Also, you can set up the master branch to show this history too.  Simply do::

     git fetch origin 'refs/replace/*:refs/replace/*'

You can use git to see the biggest developers::

     git shortlog -ns

This will show the top developers from the last year::

     git shortlog -ns --since="1 year"
