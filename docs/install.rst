.. _installation:

Installation
------------

The Diofant CAS can be installed on virtually any computer with Python
3.4 or above.  Diofant requires `setuptools`_.  You can install latest
release with pip::

    $ pip install --pre diofant

.. note::

    You could use `pyvenv`_ to create isolated Python environment first,
    instead of installing everything system-wide.

Source
======

If you are a developer or like to get the latest updates as they come, be
sure to install from git::

    $ git clone git://github.com/diofant/diofant.git
    $ cd diofant
    $ pip install -e .

or to install also extra dependencies::

    $ pip install -e .[gmpy,plot]

To update to the latest version, go into your repository and execute::

    $ git pull origin master

You can see old Diofant's history (from Hg and SVN repos) in the
branch sympy-svn-history.  To see this history as part of
master's, simply do::

    $ git fetch origin 'refs/replace/*:refs/replace/*'

Run Diofant
===========

After installation, it is best to verify that your freshly-installed Diofant
works.  To do this, start up the Python interpreter and import the
Diofant libraries::

    >>> from diofant import *

From here, execute some simple Diofant statements like the ones below::

    >>> x = Symbol('x')
    >>> limit(sin(x)/x, x, 0)
    1
    >>> integrate(1/x, x)
    log(x)

For a starter guide on using Diofant effectively, refer to the :ref:`tutorial`.

Questions
=========

If you think there's a bug or you would like to request a feature, please
:ref:`open an issue ticket <reporting-issues>`.

.. _setuptools: https://setuptools.readthedocs.io/en/latest/
.. _pyvenv: https://docs.python.org/3/library/venv.html
