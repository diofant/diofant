.. _installation:

Installation
------------

The Diofant CAS can be installed on virtually any computer with Python
3.4 or above.  Diofant requires `setuptools`_.  You can install latest
release with pip::

    $ pip install --pre diofant

or to install also extra dependencies::

    $ pip install --pre diofant[gmpy,plot]

.. tip::

    You could use `pyvenv`_ to create isolated Python environment first,
    instead of installing everything system-wide.

.. _installation-src:

From sources
============

If you are a developer or like to get the latest updates as they come, be
sure to install from git::

    $ git clone git://github.com/diofant/diofant.git
    $ cd diofant
    $ pip install -e .[develop,docs]

Run Diofant
===========

To verify that your freshly-installed Diofant works, please start up the
Python interpreter and execute some simple statements like the ones below::

    >>> from diofant.abc import x
    >>> ((1 + x)**(1/x)).limit(x, 0)
    E

However, we recommend using `IPython`_ for working interactively.  Use input
transformations from the module :mod:`diofant.interactive`, that could reduce
boilerplate while interacting with Diofant due to the Python language syntax.

For a starter guide on using Diofant, refer to the :ref:`tutorial`.

Questions
=========

If you think there's a bug or you would like to request a feature, please
:ref:`open an issue ticket <reporting-issues>`.

.. _setuptools: https://setuptools.readthedocs.io/en/latest/
.. _pyvenv: https://docs.python.org/3/library/venv.html
.. _IPython: http://ipython.readthedocs.io/en/stable/
