.. _installation:

Installation
------------

The Diofant can be installed on any computer with Python 3.5
or above.  You can install latest release with pip::

    $ pip install diofant

or to install also extra dependencies::

    $ pip install diofant[gmpy,plot]

.. tip::

    Use `venv`_ to create isolated Python environment first,
    instead of installing everything system-wide.

.. _installation-src:

From Sources
============

If you are a developer or like to get the latest updates as they come,
be sure to install from git::

    $ git clone git://github.com/diofant/diofant.git
    $ cd diofant
    $ pip install -e .[develop,docs]

.. note::

    Diofant requires `setuptools`_ for installation from sources.

Run Diofant
===========

To verify that your freshly-installed Diofant works, please start up
the Python interpreter and execute some simple statements like the
ones below::

    >>> from diofant.abc import x
    >>> ((1 + x)**(1/x)).limit(x, 0)
    E

.. tip::

    Use `IPython`_ for interactive work.  Please refer to the
    documentation of module :mod:`diofant.interactive` for details
    of available configuration settings.

For a starter guide on using Diofant, refer to the :ref:`tutorial`.

Feedback
========

If you think there's a bug, you have a question or you would like to
request a feature, please :ref:`open an issue ticket
<reporting-issues>`.

.. _setuptools: https://setuptools.readthedocs.io/en/latest/
.. _venv: https://docs.python.org/3/library/venv.html
.. _IPython: http://ipython.readthedocs.io/en/stable/
