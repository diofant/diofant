.. _installation:

Installation
------------

.. warning::

    For the development version, please do first::

        pip install git+https://github.com/fredrik-johansson/mpmath.git@1bb25a4


The Diofant can be installed on any computer with Python 3.9
or above.  You can install latest release with pip::

    pip install diofant

or to install also extra dependencies::

    pip install diofant[gmpy,plot]

.. tip::

    Use :mod:`venv` to create isolated Python environment first,
    instead of installing everything system-wide.

To use :ref:`Unicode pretty printing <d-pretty-printer>` --- configure your
system to have good TTF fonts.  The `DejaVu Sans Mono
<https://dejavu-fonts.github.io/>`_ seems to be an acceptable choice.  On
Debian you can install this `font package
<https://packages.debian.org/sid/fonts-dejavu>`_ with::

    apt install fonts-dejavu

.. _installation-src:

From Sources
============

If you are a developer or like to get the latest updates as they come,
be sure to install from git::

    git clone git://github.com/diofant/diofant.git
    cd diofant
    pip install -e .[develop,docs]

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
    documentation of module :mod:`~diofant.interactive` for details
    of available configuration settings.

For a starter guide on using Diofant, refer to the :ref:`tutorial`.

Also, you may want to run full set of unit tests to make
sure everything works::

    pytest --pyargs diofant

`pytest`_ and some other packages are required for testing, so be sure to
install the Diofant first with extra dependecies::

    pip install diofant[tests]

Feedback
========

If you think there's a bug, you have a question or you would like to
request a feature, please :ref:`open an issue ticket
<reporting-issues>`.  General questions and comments can be `sent
<mailto:diofant@googlegroups.com>`_ to the `Diofant mailing list`_.

.. _setuptools: https://setuptools.readthedocs.io/en/latest/
.. _IPython: https://ipython.readthedocs.io/en/stable/
.. _pytest: https://docs.pytest.org/en/latest/
.. _Diofant mailing list: https://groups.google.com/forum/#!forum/diofant
