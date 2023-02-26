.. _installation:

Installation
------------

The Diofant can be installed on any computer with Python 3.10
or above.  You can install latest release with pip::

    pip install diofant

or to install also extra dependencies::

    pip install diofant[gmpy,interactive]

To use :ref:`Unicode pretty printing <d-pretty-printer>` --- configure your
system to have good TTF fonts.  The `DejaVu Sans Mono
<https://dejavu-fonts.github.io/>`_ seems to be an acceptable choice.  On
Debian you can install this `font package
<https://packages.debian.org/sid/fonts-dejavu>`_ with::

    apt install fonts-dejavu

.. _installation-src:

From Sources
============

If you are a developer or like to get the latest updates as they come, be sure
to install from the git repository and include required extra dependencies::

    git clone git://github.com/diofant/diofant.git
    cd diofant
    pip install -e .[develop,docs,tests]

Run Diofant
===========

To verify that your freshly-installed Diofant works, please start up
the Python interpreter::

    python

and execute some simple statements like the ones below::

    >>> from diofant.abc import x
    >>> ((1 + x)**(1/x)).limit(x, 0)
    E

.. tip::

    :ref:`Run Diofant as a module <cli>` for interactive work::

        python -m diofant

For a starter guide on using Diofant, refer to the :ref:`tutorial`.

Also, you may want to run full set of unit tests to make
sure everything works::

    pytest --pyargs diofant

`pytest`_ and some other packages are required for testing, so be sure to
install the Diofant first with the optional "tests" list of dependencies::

    pip install diofant[tests]

Feedback
========

If you think there's a bug, you have a question or you would like to
request a feature, please :ref:`open an issue ticket
<reporting-issues>`.  General questions and comments can be `sent
<mailto:diofant@googlegroups.com>`_ to the `Diofant mailing list`_.

.. _pytest: https://docs.pytest.org/en/latest/
.. _Diofant mailing list: https://groups.google.com/forum/#!forum/diofant
