Diofant's documentation
=======================

.. image:: https://mybinder.org/badge_logo.svg
    :target: https://mybinder.org/v2/gh/diofant/diofant/master

`Diofant <https://diofant.readthedocs.io/en/latest/>`_ is a Python
library for symbolic mathematics.  If you are new to Diofant, start
with the :ref:`Tutorial <tutorial>`.

This is the central page for all of Diofant's documentation.

.. note::

   Documentation examples assume (unless otherwise clearly stated)
   that these statements were executed in the beginning of the
   interactive session:

      >>> from diofant import *
      >>> a, b, c, d, x, y, z, t = symbols('a b c d x y z t')
      >>> k, m, n = symbols('k m n', integer=True)
      >>> f, g, h = symbols('f g h', cls=Function)

.. toctree::
   :maxdepth: 2

   install
   tutorial/index
   modules/index
   internals/index
   guide
   aboutus
   release/index
   sources
