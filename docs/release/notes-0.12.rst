============
Diofant 0.12
============

Not Released Yet

New features
============

Major changes
=============

Compatibility breaks
====================

* Removed ``vring()`` and ``vfield()`` functions, see :pull:`1016`.

Minor changes
=============

* Module :mod:`~diofant.polys.sqfreetools` was ported to use sparse polynomial representation, see :pull:`1009`.
* Module :mod:`~diofant.polys.factortools` was ported to use sparse polynomial representation, see :pull:`1015`, :pull:`1018`, :pull:`1019`, :pull:`1020` and :pull:`1021`.
* Special case univariate polynomials with :class:`~diofant.polys.univar.UnivarPolynomialRing` and :class:`diofant.polys.univar.UnivarPolyElement`, see :pull:`1024`.

Developer changes
=================

* Depend on `flake8-sfs <https://github.com/peterjc/flake8-sfs>`_, see :pull:`983`.

Issues closed
=============

See the `release milestone <https://github.com/diofant/diofant/milestone/6?closed=1>`_
for complete list of issues and pull requests involved in this release.

These Sympy issues also were addressed:
