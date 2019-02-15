============
Diofant 0.11
============

Not Released Yet

New features
============

Major changes
=============

Compatibility breaks
====================

* Removed support for Python 3.5 and 3.6, see :pull:`775`.
* ``is_monomial`` attribute of :class:`~diofant.polys.polytools.Poly` renamed to :attr:`~diofant.polys.polytools.Poly.is_term`, see :pull:`780`.
* Removed ``log()`` helper from :class:`~diofant.domains.RationalField`, see :pull:`787`.

Minor changes
=============

* Support truncation for elements of :class:`~diofant.domains.RealAlgebraicField` to :class:`int`, see :pull:`788`.

Developer changes
=================

* Depend on `sphinxcontrib-bibtex <https://sphinxcontrib-bibtex.readthedocs.io/en/latest/>`_ to track the bibliography, see :pull:`766`.

Issues closed
=============

See the `release milestone <https://github.com/diofant/diofant/milestone/4?closed=1>`_
for complete list of issues and pull requests involved in this release.

These Sympy issues also were addressed:

* :sympyissue:`15943` Wrong result from summation
* :sympyissue:`12163` matematica code printer does not handle floats and derivatives correctly
* :sympyissue:`11642` Geometric sum doesn't evaluate with float base
* :sympyissue:`15984` Value error in limit
