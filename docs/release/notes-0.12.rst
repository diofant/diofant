============
Diofant 0.12
============

Not Released Yet

New features
============

* Support modular exponentiation of :class:`~diofant.polys.rings.PolyElement`'s, see :pull:`1032`.

Major changes
=============

Compatibility breaks
====================

* Removed ``vring()`` and ``vfield()`` functions, see :pull:`1016`.
* Drop support for ``from_list()`` initialization for multivariate polynomials, see :pull:`1035`.
* Drop ``to_dense()``, ``tail_degrees()`` and ``almosteq`` methods and ``is_monic``, ``is_primitive`` attributes of :class:`~diofant.polys.rings.PolyElement`, see :pull:`1035` and :pull:`1036`.
* Drop ``is_monic``, ``is_primitive``, ``zero``, ``one`` and ``unit`` attributes of :class:`~diofant.polys.polytools.Poly`, see :pull:`1036` and :pull:`1039`.
* Drop ``sring()`` and ``poly_from_expr()`` functions, see :pull:`1037`.
* Functions and classes of the :mod:`~diofant.polys.polytools` module do not support anymore iterables as polynomial generator, see :pull:`1039`.
* Drop unused function ``dispersion()``, see :pull:`1051`.

Minor changes
=============

* Module :mod:`~diofant.polys.sqfreetools` was ported to use sparse polynomial representation, see :pull:`1009`.
* Module :mod:`~diofant.polys.factortools` was ported to use sparse polynomial representation, see :pull:`1015`, :pull:`1018`, :pull:`1019`, :pull:`1020` and :pull:`1021`.
* Special case univariate polynomials with :class:`~diofant.polys.univar.UnivarPolynomialRing` and :class:`diofant.polys.univar.UnivarPolyElement`, see :pull:`1024`.
* Module :mod:`~diofant.polys.rootisolation` was ported to use sparse polynomial representation, see :pull:`1030`, :pull:`1031` and :pull:`1035`.
* Implement :attr:`~diofant.domains.finitefield.ModularInteger.is_primitive`, see :pull:`1035`.

Developer changes
=================

* Depend on `flake8-sfs <https://github.com/peterjc/flake8-sfs>`_, see :pull:`983`.
* Depend on `mypy <http://mypy-lang.org/>`_, see :pull:`1022`.

Issues closed
=============

See the `release milestone <https://github.com/diofant/diofant/milestone/6?closed=1>`_
for complete list of issues and pull requests involved in this release.

These Sympy issues also were addressed:

* :sympyissue:`19630` ``rsolve`` gives None for linear homogeneous recurrence relation
* :sympyissue:`19076` modular exponentiation of poly
* :sympyissue:`19670` Poly(E**100000000) is slow to create
* :sympyissue:`19755` poly gives coercion error when integers and rationals are mixed
* :sympyissue:`19760` minimal_polynomial using Groebner basis can give wrong result
* :sympyissue:`19770` Limit involving cosine
* :sympyissue:`19766` Incorrect limit
* :sympyissue:`19774` evalf() doesn't evaluate terms in an exponential
