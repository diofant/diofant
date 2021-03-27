============
Diofant 0.13
============

Not Released Yet

New features
============

* Support square-free factorization of multivariate polynomials over finite fields (with adaptation of Musser's algorithm), see :pull:`1132`.

Major changes
=============

Compatibility breaks
====================

* Removed ``n()`` method from :class:`~diofant.core.evalf.EvalfMixin`, see :pull:`1114`.
* Former submodule ``diofant.polys.polyconfig`` now is :mod:`diofant.config`, see :pull:`1115`.
* Drop support for ``DIOFANT_DEBUG`` environment variable, see :pull:`1115`.
* Renamed ``Ring`` as :class:`~diofant.domains.ring.CommutativeRing`, see :pull:`1123`.
* Removed support for Python 3.7 and 3.8, see :pull:`1118` and :pull:`1124`.
* ``FiniteRing`` renamed to :class:`~diofant.domains.IntegerModRing`, see :pull:`1124`.
* Removed ``igcd()``, ``ilcm()`` and ``prod()`` functions, see :pull:`1125`.
* Changed the :class:`~diofant.core.function.Derivative` (and similary :func:`~diofant.core.function.diff`) syntax to ``Derivative(foo, (x, 2))`` from ``Derivative(foo, x, 2)``, see :pull:`1131`.

Minor changes
=============

* Protect hashed :class:`~diofant.polys.rings.PolyElement`'s from modifications, see :pull:`1033`.

Developer changes
=================

* Turn on type checking for the whole codebase, see :pull:`1114`.
* Don't include regression tests in the coverage statistics, see :pull:`1060`.

Issues closed
=============

See the `release milestone <https://github.com/diofant/diofant/milestone/7?closed=1>`_
for complete list of issues and pull requests involved in this release.

These Sympy issues also were addressed:

* :sympyissue:`20861`: reduce_inequalities() gives impossible answer
* :sympyissue:`20874`: Port the PRS algorithm to the sparse polynomial implementation
* :sympyissue:`20902`: Incorrect inequality solving: False returned instead of answer
* :sympyissue:`20941`: Fails to Solve Definite Integral
* :sympyissue:`20973`: cancel raises PolynomialError for exp(1+O(x))
* :sympyissue:`20985`: TypeErrors appearing for simple plynomial manipulations (did not happen in v1.6.1)
* :sympyissue:`21031`: Limit of "limit (((1+x)**(1/x)-(1+2*x)**(1/(2*x)))/asin (x),x,0)" is wrong with v1.7.1
* :sympyissue:`21034`: (Integration) regressions?
* :sympyissue:`21038`: Incorrect computation of a basic limit, regression from 1.6.2 to 1.7.1
* :sympyissue:`21041`: integrate error
* :sympyissue:`21063`: Wrong value of improper integral when using unevaluated -oo as boundary
* :sympyissue:`21075`: Order term being added to exact expansion
* :sympyissue:`21091`: Invalid comparison of non-real when using integrate()
* :sympyissue:`19590`: Poly.diff() doesn't support higher order derivatives
* :sympyissue:`21121`: Same symbols created in different processes are not resolved as being equal
* :sympyissue:`21107`: S.Infinity.is_nonzero returns False
* :sympyissue:`21132`: Integral with parametres: wrong and too long result
