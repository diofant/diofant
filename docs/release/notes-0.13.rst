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
* Removed ``prem()`` function, see :pull:`1140`.

Minor changes
=============

* Protect hashed :class:`~diofant.polys.rings.PolyElement`'s from modifications, see :pull:`1033`.
* Add gaussian rationals as an exact domain, associated with :class:`~diofant.domains.ComplexField`, see :pull:`1138`.

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
* :sympyissue:`21180`: Bug: sympy.factor doesn't work for Poly !!!
* :sympyissue:`21167`: Empty list of solutions returned for equation with cubic roots
* :sympyissue:`21029`: Continuous limits involving division by x
* :sympyissue:`20697`: Series is not simplified to final answer in output in sympy 1.7.1
* :sympyissue:`20578`: A strange behavior of limit function
* :sympyissue:`20444`: Leading Term with log
* :sympyissue:`19453`: Limit changes from simplification of original expression
* :sympyissue:`19442`: Non-existent bi-directional limit gives ValueError
* :sympyissue:`11667`: limit(1/x, x, 0) == oo ??
* :sympyissue:`21202`: laplace_transform(cosh(2*x), x, s) raises RecursionError
* :sympyissue:`21227`: Nested logarithms add unnecessary order term to series expansions
* :sympyissue:`21263`: Solutions of cubic equation
* :sympyissue:`21334`: RecursionError while calculating leading term
* :sympyissue:`21342`: 1/(exp(it) - 2) integrates wrong
* :sympyissue:`21319`: Primitive part of zero polynomial
* :sympyissue:`21341`: Issues with continued fraction for real roots of cubic polynomials
* :sympyissue:`21024`: sympy.polys.polyerrors.CoercionFailed integration regressions?
