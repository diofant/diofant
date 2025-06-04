============
Diofant 0.16
============

Not Released Yet

New features
============

Major changes
=============

Compatibility breaks
====================

* Reorder tuples, returned by :meth:`~diofant.polys.polytools.Poly.gcdex` and :meth:`~diofant.polys.polytools.Poly.half_gcdex`, see :pull:`1454`.
* Removed ``Expr.lseries()`` (``n=None`` option for :meth:`~diofant.core.expr.Expr.series`), see :pull:`1418`.

Minor changes
=============

Developer changes
=================

Issues closed
=============

See the `release milestone <https://github.com/diofant/diofant/milestone/10?closed=1>`_
for complete list of issues and pull requests involved in this release.

These Sympy issues also were addressed:

* :sympyissue:`26990`: Getting an AtrributeError when finding the limit of a complex equation with operators like sinh() and tanh()
* :sympyissue:`27001`: Possible bug with sympy.solve
* :sympyissue:`27048`: Order term does not simplify with terms containing log
* :sympyissue:`27050`: Wrong result of definite integral
* :sympyissue:`26707`: Wrong translation of lerchphi into Mathematica
* :sympyissue:`27074`: Sum ignores an undetermined value
* :sympyissue:`27108`: Wrong integration involving Dirac Delta
* :sympyissue:`27195`: Pow of Add should not extract Floats
* :sympyissue:`27234`: Integral of Abs(cos(x + y)) typeerror
* :sympyissue:`27236`: bug with limit() function
* :sympyissue:`27238`: perfect_power(-64) fails to find (-4)**3
* :sympyissue:`27256`: Geometric series with free symbols
* :sympyissue:`27298`: Wrong result when integrating Legendre Polynomial (missing case distinction)
* :sympyissue:`27300`: Wrong result for an integral over complex exponential with a Diracdelta function
* :sympyissue:`27551`: Invalid limit
* :sympyissue:`27624`: sympy takes so long to solve a solvable system of polynomial equations symbolically
* :sympyissue:`27675`: A simple example of a wrong definite integral
* :sympyissue:`27683`: RecursionError: maximum recursion depth exceeded for SympifyError(a) when exchange an Integer into a Float
* :sympyissue:`27712`: Infinite hang and resource consumption when solving a specially crafted equation due to call to Integer
* :sympyissue:`27794`: PolyRing.index doesn't appear to follow Python's list indexing convention
* :sympyissue:`27786`: A sum of positive elements returns 0
* :sympyissue:`27798`: Bug of domain.unify
* :sympyissue:`27819`: PolyRing: Issue in index Method of PolyRing When Using String as Generator
* :sympyissue:`27874`: What should the extended Euclidean algorithm return when all inputs are zero?
* :sympyissue:`27901`: rsolve raises AttributeError for some of the hypergeometric univariate functions
* :sympyissue:`14120`: GeneratorsError for primitive_element([Poly(x**2 - 2)], x)
* :sympyissue:`28006`: Mul(0, x, evaluate=False).is_zero gives False
* :sympyissue:`28033`: Incorrect result when calculating a definite integral as the limit of a Rieman sum
* :sympyissue:`28089`: LP solver may loop or give wrong answer on infeasible LPs
