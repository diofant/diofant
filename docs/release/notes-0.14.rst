============
Diofant 0.14
============

Not Released Yet

New features
============

* Support calculating limits with :class:`~diofant.functions.elementary.piecewise.Piecewise` functions, see :pull:`1214`.

Major changes
=============

* Use recursive (former ``poly()`` method, without using :func:`~diofant.core.function.expand`) algorithm of creating polynomials from expressions, see :pull:`1047`.

Compatibility breaks
====================

* Removed support for CPython 3.9, see :pull:`1192`.
* Removed ``to_mpi()`` method of :class:`~diofant.sets.sets.Interval`, see :pull:`1194`.
* Removed ``poly()`` function, use :meth:`~diofant.core.expr.Expr.as_poly` method to create a :class:`~diofant.polys.polytools.Poly` instance from :class:`~diofant.core.expr.Expr`, see :pull:`1047`.
* Removed functions ``bool_map()``, ``POSform()`` and ``SOPform()``, see :commit:`04ea41a220` and :commit:`be319badf5`.

Minor changes
=============

* Support unevaluated :class:`~diofant.polys.rootoftools.RootOf`'s over finite fields, see :pull:`1209`.
* Provide default clause (condition :class:`~diofant.logic.boolalg.BooleanTrue`) for :class:`~diofant.functions.elementary.piecewise.Piecewise`, see :pull:`1215`.

Developer changes
=================

Issues closed
=============

See the `release milestone <https://github.com/diofant/diofant/milestone/8?closed=1>`_
for complete list of issues and pull requests involved in this release.

These Sympy issues also were addressed:

* :sympyissue:`22487`: [integrals] Wrong result for Integral((cos(x**2)-cos(x))/x**2, (x, -oo, oo))
* :sympyissue:`22493`: Series expansion introduces new variables
* :sympyissue:`22558`: Error in ODE-Solver-Documentation
* :sympyissue:`22837`: Solve simplest algebraic equations with dummy parameter
* :sympyissue:`22836`: Series: Possible improvements for Order of expressions involving factorials
* :sympyissue:`22788`: RecursionError for unevluated expression in latex
* :sympyissue:`22863`: Hangs: integrate((3*x**3-x**2+2*x-4)/sqrt(x**2-3*x+2), (x, 0, 1))
* :sympyissue:`22862`: Problem with separable differential equation
* :sympyissue:`22893`: 'limit' in combination with 'positive=True' gives wrong result
* :sympyissue:`22878`: RecursionError in trigsimp
* :sympyissue:`22982`: limit((log(E + 1/x) - 1)**(1 - sqrt(E + 1/x)), x, oo) returns 0 instead of oo
* :sympyissue:`22986`: limit(acosh(1 + 1/x)*sqrt(x), x, oo) is evaluated incorrectly.
* :sympyissue:`14433`: x not in QQ.frac_field(1/x)
* :sympyissue:`23069`: integrate(r**4*sqrt(1 - r**2), (r, 0, 1)) gives incorrect result
* :sympyissue:`19639`: TypeError in integrate
* :sympyissue:`23086`: Incorrect result of simplify
* :sympyissue:`23156`: sympy.Sum() bug when summing up reciprocal of gamma
* :sympyissue:`23174`: Problem with gf_edf_zassenhaus()
* :sympyissue:`21409`: Printing of polynomial over FF
* :sympyissue:`22673`: Roots of a polynomial over a finite field computed regardless of specified polynomial domain
* :sympyissue:`12531`: cancel does not return expanded form
* :sympyissue:`6322`: degree((x+1)**10000) takes too long
* :sympyissue:`22583`: is_polynomial right for wrong reasons (and sometimes wrong)
* :sympyissue:`23202`: Dropping "all" __ne__ methods?
* :sympyissue:`23223`: Wrong integration results of trigonometric functions
* :sympyissue:`23224`: Python code printer not respecting tuple with one element
* :sympyissue:`23231`: Sympy giving the wrong solution
* :sympyissue:`14387`: Tutorial on limits creates impression that they are two-sided by default
* :sympyissue:`8166`: Limit assumes function is continuous?
* :sympyissue:`14502`: Problem with limit including factorial.
* :sympyissue:`18492`: Limit of Piecewise function - NotImplementedError: Don't know how to calculate the mrv
* :sympyissue:`23266`: Regression(?) in 1.10 for limits
* :sympyissue:`7391`: Limits for expressions with undetermined functions give wrong results
* :sympyissue:`23287`: Regression in is_integer for Mul of Pow
* :sympyissue:`11496`: Wrong result in limit calculation of limit(erfc(ln(1/x)),x,oo)?
* :sympyissue:`3663`: series expansion of acosh and acoth
