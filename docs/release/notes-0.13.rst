============
Diofant 0.13
============

7 Nov 2021

New features
============

* Support square-free factorization of multivariate polynomials over finite fields (with adaptation of Musser's algorithm), see :pull:`1132`.

Major changes
=============

* Support calling from the command-line as ``python -m diofant``, see :pull:`853`.  Thanks to André Roberge.

Compatibility breaks
====================

* Removed ``n()`` method from :class:`~diofant.core.evalf.EvalfMixin`, see :pull:`1114`.
* Former submodule ``diofant.polys.polyconfig`` now is :mod:`diofant.config`, see :pull:`1115`.
* Drop support for ``DIOFANT_DEBUG`` environment variable, see :pull:`1115`.
* Drop support for CPython 3.7 and 3.8, see :pull:`1118` and :commit:`5cae972`.
* Renamed ``Ring`` as :class:`~diofant.domains.ring.CommutativeRing`, see :pull:`1123`.
* Removed support for Python 3.7 and 3.8, see :pull:`1118` and :pull:`1124`.
* ``FiniteRing`` renamed to :class:`~diofant.domains.IntegerModRing`, see :pull:`1124`.
* Removed ``igcd()``, ``ilcm()`` and ``prod()`` functions, see :pull:`1125`.
* Changed the :class:`~diofant.core.function.Derivative` (and similary :func:`~diofant.core.function.diff`) syntax to ``Derivative(foo, (x, 2))`` from ``Derivative(foo, x, 2)``, see :pull:`1131`.
* Removed ``prem()`` function, see :pull:`1140`.
* Removed ``lseries()`` method of :class:`~diofant.core.expr.Expr`, use :meth:`~diofant.core.expr.Expr.series` with ``n=None``, see :pull:`1146`.

Minor changes
=============

* Protect hashed :class:`~diofant.polys.rings.PolyElement`'s from modifications, see :pull:`1033`.
* Add gaussian rationals as an exact domain, associated with :class:`~diofant.domains.ComplexField`, see :pull:`1138`.
* Support :class:`~diofant.functions.elementary.trigonometric.tan` in :func:`~diofant.polys.numberfields.minimal_polynomial`, see :pull:`1159`.
* 100% test coverage for ``plotting`` module, see :pull:`1175`.
* Support CPython 3.10, see :pull:`1162`.

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
* :sympyissue:`21396`: Pow.as_base_exp inconsistent with I.as_base_exp
* :sympyissue:`21410`: Polynomial power raises KeyError
* :sympyissue:`21437`: log(Abs)
* :sympyissue:`21460`: Polynomial GCD result is different for dense trivial polynomial
* :sympyissue:`21466`: Regression for match for differential binomial expression
* :sympyissue:`21166`: Wrong integration result involving square root of absolute value
* :sympyissue:`21486`: expand_func(besselj(oo, x)) -> RecursionError
* :sympyissue:`21530`: Incorrect limit
* :sympyissue:`21549`: Bug: integrate(x*sqrt(abs(x)),(x,-1,0)) returns wrong result
* :sympyissue:`21557`: Summation of geometric series with non-real exponent does not evaluate
* :sympyissue:`21550`: Bug: limit returns wrong result for rational function
* :sympyissue:`21177`: Incorrect residue for cot(pi*x)/(x**2 - 3*x + 3)
* :sympyissue:`21245`: laurent series Fibonacci generating fuction
* :sympyissue:`11833`: error in limit involving exp, sinh and an assumption (maybe related to caching)
* :sympyissue:`9127`: ntheory.AskEvenHandler.Mul is order-dependent
* :sympyissue:`21606`: Notimplemented in simple limit
* :sympyissue:`21641`: Simplify hangs
* :sympyissue:`21651`: doit() method *sometimes* ignores floor and ceiling within Sum
* :sympyissue:`20461`: Eq(Product(4*n**2/(4*n**2 - 1), (n, 1, oo)), pi/2) incorrectly gives False
* :sympyissue:`13029`: with gens, time taken for sqf increases orders of magnitude faster than factor as input size increases
* :sympyissue:`21711`: odd result for integrate(sqrt(1 - (x-1)*(x-1)), (x, 0, 1))
* :sympyissue:`21721`: Bug in integration solver
* :sympyissue:`21716`: isympy -c python tab triggered auto completion not working
* :sympyissue:`21741`: integrate() does not work with multivariable function that is solved by simple substitution. DomainError: there is no ring associated with CC
* :sympyissue:`21756`: Incorrect limit with ratio of complex exponentials
* :sympyissue:`21760`: Poly div is slow
* :sympyissue:`21761`: sympy.polys.polyerrors.NotAlgebraic Exception
* :sympyissue:`21430`: minpoly raises 'NotAlgebraic' for tan(13*pi/45)
* :sympyissue:`21766`: solve breaks on certain repeated inputs
* :sympyissue:`21773`: TypeError multiplying Subs expressions
* :sympyissue:`21785`: Limit gives TypeError from as_leading_term
* :sympyissue:`21812`: LambertW displaying in jupyter lab
* :sympyissue:`21814`: Printing of unevaluated Mul needs brackets
* :sympyissue:`21176`: Incorrect residue of x**2*cot(pi*x)/(x**4 + 1)
* :sympyissue:`21852`: simple quadratic not solving
* :sympyissue:`21859`: AttributeError: 'mpz' object has no attribute 'denominator' with sp.series()
* :sympyissue:`21882`: Incorrect solutions given by solve
* :sympyissue:`21890`: RecursionError and TypeError in nonlinsolve
* :sympyissue:`21888`: TypeError raised for evalf containing summations
* :sympyissue:`5822`: What should summation() do with non-integer limits?
* :sympyissue:`19745`: Weird value of a sum
* :sympyissue:`9358`: summation: Wrong out for non-integral range
* :sympyissue:`21905`: raise NotImplementedError("Equation not in exact domain. Try converting to rational") Error
* :sympyissue:`21938`: Series raises an error at infinity for an example which can be solved by aseries
* :sympyissue:`21984`: ValueError: list.remove(x): x not in list occurs in nonlinsolve
* :sympyissue:`21999`: detection of infinite solution request
* :sympyissue:`22020`: Comparing two operations that contain log sometimes leads to TypeError exception
* :sympyissue:`22051`: Nonlinsolve incorrect result
* :sympyissue:`22058`: Regression in solveset for quadratic with symbolic coefficients
* :sympyissue:`22073`: Interval with oo
* :sympyissue:`22093`: sympy.polys.polyerrors.HeuristicGCDFailed: no luck
* :sympyissue:`22155`: Problem with solving simple separable ODE
* :sympyissue:`22220`: Bug in the evaluation of a log limit
* :sympyissue:`22248`: solve running forever
* :sympyissue:`22294`: Bernoulli differential equation
* :sympyissue:`22322`: 'abs' is not parsed correctly
* :sympyissue:`22334`: Wrong answer returned while calculating limit for different arrangements of the same expression
* :sympyissue:`22400`: Minpoly doesn't terminate
* :sympyissue:`22435`: sympy integration error
