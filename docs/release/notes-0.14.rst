============
Diofant 0.14
============

Not Released Yet

New features
============

* Support calculating limits with :class:`~diofant.functions.elementary.piecewise.Piecewise` functions and boolean expressions, see :pull:`1214` and :pull:`1218`.
* Support directional limits on the complex plane, see :pull:`1232`.

Major changes
=============

* Use recursive (former ``poly()`` method, without using :func:`~diofant.core.function.expand`) algorithm of creating polynomials from expressions, see :pull:`1047`.

Compatibility breaks
====================

* Removed support for CPython 3.9, see :pull:`1192`.
* Removed ``to_mpi()`` method of :class:`~diofant.sets.sets.Interval`, see :pull:`1194`.
* Removed ``poly()`` function, use :meth:`~diofant.core.expr.Expr.as_poly` method to create a :class:`~diofant.polys.polytools.Poly` instance from :class:`~diofant.core.expr.Expr`, see :pull:`1047`.
* Removed functions ``bool_map()``, ``POSform()`` and ``SOPform()``, see :commit:`04ea41a220` and :commit:`be319badf5`.
* Changed semantics of the ``dir`` kwarg for the :class:`~diofant.calculus.limits.Limit`, now '+' is -1, '-' is 1 and 'real' is :class:`~diofant.sets.fancysets.Reals`, see :pull:`1234` and :pull:`1235`.
* Removed ``diofant.calculus.euler`` and ``diofant.calculus.finite_diff`` modules, see :pull:`1271`.
* Removed ``diofant.vector`` module, see :pull:`1274`.
* Removed ``diofant.diffgeom`` module, see :pull:`1281`.
* Removed ``diofant.stats`` module, see :pull:`1276`.
* Removed ``diofant.geometry`` module and ``line_integrate`` function, see :pull:`1283`.
* Removed ``diofant.plotting`` module, see :pull:`1284`.
* Removed unused ``prefixes``, ``postfixes``, ``capture`` and ``variations`` functions, see :pull:`1282` and :pull:`1290`.

Minor changes
=============

* Support unevaluated :class:`~diofant.polys.rootoftools.RootOf`'s over finite fields, see :pull:`1209`.
* Provide default clause (condition :class:`~diofant.logic.boolalg.BooleanTrue`) for :class:`~diofant.functions.elementary.piecewise.Piecewise`, see :pull:`1215`.

Developer changes
=================

* Drop dependency on the `flake8-rst <https://github.com/flake8-docs/flake8-rst>`_ and depend on the `flake518 <https://github.com/carstencodes/flake518>`_ instead, see :pull:`1268`.

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
* :sympyissue:`23299`: Sympy is unable to integrate this
* :sympyissue:`23319`: testing limit of n*tan(pi/n) results in incorrect answer in 1.7rc1+
* :sympyissue:`5539`: Equal Integrals compare different when using different variables
* :sympyissue:`23425`: PolynomialError when I try to call classify_ode
* :sympyissue:`23432`: Series expansion around float fails with NotImplementedError
* :sympyissue:`8433`: limit involving error function returns bad result
* :sympyissue:`13750`: erf has wrong limit in -oo
* :sympyissue:`23497`: binomial(-1, -1) returns 0, should return 1
* :sympyissue:`23562`: In new version of sympy, dsolve does not give a solution when another derivative is involved
* :sympyissue:`23585`: FiniteSet documentation inconsistent with usage in sympy
* :sympyissue:`23596`: Integral of real function has complex result
* :sympyissue:`23605`: Inefficiency in the Integrator with a Rational Expression
* :sympyissue:`23637`: Missing solutions from polynomial system (various solvers)
* :sympyissue:`23479`: Sparse poly gcd fails with HeuristicGCDFailed('no luck')
* :sympyissue:`22605`: Incorrect result from minpoly(cos(pi/9))
* :sympyissue:`23677`: minimal_polynomial fails for very complicated algebraic number
* :sympyissue:`23836`: Incorrect results for limits of Piecewise at discontinuity
* :sympyissue:`23845`: Gruntz should have been free of _w, value error, recursion error
* :sympyissue:`23855`: linsolve gives odd result if symbols are duplicated
* :sympyissue:`24067`: incorrect limit in simple parametric rational polynomial
* :sympyissue:`24127`: Error on all limits with Piecewise
* :sympyissue:`23702`: Cannot specify ODE initial conditions as just f(0)
* :sympyissue:`23707`: AttributeError in integral
* :sympyissue:`24210`: Error on limits regarding terms like (1+u)^v.
* :sympyissue:`24225`: Multivariable limit should be undefined, but gives unity.
* :sympyissue:`24266`: Changed behaviour of series() involving exp, I
* :sympyissue:`24331`: Limit of log(z) as z goes to 0 with z complex returns '-oo' instead of 'zoo'
* :sympyissue:`23766`: Factor hangs on exponential functions with base e
* :sympyissue:`24360`: Remove usage of numpy.distutils in autowrap module
* :sympyissue:`24346`: factor with extension=True fails for rational expression
* :sympyissue:`20913`: Poly(x + 9671406556917067856609794, x).real_roots() is slow
* :sympyissue:`24386`: sympy.limit yields wrong limit in sigmoidal expression
* :sympyissue:`24390`: Incorrectly evaluated expression
* :sympyissue:`24461`: sympy.polys.polyerrors.HeuristicGCDFailed: no luck -- when multiplying two Polys
* :sympyissue:`24543`: Rational calc value error
* :sympyissue:`6326`: PolynomialRing should not derive from CharacteristicZero
