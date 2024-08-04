============
Diofant 0.15
============

Not Released Yet

New features
============

* New configuration option (``MAX_INTEGER_NBITS``) to control the maximal size of evaluated integers, see :pull:`1327`.
* Added :func:`~diofant.polys.polytools.eliminate` to eliminate symbols from the equations, see :pull:`1331`.

Major changes
=============

Compatibility breaks
====================

* Removed ``itermonomials()`` and ``topological_sort()`` functions, see :pull:`1321` and :pull:`1322`.
* Removed ``Float.num`` property, use :func:`mpmath.mpmathify`, see :pull:`1323`.
* Removed support for CPython 3.10, see :pull:`1344`.
* Removed ``rcall()`` method of :class:`~diofant.core.basic.Basic`, see :pull:`1346`.
* Removed ``method`` argument of :func:`~diofant.functions.special.bessel.jn_zeros`, see :pull:`1352`.
* Removed ``sstrrepr()`` function, see :pull:`1362`.
* Removed support for ASCII pretty-printing and ``pprint_use_unicode()`` function, see :pull:`1369`.
* Removed ``bottom_up()``, ``has_variety()`` and ``has_dups()`` functions, see :pull:`1380`.
* Removed ``diofant.tensor.tensor`` module, see :pull:`1380`.
* Removed ``symarray()`` function, see :pull:`1383`.
* Removed ``base`` kwarg for :func:`~diofant.core.power.log`, see :pull:`1407`.
* Removed ``ln`` alias for :func:`~diofant.core.power.log`, see :pull:`1410`.

Minor changes
=============

* Support CPython 3.12, see :pull:`1325`.
* Support :meth:`~diofant.core.expr.Expr.as_real_imag` for :class:`~diofant.polys.rootoftools.RootOf`, see :pull:`1415`.
* Support CPython 3.13, see :pull:`1425`.

Developer changes
=================

Issues closed
=============

See the `release milestone <https://github.com/diofant/diofant/milestone/9?closed=1>`_
for complete list of issues and pull requests involved in this release.

These Sympy issues also were addressed:

* :sympyissue:`25142`: incorrect simplification of a complex relational
* :sympyissue:`19813`: logcombine hangs
* :sympyissue:`22450`: Rational raised to the big power hangs
* :sympyissue:`25165`: Series expansion not working
* :sympyissue:`25197`: Simple exponential integral error in an otherwise case
* :sympyissue:`23399`: Simplifying equation with function seemingly gets stuck
* :sympyissue:`20427`: Result from clear_denoms() prints like zero poly but behaves wierdly (due to unstripped DMP)
* :sympyissue:`2720` eliminate()
* :sympyissue:`16951`: integrate(sqrt(2*m*(E - x)), x)
* :sympyissue:`25341`: CoercionFailed on eq: 2*sqrt(x)/(x + 1)**2 - 1/(sqrt(x)*(x + 1)) - 1/(4*x**(3/2)))/(x + 1) = 0
* :sympyissue:`20327`: Finite Field coercion fails from Rational type
* :sympyissue:`25406`: Resultant of Polynomials Returns Wrong Output
* :sympyissue:`25451`: Incorrect simplification when mixing basic logical operators and equality
* :sympyissue:`25496`: Privileging expr.__class__ over expr.func for reconstruction
* :sympyissue:`25521`: integrate raises HeuristicGCDFailed
* :sympyissue:`25520`: RecursionError in inverse_laplace_transform
* :sympyissue:`25399`: Cannot use typing.Generic[T] with Symbol
* :sympyissue:`25582`: Incorrect limit for atan
* :sympyissue:`25592`: factor_list sometimes generates PolificationFailed errors with algebraic extensions
* :sympyissue:`25590`: simplify produces wrong answer with non-commuting symbols
* :sympyissue:`25572`: simplify reorders noncommutative factors
* :sympyissue:`25603`: Simplifying And boolean operation removes a condition
* :sympyissue:`25612`: Lack of is_real attribute for Mul class
* :sympyissue:`25624`: lcm(-1,1) and lcm(Poly(-1,x), Poly(1,x)) gives different output
* :sympyissue:`25627`: solve does not take positive=True into account
* :sympyissue:`25681`: Issues with limits while using abs function
* :sympyissue:`25682`: Branches for series expansions involving the abs function is not handled correctly
* :sympyissue:`25679`: hypersimp does not work correctly
* :sympyissue:`25698`: n=6000002; int(n*(log(n) + log(log(n)))) takes more than 200 s to compute on sympy.live.org
* :sympyissue:`25701`: TypeError on Eq(2*sign(x + 3)/(5*Abs(x + 3)**(3/5)), 0)
* :sympyissue:`25723`: GCD missing polynomial factor
* :sympyissue:`25738`: Incorrect result of reduce_inequalities involving pi and abs
* :sympyissue:`25697`: can not reduce log inequalities
* :sympyissue:`25806`: Integrate a simple function
* :sympyissue:`25833`: Limit at infinity of arctan(expression that goes to infinity) erroneously gives NaN or it doesn't compute.
* :sympyissue:`25882`: IndexError when run classify_ode
* :sympyissue:`25885`: Wrong result for a limit
* :sympyissue:`25886`: CeortionError in integrate()
* :sympyissue:`25896`: ratint(e,x).diff().equals(e) is not always True (terms lost)
* :sympyissue:`25899`: surprising error message with Poly('a-a')
* :sympyissue:`23843`: Asymptotic series for atan/acot functions shifted to their branch cuts gives wrong answer
* :sympyissue:`25965`: ceiling(CRootOf()) not implemented, leads to exception in Range
* :sympyissue:`25983`: Incorrect result of reduce_inequalities
* :sympyissue:`25991`: Inconsistencies in as_leading_term, Series Expansion, and Limit Computations for Expressions Involving Square Roots
* :sympyissue:`26071`: Definite integral error
* :sympyissue:`26119`: Lambdify crashes on empty tuple
* :sympyissue:`26178`: Wrong result of sqf_list for PolyElement with excluded symbols
* :sympyissue:`26250`: Incorrect limit involving elliptic functions
* :sympyissue:`25786`: Wrong result for a simple integral
* :sympyissue:`26343`: TypeError: Invalid NaN Comparison using dsolve for ODE with ics={v(0) : 0}
* :sympyissue:`26313`: Error result for limit of a piecewise
* :sympyissue:`26477`: Error in integral result using hyper
* :sympyissue:`26497`: factor produces wrong output
* :sympyissue:`26501`: TypeError: '>' not supported between instances of 'Poly' and 'int' calling integrate in sympy 1.12
* :sympyissue:`26503`: TypeError: Invalid NaN comparison calling integrate in sympy 1.12
* :sympyissue:`26504`: IndexError: Index out of range: calling integrate in sympy 1.12
* :sympyissue:`26506`: RecursionError: maximum recursion depth exceeded in comparison calling integrate in sympy 1.12
* :sympyissue:`26513`: Wrong limit result for Abs((-n/(n+1))**n)
* :sympyissue:`26502`: lots of PolynomialError contains an element of the set of generators exceptions calling integrate in sympy 1.12
* :sympyissue:`14069`: Condition for TODO in zeta_functions.py is now satisfied
* :sympyissue:`25931`: Possible improvements in gruntz
* :sympyissue:`26566`: wrong result integrating sin(a*(x+pi))**2
* :sympyissue:`26571`: integrate((x**8+1)**(-1/2),x) raises _CoeffExpValueError
* :sympyissue:`26577`: factor should use the same domain for numerator and denominator
* :sympyissue:`26682`: Incorrect behaviour for solve_poly_system
* :sympyissue:`26787`: Poly.count_roots wrong for algebraic fields because is_negative is wrong
* :sympyissue:`14738`: dup_count_complex_roots() can't handle degenerate cases
* :sympyissue:`26856`: very simple example of series bug
* :sympyissue:`26873`: Saving a long expression to load in another python notebook/file
* :sympyissue:`26903`: minpoly returns a not irreducible polynomial
