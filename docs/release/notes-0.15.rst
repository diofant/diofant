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

Minor changes
=============

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
