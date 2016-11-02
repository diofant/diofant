===========
Diofant 0.8
===========

Not Released Yet.

New features
============

* MrvAsympt algorithm to find asymptotic expansion, see :func:`~diofant.core.expr.Expr.aseries` method and `#6 <https://github.com/diofant/diofant/pull/6>`_.  Thanks to Avichal Dayal.
* :func:`~diofant.concrete.summations.Sum.findrecur` method to find recurrence relations (with Sister Celine's algorithm), see `#15 <https://github.com/diofant/diofant/pull/15>`_.
* Support for Pow/log branch-cuts in limits, see `#140 <https://github.com/diofant/diofant/pull/140>`_.
* Added basic optimization package, see :func:`~diofant.calculus.optimization.minimize` and `#108 <https://github.com/diofant/diofant/pull/108>`_.
* Cartesian product of iterables using Cantor pairing, see :func:`~diofant.utilities.iterables.cantor_product` and `#276 <https://github.com/diofant/diofant/pull/276>`_.
* :class:`~diofant.sets.fancysets.Rationals` set, `#255 <https://github.com/diofant/diofant/pull/255>`_.
* New simple and robust solver for systems of linear ODEs, see `#286 <https://github.com/diofant/diofant/pull/286>`_.  Thanks to Colin B. Macdonald.
* Added mutable/immutable N-dim arrays, sparse and dense, see `#275 <https://github.com/diofant/diofant/pull/275>`_.
* :func:`~diofant.solvers.ode.dsolve` now support initial conditions for ODEs, see `#307 <https://github.com/diofant/diofant/pull/307>`_.  Thanks to Aaron Meurer.

Major changes
=============

* Depend on setuptools, see `#44 <https://github.com/diofant/diofant/pull/44>`_.
* The Gruntz algorithm reimplemented correctly, see `#68 <https://github.com/diofant/diofant/pull/68>`_.
* Replaced ``exp(x)`` with ``E**x`` internally, see `#79 <https://github.com/diofant/diofant/pull/79>`_.
* Used :func:`~diofant.printing.repr.srepr` instead of :func:`~diofant.printing.str.sstr` for :meth:`~object.__repr__` printing, see `#39 <https://github.com/diofant/diofant/pull/39>`_.
* Major cleanup for series methods, see `#187 <https://github.com/diofant/diofant/pull/187>`_.
* Depend on cachetools to implement caching, see `#72 <https://github.com/diofant/diofant/pull/72>`_ and `#209 <https://github.com/diofant/diofant/pull/209>`_.
* Assumption system (old) was validated (`#316 <https://github.com/diofant/diofant/pull/316>`_ and `#334 <https://github.com/diofant/diofant/pull/334>`_) and improved:

    * 0 now is imaginary, see `#8 <https://github.com/diofant/diofant/pull/8>`_
    * extended_real fact added, reals are finite now, see `#36 <https://github.com/diofant/diofant/pull/36>`_
    * complex are finite now, see `#42 <https://github.com/diofant/diofant/pull/42>`_.
    * added docstrings for assumption properties, see `#354 <https://github.com/diofant/diofant/pull/554>`_.

Backwards-incompatible changes
==============================

* Removed physics submodule, see `#23 <https://github.com/diofant/diofant/pull/23>`_.
* Removed galgebra submodule, see `#45 <https://github.com/diofant/diofant/pull/45>`_.
* Removed pyglet plotting, see `#50 <https://github.com/diofant/diofant/pull/50>`_.
* Removed TextBackend from plotting, see `#67 <https://github.com/diofant/diofant/pull/67>`_.
* Removed SageMath support, see `#84 <https://github.com/diofant/diofant/pull/84>`_.
* Removed unify submodule, see `#88 <https://github.com/diofant/diofant/pull/88>`_.
* Removed crypto submodule, see `#102 <https://github.com/diofant/diofant/pull/102>`_.
* Removed print_gtk, see `#114 <https://github.com/diofant/diofant/pull/114>`_.
* Unbundle strategies module, see `#103 <https://github.com/diofant/diofant/pull/103>`_.
* Removed "old" argument for match/matches, see `#141 <https://github.com/diofant/diofant/pull/141>`_.
* Removed when_multiple kwarg in Piecewise, see `#156 <https://github.com/diofant/diofant/pull/156>`_.
* Support for Python 2 was removed, see `#160 <https://github.com/diofant/diofant/pull/160>`_.
* Removed core.py, see `#60 <https://github.com/diofant/diofant/pull/60>`_ and `#164 <https://github.com/diofant/diofant/pull/164>`_.
* Removed S(foo) syntax, see `#115 <https://github.com/diofant/diofant/pull/115>`_.
* Removed (new) assumptions submodule, see `#122 <https://github.com/diofant/diofant/pull/122>`_.
* Removed undocumented Symbol.__call__, see `#201 <https://github.com/diofant/diofant/pull/201>`_
* Removed categories and liealgebras submodules, see `#280 <https://github.com/diofant/diofant/pull/280>`_.
* Rename module sympy -> diofant, see `#315 <https://github.com/diofant/diofant/pull/315>`_.
* Use gmpy2, drop gmpy support, see `#292 <https://github.com/diofant/diofant/pull/292>`_.
* Removed redundant dom properties in polys, see `#308 <https://github.com/diofant/diofant/pull/308>`_.
* Removed manualintegrate function, see `#279 <https://github.com/diofant/diofant/pull/279>`_.

Minor changes
=============

* Add support for bidirectional limits, see `#10 <https://github.com/diofant/diofant/pull/10>`_.
* Reimplement cot, see `#113 <https://github.com/diofant/diofant/pull/113>`_.
* A better implementation of singularities(), see `#147 <https://github.com/diofant/diofant/pull/147>`_.
* Fix "flip" of arguments in relational expressions, see `#30 <https://github.com/diofant/diofant/pull/30>`_.
* Make Gosper code use new dispersion algorithm, see `#205 <https://github.com/diofant/diofant/pull/205>`_.  Thanks to Raoul Bourquin.
* Consolidate code for solving linear systems, see `#253 <https://github.com/diofant/diofant/pull/253>`_.
* Ugly hacks for automatic symbols and wrapping int's now replaced with AST transformers, see `#278 <https://github.com/diofant/diofant/pull/278>`_ and `#167 <https://github.com/diofant/diofant/pull/167>`_.
* Build correct inhomogeneous solution in rsolve_hyper, see `#298 <https://github.com/diofant/diofant/pull/298>`_.
* Evaluate matrix powers for non-diagonalizable matrices, see `#275 <https://github.com/diofant/diofant/pull/275>`_.
* Support non-orthogonal Jordan blocks, see `#275 <https://github.com/diofant/diofant/pull/275>`_.
* Make risch_integrate(x**x, x) work, see `#275 <https://github.com/diofant/diofant/pull/275>`_.
* Support CPython 3.6, see `#337 <https://github.com/diofant/diofant/pull/337>`_ and `#356 <https://github.com/diofant/diofant/pull/356>`_.

Developer changes
=================

* Unbundle numpydoc, see `#26 <https://github.com/diofant/diofant/pull/26>`_.
* Deprecate AUTHORS file, all credits go to the aboutus.rst, see `#87 <https://github.com/diofant/diofant/pull/87>`_.
* Use python's :func:`~tokenize.tokenize`, see `#120 <https://github.com/diofant/diofant/pull/120>`_.
* Drop using bundled pytest fork, depend on pytest for testing, see `#38 <https://github.com/diofant/diofant/pull/38>`_, `#152 <https://github.com/diofant/diofant/pull/152>`_, `#91 <https://github.com/diofant/diofant/pull/91>`_, `#48 <https://github.com/diofant/diofant/pull/48>`_, `#90 <https://github.com/diofant/diofant/pull/90>`_, `#96 <https://github.com/diofant/diofant/pull/96>`_ and `#99 <https://github.com/diofant/diofant/pull/99>`_.
* Adopt No Code Of Conduct, see `#212 <https://github.com/diofant/diofant/pull/212>`_.
* Measure code coverage, enable codecov.io reports.  See `#217 <https://github.com/diofant/diofant/pull/217>`_.
* Adopt pep8 (`#2 <https://github.com/diofant/diofant/pull/2>`_) and then flake8 (`#214 <https://github.com/diofant/diofant/pull/214>`_) for code quality testing.
* Add regression tests with DIOFANT_USE_CACHE=False `#323 <https://github.com/diofant/diofant/pull/323>`_.
* Add interface tests, see `#219 <https://github.com/diofant/diofant/pull/219>`_ and `#307 <https://github.com/diofant/diofant/pull/307>`_.
* Test for no DeprecationWarning in the codebase, see `#356 <https://github.com/diofant/diofant/pull/356>`_.

Issues closed
=============

* `#3 <https://github.com/diofant/diofant/issues/3>`_ Set up documentation on the readthedocs
* `#20 <https://github.com/diofant/diofant/issues/20>`_ Add CONTRIBUTING.rst
* `#24 <https://github.com/diofant/diofant/issues/24>`_ Remove support for some python versions
* `#46 <https://github.com/diofant/diofant/issues/46>`_ Use rtd theme locally
* `#55 <https://github.com/diofant/diofant/issues/55>`_ limit((x+exp(x))/(x-1), x, -oo) should be 1
* `#56 <https://github.com/diofant/diofant/issues/56>`_ gruntz((ln(x)-1)**(1-sqrt(x)), x, E) should be oo
* `sympy/sympy#9351 <https://github.com/sympy/sympy/issues/9351>`_ order-1 series wrong with non-zero expansion point
* `#16 <https://github.com/diofant/diofant/issues/16>`_ solveset(sinh(x)) doesn't returns all solutions
* `#22 <https://github.com/diofant/diofant/issues/22>`_ Use py.test for testing
* `sympy/sympy#9034 <https://github.com/sympy/sympy/issues/9034>`_ Unicode printing problem with mixture of logs and powers
* `sympy/sympy#7927 <https://github.com/sympy/sympy/issues/7927>`_ pretty print incorrect result with powers of sin
* `sympy/sympy#9283 <https://github.com/sympy/sympy/issues/9283>`_ KroneckerDelta(p, 0) raises IndexError
* `sympy/sympy#9274 <https://github.com/sympy/sympy/issues/9274>`_ Wrong Jordan form: complex eigenvalues w/ geo. mult. > alg. mult.
* `sympy/sympy#9398 <https://github.com/sympy/sympy/issues/9398>`_ Simplify of small imaginary number yields 0
* `sympy/sympy#7259 <https://github.com/sympy/sympy/issues/7259>`_ LambertW has no series expansion at x=0 (nan)
* `#21 <https://github.com/diofant/diofant/issues/21>`_ Remove unsupported and obsoleted modules
* `#124 <https://github.com/diofant/diofant/issues/124>`_ exp(n*x).subs({exp(x): x}) doesn't work for integer symbol n
* `sympy/sympy#9832 <https://github.com/sympy/sympy/issues/9832>`_ ``x**2 < oo`` returns ``True`` but ``x < oo`` un-evaluated for real ``x``
* `sympy/sympy#9053 <https://github.com/sympy/sympy/issues/9053>`_ ``MatMul(2, Matrix(...)).doit()`` doesn't do it
* `sympy/sympy#9052 <https://github.com/sympy/sympy/issues/9052>`_ ``trace(2*A) != 2*Trace(A)`` because LHS still has an MatMul
* `sympy/sympy#9053 <https://github.com/sympy/sympy/issues/9053>`_ ``MatMul(2, Matrix(...)).doit()`` doesn't do it
* `sympy/sympy#9052 <https://github.com/sympy/sympy/issues/9052>`_ ``trace(2*A) != 2*Trace(A)`` because LHS still has an MatMul
* `sympy/sympy#9533 <https://github.com/sympy/sympy/issues/9533>`_ Logical operators in octave_code
* `sympy/sympy#9545 <https://github.com/sympy/sympy/issues/9545>`_ ``Mod(zoo, 0)`` causes RunTime Error
* `sympy/sympy#9652 <https://github.com/sympy/sympy/issues/9652>`_ Fail in plot_implicit test on OSX 10.8.5
* `sympy/sympy#8432 <https://github.com/sympy/sympy/issues/8432>`_ Tests fail, seems like Cython is not configured to compile with numpy correctly
* `sympy/sympy#9542 <https://github.com/sympy/sympy/issues/9542>`_ codegen octave global vars should print "global foo" at top of function
* `sympy/sympy#9326 <https://github.com/sympy/sympy/issues/9326>`_ Bug with Dummy
* `sympy/sympy#9413 <https://github.com/sympy/sympy/issues/9413>`_ Circularity in assumptions of products
* `sympy/sympy#8840 <https://github.com/sympy/sympy/issues/8840>`_ sympy fails to construct (1 + x)*x with disabled cache
* `sympy/sympy#4898 <https://github.com/sympy/sympy/issues/4898>`_ Replace exp(x) with E**x internally
* `#138 <https://github.com/diofant/diofant/issues/138>`_ Wrong polylog.eval for z=-1
* `sympy/sympy#10195 <https://github.com/sympy/sympy/issues/10195>`_ Simplification bug on alternating series.
* `#143 <https://github.com/diofant/diofant/issues/143>`_ powsimp((-1)**(odd/2)) != ImaginaryUnit
* `sympy/sympy#10196 <https://github.com/sympy/sympy/issues/10196>`_ reduce_inequalities error
* `sympy/sympy#10198 <https://github.com/sympy/sympy/issues/10198>`_ solving abs with negative powers
* `sympy/sympy#7917 <https://github.com/sympy/sympy/issues/7917>`_ Implement cot as a ReciprocalTrigonometricFunction
* `sympy/sympy#8649 <https://github.com/sympy/sympy/issues/8649>`_ If t is transcendental, t**n is determined (wrongly) to be non-integer
* `#74 <https://github.com/diofant/diofant/issues/74>`_ Trivial limit's of sign fails
* `#31 <https://github.com/diofant/diofant/issues/31>`_ Wrong automatical cancelation of expr with O terms
* `sympy/sympy#10258 <https://github.com/sympy/sympy/issues/10258>`_ Relational involving Piecewise evaluates incorrectly as True
* `sympy/sympy#10205 <https://github.com/sympy/sympy/issues/10205>`_ 10203: handle Eq and Ne with _solve_inequality
* `sympy/sympy#10268 <https://github.com/sympy/sympy/issues/10268>`_ solving inequality involving exp fails for large values
* `sympy/sympy#10237 <https://github.com/sympy/sympy/issues/10237>`_ improper inequality reduction
* `sympy/sympy#10255 <https://github.com/sympy/sympy/issues/10255>`_ solving a Relational involving Piecewise fails
* `sympy/sympy#10290 <https://github.com/sympy/sympy/issues/10290>`_ Computing series where the free variable is not just a symbol is broken
* `sympy/sympy#10304 <https://github.com/sympy/sympy/issues/10304>`_ Equality involving expression with known real part and 0 should evaluate
* `#148 <https://github.com/diofant/diofant/issues/148>`_ Drop py2 support?
* `sympy/sympy#9471 <https://github.com/sympy/sympy/issues/9471>`_ Wrong limit with log and constant in exponent
* `sympy/sympy#9449 <https://github.com/sympy/sympy/issues/9449>`_ limit fails with "maximum recursion depth exceeded" / Python crash
* `sympy/sympy#8462 <https://github.com/sympy/sympy/issues/8462>`_ Trivial bounds on binomial coefficients
* `sympy/sympy#9917 <https://github.com/sympy/sympy/issues/9917>`_ O(n*sin(n) + 1, (n, oo)) returns O(n*sin(n), (n, oo))
* `sympy/sympy#7383 <https://github.com/sympy/sympy/issues/7383>`_ Integration error
* `sympy/sympy#7098 <https://github.com/sympy/sympy/issues/7098>`_ Incorrect expression resulting from integral evaluation
* `sympy/sympy#10323 <https://github.com/sympy/sympy/issues/10323>`_ bad ceiling(sqrt(big integer))
* `sympy/sympy#10326 <https://github.com/sympy/sympy/issues/10326>`_ Interval(-oo, oo) contains oo
* `sympy/sympy#10095 <https://github.com/sympy/sympy/issues/10095>`_ simplify((1/(2*E))**oo) returns `nan`
* `sympy/sympy#4187 <https://github.com/sympy/sympy/issues/4187>`_ integrate(log(x)*exp(x), (x, 0, oo)) should return -EulerGamma
* `sympy/sympy#10383 <https://github.com/sympy/sympy/issues/10383>`_ det of empty matrix is 1
* `sympy/sympy#10382 <https://github.com/sympy/sympy/issues/10382>`_ limit(fibonacci(n + 1)/fibonacci(n), n, oo) does not give GoldenRatio
* `sympy/sympy#10388 <https://github.com/sympy/sympy/issues/10388>`_ ``factorial2`` runs into ``RunTimeError`` for non-integer
* `sympy/sympy#10391 <https://github.com/sympy/sympy/issues/10391>`_ solve((2*x + 8)*exp(-6*x), x) can't find any solution
* `#32 <https://github.com/diofant/diofant/issues/32>`_ repr printing oddness
* `sympy/sympy#8241 <https://github.com/sympy/sympy/issues/8241>`_ Wrong assumption/result in a parametric limit
* `sympy/sympy#3539 <https://github.com/sympy/sympy/issues/3539>`_ Symbol.__call__ should not create a Function
* `#203 <https://github.com/diofant/diofant/issues/203>`_ Wrong hyperexpand(hyper((-6, -7, -5), (-6, -6), 1))
* `sympy/sympy#7216 <https://github.com/sympy/sympy/issues/7216>`_ Limits involving branch cuts of elementary functions not handled
* `#19 <https://github.com/diofant/diofant/issues/19>`_ Remove obsoleted/redundant docs
* `sympy/sympy#10503 <https://github.com/sympy/sympy/issues/10503>`_ Series return an incorrect result
* `#210 <https://github.com/diofant/diofant/issues/210>`_ Incorrect nseries for cos(x**6)
* `sympy/sympy#10567 <https://github.com/sympy/sympy/issues/10567>`_ Integral(v,t).doit() differs from integrate(v,t)
* `sympy/sympy#9075 <https://github.com/sympy/sympy/issues/9075>`_ sympy.limit yields incorrect result
* `sympy/sympy#10610 <https://github.com/sympy/sympy/issues/10610>`_ limit(3**n*3**(-n - 1)*(n + 1)**2/n**2, n, oo) is wrong
* `#238 <https://github.com/diofant/diofant/issues/238>`_ Wrong coeff in \*_factor_list with RR domain
* `#236 <https://github.com/diofant/diofant/issues/236>`_ simplify(summation(n/((n+2)*(n+4)*(n+8)), (n, 1, oo))) returns 521/25200
* `sympy/sympy#4173 <https://github.com/sympy/sympy/issues/4173>`_ implement maximize([x**(1/x), x>0], x)
* `sympy/sympy#10803 <https://github.com/sympy/sympy/issues/10803>`_ Bad pretty printing of power of Limit
* `sympy/sympy#10836 <https://github.com/sympy/sympy/issues/10836>`_ Latex generation error for .series expansion for \rightarrow term
* `#241 <https://github.com/diofant/diofant/issues/241>`_ Wrong hyperexpand(hyper((2, 3, 5, 9, 1), (1, 4, 6, 10), 1))
* `#172 <https://github.com/diofant/diofant/issues/172>`_ limit(sin(x)**15,x,0,'-') is slow
* `sympy/sympy#9558 <https://github.com/sympy/sympy/issues/9558>`_ Bug with limit
* `#251 <https://github.com/diofant/diofant/issues/251>`_ Random MemoryError in test_gruntz_eval_special
* `sympy/sympy#4949 <https://github.com/sympy/sympy/issues/4949>`_ solve_linear_system contains duplicate rref algorithm
* `#213 <https://github.com/diofant/diofant/issues/213>`_ Consolidate all code for solving linear systems
* `sympy/sympy#5952 <https://github.com/sympy/sympy/issues/5952>`_ Standard sets (ZZ, QQ, RR, etc.) for the sets module
* `sympy/sympy#9608 <https://github.com/sympy/sympy/issues/9608>`_ Partition can't be ordered
* `sympy/sympy#10961 <https://github.com/sympy/sympy/issues/10961>`_ fractional order Laguerre gives wrong result
* `sympy/sympy#10976 <https://github.com/sympy/sympy/issues/10976>`_ incorrect answer for limit involving erf
* `sympy/sympy#10995 <https://github.com/sympy/sympy/issues/10995>`_ acot(-x) evaluation
* `sympy/sympy#11011 <https://github.com/sympy/sympy/issues/11011>`_ Scientific notation should be delimited for LaTeX
* `#263 <https://github.com/diofant/diofant/issues/263>`_ Workaround decreased coverage due to randomness
* `sympy/sympy#11062 <https://github.com/sympy/sympy/issues/11062>`_ Error while simplifying equations containing csc and sec using trigsimp_groebner
* `sympy/sympy#10804 <https://github.com/sympy/sympy/issues/10804>`_ 1/limit(airybi(x)*root(x, 4)*exp(-2*x**(S(3)/2)/3), x, oo)**2 is wrong
* `sympy/sympy#11063 <https://github.com/sympy/sympy/issues/11063>`_ Some wrong answers from rsolve
* `#282 <https://github.com/diofant/diofant/issues/282>`_ Random test failure in master (minimize tests)
* `sympy/sympy#9480 <https://github.com/sympy/sympy/issues/9480>`_ Matrix.rank() incorrect results
* `#288 <https://github.com/diofant/diofant/issues/288>`_ Wrong rank for matrix with det = 0
* `sympy/sympy#10497 <https://github.com/sympy/sympy/issues/10497>`_ next(iter(S.Integers*S.Integers)) hangs (expected (0, 0), ...)
* `sympy/sympy#5383 <https://github.com/sympy/sympy/issues/5383>`_ Calculate limit error
* `sympy/sympy#11270 <https://github.com/sympy/sympy/issues/11270>`_ Limit erroneously reported as infinity
* `#296 <https://github.com/diofant/diofant/issues/296>`_ limit produces bad results with Floats
* `sympy/sympy#5172 <https://github.com/sympy/sympy/issues/5172>`_ limit() throws TypeError: an integer is required
* `sympy/sympy#7055 <https://github.com/sympy/sympy/issues/7055>`_ Failures in rsolve_hyper from test_rsolve_bulk()
* `sympy/sympy#11261 <https://github.com/sympy/sympy/issues/11261>`_ Recursion solver fails
* `#294 <https://github.com/diofant/diofant/issues/294>`_ Wrong rsolve(f(n)-f(n-1)-2*f(n-2)-2*n, f(n))
* `sympy/sympy#11313 <https://github.com/sympy/sympy/issues/11313>`_ Series of Derivative
* `#293 <https://github.com/diofant/diofant/issues/293>`_ classify_sysode should be modified to support mass matrix case in LODE
* `#65 <https://github.com/diofant/diofant/issues/65>`_ Docs todo
* `#215 <https://github.com/diofant/diofant/issues/215>`_ Replace test_code_quality.py with flake8/pep8 tests
* `sympy/sympy#11290 <https://github.com/sympy/sympy/issues/11290>`_ 1st_exact_Integral wrong result
* `sympy/sympy#10761 <https://github.com/sympy/sympy/issues/10761>`_ (1/(x**-2 + x**-3)).series(x, 0) gives wrong result
* `#312 <https://github.com/diofant/diofant/issues/312>`_ Mod(-x, 2*x) should be x, not -x
* `sympy/sympy#10024 <https://github.com/sympy/sympy/issues/10024>`_ Eq( Mod(x, 2*pi), 0 ) evaluates to False
* `sympy/sympy#7985 <https://github.com/sympy/sympy/issues/7985>`_ Indexed should work with subs on a container
* `sympy/sympy#9637 <https://github.com/sympy/sympy/issues/9637>`_ ``S.Reals - FiniteSet(n)`` returns ``EmptySet - FiniteSet(n)``
* `sympy/sympy#10003 <https://github.com/sympy/sympy/issues/10003>`_ P(X < -1) of ExponentialDistribution
* `sympy/sympy#10052 <https://github.com/sympy/sympy/issues/10052>`_ P(X < oo ) for any Continuous Distribution raises AttributeError
* `sympy/sympy#10063 <https://github.com/sympy/sympy/issues/10063>`_ Integer raised to Float power does not evaluate
* `sympy/sympy#10075 <https://github.com/sympy/sympy/issues/10075>`_ X.pdf(x) for Symbol x returns 0
* `sympy/sympy#9823 <https://github.com/sympy/sympy/issues/9823>`_ Matrix power of identity matrix fails
* `sympy/sympy#10156 <https://github.com/sympy/sympy/issues/10156>`_ do not use `has` to test against self.variables when factoring Sum
* `sympy/sympy#10113 <https://github.com/sympy/sympy/issues/10113>`_ imageset(lambda x: x**2/(x**2 - 4), S.Reals) returns (1, ∞)
* `sympy/sympy#10020 <https://github.com/sympy/sympy/issues/10020>`_ oo**I raises RunTimeError
* `sympy/sympy#10240 <https://github.com/sympy/sympy/issues/10240>`_ Not(And(x>2, x<3)) does not evaluate
* `sympy/sympy#8510 <https://github.com/sympy/sympy/issues/8510>`_ Differentiation of general functions
* `sympy/sympy#10220 <https://github.com/sympy/sympy/issues/10220>`_ Matrix.jordan_cells() fails
* `sympy/sympy#10092 <https://github.com/sympy/sympy/issues/10092>`_ subs into inequality involving RootOf raises GeneratorsNeeded
* `sympy/sympy#10161 <https://github.com/sympy/sympy/issues/10161>`_ factor gives an invalid expression
* `sympy/sympy#10243 <https://github.com/sympy/sympy/issues/10243>`_ Run the examples during automated testing or at release
* `sympy/sympy#10274 <https://github.com/sympy/sympy/issues/10274>`_ The helpers kwarg in autowrap method is probably broken.
* `sympy/sympy#10210 <https://github.com/sympy/sympy/issues/10210>`_ LaTex printing of Cycle
* `sympy/sympy#9539 <https://github.com/sympy/sympy/issues/9539>`_ diophantine(6\*k + 9\*n + 20\*m - x) gives TypeError: unsupported operand type(s) for \*: 'NoneType' and 'Symbol'
* `sympy/sympy#11407 <https://github.com/sympy/sympy/issues/11407>`_ Series expansion of the square root gives wrong result
* `sympy/sympy#11413 <https://github.com/sympy/sympy/issues/11413>`_ Wrong result from Matrix norm
* `sympy/sympy#11434 <https://github.com/sympy/sympy/issues/11434>`_ Matrix rank() produces wrong result
* `#135 <https://github.com/diofant/diofant/issues/135>`_ Rename project and adapt imports (sympy -> diofant)
* `#129 <https://github.com/diofant/diofant/issues/129>`_ Use gmpy2 in travis, get rid of gmpy support
* `#133 <https://github.com/diofant/diofant/issues/133>`_ Test regressions with cache on/off
* `#220 <https://github.com/diofant/diofant/issues/220>`_ Update docs/aboutus.rst with more actual info (and move this file?)
* `sympy/sympy#11526 <https://github.com/sympy/sympy/issues/11526>`_ Different result of limit after simplify
* `sympy/sympy#11553 <https://github.com/sympy/sympy/issues/11553>`_ Polynomial solve with GoldenRatio causes Traceback
* `sympy/sympy#8045 <https://github.com/sympy/sympy/issues/8045>`_ make all NaN is_* properties that are now None -> False (including is_complex)
* `#34 <https://github.com/diofant/diofant/issues/34>`_ assumptions todo
* `#203 <https://github.com/diofant/diofant/issues/203>`_ Add changelog (in sphinx docs)
* `sympy/sympy#11553 <https://github.com/sympy/sympy/issues/11553>`_ Polynomial solve with GoldenRatio causes Traceback
* `sympy/sympy#11602 <https://github.com/sympy/sympy/issues/11602>`_ Replace \dots with \ldots or \cdots
* `sympy/sympy#4720 <https://github.com/sympy/sympy/issues/4720>`_ Initial conditions in dsolve()
* `sympy/sympy#11623 <https://github.com/sympy/sympy/issues/11623>`_ Wrong groebner basis
* `sympy/sympy#10292 <https://github.com/sympy/sympy/issues/10292>`_ poly cannot generically be rebuilt from its args
* `#333 <https://github.com/diofant/diofant/issues/333>`_ Expose docs for diofant.interactive (both entry-level and api)
* `#218 <https://github.com/diofant/diofant/issues/218>`_ Remove manualintegrate?
* `sympy/sympy#6572 <https://github.com/sympy/sympy/issues/6572>`_ Remove "#doctest: +SKIP" comments on valid docstrings
* `sympy/sympy#10134 <https://github.com/sympy/sympy/issues/10134>`_ Remove "raise StopIteration"
* `#329 <https://github.com/diofant/diofant/issues/329>`_ Drop examples/
* `sympy/sympy#11672 <https://github.com/sympy/sympy/issues/11672>`_ limit(Rational(-1,2)**k, k, oo) fails
* `#338 <https://github.com/diofant/diofant/issues/338>`_ Rosetta stone for dev's
* `#351 <https://github.com/diofant/diofant/issues/351>`_ Test on CPython 3.6
* `#352 <https://github.com/diofant/diofant/issues/352>`_ Enable testing for DeprecationWarning's
* `sympy/sympy#11678 <https://github.com/sympy/sympy/issues/11678>`_ Invalid limit of floating point matrix power
* `sympy/sympy#11746 <https://github.com/sympy/sympy/issues/11746>`_ undesired (wrong) substition behavior in sympy?
* `sympy/sympy#3904 <https://github.com/sympy/sympy/issues/3904>`_ missing docstrings in core
* `#364 <https://github.com/diofant/diofant/issues/364>`_ Random test failure in combinatorics
* `sympy/sympy#3112 <https://github.com/sympy/sympy/issues/3112>`_ Asymptotic expansion
* `sympy/sympy#9173 <https://github.com/sympy/sympy/issues/9173>`_ Series/limit fails unless expression is simplified first.
* `sympy/sympy#9808 <https://github.com/sympy/sympy/issues/9808>`_ Complements with symbols should remain unevaluated
* `sympy/sympy#9341 <https://github.com/sympy/sympy/issues/9341>`_ Cancelling very long polynomial expression
* `sympy/sympy#9908 <https://github.com/sympy/sympy/issues/9908>`_ Sum(1/(n**3 - 1), (n, -oo, -2)).doit() raise UnboundLocalVariable
* `sympy/sympy#6171 <https://github.com/sympy/sympy/issues/6171>`_ Limit of a piecewise function
* `sympy/sympy#9276 <https://github.com/sympy/sympy/issues/9276>`_ ./bin/diagnose_imports: does it work at all?!
* `sympy/sympy#10201 <https://github.com/sympy/sympy/issues/10201>`_ Solution of "first order linear non-homogeneous ODE-System" is wrong
* `sympy/sympy#9057 <https://github.com/sympy/sympy/issues/9057>`_ segfault on printing Integral of phi(t)
* `sympy/sympy#11159 <https://github.com/sympy/sympy/issues/11159>`_ Substitution with E fails
* `sympy/sympy#2839 <https://github.com/sympy/sympy/issues/2839>`_ init_session(auto_symbols=True) and init_session(auto_int_to_Integer=True) do not work
* `sympy/sympy#11081 <https://github.com/sympy/sympy/issues/11081>`_ where possible, use python fractions for Rational
* `sympy/sympy#10974 <https://github.com/sympy/sympy/issues/10974>`_ solvers.py contains BOM character
* `sympy/sympy#10806 <https://github.com/sympy/sympy/issues/10806>`_ LaTeX printer: Integral not surrounded in brackets
* `sympy/sympy#10801 <https://github.com/sympy/sympy/issues/10801>`_ Make limit work with binomial
* `sympy/sympy#9549 <https://github.com/sympy/sympy/issues/9549>`_ series expansion: (x**2 + x + 1)/(x**3 + x**2) about oo gives wrong result
* `sympy/sympy#4231 <https://github.com/sympy/sympy/issues/4231>`_ add a test for complex integral from wikipedia
* `sympy/sympy#8634 <https://github.com/sympy/sympy/issues/8634>`_ limit(x**n, x, -oo) is sometimes wrong
* `sympy/sympy#8481 <https://github.com/sympy/sympy/issues/8481>`_ Wrong error raised trying to calculate limit of Poisson PMF
* `sympy/sympy#9956 <https://github.com/sympy/sympy/issues/9956>`_ Union(Interval(-oo, oo), FiniteSet(1)) not evaluated
* `sympy/sympy#9747 <https://github.com/sympy/sympy/issues/9747>`_ test_piecewise_lambdify fails locally
* `sympy/sympy#7853 <https://github.com/sympy/sympy/issues/7853>`_ Deprecation of lambdify converting `Matrix` -> `numpy.matrix`
* `sympy/sympy#9634 <https://github.com/sympy/sympy/issues/9634>`_ Repeated example in the docstring of hermite
* `sympy/sympy#8500 <https://github.com/sympy/sympy/issues/8500>`_ Using and operator vs fuzzy_and while querying assumptions
* `sympy/sympy#9192 <https://github.com/sympy/sympy/issues/9192>`_ O(y + 1) = O(1)
* `sympy/sympy#7130 <https://github.com/sympy/sympy/issues/7130>`_ Definite integral returns an answer with indefinite integrals
* `sympy/sympy#8514 <https://github.com/sympy/sympy/issues/8514>`_ Inverse Laplace transform of a simple function fails after updating from 0.7.5 to 0.7.6
* `sympy/sympy#9334 <https://github.com/sympy/sympy/issues/9334>`_ Numexpr must be string argument to lambdify
* `sympy/sympy#8229 <https://github.com/sympy/sympy/issues/8229>`_ limit((x**Rational(1,4)-2)/(sqrt(x)-4)**Rational(2, 3), x, 16) NotImplementedError
* `sympy/sympy#8061 <https://github.com/sympy/sympy/issues/8061>`_ limit(4**(acos(1/(1+x**2))**2)/log(1+x, 4), x, 0) raises NotImplementedError
* `sympy/sympy#7872 <https://github.com/sympy/sympy/issues/7872>`_ Substitution in Order fails
* `sympy/sympy#3496 <https://github.com/sympy/sympy/issues/3496>`_ limits for complex variables
* `sympy/sympy#2929 <https://github.com/sympy/sympy/issues/2929>`_ limit((x*exp(x))/(exp(x)-1), x, -oo) gives -oo
* `sympy/sympy#8203 <https://github.com/sympy/sympy/issues/8203>`_ Why is oo real?
* `sympy/sympy#7649 <https://github.com/sympy/sympy/issues/7649>`_ S.Zero.is_imaginary should be True?
* `sympy/sympy#7256 <https://github.com/sympy/sympy/issues/7256>`_ use old assumptions in code
* `sympy/sympy#6783 <https://github.com/sympy/sympy/issues/6783>`_ Get rid of confusing assumptions
* `sympy/sympy#5662 <https://github.com/sympy/sympy/issues/5662>`_ AssocOp._eval_template_is_attr is wrong or misused
* `sympy/sympy#5295 <https://github.com/sympy/sympy/issues/5295>`_ Document assumptions
* `sympy/sympy#4856 <https://github.com/sympy/sympy/issues/4856>`_ coding style
* `sympy/sympy#4555 <https://github.com/sympy/sympy/issues/4555>`_ use pyflakes to identify simple bugs in sympy and fix them
* `sympy/sympy#5773 <https://github.com/sympy/sympy/issues/5773>`_ Remove the cmp_to_key() helper function
* `sympy/sympy#5484 <https://github.com/sympy/sympy/issues/5484>`_ use sort_key instead of old comparison system
* `sympy/sympy#8825 <https://github.com/sympy/sympy/issues/8825>`_ Can't use both weakref's & cache
* `sympy/sympy#8635 <https://github.com/sympy/sympy/issues/8635>`_ limit(x**n-x**(n-k), x, oo) sometimes raises NotImplementedError
* `sympy/sympy#8157 <https://github.com/sympy/sympy/issues/8157>`_ Non-informative error raised when computing limit of cos(n*pi)
* `sympy/sympy#7872 <https://github.com/sympy/sympy/issues/7872>`_ Substitution in Order fails
* `sympy/sympy#7599 <https://github.com/sympy/sympy/issues/7599>`_ Addition of expression and order term fails
* `sympy/sympy#6179 <https://github.com/sympy/sympy/issues/6179>`_ wrong order in series
* `sympy/sympy#5415 <https://github.com/sympy/sympy/issues/5415>`_ limit involving multi-arg function (polygamma) fails
* `sympy/sympy#2865 <https://github.com/sympy/sympy/issues/2865>`_ gruntz doesn't work properly for big-O with point!=0
* `sympy/sympy#5907 <https://github.com/sympy/sympy/issues/5907>`_ integrate(1/(x**2 + a**2)**2, x) is wrong if a is real
* `sympy/sympy#11722 <https://github.com/sympy/sympy/issues/11722>`_ series() calculation up to O(t**k) returns invalid coefficients for t**k * log(t)
* `#347 <https://github.com/diofant/diofant/issues/347>`_ Search & mention more closed SymPy issues
* `sympy/sympy#8804_ <https://github.com/sympy/sympy/issues/8804>`_ series expansion of 1/x ignores order parameter

See also full `list of closed issues <https://github.com/diofant/diofant/issues?q=is%3Aissue+milestone%3A0.8.0+is%3Aclosed>`_ in the Diofant repository.

Pull requests
=============

* `#1 <https://github.com/diofant/diofant/pull/1>`_ Start the fork, adopt README.txt
* `#4 <https://github.com/diofant/diofant/pull/4>`_ Enhance setup.py
* `#2 <https://github.com/diofant/diofant/pull/2>`_ Add pep8 config, use pep8 in travis
* `#5 <https://github.com/diofant/diofant/pull/5>`_ Don't evaluate derivatives for O expressions
* `#14 <https://github.com/diofant/diofant/pull/14>`_ Set zoo.is_complex to True and zoo.is_real to False
* `#17 <https://github.com/diofant/diofant/pull/17>`_ Replace subs with xreplace (less smart) in the gruntz module
* `#18 <https://github.com/diofant/diofant/pull/18>`_ Remove C (part 1)
* `#8 <https://github.com/diofant/diofant/pull/8>`_ set zero to be imaginary (for old assumptions)
* `#10 <https://github.com/diofant/diofant/pull/10>`_ Add support for bidirectional limits (dir="real")
* `#25 <https://github.com/diofant/diofant/pull/25>`_ Travis tests against pypy3 (not pypy)
* `#26 <https://github.com/diofant/diofant/pull/26>`_ Unbundle numpydoc
* `#6 <https://github.com/diofant/diofant/pull/6>`_ MrvAsympt algorithm to find asymptotic expansion
* `#15 <https://github.com/diofant/diofant/pull/15>`_ Implement findrecur (with Sister Celine's algorithm)
* `#28 <https://github.com/diofant/diofant/pull/28>`_ Fix cross-references in the sphinx documentation
* `#27 <https://github.com/diofant/diofant/pull/27>`_ Removed support for some python versions
* `#29 <https://github.com/diofant/diofant/pull/29>`_ Removed few remaining C imports/exports, fix docs
* `#23 <https://github.com/diofant/diofant/pull/23>`_ Removed physics module
* `#12 <https://github.com/diofant/diofant/pull/12>`_ Q.positive/negative are meaningfull now for Q.extended_real
* `#13 <https://github.com/diofant/diofant/pull/13>`_ Keep trivial sums/products unevaluated
* `#35 <https://github.com/diofant/diofant/pull/35>`_ Add guidelines for contributing (CONTRIBUTING.rst)
* `#36 <https://github.com/diofant/diofant/pull/36>`_ Add extended_real fact for old assumptions.
* `#37 <https://github.com/diofant/diofant/pull/37>`_ Cleanup
* `#40 <https://github.com/diofant/diofant/pull/40>`_ Several modifications for consistency with old assumptions
* `#43 <https://github.com/diofant/diofant/pull/43>`_ Removed doc/python-comparisons.rst
* `#44 <https://github.com/diofant/diofant/pull/44>`_ Use setuptools
* `#41 <https://github.com/diofant/diofant/pull/41>`_ Add noninteger predicate for new assumptions.
* `#45 <https://github.com/diofant/diofant/pull/45>`_ Removed galgebra module
* `#47 <https://github.com/diofant/diofant/pull/47>`_ Remove deprecated stuff
* `#38 <https://github.com/diofant/diofant/pull/38>`_ Use py.test for regular tests and for slow tests
* `#50 <https://github.com/diofant/diofant/pull/50>`_ Removed pyglet plotting
* `#53 <https://github.com/diofant/diofant/pull/53>`_ Remove useless diagnose_imports.py
* `#52 <https://github.com/diofant/diofant/pull/52>`_ Reorder known_facts to be more consistent with sympy/core
* `#49 <https://github.com/diofant/diofant/pull/49>`_ Enable coveralls.io reports
* `#51 <https://github.com/diofant/diofant/pull/51>`_ Use rtd theme
* `#57 <https://github.com/diofant/diofant/pull/57>`_ Use ordered set of monoms in heurisch
* `#60 <https://github.com/diofant/diofant/pull/60>`_ Removed last traces of sympy.core.core.C from SymPy
* `#54 <https://github.com/diofant/diofant/pull/54>`_ Backport some bugfixes from SymPy
* `#62 <https://github.com/diofant/diofant/pull/62>`_ Revert "Removing Kirill from credits."
* `#59 <https://github.com/diofant/diofant/pull/59>`_ Misc fixes
* `#63 <https://github.com/diofant/diofant/pull/63>`_ Revert "Revert "Removing Kirill from credits.""
* `#64 <https://github.com/diofant/diofant/pull/64>`_ Cherry-pick'ed commits from use-py.test-doctests
* `#67 <https://github.com/diofant/diofant/pull/67>`_ Removed TextBackend
* `#70 <https://github.com/diofant/diofant/pull/70>`_ Fix skirpichev/omg#55
* `#69 <https://github.com/diofant/diofant/pull/69>`_ Cleanup of the series docs
* `#71 <https://github.com/diofant/diofant/pull/71>`_ Use set/dict literals, misc fixes
* `#72 <https://github.com/diofant/diofant/pull/72>`_ Revert back new cache stuff (cache.py restored to b4352dd)
* `#68 <https://github.com/diofant/diofant/pull/68>`_ Removed SubsSet in gruntz, use xreplace()
* `#77 <https://github.com/diofant/diofant/pull/77>`_ Fix O.contains expr.is_Add heuristics (was invalid for point != 0)
* `#73 <https://github.com/diofant/diofant/pull/73>`_ Removed "Contributions to docs" section, misc fixes
* `#84 <https://github.com/diofant/diofant/pull/84>`_ Removed sage support
* `#85 <https://github.com/diofant/diofant/pull/85>`_ Removed (broken long time ago) benchmarks support
* `#80 <https://github.com/diofant/diofant/pull/80>`_ Make Q.nonzero compatible with old assumptions
* `#87 <https://github.com/diofant/diofant/pull/87>`_ Deprecate AUTHORS file, all credits go to the aboutus.rst
* `#88 <https://github.com/diofant/diofant/pull/88>`_ Removed (unused, undocumented) unify module
* `#89 <https://github.com/diofant/diofant/pull/89>`_ Restore broken (in sympy) support for matplotlib-enabled tests
* `#91 <https://github.com/diofant/diofant/pull/91>`_ Adopt doctests for py.test + misc fixes
* `#48 <https://github.com/diofant/diofant/pull/48>`_ Enable regular doctest testing with py.test
* `#94 <https://github.com/diofant/diofant/pull/94>`_ Mark more tests as @slow
* `#92 <https://github.com/diofant/diofant/pull/92>`_ Implement helper function _zetas to make zeta tractable by the Gruntz algorithm
* `#90 <https://github.com/diofant/diofant/pull/90>`_ Use py.test to test sphinx docs
* `#96 <https://github.com/diofant/diofant/pull/96>`_ Test examples in travis, runtests.py removed
* `#97 <https://github.com/diofant/diofant/pull/97>`_ Fix infinite recursion for oo**zoo, misc fixes
* `#99 <https://github.com/diofant/diofant/pull/99>`_ Use py.test in setup.py
* `#95 <https://github.com/diofant/diofant/pull/95>`_ Try to preserve decorated signatures
* `#102 <https://github.com/diofant/diofant/pull/102>`_ Removed crypto module
* `#98 <https://github.com/diofant/diofant/pull/98>`_ New set of sympy's fixes
* `#58 <https://github.com/diofant/diofant/pull/58>`_ Improve ipython support
* `#106 <https://github.com/diofant/diofant/pull/106>`_ Travis: Migrating to container-based infrastructure
* `#105 <https://github.com/diofant/diofant/pull/105>`_ Implement nseries helper for LambertW
* `#107 <https://github.com/diofant/diofant/pull/107>`_ Removed old intcache, @cacheit used instead
* `#104 <https://github.com/diofant/diofant/pull/104>`_ Resolve pep8 errors, misc fixes
* `#109 <https://github.com/diofant/diofant/pull/109>`_ Travis: less split for slow tests
* `#100 <https://github.com/diofant/diofant/pull/100>`_ Add Developer's Guide
* `#111 <https://github.com/diofant/diofant/pull/111>`_ Pep8
* `#114 <https://github.com/diofant/diofant/pull/114>`_ Removed print_gtk & sympy/utilities/mathml/
* `#119 <https://github.com/diofant/diofant/pull/119>`_ Removed --split option for pytest
* `#121 <https://github.com/diofant/diofant/pull/121>`_ Change pep8 config defaults: select -> ignore, fix few tests
* `#120 <https://github.com/diofant/diofant/pull/120>`_ use python's tokenize()
* `#118 <https://github.com/diofant/diofant/pull/118>`_ Remove redundant examples
* `#125 <https://github.com/diofant/diofant/pull/125>`_ Fix #124
* `#103 <https://github.com/diofant/diofant/pull/103>`_ Unbundle strategies module
* `#126 <https://github.com/diofant/diofant/pull/126>`_ Misc fixes
* `#130 <https://github.com/diofant/diofant/pull/130>`_ return None -> return, misc fixes
* `#123 <https://github.com/diofant/diofant/pull/123>`_ Fixes sympy/sympy#9832
* `#132 <https://github.com/diofant/diofant/pull/132>`_ Reformat references in the polys module, misc fixes
* `#116 <https://github.com/diofant/diofant/pull/116>`_ New set of sympy's fixes
* `#78 <https://github.com/diofant/diofant/pull/78>`_ Misc no-cache fixes
* `#79 <https://github.com/diofant/diofant/pull/79>`_ Consolidate exp and Pow
* `#136 <https://github.com/diofant/diofant/pull/136>`_ Fix type, returned by Interval._contains
* `#137 <https://github.com/diofant/diofant/pull/137>`_ Fix polylog eval
* `#139 <https://github.com/diofant/diofant/pull/139>`_ Catch NotImplementedError from gruntz
* `#127 <https://github.com/diofant/diofant/pull/127>`_ Travis: use setup.py test
* `#141 <https://github.com/diofant/diofant/pull/141>`_ Removed "old" argument for match/matches
* `#144 <https://github.com/diofant/diofant/pull/144>`_ Stop brave "simplifications" of complex powers with neg bases
* `#142 <https://github.com/diofant/diofant/pull/142>`_ Add a quick exit in _reduce_inequalities if inequality == True/False
* `#146 <https://github.com/diofant/diofant/pull/146>`_ Allow negative powers of abs in the reduce_abs_inequality
* `#113 <https://github.com/diofant/diofant/pull/113>`_ Implement cot as a ReciprocalTrigonometricFunction
* `#147 <https://github.com/diofant/diofant/pull/147>`_ A better implementation of singularities()
* `#150 <https://github.com/diofant/diofant/pull/150>`_ Correct Pow._eval_is_algebraic in case exp is rational
* `#154 <https://github.com/diofant/diofant/pull/154>`_ Add sign._eval_nseries, fixes skirpichev/omg#74
* `#153 <https://github.com/diofant/diofant/pull/153>`_ Fix wrong cancelation of expr with O terms in Add/Mul.flatten
* `#152 <https://github.com/diofant/diofant/pull/152>`_ Last remnants of bundled pytest removed
* `#82 <https://github.com/diofant/diofant/pull/82>`_ Correct Abs._eval_nseries
* `#156 <https://github.com/diofant/diofant/pull/156>`_ Drop errorneous when_multiple kwargs in Piecewise
* `#145 <https://github.com/diofant/diofant/pull/145>`_ Remove _solve_inequality helper
* `#157 <https://github.com/diofant/diofant/pull/157>`_ Fix precision issues in Rel._eval_simplify
* `#151 <https://github.com/diofant/diofant/pull/151>`_ Correct logic of reduce_rational_inequalities
* `#155 <https://github.com/diofant/diofant/pull/155>`_ Support inequalities with piecewise functions
* `#101 <https://github.com/diofant/diofant/pull/101>`_ calculate_leading_term: raise an exception for zero-decision problems
* `#159 <https://github.com/diofant/diofant/pull/159>`_ Improve tutorial (pretty printing), removed support for old IPython versions
* `#158 <https://github.com/diofant/diofant/pull/158>`_ Add a quick exit for Expr.series if x is not a Symbol
* `#160 <https://github.com/diofant/diofant/pull/160>`_ Drop py2 support
* `#166 <https://github.com/diofant/diofant/pull/166>`_ Exclude xfail'ed tests from coverage run
* `#165 <https://github.com/diofant/diofant/pull/165>`_ Simplify Eq/Ne involving expression with known real part and 0
* `#168 <https://github.com/diofant/diofant/pull/168>`_ inspect.getargspec (removed in 3.6) -> getfullargspec
* `#167 <https://github.com/diofant/diofant/pull/167>`_ Replace ugly hack for wrapping int with Integer in the IPython
* `#164 <https://github.com/diofant/diofant/pull/164>`_ Drop use ordering_of_classes and core.py
* `#173 <https://github.com/diofant/diofant/pull/173>`_ Add regression tests for some SymPy's bugs
* `#175 <https://github.com/diofant/diofant/pull/175>`_ Make parallel_poly_from_expr aware of unevaluated Mul
* `#177 <https://github.com/diofant/diofant/pull/177>`_ Add a regression test for sympy/sympy#8016
* `#176 <https://github.com/diofant/diofant/pull/176>`_ Improve Piecewise._eval_interval: support cond's with Abs
* `#179 <https://github.com/diofant/diofant/pull/179>`_ Use mpmath's floor/ceil to calculate round/ceiling, drop get_integer_part()
* `#181 <https://github.com/diofant/diofant/pull/181>`_ Drop redundant ExpBase class
* `#163 <https://github.com/diofant/diofant/pull/163>`_ Make Basic.is_comparable more conservative for extended_real's
* `#184 <https://github.com/diofant/diofant/pull/184>`_ Interval now support extended_real end points, correct S.Reals
* `#42 <https://github.com/diofant/diofant/pull/42>`_ Make complex numbers - finite in old assumptions
* `#183 <https://github.com/diofant/diofant/pull/183>`_ Use more py3 idioms, misc fixes
* `#170 <https://github.com/diofant/diofant/pull/170>`_ Correct Pow.as_numer_denom for cases where base=1, 1/d or n/1
* `#187 <https://github.com/diofant/diofant/pull/187>`_ Major rewrite of ancient garbage in Pow._eval_nseries
* `#186 <https://github.com/diofant/diofant/pull/186>`_ Integral.doit: Vectorize _eval_interval calls only if antideriv has Integral
* `#188 <https://github.com/diofant/diofant/pull/188>`_ Document that det(Matrix()) == 1, misc fixes
* `#115 <https://github.com/diofant/diofant/pull/115>`_ Remove S(foo) syntax from library & tests
* `#174 <https://github.com/diofant/diofant/pull/174>`_ Add some docstrings for gruntz module
* `#189 <https://github.com/diofant/diofant/pull/189>`_ Add rewrite helpers for fibonacci
* `#134 <https://github.com/diofant/diofant/pull/134>`_ Add build_sphinx comand for setup.py
* `#190 <https://github.com/diofant/diofant/pull/190>`_ Fix RuntimeError for factorial2(noninteger)
* `#191 <https://github.com/diofant/diofant/pull/191>`_ Add quick tests to checksol: is_nonzero
* `#192 <https://github.com/diofant/diofant/pull/192>`_ Drop support for "old" order in printers, misc fixes
* `#39 <https://github.com/diofant/diofant/pull/39>`_ Use srepr instead of sstr for __repr__ printing
* `#122 <https://github.com/diofant/diofant/pull/122>`_ Remove new assumptions
* `#197 <https://github.com/diofant/diofant/pull/197>`_ Fixed str() printing of Poly with non-atomic generators
* `#30 <https://github.com/diofant/diofant/pull/30>`_ Fix "flip" of arguments in relational expressions
* `#196 <https://github.com/diofant/diofant/pull/196>`_ Impove coverage
* `#198 <https://github.com/diofant/diofant/pull/198>`_ Fix more pep8 errors, misc fixes
* `#93 <https://github.com/diofant/diofant/pull/93>`_ Complete XPOS todo in Expr.series
* `#202 <https://github.com/diofant/diofant/pull/202>`_ Correct general case in _linear_2eq_order1_type7
* `#199 <https://github.com/diofant/diofant/pull/199>`_ PEP E712
* `#204 <https://github.com/diofant/diofant/pull/204>`_ Fix #203
* `#201 <https://github.com/diofant/diofant/pull/201>`_ Remove undocumented Symbol.__call__ helper
* `#206 <https://github.com/diofant/diofant/pull/206>`_ Enable more PEP8 tests
* `#205 <https://github.com/diofant/diofant/pull/205>`_ Make Gosper code use new dispersion algorithm
* `#140 <https://github.com/diofant/diofant/pull/140>`_ Take into account branch cut for Pow/Log series
* `#207 <https://github.com/diofant/diofant/pull/207>`_ Misc fixes
* `#212 <https://github.com/diofant/diofant/pull/212>`_ Adopt No Code Of Conduct
* `#182 <https://github.com/diofant/diofant/pull/182>`_ Remove __slots__ from core
* `#211 <https://github.com/diofant/diofant/pull/211>`_ Function._eval_nseries: Drop heuristic prediction for number of terms
* `#217 <https://github.com/diofant/diofant/pull/217>`_ Use codecov instead of coveralls
* `#221 <https://github.com/diofant/diofant/pull/221>`_ Add link to aboutus.rst and note about LICENSE in README.rst
* `#219 <https://github.com/diofant/diofant/pull/219>`_ Partial fix for sympy/sympy#4064, test doit
* `#223 <https://github.com/diofant/diofant/pull/223>`_ license stuff
* `#225 <https://github.com/diofant/diofant/pull/225>`_ Optimize travis tests
* `#228 <https://github.com/diofant/diofant/pull/228>`_ Improve collect() docstring
* `#226 <https://github.com/diofant/diofant/pull/226>`_ Don't use xthreaded decorator in integrals
* `#214 <https://github.com/diofant/diofant/pull/214>`_ Use flake8, fix errors N804, N805
* `#222 <https://github.com/diofant/diofant/pull/222>`_ Improve coverage status
* `#61 <https://github.com/diofant/diofant/pull/61>`_ Removed is_Mul heuristic in Limit.doit()
* `#231 <https://github.com/diofant/diofant/pull/231>`_ Test some sympy bugs
* `#233 <https://github.com/diofant/diofant/pull/233>`_ Revert redundant return statement, introduced in ea4ff5a
* `#234 <https://github.com/diofant/diofant/pull/234>`_ Add tests
* `#209 <https://github.com/diofant/diofant/pull/209>`_ Use cachetools
* `#240 <https://github.com/diofant/diofant/pull/240>`_ Try gosper_sum before eval_sum_hyper
* `#237 <https://github.com/diofant/diofant/pull/237>`_ Remove redundant print/sstr/pprint for doctests, misc fixes
* `#108 <https://github.com/diofant/diofant/pull/108>`_ Add minimize/maximize
* `#239 <https://github.com/diofant/diofant/pull/239>`_ Correct wrong coeff for RR domain in \*_factor_list()'s
* `#232 <https://github.com/diofant/diofant/pull/232>`_ Improve coverage
* `#244 <https://github.com/diofant/diofant/pull/244>`_ Add evaluate option for LatticeOp constructor
* `#243 <https://github.com/diofant/diofant/pull/243>`_ Fix pretty printing for powers of Limit's, add regression tests
* `#245 <https://github.com/diofant/diofant/pull/245>`_ Improve coverage
* `#246 <https://github.com/diofant/diofant/pull/246>`_ Use limit in hyperexpand
* `#248 <https://github.com/diofant/diofant/pull/248>`_ Fix some printing bugs, misc fixes
* `#252 <https://github.com/diofant/diofant/pull/252>`_ is_constant should do evalf on results of substitutions 0's and 1's
* `#250 <https://github.com/diofant/diofant/pull/250>`_ Improve coverage
* `#249 <https://github.com/diofant/diofant/pull/249>`_ Fix flake8 errors
* `#253 <https://github.com/diofant/diofant/pull/253>`_ Consolidate code for solving linear systems
* `#255 <https://github.com/diofant/diofant/pull/255>`_ Add primitive implementation for Rationals set, misc fixes
* `#112 <https://github.com/diofant/diofant/pull/112>`_ Improve evaluation of Intersection's for FiniteSet with symbolic elements
* `#258 <https://github.com/diofant/diofant/pull/258>`_ Fix _rebuild in rings like for FracField, misc fixes
* `#259 <https://github.com/diofant/diofant/pull/259>`_ Improve coverage
* `#260 <https://github.com/diofant/diofant/pull/260>`_ Implement _erfs.eval helper
* `#262 <https://github.com/diofant/diofant/pull/262>`_ Improve coverage
* `#261 <https://github.com/diofant/diofant/pull/261>`_ Add notes about acot definition, misc fixes
* `#267 <https://github.com/diofant/diofant/pull/267>`_ Add codecov.yml
* `#264 <https://github.com/diofant/diofant/pull/264>`_ Improve coverage
* `#265 <https://github.com/diofant/diofant/pull/265>`_ Update docs URL: rtfd.org -> rtfd.io, misc fixes
* `#270 <https://github.com/diofant/diofant/pull/270>`_ Update project name references: omg -> diofant
* `#273 <https://github.com/diofant/diofant/pull/273>`_ Rsolve cleanup
* `#277 <https://github.com/diofant/diofant/pull/277>`_ Improve coverage
* `#274 <https://github.com/diofant/diofant/pull/274>`_ Use Fraction for Rational handling, misc fixes
* `#278 <https://github.com/diofant/diofant/pull/278>`_ Replace ugly hack for automatic symbols with ast transformations
* `#280 <https://github.com/diofant/diofant/pull/280>`_ Improve coverage, drop liealgebras and categories modules
* `#272 <https://github.com/diofant/diofant/pull/272>`_ Implement rewrite('tractable') for airyai/airybi
* `#285 <https://github.com/diofant/diofant/pull/285>`_ Improve coverage
* `#284 <https://github.com/diofant/diofant/pull/284>`_ Add regression tests, misc fixes
* `#290 <https://github.com/diofant/diofant/pull/290>`_ Improve coverage
* `#276 <https://github.com/diofant/diofant/pull/276>`_ Cartesian product of iterables using Cantor pairing
* `#291 <https://github.com/diofant/diofant/pull/291>`_ Better zero-equivalence testing in Matrix.rref
* `#289 <https://github.com/diofant/diofant/pull/289>`_ Support Derivative printing in mathematica.py, misc fixes
* `#286 <https://github.com/diofant/diofant/pull/286>`_ dsolve: expm/jordan solver
* `#295 <https://github.com/diofant/diofant/pull/295>`_ Fix getargspec -> getfullargspec, misc fixes
* `#298 <https://github.com/diofant/diofant/pull/298>`_ Build correct inhomogeneous solution in rsolve_hyper
* `#300 <https://github.com/diofant/diofant/pull/300>`_ Support Matrix printing for Mathematica, misc fixes
* `#299 <https://github.com/diofant/diofant/pull/299>`_ Fix "Unknown section" warnings from numpydoc
* `#301 <https://github.com/diofant/diofant/pull/301>`_ This should allow mass matrix in LODE
* `#304 <https://github.com/diofant/diofant/pull/304>`_ Use fractions.Fraction for PythonRational
* `#306 <https://github.com/diofant/diofant/pull/306>`_ Cleanup test_code_quality.py
* `#310 <https://github.com/diofant/diofant/pull/310>`_ Add Relational's printing for Mathematica, misc fixes
* `#313 <https://github.com/diofant/diofant/pull/313>`_ Correct ratio test in Mod.eval
* `#275 <https://github.com/diofant/diofant/pull/275>`_ New set of sympy's fixes
* `#315 <https://github.com/diofant/diofant/pull/315>`_ rename sympy -> diofant
* `#314 <https://github.com/diofant/diofant/pull/314>`_ return (a, b, c, ...) -> return a, b, c, ..., misc fixes
* `#317 <https://github.com/diofant/diofant/pull/317>`_ Cleanup Rational.__new__, reuse Fraction's, misc fixes
* `#318 <https://github.com/diofant/diofant/pull/318>`_ The Diofant's 0.8.0a1 release
* `#292 <https://github.com/diofant/diofant/pull/292>`_ Use gmpy2, drop gmpy support
* `#308 <https://github.com/diofant/diofant/pull/308>`_ Remove redundant .dom (== domain) properties in polys
* `#302 <https://github.com/diofant/diofant/pull/302>`_ Improve coverage
* `#320 <https://github.com/diofant/diofant/pull/320>`_ Version 0.8.0a2
* `#322 <https://github.com/diofant/diofant/pull/322>`_ v0.8.0a2
* `#323 <https://github.com/diofant/diofant/pull/323>`_ Add regression tests with DIOFANT_USE_CACHE=False
* `#324 <https://github.com/diofant/diofant/pull/324>`_ Update docs/aboutus.rst
* `#325 <https://github.com/diofant/diofant/pull/325>`_ Add sanity checks for meijerg parameters
* `#330 <https://github.com/diofant/diofant/pull/330>`_ Add regression test for sympy/sympy#11526
* `#327 <https://github.com/diofant/diofant/pull/327>`_ v0.8.0a3
* `#316 <https://github.com/diofant/diofant/pull/316>`_ Check & fix all assumptions helpers
* `#334 <https://github.com/diofant/diofant/pull/334>`_ Check & fix explicit assumption properties (i.e. is_real = False)
* `#305 <https://github.com/diofant/diofant/pull/305>`_ Add release notes
* `#331 <https://github.com/diofant/diofant/pull/331>`_ v0.8.0a4
* `#341 <https://github.com/diofant/diofant/pull/341>`_ Fix v0.8.0a4
* `#307 <https://github.com/diofant/diofant/pull/307>`_ Support IVP for dsolve
* `#342 <https://github.com/diofant/diofant/pull/342>`_ Stop groebner bases computation, if domain is not exact, like RR
* `#339 <https://github.com/diofant/diofant/pull/339>`_ Test args invariant
* `#335 <https://github.com/diofant/diofant/pull/335>`_ Add sphinx docs for interactive module
* `#279 <https://github.com/diofant/diofant/pull/279>`_ Removed manualintegrate()
* `#343 <https://github.com/diofant/diofant/pull/343>`_ First beta
* `#326 <https://github.com/diofant/diofant/pull/326>`_ Improve coverage
* `#344 <https://github.com/diofant/diofant/pull/344>`_ v0.8.0b2
* `#348 <https://github.com/diofant/diofant/pull/348>`_ Configure rtfd.io builds from file
* `#354 <https://github.com/diofant/diofant/pull/554>`_ Add docstrings for assumptions to the Basic class
* `#337 <https://github.com/diofant/diofant/pull/337>`_ Test Python 3.6
* `#357 <https://github.com/diofant/diofant/pull/357>`_ Readd docs/requirements.txt & rename readthedocs.yml
* `#356 <https://github.com/diofant/diofant/pull/356>`_ Test for DeprecationWarning's
* `#355 <https://github.com/diofant/diofant/pull/355>`_ Run printing setup on (interactive) session startup
* `#359 <https://github.com/diofant/diofant/pull/359>`_ Beta 3

Full `list of merged pull requests <https://github.com/diofant/diofant/pulls?utf8=%E2%9C%93&q=is%3Apr%20is%3Amerged%20milestone%3A0.8.0>`_.
