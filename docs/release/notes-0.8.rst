===========
Diofant 0.8
===========

7 Nov 2016

New features
============

* MrvAsympt algorithm to find asymptotic expansion, see :func:`~diofant.core.expr.Expr.aseries` method and :pull:`6`.  Thanks to Avichal Dayal.
* :func:`~diofant.concrete.summations.Sum.findrecur` method to find recurrence relations (with Sister Celine's algorithm), see :pull:`15`.
* Support for :class:`~diofant.core.power.Pow`/:class:`~diofant.functions.elementary.exponential.log` branch-cuts in limits, see :pull:`140`.
* Added basic optimization package, see :func:`~diofant.calculus.optimization.minimize` and :pull:`108`.
* Cartesian product of iterables using Cantor pairing, see :func:`~diofant.utilities.iterables.cantor_product` and :pull:`276`.
* :class:`~diofant.sets.fancysets.Rationals` set, :pull:`255`.
* New simple and robust solver for systems of linear ODEs, see :pull:`286`.  Thanks to Colin B. Macdonald.
* Added mutable/immutable N-dim arrays, sparse and dense, see :pull:`275`.
* :func:`~diofant.solvers.ode.dsolve` now support initial conditions for ODEs, see :pull:`307`.  Thanks to Aaron Meurer.

Major changes
=============

* Depend on setuptools, see :pull:`44`.
* :mod:`The Gruntz Algorithm <diofant.series.gruntz>` reimplemented correctly, see :pull:`68`.
* Replaced ``exp(x)`` with ``E**x`` internally, see :pull:`79`.
* Used :func:`~diofant.printing.repr.srepr` instead of :func:`~diofant.printing.str.sstr` for :meth:`~object.__repr__` printing, see :pull:`39`.
* Major cleanup for series methods, see :pull:`187`.
* Depend on cachetools to implement caching, see :pull:`72` and :pull:`209`.
* Assumption system (old) was validated (:pull:`316` and :pull:`334`) and improved:

    * 0 now is imaginary, see :pull:`8`
    * extended_real fact added, reals are finite now, see :pull:`36`
    * complex are finite now, see :pull:`42`.
    * added docstrings for assumption properties, see :pull:`354`.

Compatibility breaks
====================

* Removed physics submodule, see :pull:`23`.
* Removed galgebra submodule, see :pull:`45`.
* Removed pyglet plotting, see :pull:`50`.
* Removed TextBackend from plotting, see :pull:`67`.
* Removed SageMath support, see :pull:`84`.
* Removed unify submodule, see :pull:`88`.
* Removed crypto submodule, see :pull:`102`.
* Removed print_gtk, see :pull:`114`.
* Unbundle strategies module, see :pull:`103`.
* Removed "old" argument for match/matches, see :pull:`141`.
* Removed when_multiple kwarg in Piecewise, see :pull:`156`.
* Support for Python 2 was removed, see :pull:`160`.
* Removed core.py, see :pull:`60` and :pull:`164`.
* Removed S(foo) syntax, see :pull:`115`.
* Removed (new) assumptions submodule, see :pull:`122`.
* Removed undocumented Symbol.__call__, see :pull:`201`
* Removed categories and liealgebras submodules, see :pull:`280`.
* Rename module sympy -> diofant, see :pull:`315`.
* Use gmpy2, drop gmpy support, see :pull:`292`.
* Removed redundant dom properties in polys, see :pull:`308`.
* Removed manualintegrate function, see :pull:`279`.

Minor changes
=============

* Add support for bidirectional limits, see :pull:`10`.
* Reimplement :class:`~diofant.functions.elementary.trigonometric.cot`, see :pull:`113`.
* A better implementation of :func:`~diofant.calculus.singularities.singularities`, see :pull:`147`.
* Fix "flip" of arguments in relational expressions, see :pull:`30`.
* Make Gosper code use new dispersion algorithm, see :pull:`205`.  Thanks to Raoul Bourquin.
* Consolidate code for solving linear systems, see :pull:`253`.
* Hacks for automatic symbols and wrapping int's replaced with AST transformers, see :pull:`278` and :pull:`167`.
* Build correct inhomogeneous solution in :func:`~diofant.solvers.recurr.rsolve_hyper`, see :pull:`298`.
* Evaluate matrix powers for non-diagonalizable matrices, see :pull:`275`.
* Support non-orthogonal Jordan blocks, see :pull:`275`.
* Make risch_integrate(x**x, x) work, see :pull:`275`.
* Support CPython 3.6, see :pull:`337` and :pull:`356`.

Developer changes
=================

* Unbundle numpydoc, see :pull:`26`.
* Deprecate AUTHORS file, all credits go to the aboutus.rst, see :pull:`87`.
* Use python's :func:`~tokenize.tokenize`, see :pull:`120`.
* Drop using bundled pytest fork, depend on pytest for testing, see :pull:`38`, :pull:`152`, :pull:`91`, :pull:`48`, :pull:`90`, :pull:`96` and :pull:`99`.
* Adopt No Code Of Conduct, see :pull:`212`.
* Measure code coverage, enable codecov.io reports.  See :pull:`217`.
* Adopt pep8 (:pull:`2`) and then flake8 (:pull:`214`) for code quality testing.
* Add regression tests with DIOFANT_USE_CACHE=False :pull:`323`.
* Add interface tests, see :pull:`219` and :pull:`307`.
* Test for no DeprecationWarning in the codebase, see :pull:`356`.

Issues closed
=============

See the `release milestone <https://github.com/diofant/diofant/milestone/1?closed=1>`_
for complete list of issues and pull requests involved in this release.

These Sympy issues also were addressed:

* :sympyissue:`9351` order-1 series wrong with non-zero expansion point
* :sympyissue:`9034` Unicode printing problem with mixture of logs and powers
* :sympyissue:`7927` pretty print incorrect result with powers of sin
* :sympyissue:`9283` KroneckerDelta(p, 0) raises IndexError
* :sympyissue:`9274` Wrong Jordan form: complex eigenvalues w/ geo. mult. > alg. mult.
* :sympyissue:`9398` Simplify of small imaginary number yields 0
* :sympyissue:`7259` LambertW has no series expansion at x=0 (nan)
* :sympyissue:`9832` x**2 < oo returns True but x < oo un-evaluated for real x
* :sympyissue:`9053` MatMul(2, Matrix(...)).doit() doesn't do it
* :sympyissue:`9052` trace(2*A) != 2*Trace(A) because LHS still has an MatMul
* :sympyissue:`9533` Logical operators in octave_code
* :sympyissue:`9545` Mod(zoo, 0) causes RunTime Error
* :sympyissue:`9652` Fail in plot_implicit test on OSX 10.8.5
* :sympyissue:`8432` Tests fail, seems like Cython is not configured to compile with numpy correctly
* :sympyissue:`9542` codegen octave global vars should print "global foo" at top of function
* :sympyissue:`9326` Bug with Dummy
* :sympyissue:`9413` Circularity in assumptions of products
* :sympyissue:`8840` sympy fails to construct (1 + x)*x with disabled cache
* :sympyissue:`4898` Replace exp(x) with E**x internally
* :sympyissue:`10195` Simplification bug on alternating series.
* :sympyissue:`10196` reduce_inequalities error
* :sympyissue:`10198` solving abs with negative powers
* :sympyissue:`7917` Implement cot as a ReciprocalTrigonometricFunction
* :sympyissue:`8649` If t is transcendental, t**n is determined (wrongly) to be non-integer
* :sympyissue:`5641` Compatibility with py.test
* :sympyissue:`10258` Relational involving Piecewise evaluates incorrectly as True
* :sympyissue:`10268` solving inequality involving exp fails for large values
* :sympyissue:`10237` improper inequality reduction
* :sympyissue:`10255` solving a Relational involving Piecewise fails
* :sympyissue:`10290` Computing series where the free variable is not just a symbol is broken
* :sympyissue:`10304` Equality involving expression with known real part and 0 should evaluate
* :sympyissue:`9471` Wrong limit with log and constant in exponent
* :sympyissue:`9449` limit fails with "maximum recursion depth exceeded" / Python crash
* :sympyissue:`8462` Trivial bounds on binomial coefficients
* :sympyissue:`9917` O(n*sin(n) + 1, (n, oo)) returns O(n*sin(n), (n, oo))
* :sympyissue:`7383` Integration error
* :sympyissue:`7098` Incorrect expression resulting from integral evaluation
* :sympyissue:`10323` bad ceiling(sqrt(big integer))
* :sympyissue:`10326` Interval(-oo, oo) contains oo
* :sympyissue:`10095` simplify((1/(2*E))**oo) returns nan
* :sympyissue:`4187` integrate(log(x)*exp(x), (x, 0, oo)) should return -EulerGamma
* :sympyissue:`10383` det of empty matrix is 1
* :sympyissue:`10382` limit(fibonacci(n + 1)/fibonacci(n), n, oo) does not give GoldenRatio
* :sympyissue:`10388` factorial2 runs into RunTimeError for non-integer
* :sympyissue:`10391` solve((2*x + 8)*exp(-6*x), x) can't find any solution
* :sympyissue:`8241` Wrong assumption/result in a parametric limit
* :sympyissue:`3539` Symbol.__call__ should not create a Function
* :sympyissue:`7216` Limits involving branch cuts of elementary functions not handled
* :sympyissue:`10503` Series return an incorrect result
* :sympyissue:`10567` Integral(v,t).doit() differs from integrate(v,t)
* :sympyissue:`9075` sympy.limit yields incorrect result
* :sympyissue:`10610` limit(3**n*3**(-n - 1)*(n + 1)**2/n**2, n, oo) is wrong
* :sympyissue:`4173` implement maximize([x**(1/x), x>0], x)
* :sympyissue:`10803` Bad pretty printing of power of Limit
* :sympyissue:`10836` Latex generation error for .series expansion for \rightarrow term
* :sympyissue:`9558` Bug with limit
* :sympyissue:`4949` solve_linear_system contains duplicate rref algorithm
* :sympyissue:`5952` Standard sets (ZZ, QQ, RR, etc.) for the sets module
* :sympyissue:`9608` Partition can't be ordered
* :sympyissue:`10961` fractional order Laguerre gives wrong result
* :sympyissue:`10976` incorrect answer for limit involving erf
* :sympyissue:`10995` acot(-x) evaluation
* :sympyissue:`11011` Scientific notation should be delimited for LaTeX
* :sympyissue:`11062` Error while simplifying equations containing csc and sec using trigsimp_groebner
* :sympyissue:`10804` 1/limit(airybi(x)*root(x, 4)*exp(-2*x**(S(3)/2)/3), x, oo)**2 is wrong
* :sympyissue:`11063` Some wrong answers from rsolve
* :sympyissue:`9480` Matrix.rank() incorrect results
* :sympyissue:`10497` next(iter(S.Integers*S.Integers)) hangs (expected (0, 0), ...)
* :sympyissue:`5383` Calculate limit error
* :sympyissue:`11270` Limit erroneously reported as infinity
* :sympyissue:`5172` limit() throws TypeError: an integer is required
* :sympyissue:`7055` Failures in rsolve_hyper from test_rsolve_bulk()
* :sympyissue:`11261` Recursion solver fails
* :sympyissue:`11313` Series of Derivative
* :sympyissue:`11290` 1st_exact_Integral wrong result
* :sympyissue:`10761` (1/(x**-2 + x**-3)).series(x, 0) gives wrong result
* :sympyissue:`10024` Eq( Mod(x, 2*pi), 0 ) evaluates to False
* :sympyissue:`7985` Indexed should work with subs on a container
* :sympyissue:`9637` S.Reals - FiniteSet(n) returns EmptySet - FiniteSet(n)
* :sympyissue:`10003` P(X < -1) of ExponentialDistribution
* :sympyissue:`10052` P(X < oo ) for any Continuous Distribution raises AttributeError
* :sympyissue:`10063` Integer raised to Float power does not evaluate
* :sympyissue:`10075` X.pdf(x) for Symbol x returns 0
* :sympyissue:`9823` Matrix power of identity matrix fails
* :sympyissue:`10156` do not use has() to test against self.variables when factoring Sum
* :sympyissue:`10113` imageset(lambda x: x**2/(x**2 - 4), S.Reals) returns (1, oo)
* :sympyissue:`10020` oo**I raises RunTimeError
* :sympyissue:`10240` Not(And(x>2, x<3)) does not evaluate
* :sympyissue:`8510` Differentiation of general functions
* :sympyissue:`10220` Matrix.jordan_cells() fails
* :sympyissue:`10092` subs into inequality involving RootOf raises GeneratorsNeeded
* :sympyissue:`10161` factor gives an invalid expression
* :sympyissue:`10243` Run the examples during automated testing or at release
* :sympyissue:`10274` The helpers kwarg in autowrap method is probably broken.
* :sympyissue:`10210` LaTex printing of Cycle
* :sympyissue:`9539` diophantine(6\*k + 9\*n + 20\*m - x) gives TypeError: unsupported operand type(s) for \*: 'NoneType' and 'Symbol'
* :sympyissue:`11407` Series expansion of the square root gives wrong result
* :sympyissue:`11413` Wrong result from Matrix norm
* :sympyissue:`11434` Matrix rank() produces wrong result
* :sympyissue:`11526` Different result of limit after simplify
* :sympyissue:`11553` Polynomial solve with GoldenRatio causes Traceback
* :sympyissue:`8045` make all NaN is_* properties that are now None -> False (including is_complex)
* :sympyissue:`11602` Replace \dots with \ldots or \cdots
* :sympyissue:`4720` Initial conditions in dsolve()
* :sympyissue:`11623` Wrong groebner basis
* :sympyissue:`10292` poly cannot generically be rebuilt from its args
* :sympyissue:`6572` Remove "#doctest: +SKIP" comments on valid docstrings
* :sympyissue:`10134` Remove "raise StopIteration"
* :sympyissue:`11672` limit(Rational(-1,2)**k, k, oo) fails
* :sympyissue:`11678` Invalid limit of floating point matrix power
* :sympyissue:`11746` undesired (wrong) substition behavior in sympy?
* :sympyissue:`3904` missing docstrings in core
* :sympyissue:`3112` Asymptotic expansion
* :sympyissue:`9173` Series/limit fails unless expression is simplified first.
* :sympyissue:`9808` Complements with symbols should remain unevaluated
* :sympyissue:`9341` Cancelling very long polynomial expression
* :sympyissue:`9908` Sum(1/(n**3 - 1), (n, -oo, -2)).doit() raise UnboundLocalVariable
* :sympyissue:`6171` Limit of a piecewise function
* :sympyissue:`9276` ./bin/diagnose_imports: does it work at all?!
* :sympyissue:`10201` Solution of "first order linear non-homogeneous ODE-System" is wrong
* :sympyissue:`9057` segfault on printing Integral of phi(t)
* :sympyissue:`11159` Substitution with E fails
* :sympyissue:`2839` init_session(auto_symbols=True) and init_session(auto_int_to_Integer=True) do not work
* :sympyissue:`11081` where possible, use python fractions for Rational
* :sympyissue:`10974` solvers.py contains BOM character
* :sympyissue:`10806` LaTeX printer: Integral not surrounded in brackets
* :sympyissue:`10801` Make limit work with binomial
* :sympyissue:`9549` series expansion: (x**2 + x + 1)/(x**3 + x**2) about oo gives wrong result
* :sympyissue:`4231` add a test for complex integral from wikipedia
* :sympyissue:`8634` limit(x**n, x, -oo) is sometimes wrong
* :sympyissue:`8481` Wrong error raised trying to calculate limit of Poisson PMF
* :sympyissue:`9956` Union(Interval(-oo, oo), FiniteSet(1)) not evaluated
* :sympyissue:`9747` test_piecewise_lambdify fails locally
* :sympyissue:`7853` Deprecation of lambdify converting Matrix -> numpy.matrix
* :sympyissue:`9634` Repeated example in the docstring of hermite
* :sympyissue:`8500` Using and operator vs fuzzy_and while querying assumptions
* :sympyissue:`9192` O(y + 1) = O(1)
* :sympyissue:`7130` Definite integral returns an answer with indefinite integrals
* :sympyissue:`8514` Inverse Laplace transform of a simple function fails after updating from 0.7.5 to 0.7.6
* :sympyissue:`9334` Numexpr must be string argument to lambdify
* :sympyissue:`8229` limit((x**Rational(1,4)-2)/(sqrt(x)-4)**Rational(2, 3), x, 16) NotImplementedError
* :sympyissue:`8061` limit(4**(acos(1/(1+x**2))**2)/log(1+x, 4), x, 0) raises NotImplementedError
* :sympyissue:`7872` Substitution in Order fails
* :sympyissue:`3496` limits for complex variables
* :sympyissue:`2929` limit((x*exp(x))/(exp(x)-1), x, -oo) gives -oo
* :sympyissue:`8203` Why is oo real?
* :sympyissue:`7649` S.Zero.is_imaginary should be True?
* :sympyissue:`7256` use old assumptions in code
* :sympyissue:`6783` Get rid of confusing assumptions
* :sympyissue:`5662` AssocOp._eval_template_is_attr is wrong or misused
* :sympyissue:`5295` Document assumptions
* :sympyissue:`4856` coding style
* :sympyissue:`4555` use pyflakes to identify simple bugs in sympy and fix them
* :sympyissue:`5773` Remove the cmp_to_key() helper function
* :sympyissue:`5484` use sort_key instead of old comparison system
* :sympyissue:`8825` Can't use both weakref's & cache
* :sympyissue:`8635` limit(x**n-x**(n-k), x, oo) sometimes raises NotImplementedError
* :sympyissue:`8157` Non-informative error raised when computing limit of cos(n*pi)
* :sympyissue:`7599` Addition of expression and order term fails
* :sympyissue:`6179` wrong order in series
* :sympyissue:`5415` limit involving multi-arg function (polygamma) fails
* :sympyissue:`2865` gruntz doesn't work properly for big-O with point!=0
* :sympyissue:`5907` integrate(1/(x**2 + a**2)**2, x) is wrong if a is real
* :sympyissue:`11722` series() calculation up to O(t**k) returns invalid coefficients for t**k * log(t)
* :sympyissue:`8804` series expansion of 1/x ignores order parameter
* :sympyissue:`10728` Dummy(commutative=False).is_zero -> False
