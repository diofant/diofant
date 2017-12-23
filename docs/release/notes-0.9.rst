===========
Diofant 0.9
===========

Not Released Yet

New features
============

* Polynomial solvers now express all available solutions with :class:`~diofant.polys.rootoftools.RootOf`, see :pull:`400`.  Only zero-dimensional systems are supported, however.
* Support solving linear programming problems, see :pull:`283`.

Major changes
=============

* Assumptions (old) moved from :class:`~diofant.core.basic.Basic` to :class:`~diofant.core.expr.Expr`, see :pull:`311`.
* :func:`~diofant.solvers.solvers.solve` now return :class:`list` of :class:`dict`'s, see :pull:`473`.
* ``diofant.polys.domains`` module is now top-level module :mod:`~diofant.domains`, see :pull:`487`.

Compatibility breaks
====================

* Removed ``assumption0`` property, see :pull:`382`.
* :func:`~diofant.core.assumptions.check_assumptions` was moved to :mod:`~diofant.core.assumptions`, see :pull:`387`.
* ``nsolve()`` function was removed, see :pull:`387`.
* :attr:`~diofant.core.expr.Expr.is_comparable` and :meth:`~diofant.core.expr.Expr.is_hypergeometric` moved to :class:`~diofant.core.expr.Expr`, see :pull:`391`.
* Removed ``solve_triangulated()`` and ``solve_biquadratic()`` functions, :func:`~diofant.solvers.polysys.solve_poly_system` now use :class:`dict` as output, see :pull:`389` and :pull:`448`.
* Dropped support for solving undetermined coefficients in :func:`~diofant.solvers.solvers.solve`, see :pull:`389`.
* Drop ``intersect()`` alias for :meth:`~diofant.sets.sets.Set.intersection`, see :pull:`396`.
* Drop ``interactive_traversal()``, see :pull:`395`.
* Drop ``xring()`` and ``xfield()``, see :pull:`403`.
* Drop JS printer and ``TableForm`` class, see :pull:`403`.
* Removed agca submodule of :mod:`~diofant.polys`, see :pull:`404`.
* Removed ``pager_print()`` and ``print_fcode()``, see :pull:`411`.
* "Increase" precision of Floats with :meth:`~diofant.core.evalf.EvalfMixin.evalf` now disallowed, see :pull:`380`.
* Removed ``experimental_lambdify()`` and ``intervalmath`` module from plotting package, see :pull:`384`.
* Removed :func:`~diofant.solvers.solvers.solve` flags ``set``, ``manual``, ``minimal``, ``implicit``, ``particular``, ``quick``, ``exclude``, ``force`` and ``numerical`` see :pull:`426`, :pull:`542` and :pull:`554`.
* Removed support for inequalities in :func:`~diofant.solvers.solvers.solve`, please use :func:`~diofant.solvers.inequalities.reduce_inequalities` instead, see :pull:`426`.
* Removed ``get_domain()`` method of :class:`~diofant.polys.polytools.Poly`, use :attr:`~diofant.polys.polytools.Poly.domain` property instead, see :pull:`479`.
* Renamed 'prec' argument of Float to 'dps', see :pull:`510`.
* Removed 'group' option of :meth:`~diofant.core.basic.Basic.find`, which now return a :class:`dict`.
* Support for Python 3.4 was removed, see :pull:`543`.
* Second argument of :func:`~diofant.solvers.solvers.checksol` must be a :class:`dict`.  See :pull:`549`.
* Removed ``solve_undetermined_coeffs()`` function, see :pull:`554`.
* Make ``matches()`` method for :class:`~diofant.core.basic.Basic` - private, see :pull:`557`.
* Removed :meth:`~diofant.core.basic.Basic.replace` flags ``simultaneous`` and ``map``, see :pull:`557`.
* Make ``strict`` flag - default for :meth:`~diofant.core.evalf.EvalfMixin.evalf`, see :pull:`537`.
* Removed ``I`` property of the :class:`~diofant.matrices.expressions.MatrixExpr`, see :pull:`577`.

Minor changes
=============

* New integration heuristics for integrals with :class:`~diofant.functions.elementary.complexes.Abs`, see :pull:`321`.
* Support unevaluated :class:`~diofant.polys.rootoftools.RootOf`, see :pull:`400`.
* Sorting of symbolic quadratic roots now same as in :class:`~diofant.polys.rootoftools.RootOf` for numerical coefficients, see :pull:`400`.
* Support simple first-order DAE with :func:`~diofant.solvers.ode.dsolve` helper :func:`~diofant.solvers.ode.ode_lie_group`, see :pull:`413`.
* Add support for limits of relational expressions, see :pull:`414`.
* Support rewriting :class:`~diofant.functions.elementary.miscellaneous.Min` and :class:`~diofant.functions.elementary.miscellaneous.Max` as :class:`~diofant.functions.elementary.piecewise.Piecewise`, this allow solving more piecewise equations, see :pull:`426`.
* :func:`~diofant.polys.numberfields.minimal_polynomial` fixed to support generic :class:`~diofant.core.numbers.AlgebraicNumber`'s, see :pull:`433` and :pull:`438`.
* :class:`~diofant.core.numbers.AlgebraicNumber` now support arithmetic operations and exponentiation with integer exponents, see :pull:`428` and :pull:`485`.
* Add AST transformation :class:`~diofant.interactive.session.IntegerDivisionWrapper` to wrap integer division, see :pull:`519`.
* Add AST transformation :class:`~diofant.interactive.session.FloatRationalizer` to wrap :class:`float`'s, see :pull:`538`.
* Support rewrite :class:`~diofant.polys.rootoftools.RootOf` via radicals, see :pull:`563`.

Developer changes
=================

* Enabled docstring testing with flake8, see :pull:`408`.
* Use only relative imports in the codebase, see :pull:`421`.
* Enabled flake8-comprehensions plugin, see :pull:`420`.
* Sort imports with `isort <https://github.com/timothycrosley/isort>`_, see :pull:`520`.
* Depend on `hypothesis <https://hypothesis.readthedocs.io/en/latest/>`_, see :pull:`547`.
* Depend on `pytest-xdist <https://github.com/pytest-dev/pytest-xdist>`_, see :pull:`551`.

Issues closed
=============

See the `release milestone <https://github.com/diofant/diofant/milestone/2?closed=1>`_
for complete list of issues and pull requests involved in this release.

These Sympy issues also were addressed:

* :sympyissue:`11879` Strange output from common limit used in elementary calculus
* :sympyissue:`11884` Addition with Order gives wrong result
* :sympyissue:`11045` integrate(1/(x*sqrt(x**2-1)), (x, 1, 2)) Sympy latest version AttributeError: 'Or' object has no attribute 'lts'
* :sympyissue:`7165` integrate(abs(y - x**2), (y,0,2)) raises ValueError: gamma function pole
* :sympyissue:`8733` integrate(abs(x+1), (x, 0, 1)) raises gamma function pole error
* :sympyissue:`8430` integrate(abs(x), (x, 0, 1)) does not simplify
* :sympyissue:`12005` Subs._eval_derivative doubles derivatives
* :sympyissue:`11799` Something wrong with the Riemann tensor?
* :sympyissue:`12018` solution not found by Sum and gosper_sum
* :sympyissue:`5649` Bug with AlgebraicNumber.__eq__
* :sympyissue:`11538` Bug in solve maybe
* :sympyissue:`12081` integrate(x**(-S(3)/2)*exp(-x), (x, 0, oo)) encounters Runtime Error
* :sympyissue:`7214` Move old assumptions from Basic to Expr
* :sympyissue:`4678` Have solve() return RootOf when it can't solve equations
* :sympyissue:`7789` Poly(...).all_roots fails for general quadratic equation
* :sympyissue:`8255` roots_quadratic should return roots in same order as Poly.all_roots(radicals=False)
* :sympyissue:`7138` How to solve system of differential equations with symbolic solution?
* :sympyissue:`7457` TypeError when using both multiprocessing and gmpy
* :sympyissue:`12115` Cannot access imported submodules in sympy.core
* :sympyissue:`4315` series expansion of piecewise fails
* :sympyissue:`6807` atoms does not work correctly in the otherwise case of Piecewise
* :sympyissue:`12114` solve() leads to ZeroDivisionError: polynomial division
* :sympyissue:`5169` All elements of .args should be Basic
* :sympyissue:`6249` Problems with MatrixSymbol and simplifying functions
* :sympyissue:`6426` test_args.py should also test rebuilability
* :sympyissue:`11461` NameError: name 'Ne' is not defined plotting real_root((log(x/(x-2))), 3)
* :sympyissue:`10925` plot doesn't work with Piecewise
* :sympyissue:`12180` Confusing output from sympy.solve
* :sympyissue:`5786` factor(extension=[I]) gives wrong results
* :sympyissue:`9607` factor - incorrect result
* :sympyissue:`8754` Problem factoring trivial polynomial
* :sympyissue:`8697` rsolve fails to find solutions to some higer order recurrence relations
* :sympyissue:`8694` Match fail
* :sympyissue:`8710` geometry's encloses method fails for non-polygons
* :sympyissue:`10337` bad Boolean args not rejected
* :sympyissue:`9447` sets.Complement fails on certain Unions
* :sympyissue:`10305` Complement Of Universal Subsets
* :sympyissue:`10413` ascii pprint of ProductSet uses non-ascii multiplication symbol
* :sympyissue:`10414` pprint(Union, use_unicode=False) raises error (but str(Union) works)
* :sympyissue:`10375` lambdify on sympy.Min does not work with NumPy
* :sympyissue:`10433`  Dict does not accept collections.defaultdict
* :sympyissue:`9044` pretty printing: Trace could be improved (and LaTeX)
* :sympyissue:`10445` Improper integral does not evaluate
* :sympyissue:`10379` dsolve() converts floats to integers/rationals
* :sympyissue:`10633` Eq(True, False) doesn't evaluate
* :sympyissue:`7163` integrate((sign(x - 1) - sign(x - 2))*cos(x), x) raises TypeError: doit() got an unexpected keyword argument 'manual'
* :sympyissue:`11881` ZeroDivisionError: pole in hypergeometric series random test failure
* :sympyissue:`11801` Exception when printing Symbol('')
* :sympyissue:`11911` typo in docs of printing
* :sympyissue:`10489` Mathematical Symbol does not seem to serialize correctly LaTeX printer
* :sympyissue:`10336` nsimplify problems with oo and inf
* :sympyissue:`12345` nonlinsolve (solve_biquadratic) gives no solution with radical
* :sympyissue:`12375` sympy.series() is broken?
* :sympyissue:`5514` Poly(x, x) * I != I * Poly(x, x)
* :sympyissue:`12398` Limits With abs in certain cases remains unevaluated
* :sympyissue:`12400` polytool.poly() can't raise polynomial to negative power?
* :sympyissue:`12221` Issue with definite piecewise integration
* :sympyissue:`12522` BooleanTrue and Boolean False should have simplify method
* :sympyissue:`12555` limit((3**x + 2 * x**10) / (x**10 + E**x), x, -oo) gives 0 instead of 2
* :sympyissue:`12569` problem with polygamma or im
* :sympyissue:`12578` Taylor expansion wrong (likely because of wrong substitution at point of evaluation?)
* :sympyissue:`12582` Can't solve integrate(abs(x**2-3*x), (x, -15, 15))
* :sympyissue:`12747` Missing constant coefficient in Taylor series of degree 1
* :sympyissue:`12769` Slow limit() calculation?!
* :sympyissue:`12942` Remove x**1.0 == x hack from core
* :sympyissue:`12238` match can take a long time (possibly forever)
* :sympyissue:`4269` ordering of classes
* :sympyissue:`13081` Some comparisons between rational and irrational numbers are incorrect
* :sympyissue:`13078` Return NotImplemented, not False, upon rich comparison with unknown type
* :sympyissue:`13098` sympy.floor() sometimes returns the wrong answer
* :sympyissue:`13312` SymPy does not evaluate integrals of exponentials with symbolic parameter and limit
* :sympyissue:`13111` Don't use "is" to compare classes
* :sympyissue:`10488` integrate(x/(a*x+b), x) gives wrong answer
* :sympyissue:`9706` Interval(-oo, 0).closure hangs
* :sympyissue:`10740` Add a test for Interval(..) in Interval(..) == False
* :sympyissue:`10592` zeta(0, n) where n is negative is wrong
* :sympyissue:`7858` Nth root mod giving wrong solutions
* :sympyissue:`5412` N(oo*I) returns wrong result
* :sympyissue:`10710` Any dict-like object in expr.subs
* :sympyissue:`10810` Implemented function gives ValueError when constructing float expression in sympy 1.0
* :sympyissue:`10867` Getting KeyError while solving ode : dsolve(Eq(g(x).diff(x).diff(x) , (x-2)**2 +(x-3)**3), g(x))
* :sympyissue:`10782` condition_number() for empty matrices giving ValueError
* :sympyissue:`10719` eigenvals of empty matrix raises IndexError
* :sympyissue:`10680` unable to get unevaluated Integral object for  integrate ( x**log (x**log (x**log(x) ) ) , x) .
* :sympyissue:`10701` Is the empty matrix nilpotent? IndexError: Index out of range: a[0]
* :sympyissue:`10770` Adding a row or a column to an empty matrix
* :sympyissue:`10773` sympify evaluates Div Operation in case of Unary Operator when evaluate = False
* :sympyissue:`13332` limit(): AttributeError: 'NoneType' object has no attribute 'rewrite'
* :sympyissue:`13382` Incorrect Result for limit(n*(((n+1)**2+1)/((n)**2+1)-1), n ,oo)
* :sympyissue:`13403` Incorrect Result for limit(n*(-1 + (n + log(n + 1) + 1)/(n + log(n))), n ,oo)
* :sympyissue:`13416` Incorrect Result for limit((-n**3*log(n)**3 + (n - 1)*(n + 1)**2*log(n + 1)**3)/(n**2*log(n)**3), n ,oo)
* :sympyissue:`13462` Bug in sympy.limit()
* :sympyissue:`13501` Incorrect integral of a rational function with a symbolic coefficient
* :sympyissue:`13536` TypeError for integration from infinity to a positive value
* :sympyissue:`13545` Poly loses modulus after arithemetic
* :sympyissue:`13460` Integration of certain cubic rational functions is incorrect
* :sympyissue:`13071` meijerg.is_number is wrong
* :sympyissue:`13575` limit(acos(erfi(x)), x, 1) causes recursion error
* :sympyissue:`13629` bug in rsolve
* :sympyissue:`13645` sympy hangs on evaluating expression
* :sympyissue:`11378` S.Reals should be accessible as just "Reals"
* :sympyissue:`10999` diop: holzer error
* :sympyissue:`11000` diop: power_representation
* :sympyissue:`11026` diophantine(x**3+y**3-2) -> KeyError instead of {(1, 1)}
* :sympyissue:`8943` diophantine misses trivial solution
* :sympyissue:`11016` diop: sum of squares needs to try more options to satisfy conditions
* :sympyissue:`9538` diophantine() doesn't let you specify the variable order
* :sympyissue:`11049` diop: recursion error
* :sympyissue:`11021` diop: power_representation(4**5, 3, 1) -> (4,)
* :sympyissue:`11050` diop: partition(n, k) gives redundant result
* :issue:`369` Issue labeling policy
