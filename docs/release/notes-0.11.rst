============
Diofant 0.11
============

Not Released Yet

New features
============

* Added :func:`~diofant.ntheory.residue_ntheory.discrete_log` to compute discrete logarithms, see :pull:`785`.  Thanks to Gabriel Orisaka.

Major changes
=============

* :class:`~diofant.polys.polytools.Poly` now use sparse polynomial representation (via :class:`~diofant.polys.rings.PolyElement`) internally, see :pull:`795`.

Compatibility breaks
====================

* Removed support for Python 3.5 and 3.6, see :pull:`775`.
* ``is_monomial`` attribute of :class:`~diofant.polys.polytools.Poly` renamed to :attr:`~diofant.polys.polytools.Poly.is_term`, see :pull:`780`.
* Removed ``log()`` helper from :class:`~diofant.domains.RationalField`, see :pull:`787`.
* Removed ``seterr()`` function, see :pull:`794`.
* Removed ``DMP`` class, see :pull:`795`.
* Removed ``ring_series`` module, see :pull:`820`.
* :class:`~diofant.core.relational.Equality` doesn't support single-argument call, see :pull:`828`.
* Removed ``is_nonnegative()`` and ``is_nonpositive()`` methods of :class:`~diofant.domains.domain.Domain` subclasses, see :pull:`834`.
* Former ``fast=True`` option is now a default for :meth:`~diofant.polys.polytools.Poly.intervals` and :meth:`~diofant.polys.polytools.Poly.refine_root`, see :pull:`834`.
* Change order of keyword arguments for :meth:`~diofant.polys.rings.PolyElement.integrate`, see :pull:`834`.
* Removed support for ``dps=''`` in :class:`~diofant.core.numbers.Float`.  Significant digits automatically counted for :class:`int` and :class:`str` inputs, see :pull:`797`.

Minor changes
=============

* Support truncation for elements of :class:`~diofant.domains.RealAlgebraicField` to :class:`int`, see :pull:`788`.
* :class:`~diofant.matrices.Matrix`'s and :class:`~diofant.tensor.array.Array`'s support symbolic indexes, see :pull:`785`.  Thanks to Francesco Bonazzi.
* Added ``AA_FACTOR_METHOD`` configuration option to specify factorization algorithm for polynomials with algebraic coefficients, see :pull:`844`.

Developer changes
=================

* Depend on `sphinxcontrib-bibtex <https://sphinxcontrib-bibtex.readthedocs.io/en/latest/>`_ to track the bibliography, see :pull:`766`.

Issues closed
=============

See the `release milestone <https://github.com/diofant/diofant/milestone/4?closed=1>`_
for complete list of issues and pull requests involved in this release.

These Sympy issues also were addressed:

* :sympyissue:`15943` Wrong result from summation
* :sympyissue:`12163` matematica code printer does not handle floats and derivatives correctly
* :sympyissue:`11642` Geometric sum doesn't evaluate with float base
* :sympyissue:`15984` Value error in limit
* :sympyissue:`7337` Wrong integration result
* :sympyissue:`11600` re and im should work for matrix expressions
* :sympyissue:`16038` solve_poly_system works with integers but not floats
* :sympyissue:`15553` rsolve can not solve this kind of recurrences
* :sympyissue:`11581` conjugate of real expression should not change expression
* :sympyissue:`11976` Typo in ellipse.py
* :sympyissue:`11275` LaTeX printer inconsistent with pretty printer
* :sympyissue:`11841` Function('gamma') pretty prints as Γ
* :sympyissue:`11926` ccode does not accept user_functions for Max and Min
* :sympyissue:`11855` DiracDelta function is zero for nonzero arguments
* :sympyissue:`11955` diophantine gives wrong solution for -4*x**2+4*x*y-y**2+2*x-3
* :sympyissue:`11502` Discrete logarithms
* :sympyissue:`11435` str printing of logic expressions should use operators
* :sympyissue:`12200` coeff docstring is wrong
* :sympyissue:`9123` apart drops term
* :sympyissue:`12177` Wrong result with apart Wrong Result
* :sympyissue:`8129` The probability function does not handle expressions like b>=b
* :sympyissue:`9983` Product(1 + 1/n**(S(2)/3), (n, 1, oo)).doit() raise RunTimeError
* :sympyissue:`11726` pde_separate does not allow expressions as input
* :sympyissue:`11981` powsimp() fails with noncommutative variables
* :sympyissue:`12092` evalf does not call _imp_ recursively
* :sympyissue:`10472` pprint should align the middle of the matrix to the baseline?
* :sympyissue:`11959` diophantine gives wrong solution for -4*x**2+4*x*y-y**2+2*x-3
* :sympyissue:`11944` matrix vstack/hstack can fail with immutable matrix as first argument
* :sympyissue:`11732` Fails operators between Interval and some S.Sets
* :sympyissue:`12178` Empty intersection should be UniversalSet
* :sympyissue:`10681` TypeError: 'Float' object cannot be interpreted as an integer from integrate(r**2*(R**2-r**2)**0.5, r)
* :sympyissue:`11078` TypeError: 'Float' object cannot be interpreted as an integer from integrate((6-x*x)**(1.5))
* :sympyissue:`11877` integrate(log(0.5-x), (x, 0, 0.5)) wrongly produces imaginary part
* :sympyissue:`7337` Wrong integration result
* :sympyissue:`10211` integrate((1/sqrt(((y-x)**2 + h**2))**3), (x,0,w), (y,0,w)) is wrong
* :sympyissue:`11806` Incorrectly evaluating integral
* :sympyissue:`12325` string formatting error in dmp_integrate_in
* :sympyissue:`16222` Poly(E**100000000) is slow to create
* :sympyissue:`15413` rootof fails for polynomial with irrational coefficients
* :sympyissue:`16432` a.is_even does not imply a.is_finite
* :sympyissue:`16431` a.is_zero is False does not imply a.is_nonzero is True
* :sympyissue:`16530` (1/x).is_real should be None if x can be zero
* :sympyissue:`16562` Eq with 1 argument is allowed?
* :sympyissue:`16589` roots gives incorrect result
* :sympyissue:`16714` Limit ((n**(n+1) + (n+1)**n) / n**(n+1))**n recursion error
* :sympyissue:`16774` square proportion match has no result
