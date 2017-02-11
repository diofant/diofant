===========
Diofant 0.9
===========

Not Released Yet

New features
============

* Polynomial solvers now express all available solutions with :class:`~diofant.polys.rootoftools.RootOf`, see `#400 <https://github.com/diofant/diofant/pull/400>`_.  Only zero-dimentional systems are supported, however.

Major changes
=============

* Assumptions (old) moved from :class:`~diofant.core.basic.Basic` to :class:`~diofant.core.expr.Expr`, see `#311 <https://github.com/diofant/diofant/pull/311>`_.

Backwards-incompatible changes
==============================

* Removed ``assumption0`` property, see  `#382 <https://github.com/diofant/diofant/pull/382>`_.
* :func:`~diofant.core.assumptions.check_assumptions` was moved to :mod:`~diofant.core.assumptions`, see `#387 <https://github.com/diofant/diofant/pull/387>`_.
* ``nsolve()`` function was removed, see `#387 <https://github.com/diofant/diofant/pull/387>`_.
* :attr:`~diofant.core.expr.Expr.is_comparable` and :func:`~diofant.core.expr.Expr.is_hypergeometric` moved to :class:`~diofant.core.expr.Expr`, see `#391 <https://github.com/diofant/diofant/pull/391>`_.
* Removed ``solve_triangulated()`` function, :func:`~diofant.solvers.polysys.solve_biquadratic` and :func:`~diofant.solvers.polysys.solve_poly_system` now use :class:`dict` as output, see `#389 <https://github.com/diofant/diofant/pull/389>`_.
* Dropped support for solving undetermined coefficients in :func:`~diofant.solvers.solvers.solve`, see `#389 <https://github.com/diofant/diofant/pull/389>`_.
* Drop ``intersect()`` alias for :meth:`~diofant.sets.sets.Set.intersection`, see `#396 <https://github.com/diofant/diofant/pull/396>`_.
* Drop ``interactive_traversal()``, see `#395 <https://github.com/diofant/diofant/pull/395>`_.
* Drop ``xring()`` and ``xfield()``, see `#403 <https://github.com/diofant/diofant/pull/403>`_.
* Drop JS printer and ``TableForm`` class, see `#403 <https://github.com/diofant/diofant/pull/403>`_.
* Removed agca submodule of :mod:`~diofant.polys`, see `#404 <https://github.com/diofant/diofant/pull/404>`_.
* Removed ``pager_print()`` and ``print_fcode()``, see `#411 <https://github.com/diofant/diofant/pull/411>`_.
* "Increase" precision of Floats with :meth:`~diofant.core.evalf.EvalfMixin.evalf` now disallowed, see `#380 <https://github.com/diofant/diofant/pull/380>`_.

Minor changes
=============

* New integration heuristics for integrals with :class:`~diofant.functions.elementary.complexes.Abs`, see `#321 <https://github.com/diofant/diofant/pull/321>`_.
* Support unevaluated :class:`~diofant.polys.rootoftools.RootOf`, see `#400 <https://github.com/diofant/diofant/pull/400>`_.
* Sorting of symbolic quadratic roots now same as in :class:`~diofant.polys.rootoftools.RootOf` for numerical coefficients, see `#400 <https://github.com/diofant/diofant/pull/400>`_.
* Support simple first-order DAE with :func:`~diofant.solvers.ode.dsolve` helper :func:`~diofant.solvers.ode.ode_lie_group`, see `#413 <https://github.com/diofant/diofant/pull/413>`_.
* Add support for limits of relational expressions, see `#414 <https://github.com/diofant/diofant/pull/414>`_.

Developer changes
=================

* Enabled docstring testing with flake8, see `#408 <https://github.com/diofant/diofant/pull/408>`_.
* Use only relative imports in the codebase, see `#421 <https://github.com/diofant/diofant/pull/421>`_.
* Enabled flake8-comprehensions plugin, see `#420 <https://github.com/diofant/diofant/pull/420>`_.

Issues closed
=============

* `#376 <https://github.com/diofant/diofant/issues/376>`_ problem with derivative and chain rule
* `#377 <https://github.com/diofant/diofant/issues/377>`_ Substitution of unevaluated Derivatives doesn't ignore bounded symbols
* `sympy/sympy#11879 <https://github.com/sympy/sympy/issues/11879>`_ Strange output from common limit used in elementary calculus
* `sympy/sympy#11884 <https://github.com/sympy/sympy/issues/11884>`_ Addition with Order gives wrong result
* `#370 <https://github.com/diofant/diofant/issues/370>`_ Use git hook for flake8
* `sympy/sympy#11045 <https://github.com/sympy/sympy/issues/11045>`_ integrate(1/(x*sqrt(x**2-1)), (x, 1, 2)) Sympy latest version AttributeError: 'Or' object has no attribute 'lts'
* `sympy/sympy#7165 <https://github.com/sympy/sympy/issues/7165>`_ integrate(abs(y - x**2), (y,0,2)) raises ValueError: gamma function pole
* `sympy/sympy#8733 <https://github.com/sympy/sympy/issues/8733>`_ integrate(abs(x+1), (x, 0, 1)) raises gamma function pole error
* `sympy/sympy#8430 <https://github.com/sympy/sympy/issues/8430>`_ integrate(abs(x), (x, 0, 1)) does not simplify
* `sympy/sympy#12005 <https://github.com/sympy/sympy/issues/12005>`_ Subs._eval_derivative doubles derivatives
* `sympy/sympy#11799 <https://github.com/sympy/sympy/issues/11799>`_ Something wrong with the Riemann tensor?
* `sympy/sympy#12018 <https://github.com/sympy/sympy/issues/12018>`_ solution not found by Sum and gosper_sum
* `sympy/sympy#5649 <https://github.com/sympy/sympy/issues/5649>`_ Bug with AlgebraicNumber.__eq__
* `sympy/sympy#11538 <https://github.com/sympy/sympy/issues/11538>`_ Bug in solve maybe
* `sympy/sympy#12081 <https://github.com/sympy/sympy/issues/12081>`_ integrate(x**(-S(3)/2)*exp(-x), (x, 0, oo)) encounters Runtime Error
* `sympy/sympy#7214 <https://github.com/sympy/sympy/issues/7214>`_ Move old assumptions from Basic to Expr
* `sympy/sympy#4678 <https://github.com/sympy/sympy/issues/4678>`_ Have solve() return RootOf when it can't solve equations
* `sympy/sympy#7789 <https://github.com/sympy/sympy/issues/7789>`_ Poly(...).all_roots fails for general quadratic equation
* `sympy/sympy#8255 <https://github.com/sympy/sympy/issues/8255>`_ roots_quadratic should return roots in same order as Poly.all_roots(radicals=False)
* `sympy/sympy#7138 <https://github.com/sympy/sympy/issues/7138>`_ How to solve system of differential equations with symbolic solution?
* `#161 <https://github.com/diofant/diofant/issues/161>`_ Evalf can increase prec for Floats!
* `sympy/sympy#7457 <https://github.com/sympy/sympy/issues/7457>`_ TypeError when using both multiprocessing and gmpy
* `#309 <https://github.com/diofant/diofant/issues/309>`_ Missing solution for trivial ODE f(t).diff(t)**2 - 1
* `sympy/sympy#12115 <https://github.com/sympy/sympy/issues/12115>`_ Cannot access imported submodules in `sympy.core`
* `sympy/sympy#4315 <https://github.com/sympy/sympy/issues/4315>`_ series expansion of piecewise fails
* `sympy/sympy#6807 <https://github.com/sympy/sympy/issues/6807>`_ atoms does not work correctly in the otherwise case of Piecewise
* `sympy/sympy#12114 <https://github.com/sympy/sympy/issues/12114>`_ solve() leads to ZeroDivisionError: polynomial division
* `#423 <https://github.com/diofant/diofant/issues/423>`_ Problem with expr match by template (a1*x + b1)/(c1*x + d1) + (a2*x + b2)/(c2*x + d2)
* `#66 <https://github.com/diofant/diofant/issues/66>`_ polys todo
* `sympy/sympy#5169 <https://github.com/sympy/sympy/issues/5169>`_ All elements of .args should be Basic
* `sympy/sympy#6249 <https://github.com/sympy/sympy/issues/6249>`_ Problems with MatrixSymbol and simplifying functions
* `sympy/sympy#6426 <https://github.com/sympy/sympy/issues/6426>`_ test_args.py should also test rebuilability

.. last pr: #415

See also full `list of closed issues
<https://github.com/diofant/diofant/issues?q=is%3Aissue+milestone%3A0.9.0+is%3Aclosed>`_
and full `list of merged pull requests
<https://github.com/diofant/diofant/pulls?utf8=%E2%9C%93&q=is%3Apr%20is%3Amerged%20milestone%3A0.9.0>`_
in the Diofant repository.
