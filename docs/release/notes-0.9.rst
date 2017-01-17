===========
Diofant 0.9
===========

Not Released Yet

New features
============

Major changes
=============

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

Minor changes
=============

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

.. last pr: #403

See also full `list of closed issues
<https://github.com/diofant/diofant/issues?q=is%3Aissue+milestone%3A0.9.0+is%3Aclosed>`_
and full `list of merged pull requests
<https://github.com/diofant/diofant/pulls?utf8=%E2%9C%93&q=is%3Apr%20is%3Amerged%20milestone%3A0.9.0>`_
in the Diofant repository.
