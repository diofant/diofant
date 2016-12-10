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

Minor changes
=============

Issues closed
=============

* `#376 <https://github.com/diofant/diofant/issues/376>`_ problem with derivative and chain rule
* `#377 <https://github.com/diofant/diofant/issues/377>`_ Substitution of unevaluated Derivatives doesn't ignore bounded symbols
* `sympy/sympy#11879 <https://github.com/sympy/sympy/issues/11879>`_ Strange output from common limit used in elementary calculus
* `sympy/sympy#11884 <https://github.com/sympy/sympy/issues/11884>`_ Addition with Order gives wrong result
* `#370 <https://github.com/diofant/diofant/issues/370>`_ Use git hook for flake8

.. last pr: #389

See also full `list of closed issues
<https://github.com/diofant/diofant/issues?q=is%3Aissue+milestone%3A0.9.0+is%3Aclosed>`_
and full `list of merged pull requests
<https://github.com/diofant/diofant/pulls?utf8=%E2%9C%93&q=is%3Apr%20is%3Amerged%20milestone%3A0.9.0>`_
in the Diofant repository.
