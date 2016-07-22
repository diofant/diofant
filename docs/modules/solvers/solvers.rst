Solvers
==========

.. module:: diofant.solvers

The *solvers* module in Diofant implements methods for solving equations.

Algebraic equations
--------------------

Use :func:`~diofant.solvers.solvers.solve` to solve algebraic equations. We suppose all equations are equaled to 0,
so solving x**2 == 1 translates into the following code::

    >>> from diofant.solvers import solve
    >>> from diofant import Symbol
    >>> x = Symbol('x')
    >>> solve(x**2 - 1, x)
    [-1, 1]

The first argument for :func:`~diofant.solvers.solvers.solve` is an equation (equaled to zero) and the second argument
is the symbol that we want to solve the equation for.

.. autofunction:: diofant.solvers.solvers.solve

.. autofunction:: diofant.solvers.solvers.solve_linear

.. autofunction:: diofant.solvers.solvers.solve_linear_system

.. autofunction:: diofant.solvers.solvers.solve_undetermined_coeffs

.. autofunction:: diofant.solvers.solvers.nsolve

.. autofunction:: diofant.solvers.solvers.check_assumptions

.. autofunction:: diofant.solvers.solvers.checksol

.. autofunction:: diofant.solvers.solvers.unrad

Ordinary Differential equations (ODEs)
--------------------------------------

See :ref:`ode-docs`.

Partial Differential Equations (PDEs)
-------------------------------------

See :ref:`pde-docs`.

Deutils (Utilities for solving ODE's and PDE's)
-----------------------------------------------

.. autofunction:: diofant.solvers.deutils.ode_order

Recurrence Equations
--------------------

.. module:: diofant.solvers.recurr

.. autofunction:: rsolve

.. autofunction:: rsolve_poly

.. autofunction:: rsolve_ratio

.. autofunction:: rsolve_hyper

Systems of Polynomial Equations
-------------------------------

.. autofunction:: diofant.solvers.polysys.solve_poly_system

.. autofunction:: diofant.solvers.polysys.solve_triangulated

Diophantine Equations (DEs)
---------------------------

See :ref:`diophantine-docs`

Inequalities
------------

See :ref:`inequality-docs`
