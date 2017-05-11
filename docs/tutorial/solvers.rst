=========
 Solvers
=========

..
    >>> from diofant import *
    >>> x, y, z = symbols('x y z')
    >>> init_printing(pretty_print=True, use_unicode=True)

This section covers equations solving.

.. note::

    Any expression in input, that not in an
    :class:`~diofant.core.relational.Eq` is automatically assumed to
    be equal to 0 by the solving functions.

Algebraic Equations
===================

The main function for solving algebraic equations is
:func:`~diofant.solvers.solvers.solve`.

When solving a single equation, the output is a list of the solutions.

    >>> solve(x**2 - x, x)
    [{x: 0}, {x: 1}]

If no solutions are found, an empty list is returned.

    >>> solve(exp(x), x)
    []

:func:`~diofant.solvers.solvers.solve` can also solve systems of equations.

    >>> solve([x - y + 2, x + y - 3], [x, y])
    [{x: 1/2, y: 5/2}]
    >>> solve([x*y - 7, x + y - 6], [x, y])
    ⎡⎧       ___           ___    ⎫  ⎧     ___             ___    ⎫⎤
    ⎢⎨x: - ╲╱ 2  + 3, y: ╲╱ 2  + 3⎬, ⎨x: ╲╱ 2  + 3, y: - ╲╱ 2  + 3⎬⎥
    ⎣⎩                            ⎭  ⎩                            ⎭⎦

:func:`~diofant.solvers.solvers.solve` reports each solution only once.

    >>> solve(x**3 - 6*x**2 + 9*x, x)
    [{x: 0}, {x: 3}]

To get the solutions of a polynomial including multiplicity use
:func:`~diofant.polys.polyroots.roots`.

    >>> roots(x**3 - 6*x**2 + 9*x, x)
    {0: 1, 3: 2}

Differential Equations
======================

To solve differential equations, use
:func:`~diofant.solvers.ode.dsolve`.  First, create an undefined
function by passing ``cls=Function`` to the
:func:`~diofant.core.symbol.symbols` function.

    >>> f, g = symbols('f g', cls=Function)

``f`` and ``g`` are now undefined functions.  We can call ``f(x)``,
and it will represent an unknown function application.  Derivatives of
``f(x)`` are unevaluated.

    >>> f(x).diff(x)
    d
    ──(f(x))
    dx

To represent the differential equation `f''(x) - 2f'(x) + f(x) =
\sin(x)`, we would thus use

    >>> Eq(f(x).diff(x, x) - 2*f(x).diff(x) + f(x), sin(x))
                          2
             d           d
    f(x) - 2⋅──(f(x)) + ───(f(x)) = sin(x)
             dx           2
                        dx

To solve the ODE, pass it and the function to solve for to
:func:`~diofant.solvers.ode.dsolve`.

    >>> dsolve(_, f(x))
            x               cos(x)
    f(x) = ℯ ⋅(C₁ + C₂⋅x) + ──────
                              2

:func:`~diofant.solvers.ode.dsolve` returns an instance of
:class:`~diofant.core.relational.Eq`.  This is because in general,
solutions to differential equations cannot be solved explicitly for
the function.

    >>> dsolve(f(x).diff(x)*(1 - sin(f(x))), f(x))
    f(x) + cos(f(x)) = C₁

The arbitrary constants in the solutions from dsolve are symbols of
the form ``C1``, ``C2``, ``C3``, and so on.

:func:`~diofant.solvers.ode.dsolve` can also solve systems of
equations, like :func:`~diofant.solvers.solvers.solve`.

    >>> dsolve([f(x).diff(x) - g(x), g(x).diff(x) - f(x)], [f(x), g(x)])
    ⎡          ⎛ x    -x⎞      ⎛ x    -x⎞            ⎛ x    -x⎞      ⎛ x    -x⎞⎤
    ⎢          ⎜ℯ    ℯ  ⎟      ⎜ℯ    ℯ  ⎟            ⎜ℯ    ℯ  ⎟      ⎜ℯ    ℯ  ⎟⎥
    ⎢f(x) = C₁⋅⎜── + ───⎟ + C₂⋅⎜── - ───⎟, g(x) = C₁⋅⎜── - ───⎟ + C₂⋅⎜── + ───⎟⎥
    ⎣          ⎝2     2 ⎠      ⎝2     2 ⎠            ⎝2     2 ⎠      ⎝2     2 ⎠⎦
