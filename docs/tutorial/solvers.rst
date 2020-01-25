=========
 Solvers
=========

..
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

    >>> solve(x**2 - x)
    [{x: 0}, {x: 1}]

If no solutions are found, an empty list is returned.

    >>> solve(exp(x))
    []

:func:`~diofant.solvers.solvers.solve` can also solve systems of equations.

    >>> solve([x - y + 2, x + y - 3])
    [{x: 1/2, y: 5/2}]
    >>> solve([x*y - 7, x + y - 6])
    ⎡⎧       ___           ___    ⎫  ⎧     ___             ___    ⎫⎤
    ⎢⎨x: - ╲╱ 2  + 3, y: ╲╱ 2  + 3⎬, ⎨x: ╲╱ 2  + 3, y: - ╲╱ 2  + 3⎬⎥
    ⎣⎩                            ⎭  ⎩                            ⎭⎦

Each solution reported only once:

    >>> solve(x**3 - 6*x**2 + 9*x)
    [{x: 0}, {x: 3}]

To get the solutions of a polynomial including multiplicity use
:func:`~diofant.polys.polyroots.roots`.

    >>> roots(x**3 - 6*x**2 + 9*x)
    {0: 1, 3: 2}

Recurrence Equations
====================

To solve recurrence equations, use
:func:`~diofant.solvers.recurr.rsolve`.  First, create an undefined
function by passing ``cls=Function`` to the
:func:`~diofant.core.symbol.symbols` function.

    >>> f = symbols('f', cls=Function)

We can call ``f(x)``, and it will represent an unknown function application.

.. note::

   From here on in this tutorial we assume that these statements were
   executed:

      >>> from diofant import *
      >>> a, b, c, d, t, x, y, z = symbols('a:d t x:z')
      >>> k, m, n = symbols('k m n', integer=True)
      >>> f, g, h = symbols('f:h', cls=Function)
      >>> init_printing(pretty_print=True, use_unicode=True)

As for algebraic equations, the output is a list of :class:`dict`'s

    >>> rsolve(f(n + 1) - 3*f(n) - 1)
    ⎡⎧        n      1⎫⎤
    ⎢⎨f: n ↦ 3 ⋅C₀ - ─⎬⎥
    ⎣⎩               2⎭⎦

The arbitrary constants in the solutions are symbols of the
form ``C0``, ``C1``, and so on.

Differential Equations
======================

To solve the differential equation

    >>> Eq(f(x).diff(x, x) - 2*f(x).diff(x) + f(x), sin(x))
                          2
             d           d
    f(x) - 2⋅──(f(x)) + ───(f(x)) = sin(x)
             dx           2
                        dx

.. note::

    Derivatives of the unknown function ``f(x)`` are unevaluated.

we would use

    >>> dsolve(_)
            x               cos(x)
    f(x) = ℯ ⋅(C₁ + C₂⋅x) + ──────
                              2

:func:`~diofant.solvers.ode.dsolve` can also solve systems of
equations, like :func:`~diofant.solvers.solvers.solve`.

    >>> dsolve([f(x).diff(x) - g(x), g(x).diff(x) - f(x)])
    ⎡        x       -x             x       -x   ⎤
    ⎣f(x) = ℯ ⋅C₂ - ℯ  ⋅C₁, g(x) = ℯ ⋅C₂ + ℯ  ⋅C₁⎦
