from itertools import combinations_with_replacement

from ..core import Derivative, Eq, Function, S, Symbol, diff, sympify
from ..core.compatibility import iterable


def euler_equations(L, funcs=(), vars=()):
    r"""
    Find the Euler-Lagrange equations [1]_ for a given Lagrangian.

    Parameters
    ==========

    L : Expr
        The Lagrangian that should be a function of the functions listed
        in the second argument and their derivatives.

        For example, in the case of two functions `f(x,y)`, `g(x,y)` and
        two independent variables `x`, `y` the Lagrangian would have the form:

            .. math:: L\left(f(x,y),g(x,y),\frac{\partial f(x,y)}{\partial x},
                      \frac{\partial f(x,y)}{\partial y},
                      \frac{\partial g(x,y)}{\partial x},
                      \frac{\partial g(x,y)}{\partial y},x,y\right)

        In many cases it is not necessary to provide anything, except the
        Lagrangian, it will be auto-detected (and an error raised if this
        couldn't be done).

    funcs : Function or an iterable of Functions
        The functions that the Lagrangian depends on. The Euler equations
        are differential equations for each of these functions.

    vars : Symbol or an iterable of Symbols
        The Symbols that are the independent variables of the functions.

    Returns
    =======

    eqns : list of Eq
        The list of differential equations, one for each function.

    Examples
    ========

    >>> from diofant import Symbol, Function

    >>> x = Function('x')
    >>> t = Symbol('t')
    >>> L = (x(t).diff(t))**2/2 - x(t)**2/2
    >>> euler_equations(L, x(t), t)
    [Eq(-x(t) - Derivative(x(t), t, t), 0)]
    >>> u = Function('u')
    >>> x = Symbol('x')
    >>> L = (u(t, x).diff(t))**2/2 - (u(t, x).diff(x))**2/2
    >>> euler_equations(L, u(t, x), [t, x])
    [Eq(-Derivative(u(t, x), t, t) + Derivative(u(t, x), x, x), 0)]

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Euler%E2%80%93Lagrange_equation
    """
    L = sympify(L)
    funcs = tuple(funcs) if iterable(funcs) else (funcs,)

    if not funcs:
        funcs = tuple(L.atoms(Function))
    else:
        for f in funcs:
            if not isinstance(f, Function):
                raise TypeError('Function expected, got: %s' % f)

    vars = tuple(vars) if iterable(vars) else (vars,)

    if not vars:
        vars = funcs[0].args
    else:
        vars = tuple(sympify(var) for var in vars)

    if not all(isinstance(v, Symbol) for v in vars):
        raise TypeError('Variables are not symbols, got %s' % vars)

    for f in funcs:
        if not vars == f.args:
            raise ValueError("Variables %s don't match args: %s" % (vars, f))

    order = max(len(d.variables) for d in L.atoms(Derivative)
                if d.expr in funcs)

    eqns = []
    for f in funcs:
        eq = diff(L, f)
        for i in range(1, order + 1):
            for p in combinations_with_replacement(vars, i):
                eq += S.NegativeOne**i * diff(L, diff(f, *p), *p)
        eqns.append(Eq(eq))

    return eqns
