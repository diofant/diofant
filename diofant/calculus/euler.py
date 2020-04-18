from itertools import combinations_with_replacement

from ..core import Derivative, Eq, Function, Symbol, diff, sympify
from ..core.compatibility import iterable


def euler_equations(L, funcs=(), vars=()):
    r"""
    Find the Euler-Lagrange equations for a given Lagrangian.

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

    >>> L = (f(t).diff(t))**2/2 - f(t)**2/2
    >>> euler_equations(L, f(t), t)
    [Eq(-f(t) - Derivative(f(t), t, t), 0)]
    >>> L = (f(t, x).diff(t))**2/2 - (f(t, x).diff(x))**2/2
    >>> euler_equations(L, f(t, x), [t, x])
    [Eq(-Derivative(f(t, x), t, t) + Derivative(f(t, x), x, x), 0)]

    References
    ==========

    * https://en.wikipedia.org/wiki/Euler%E2%80%93Lagrange_equation

    """
    L = sympify(L)
    funcs = tuple(funcs) if iterable(funcs) else (funcs,)

    if not funcs:
        funcs = tuple(L.atoms(Function))
    else:
        for f in funcs:
            if not isinstance(f, Function):
                raise TypeError(f'Function expected, got: {f}')

    vars = tuple(vars) if iterable(vars) else (vars,)

    if not vars:
        vars = funcs[0].args
    else:
        vars = tuple(sympify(var) for var in vars)

    if not all(isinstance(v, Symbol) for v in vars):
        raise TypeError(f'Variables are not symbols, got {vars}')

    for f in funcs:
        if not vars == f.args:
            raise ValueError(f"Variables {vars} don't match args: {f}")

    order = max(len(d.variables) for d in L.atoms(Derivative)
                if d.expr in funcs)

    eqns = []
    for f in funcs:
        eq = diff(L, f)
        for i in range(1, order + 1):
            for p in combinations_with_replacement(vars, i):
                eq += (-1)**i * diff(L, diff(f, *p), *p)
        eqns.append(Eq(eq, 0))

    return eqns
