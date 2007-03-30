"""
This module contain solvers for all kinds of equations,
both algebraic (solve) and differential (dsolve).
"""

from sympy import Basic, Symbol, Number, Mul, log, Add, Derivative, \
        sin, cos, integrate

def solve(eq, vars):
    """
    Solves any (supported) kind of equation (not differential).

    Examples
    ================
      >>> from sympy import Symbol
      >>> x, y = Symbol('x'), Symbol('y')
      >>> solve(2*x-3, [x])
      3/2

    """

    #currently only solve for one function
    if len(vars) == 1:
        x = vars[0]
        a,b,c = [Symbol(s, dummy = True) for s in ["a","b","c"]]

        r = eq.match(a*x + b, [a,b])
        if r and wo(r,x): return solve_linear(r[a], r[b])

        r = eq.match(a*x**2 + c, [a,c])
        if r and wo(r,x): return solve_quadratic(r[a], 0, r[c])

        r = eq.match(a*x**2 + b*x + c, [a,b,c])
        if r and wo(r,x): return solve_quadratic(r[a], r[b], r[c])

    raise "Sorry, can't solve it (yet)."

def solve_linear(a, b):
    return -b/a

def solve_quadratic(a, b, c):
    D = b**2-4*a*c
    if D == 0:
        return [-b/(2*a)]
    else:
        return [
                (-b+D.sqrt())/(2*a),
                (-b-D.sqrt())/(2*a)
               ]

def dsolve(eq, funcs):
    """
    Solves any (supported) kind of differential equation.

    """

    #currently only solve for one function
    if len(funcs) == 1:
        f = funcs[0]
        x = f[0]
        a,b,c = [Symbol(s, dummy = True) for s in ["a","b","c"]]

        r = eq.match(a*Derivative(f,x) + b, [a,b])
        if r and wo(r,f): return solve_ODE_first_order(r[a], r[b], f, x)

    raise "Sorry, can't solve it (yet)."

def solve_ODE_first_order(a, b, f, x):
    return integrate(-b/a, x) + Symbol("C1")


def wo(di, x):
    """Are all items in the dictionary "di" without "x"?"""
    for d in di:
        if di[d].has(x):
            return False
    return True
