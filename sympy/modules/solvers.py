from sympy.core import Basic, Symbol, Number, Mul, Pow, log, Add
from sympy.modules import cos, sin

def solve(eq, vars):
    """
    Solves any (supported) kind of equations (not differential).

    Examples
    ================
      >>> from sympy import Symbol
      >>> x, y = Symbol('x'), Symbol('y')
      >>> solve(2*x-3, [x])
      3/2

    """
    #currently solve only for one variable
    if len(vars) == 1:
        x = vars[0]
        a,b,c = [Symbol(s, dummy = True) for s in ["a","b","c"]]
        r = eq.match(a*x + b, [a,b])
        if r: return solve_linear(r[a], r[b])
        r = eq.match(a*x**2 + b*x +c, [a,b,c])
        if r: return solve_quadratic(r[a], r[b], r[c])
    raise "sorry"

def solve_linear(a, b):
    return -b/a

def solve_quadratic(a, b, c):
    return a+b
