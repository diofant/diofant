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
    def wo(di, x):
        """Are all items in the dictionary "di" without "x"?"""
        for d in di:
            if di[d].has(x):
                return False
        return True

    #currently solve only for one variable
    if len(vars) == 1:
        x = vars[0]
        a,b,c = [Symbol(s, dummy = True) for s in ["a","b","c"]]

        r = eq.match(a*x + b, [a,b])
        if r and wo(r,x): return solve_linear(r[a], r[b])

        r = eq.match(a*x**2 + c, [a,c])
        if r and wo(r,x): return solve_quadratic(r[a], 0, r[c])

        r = eq.match(a*x**2 + b*x + c, [a,b,c])
        if r and wo(r,x): return solve_quadratic(r[a], r[b], r[c])

    raise "sorry"

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
