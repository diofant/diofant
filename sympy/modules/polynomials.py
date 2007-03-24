"""Module with some routines for polynomials"""

from sympy.core import Pow, Add, Mul, Rational, Number, Symbol

class PolynomialException(Exception):
    pass

def ispoly(p,x):
    try:
        get_poly(p,x)
    except PolynomialException:
        return False
    return True

def fact(n):
    "Returns n!"
    if n == 0: return 1
    else: return fact(n-1)*n

def coeff(poly, x, n):
    """Returns the coefficient of x**n in the polynomial"""
    assert ispoly(poly,x)
    return poly.diffn(x,n).subs(x,0)/fact(n)

def get_poly(p, x):
    if not p.has(x):
        return [(p,0)]
    if p==x:
        return [(1,1)]
    if isinstance(p, Pow):
        if isinstance(p.exp, Rational) and p.exp.isinteger():
            n = int(p.exp)
            if n>0 and isinstance(p.base, Symbol):
                return [(1,n)]
    if isinstance(p,Add):
        a,b = p.getab()
        r = get_poly(a,x) + get_poly(b,x)
        r.sort(key = lambda x: x[1])
        return r
    if isinstance(p,Mul):
        a,b = p.getab()
        if isinstance(a, Number):
            c,n = get_poly(b,x)[0]
            return [(a,n)]
        a, b = get_poly(a,x), get_poly(b,x)
        assert len(a) == 1
        assert len(b) == 1
        a, b = a[0], b[0]
        r = (a[0]*b[0], a[1]+b[1])
        return [r]
    raise PolynomialException("p is not a polynomial")
