"""Module with some routines for polynomials"""

from sympy.core import Pow, Add, Mul, Rational, Number, Symbol

class PolynomialException(Exception):
    pass

def ispoly(p,x):
    """Is 'p' a polynomial in 'x'? Return True or False"""
    try:
        #basically, the polynomial is whatever we can convert using
        #the get_poly(). See it's docstring for more info.
        get_poly(p,x)
    except PolynomialException:
        return False
    return True

def fact(n):
    """Returns n!"""
    if n == 0: return 1
    else: return fact(n-1)*n

def coeff(poly, x, n):
    """Returns the coefficient of x**n in the polynomial"""
    assert ispoly(poly,x)
    return poly.diffn(x,n).subs(x,0)/fact(n)

def get_poly(p, x):
    """Returns a python list representing the polynomial 'p(x)'.
    
    'p' is a polynomial, for example: x**2 + 3*x*y.sqrt() - 8
    'x' is the variable of the polynomial, for example: x
    get_poly returns a python list of the form [(coeff0,n0), (coeff1,n1), ...]
    where p = coeff0*x**n0 + coeff1*x**n1 + ...
    and n0, n1, ... are sorted from the lower exponents up.

    Example:
    >>> from sympy import *
    >>> from sympy.modules.polynomials import get_poly
    >>> x = Symbol("x")
    >>> y = Symbol("y")
    >>> get_poly(x**2 + 3*x**7*y.sqrt() - 8, x)
    [(-8, 0), (1, 2), (3*y**(1/2), 7)]

    
    """
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
            return [(a*c,n)]
        a, b = get_poly(a,x), get_poly(b,x)
        assert len(a) == 1
        assert len(b) == 1
        a, b = a[0], b[0]
        r = (a[0]*b[0], a[1]+b[1])
        return [r]
    raise PolynomialException("p is not a polynomial")
