"""Module with some routines for polynomials"""

from sympy.core import Pow, Add, Mul, Rational, Number, Symbol, Basic

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
    p = Basic.sympify(p)
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

def poly(p, x):
    """Returns a sympy polynomial from the representation "p" returned by
    get_poly().
    """
    r = 0
    for t in p:
        r+=t[0]*x**t[1]
    return r

def rep(n, base):
    """Returns a representation of the integer 'n' in the base 'base'."""
    r = []
    while n!=0:
        r.append(n % base)
        n = (n - r[-1])/base
    return tuple(r)

def gcd(a, b, x):
    """Calculates a greatest common divisor of two polynomials.

    Currently using a heuristics algorithm.
    """

    def getcandidate(a, b, x, x0):
        n1 = a.subs(x, x0)
        n2 = b.subs(x, x0)
        n3 = n1.gcd(int(n1),int(n2))
        c = []
        for n, t in enumerate(rep(n3, x0)):
            if t != 0:
                c.append((t,n))
        return poly(c, x)

    c = getcandidate(a, b, x, 100)
    if div(a, c, x)[1] == 0 and div(b, c, x)[1] == 0: return c
    c = getcandidate(a, b, x, 101)
    if div(a, c, x)[1] == 0 and div(b, c, x)[1] == 0: return c

    raise PolynomialException("Can't calculate gcd for these polynomials")

def sqf(p, x):
    """Calculates the square free decomposition of 'p'.
    """
    g = gcd(p, p.diff(x), x)
    if g == 1: return p
    a, b = div(p, g, x)
    assert b == 0
    return sqf(a, x) * g

def div(f, g, x):
    """Expresses f = g*q + r, returns (q,r)

    All coefficients of 'f' and 'g' are assumed to be integers,
    and coefficients in 'q' and 'r' are then guaranteed to be integers.
    """
    fp = get_poly(f, x)
    gp = get_poly(g, x)
    q = 0
    while fp[-1][1] >= gp[-1][1] and fp[-1][0]!=0:
        s1 = poly([fp[-1]], x) / poly([gp[-1]], x)
        if isinstance(s1, Mul):
            a,b = s1.getab()
            if isinstance(a, Number) and not a.isinteger():
                #the coefficient is rational but not real, let's
                #put it in the remainder and we are done
                return q, f
        f = (f - g*s1).expand()
        fp = get_poly(f, x)
        q+=s1
    return q, f
