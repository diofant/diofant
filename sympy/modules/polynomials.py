"""Module with some routines for polynomials"""

from sympy.core import Pow, Add, Mul, Rational, Number, Symbol, Basic
from sympy.core.functions import diff

class PolynomialException(Exception):
    pass

def ispoly(p, var=None):
    """
    Usage
    =====
      ispoly(p, var) -> Returns True if p is a polynomial in variable var. 
                        Returns False otherwise.
        
    Notes
    =====
        You can check wether it's a polynomial in several variables at once giving a 
        tuple of symbols second argument (like ispoly(x**2 + y + 1, (x,y)) ).See
        examples for more info.
    
    Examples
    ========
        >>> from sympy import *
        >>> from sympy.modules.polynomials import *
        >>> x = Symbol('x')
        >>> ispoly(x**2+x+1, x)
        True
        >>> y = Symbol('y')
        >>> ispoly(x**2 + y + 1, (x,y)) #polynomial in variables x and y
        True
        >>> ispoly(x**2 + exp(y) + 1, (x,y))
        False
        
    See also
    ========
       L{get_poly}, L{coeff}
       
    """
    if isinstance(var, Basic):
        # if the var argument is not a tuple or list
        var = [var] # so we can iterate over it
    try:
        #basically, the polynomial is whatever we can convert using
        #the get_poly(). See it's docstring for more info.
        if var is None:
            var = p.atoms(type=Symbol)[0] # make it work even if the user doesen't issue a variable
            print "\t*** Warning. You have not issued a variable as argument."
            print "\t*** Please see the interactive help on this function for more info"
            print "\t*** Using %s as variable" % str(var)
        for v in var:
            get_poly(p, v)
    except PolynomialException:
        return False
    except IndexError:
        # if p.atoms() is empty
        raise TypeError("Wrong arguments")
    return True

def fact(n):
    """Returns n!"""
    if n == 0: return 1
    else: return fact(n-1)*n

def coeff(poly, x, n):
    """Returns the coefficient of x**n in the polynomial"""
    assert ispoly(poly,x)
    return diff(poly, x,n).subs(x,0)/fact(n)

def get_poly(p, x):
    """Returns a python list representing the polynomial 'p(x)'.
    
    'p' is a polynomial, for example: x**2 + 3*x*sqrt(y) - 8
    'x' is the variable of the polynomial, for example: x
    get_poly returns a python list of the form [(coeff0,n0), (coeff1,n1), ...]
    where p = coeff0*x**n0 + coeff1*x**n1 + ...
    and n0, n1, ... are sorted from the lower exponents up.

    Example:
    >>> from sympy import *
    >>> from sympy.modules.polynomials import get_poly
    >>> x = Symbol("x")
    >>> y = Symbol("y")
    >>> get_poly(x**2 + 3*x**7*sqrt(y) - 8, x)
    [(-8, 0), (1, 2), (3*y**(1/2), 7)]

    
    """
    p = Basic.sympify(p)
    if not p.has(x):
        return [(p,0)]
    if p==x:
        return [(1,1)]
    if isinstance(p, Pow):
        if isinstance(p.exp, Rational) and p.exp.is_integer:
            n = int(p.exp)
            if n>0:
                if isinstance(p.base, Symbol):
                    return [(1,n)]
                else:
                    # FIXME The return value isn't correct, but at least it
                    #       doesn't break is_poly
                    get_poly(p.base, x)
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

    #try some values of x0. If you find polynomials for which gcd doesn't
    #work, just find a number of x0, that works and add it to the end
    #of this list:
    for x0 in [100, 101]:
        c = getcandidate(a, b, x, x0)
        if div(a, c, x)[1] == 0 and div(b, c, x)[1] == 0: 
            return c

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
            if isinstance(a, Number) and not a.is_integer:
                #the coefficient is rational but not integer, let's
                #put it in the remainder and we are done
                return q, f
        f = (f - g*s1).expand()
        fp = get_poly(f, x)
        q+=s1
    return q, f
    
def collect(expr, syms):
    """Collect additive terms with respect to a list of variables in a linear
       multivariate polynomial. This function assumes the input expression is
       in an expanded form and will return None if this is not a linear
       polynomial or else a pair of the following form:

          ( { variable : coefficient }, free term )

       Example:
       >>> from sympy import *
       >>> from sympy.modules.polynomials import collect
       >>> x, y, z = Symbol('x'), Symbol('y'), Symbol('z')
       >>> collect(2*x + sin(z)*x + cos(z)*y + 1, [x, y])
       ({x: 2+sin(z), y: cos(z)}, 1)       

    """

    if isinstance(expr, (Add, Mul)):
        content, tail = {}, 0

        if isinstance(expr, Mul):
            expr = [expr]

        for term in expr:
            coeff = 1

            while isinstance(term, Mul):
                a, term = term.getab()

                if isinstance(a, Symbol) and a in syms:
                    if (term.has_any(syms)):
                        return None
                    else:
                        coeff *= term

                        if a in content:
                            content[a] += coeff
                        else:
                            content[a] = coeff

                        break
                else:
                    if (a.has_any(syms)):
                        return None
                    else:
                        coeff *= a
            else:
                if isinstance(term, Symbol) and term in syms:
                    if term in content:
                        content[term] += coeff
                    else:
                        content[term] = coeff
                else:
                    tail += coeff * term

        return (content, tail)
    elif isinstance(expr, Symbol) and expr in syms:
        return ({expr : 1}, 0)
    elif isinstance(expr, Basic) and not expr.has_any(syms):
        return ({}, expr)
    else:
        return None

