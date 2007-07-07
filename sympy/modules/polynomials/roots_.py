"""Algorithms to determine the roots of polynomials"""

from sympy import exp, I, pi, sqrt 

from sympy.modules.polynomials.base import *
from sympy.modules.polynomials import div_

def cubic(f):
    """Returns the real or complex roots of a cubic polynomial.

    Works for univariate instances of Polynomial only."""
    assert f.cl[0][1] == 3

    # Get monic polynomial, for p-q formula
    f.cl = map(lambda t:[t[0]/f.cl[0][0]] + t[1:], f.cl)

    a = f.nth_coeff(2)
    b = f.nth_coeff(1)
    c = f.nth_coeff(0)

    # Substitute variable to get depressed cubic: t**3 + p * t + q
    p = b + -a**2/3
    q = c + Rational(2, 27)*a**3 - a*b/9

    if p == 0: # Solve special cases:
        if q == 0:
            return [-a/3]
        else:
            u1 = q**Rational(1, 3)
    else: 
        u1 = (q/2 + sqrt(q**2/4 + p**3/27))**Rational(1, 3)

    u2 = u1*(Rational(-1, 2) + I*sqrt(3)/2)
    u3 = u1*(Rational(-1, 2) - I*sqrt(3)/2)
    result = map(lambda u: (p/(u*3) - u - a/3).expand(), [u1, u2, u3])

    if f.coeff in ['cplx', 'sym']:
        return result
    else:
        return filter(lambda r: r.is_real , result)

def n_poly(f):
    """Checks if the polynomial can be simplifed by substituting the
    variable by a power.

    Returns a list of real or complex roots, or None, if no solution
    was found. Works for univariate instances of Polynomial only.
    """
    def roots_of_unity(n):
        result = []
        for i in range(0,n):
            result.append((exp(2*i*pi*I/n)).evalc())
        return result

    exponents = map(lambda t:int(t[1]), f.cl)
    g = reduce(Rational(0).gcd, exponents)
    if g == 1:
        return None
    n = int(f.cl[0][1]/g)
    if not n in [1, 2, 3]: # Cases where solution can be computed
        return None

    ff = f.copy()
    ff.cl = map(lambda t:[t[0], t[1]/g], ff.cl)
    # TODO: use uv() here?
    if n == 1:
        sol = [-(ff.cl[1][0]/ff.cl[0][0])]
    if n == 2:
        sol = quadratic(ff)
    elif n == 3:
        sol = cubic(ff)
    return [(zeta*s**Rational(1,g)).expand()
            for s in sol for zeta in roots_of_unity(g)]
        
def quadratic(f):
    """Returns the real or complex roots of a quadratic polynomial.

    Works for univariate instances of Polynomial only."""

    assert f.cl[0][1] == 2

    # Get monic polynomial, for p-q formula
    f.cl = map(lambda t:[t[0]/f.cl[0][0]] + t[1:], f.cl)

    # Solve special cases:
    if len(f.cl) == 1:
        return [Rational(0)]
    if len(f.cl) == 2:
        if f.cl[1][1] == 1: # No constant term
            return [Rational(0), -(f.cl[1][0])]
        else: # No linear term
            q = -(f.cl[1][0])
            if q > 0:
                return [-sqrt(q), sqrt(q)]
            else:
                if f.coeff in ['cplx', 'sym']:
                    return [-sqrt(q), sqrt(q)]
                else:
                    return []

    p = f.cl[1][0]
    q = f.cl[2][0]
    discr = p**2 - 4*q
    if (not discr.is_number) or discr > 0:
        return [-p/2 + sqrt(discr)/2,
                -p/2 - sqrt(discr)/2]
    elif discr == 0:
        return [-p/2]
    else: # discr < 0
        if f.coeff in ['cplx', 'sym']:
            return [-p/2 + I*sqrt(-discr)/2,
                    -p/2 - I*sqrt(-discr)/2]
        else:
            return []

# TODO: Implement function to find roots of quartic polynomials?

def rat_roots(f):
    """Returns a list of rational roots of an integer Polynomial.

    For an polynomial an*x**n + ... + a0, all rational roots are of
    the form p/q, where p and q are integer factors of a0 and an.
    """
    assert isinstance(f, Polynomial) and f.coeff == 'int'
    an_divs = integer_divisors(int(f.cl[0][0]))
    a0_divs = integer_divisors(int(f.cl[-1][0]))
    result = []
    for p in a0_divs:
        for q in an_divs:
            if f(Rational(p, q)) == 0:
                result.append(Rational(p, q))
            if f(Rational(-p, q)) == 0:
                result.append(Rational(-p, q))
    return result

def real_roots(s, a=None, b=None):
    """Returns the number of unique real roots of f in the interval (a, b].

    Assumes a sturm sequence of an univariate, square-free instance of
    Polynomial with real coefficients. The boundaries a and b can be omitted
    to check the whole real line.
    """
    def sign_changes(list):
        counter = 0
        current = list[0]
        for el in list:
            if (current < 0 and el >= 0) or \
               (current > 0 and el <= 0):
                counter += 1
            current = el
        return counter
    
    if a == None: # a = -oo
        sa = sign_changes(map(lambda p:p.cl[0][0]*(-1)**p.cl[0][1], s))
    else:
        sa = sign_changes(map(lambda p:p.basic.subs(p.var[0], a), s))        
    if b == None: # b = oo
        sb = sign_changes(map(lambda p:p.cl[0][0], s))
    else:
        sb = sign_changes(map(lambda p:p.basic.subs(p.var[0], b), s))
    return sa - sb
    
def sturm(f):
    """Compute the Sturm sequence of given polynomial

    Assumes an univariate instance of Polynomial with real coefficients.
    """
    seq = [f]
    seq.append(Polynomial(f.basic.diff(f.var[0]), f.var, f.order, f.coeff))
    while seq[-1].cl[0][0] != Rational(0):
        seq.append(-(div_.mv(seq[-2], seq[-1])[-1]))
    return seq[:-1]

def uv(f, verbose=False):
    """Compute the roots of an univariate polynomial.
    """
    n = f.cl[0][1] # degree
    if n == 0: # constant
        return []
    elif n == 1:
        return [-(f.cl[1][0] / f.cl[0][0])]
    elif n == 2:
        return quadratic(f)
    elif n == 3:
        return cubic(f)
    
    res = n_poly(f)
    if res != None:
        return res
    else:
        return []
