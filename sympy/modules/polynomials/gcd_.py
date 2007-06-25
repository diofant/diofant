"""Greatest common divisor for the Polynomial class"""

from sympy.modules.polynomials.base import *

def uv(f, g):
    """Euclidean algorithm for univariate polynomials.

    Coefficients assumed to be in a field.
    """
    from sympy.modules.polynomials import div_

    while True:
        if g.cl[0][0] != 0:
            g.cl = map(lambda t: [t[0]/g.cl[0][0]] + t[1:], g.cl)
            f, g = g, div_.mv(f, g)[-1]
        else:
            break
    return f

def mv(f, g):
    """Computes the gcd of 2 polynomials by dividing their product by their lcm.

    It is assumed that f and g are instances of the Polynomial class with
    matching variables and orders.
    """
    from sympy.modules.polynomials import lcm_
    from sympy.modules.polynomials import div_

    lcm = lcm_.mv(f, g)
    q, r = div_.mv(f*g, lcm)
    assert r == Rational(0)
    q = q[0] # q is a list!
    q.cl = map(lambda t:[t[0]/q.cl[0][0]] + t[1:], q.cl)
    return q
