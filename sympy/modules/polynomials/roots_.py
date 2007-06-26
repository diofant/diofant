"""Algorithms to determine the roots of polynomials"""

#from sympy.core.numbers import sign

from sympy.modules.polynomials.base import *
from sympy.modules.polynomials import div_

def sturm(f):
    """Compute the Sturm sequence of given polynomial

    Assumes an univariate instance of Polynomial with real coefficients.
    """
    seq = [f]
    seq.append(Polynomial(f.basic.diff(f.var[0]), f.var, f.order, f.coeff))
    while seq[-1].cl[0][0] != Rational(0):
        seq.append(-(div_.mv(seq[-2], seq[-1])[-1]))
    return seq[:-1]

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
    
