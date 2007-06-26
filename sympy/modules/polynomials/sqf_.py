"""Algorithms for square-free decomposition of polynomials"""

from sympy.modules.polynomials.base import *
from sympy.modules.polynomials import gcd_
from sympy.modules.polynomials import div_

def uv(f):
    """Returns a decomposition of f in a1 * a2**2 * ... * an**n.

    Here, the ai are pairwise prime and square-free polynomials, returned
    in a list. f is assumed to be a univariate instance of Polynomial.

    """
    f = [f]
    while f[-1].cl[0][1] != 0:
        f.append(gcd_.uv(f[-1], f[-1].diff(f[-1].var[0])))
    g = []
    for i in range(1, len(f)):
        g.append(div_.mv(f[i-1], f[i])[0][0])
    a = []
    for i in range(0, len(g)-1):
        a.append(div_.mv(g[i], g[i+1])[0][0])
    a.append(g[-1])
    return a

def uv_part(f):
    """Returns the square-free part of f.

    f is assumed to be a univariate instance of Polynomial.

    """
    ff = gcd_.uv(f, f.diff(f.var[0]))
    return div_.mv(f, ff)[0][0]
