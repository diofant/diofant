import sys
sys.path.append(".")

import py

from sympy import *
from sympy.modules.specfun.orthogonal_polynomials import *

x = Symbol('x')

def test_legendre():
    assert legendre(0, x) == 1
    assert legendre(1, x) == x
    assert legendre(2, x) == simplify((3*x**2-1)/2)
    assert legendre(3, x) == simplify((5*x**3-3*x)/2)
    assert legendre(10, -1) == 1
    assert legendre(11, -1) == -1
    assert legendre(10, 1) == 1
    assert legendre(11, 1) == 1
    assert legendre(10, 0) != 0
    assert legendre(11, 0) == 0
    for n in range(1, 5):
        for k in range(n):
            z = legendre_zero(n, k)
            assert legendre(n, z) == 0
            assert abs(legendre(n, z.evalf())) < 1e-8
            assert abs(legendre(n+1, z.evalf())) > 1e-8

test_legendre()