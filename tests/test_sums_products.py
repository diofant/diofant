import sys
sys.path.append(".")

import py

from sympy import *
from sympy.modules.sums_products import *

n = Symbol('n')
a = Symbol('a')
b = Symbol('b')

def test_str():
    assert str(Sum(cos(3*n), (n, a, b))) == "sum_{n=a}^{b} cos(3*n)"

def test_arithmetic_sums():
    assert Sum(1, (n, a, b)) == b-a+1
    assert Sum(1, (n, 1, 10)) == 10
    assert Sum(2*n, (n, 0, 10**10)) == 100000000010000000000

def test_polynomial_sums():
    assert Sum(n**2, (n, 3, 8)) == 199
    assert Sum(n, (n, a, b)) == \
        ((a+b)*(b-a+1)/2).expand()
    assert Sum(n**2, (n, 1, b)) == \
        ((2*b**3+3*b**2+b)/6).expand()
    assert Sum(n**3, (n, 1, b)) == \
        ((b**4+2*b**3+b**2)/4).expand()
    assert Sum(n**6, (n, 1, b)) == \
        ((6*b**7+21*b**6+21*b**5-7*b**3+b)/42).expand()

def test_geometric_sums():
    assert Sum(pi**n, (n, 0, b)) == (1-pi**(b+1)) / (1-pi)
    assert Sum(2 * 3**n, (n, 0, b)) == 3**(b+1) - 1
    assert Sum(Rational(1,2)**n, (n, 1, oo)) == 1
    assert Sum(2**n, (n, 0, b)) == 2**(b+1) - 1

def test_composite_sums():
    f = Rational(1,2)*(7 - 6*n + Rational(1,7)*n**3)
    s = Sum(f, (n, a, b))
    assert not isinstance(s, Sum)
    A = 0
    for i in range(-3, 5):
        A += f.subs(n, i)
    B = s.subs(a,-3).subs(b,4)
    assert A == B

def test_finite_sums():
    assert Sum(cos(n), (n, -2, 1)) == cos(-2)+cos(-1)+cos(0)+cos(1)
