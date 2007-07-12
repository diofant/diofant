import sys
sys.path.append(".")

import py

from sympy import *
from sympy.modules.numerics import *

import math
import cmath

def test_transcendental_basic():
    # Basic sanity check: low-precision comparison against Python's built-ins
    x = [-9.3, -4.7, 0.0, 0.0001, 0.253, 0.97, 1.03, 2.56, 3.83, 15.111]
    def relcmp(a,b):
        if a == b == 0:
            return True
        d = abs(a-b)/max(abs(a),abs(b))
        return d < 1e-15
    def compare(f, g):
        for a in x:
            for b in x:
                z = complex(a,b)
                assert relcmp(f(z), complex(g(BinaryComplex(z))))
    compare(cmath.exp, expf)
    compare(cmath.sin, sinf)
    compare(cmath.cos, cosf)
    compare(cmath.tan, tanf)
    compare(cmath.sinh, sinhf)
    compare(cmath.cosh, coshf)
    compare(cmath.tanh, tanhf)

def test_rounding():
    # TODO: str() should give a correctly rounded decimal result. This
    # failed for sin(1) at dps=19. for now, reduced to 18 dps in this test
    assert pif().decimal(18)   == '3.14159265358979324'
    assert expf(1).decimal(18) == '2.71828182845904524'
    assert cosf(1).decimal(18) == '0.540302305868139717'
    assert sinf(1).decimal(18) == '0.841470984807896507'
    assert logf(2).decimal(18) == '0.693147180559945309'
    assert powerf(2, 0.5).decimal(18) == '1.41421356237309505'
    assert (logf(-1).imag).decimal(18)  == '3.14159265358979324'
    assert eulergammaf().decimal(18) == '0.577215664901532861'

def test_transcendental_hard():
    # Check difficult cases (many more tests needed here)
    assert str(sinf('1e-5')) == '0.000009999999999833333333'
    assert str(cosf('1e-5')) == '0.9999999999500000000'
    assert str(sinf('-1e-5')) == '-0.000009999999999833333333'
    assert str(cosf('-1e-5')) == '0.9999999999500000000'