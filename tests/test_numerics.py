import sys
sys.path.append(".")

import py

from sympy import *
from sympy.modules.numerics import *

def test_transcendental():
    # All should be displayed correctly rounded to 19 decimal places
    assert str(pif())   == '3.141592653589793238'
    assert str(expf(1)) == '2.718281828459045235'
    assert str(cosf(1)) == '0.5403023058681397174'
    assert str(sinf(1)) == '0.8414709848078965067'
    assert str(logf(2)) == '0.6931471805599453094'
    sq2 = str(powerf(2, '0.5'))
    assert sq2          == '1.414213562373095049'
    pi_again = str(logf(-1).imag)
    assert pi_again     == '3.141592653589793238'
