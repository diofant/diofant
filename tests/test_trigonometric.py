import sys
sys.path.append(".")

from sympy import *
from decimal import Decimal

def test_cos():
    assert cos(pi) == -1
    assert float(cos(1)) == 0.54030230586813977
    assert cos(1).evalf() == Decimal("0.540302305868171952183697906453")
    assert float(cos(1) + cos(2)) == 0.12415546932099736
    assert float(cos(1)*cos(2)*cos(3)) == 0.22259495730990297

def test_sin():
    assert sin(pi) == 0
    assert float(sin(1)) == 0.8414709848078965
