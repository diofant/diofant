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
    assert sin(0) == 0
    assert sin(pi) == 0
    assert sin(-pi) == 0
    assert sin(2*pi) == 0
    assert sin(-3*10**73*pi) == 0
    assert sin(7*10**103*pi) == 0

    assert sin(pi/2) == 1
    assert sin(-pi/2) == -1
    assert sin(5*pi/2) == 1
    assert sin(7*pi/2) == -1

    assert sin(pi/3) == Rational(1, 2)*sqrt(3)
    assert sin(-2*pi/3) == -Rational(1, 2)*sqrt(3)

    assert sin(pi/4) == Rational(1, 2)*sqrt(2)
    assert sin(-pi/4) == -Rational(1, 2)*sqrt(2)
    assert sin(17*pi/4) == Rational(1, 2)*sqrt(2)
    assert sin(-3*pi/4) == -Rational(1, 2)*sqrt(2)

    assert sin(pi/6) == Rational(1, 2)
    assert sin(-pi/6) == -Rational(1, 2)
    assert sin(7*pi/6) == -Rational(1, 2)
    assert sin(-5*pi/6) == -Rational(1, 2)

    assert float(sin(1)) == 0.8414709848078965
