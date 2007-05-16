import sys
sys.path.append(".")

from sympy import *
from decimal import Decimal

def test_cos():
    assert cos(0) == 1
    assert cos(pi) == -1
    assert cos(-pi) == -1
    assert cos(2*pi) == 1
    assert cos(3*pi) == -1
    assert cos(10**157*pi) == 1
    assert cos((10**157+1)*pi) == -1

    assert cos(pi/2) == 0
    assert cos((10**157+123)*pi/2) == 0

    assert cos(pi/3) == Rational(1, 2)
    assert cos(2*pi/3) == -Rational(1, 2)
    assert cos(-2*pi/3) == -Rational(1, 2)
    assert cos(4*pi/3) == -Rational(1, 2)
    assert cos(5*pi/3) == Rational(1, 2)
    assert cos(7*pi/3) == Rational(1, 2)
    assert cos(8*pi/3) == -Rational(1, 2)

    assert cos(pi/4) == Rational(1, 2)*sqrt(2)
    assert cos(-pi/4) == Rational(1, 2)*sqrt(2)
    assert cos(3*pi/4) == -Rational(1, 2)*sqrt(2)
    assert cos(5*pi/4) == -Rational(1, 2)*sqrt(2)
    assert cos(7*pi/4) == Rational(1, 2)*sqrt(2)

    assert cos(pi/6) == Rational(1, 2)*sqrt(3)
    assert cos(5*pi/6) == -Rational(1, 2)*sqrt(3)

    assert float(cos(1)) == 0.54030230586813977
    assert cos(1).evalf(precision=28) == Decimal("0.540302305868171952183697906453")
    assert float(cos(1) + cos(2)) == 0.12415546932099736
    assert abs(float(cos(1)*cos(2)*cos(3)) - 0.22259495730990297) < 1e-15 

    x = Symbol('x')
    y = Symbol('y')
    z = Symbol('z')
    assert cos(-x).expand() == cos(x)
    assert cos(-x-y).expand() == cos(x)*cos(y) - sin(x)*sin(y)
    assert cos(x+y-z).expand() == cos(x)*cos(y)*cos(z) + cos(x)*sin(y)*sin(z) - sin(x)*sin(y)*cos(z) + sin(x)*cos(y)*sin(z)

    assert cos(2*x).expand() == 2*cos(x)**2 - 1
    assert cos(-3*x).expand() == 4*cos(x)**3 - 3*cos(x)

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

    x = Symbol('x')
    y = Symbol('y')
    z = Symbol('z')
    assert sin(-x).expand() == -sin(x)
    assert sin(-x-y).expand() == -sin(x)*cos(y) - cos(x)*sin(y)
    assert sin(x+y+z).expand() == sin(x)*cos(y)*cos(z) - sin(x)*sin(y)*sin(z) + cos(x)*sin(y)*cos(z) + cos(x)*cos(y)*sin(z)

    assert sin(2*x).expand() == 2*sin(x)*cos(x)
    assert sin(-3*x).expand() == -4*sin(x)*cos(x)**2 + sin(x)

def test_tan():
    x = Symbol('x')
    y = Symbol('y')
    z = Symbol('z')
    assert tan(-x).expand() == -tan(x)
    assert tan(-x-y).expand() == (-tan(x) - tan(y)) / (1 - tan(x)*tan(y))
    assert tan(x+y+z).expand() == (tan(x) + tan(y) + tan(z) - tan(x)*tan(y)*tan(z)) / (1 - tan(x)*tan(y) - tan(x)*tan(z) - tan(y)*tan(z)) 

    assert tan(2*x).expand() == 2*tan(x) / (1 - tan(x)**2)
    assert tan(-3*x).expand() == -(3*tan(x) - tan(x)**3) / (1 - 3*tan(x)**2)