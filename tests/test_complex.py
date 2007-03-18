import sys
sys.path.append(".")

from sympy import Rational, Symbol, I

def test_complex():
    a=Symbol("a")
    b=Symbol("b")
    e=(a+I*b)*(a-I*b)
    assert e.expand() == a**2+b**2
    assert e.expand() != a**2-b**2

    assert (a+I*b).conjugate() !=  a+I*b
    assert (a+I*b).conjugate() ==  a-I*b

    assert str(abs(a))=="abs(a)"

def test_abs1():
    a=Symbol("a", real=True)
    b=Symbol("b", real=True)
    assert abs(a) == a
    assert abs(-a) == a
    assert abs(-a) != -a
    assert abs(a+I*b) == (a*a+b*b).sqrt()

def test_abs2():
    a=Symbol("a", real=False)
    b=Symbol("b", real=False)
    assert abs(a) != a
    assert abs(-a) != a
    assert abs(a+I*b) != (a*a+b*b).sqrt()
