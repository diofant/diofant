import sys
sys.path.append(".")

import sympy as g
from sympy import Rational, Symbol, Real, sqrt
import py

def test_Rational():
    n1=g.Rational(1,4)
    n2=g.Rational(1,3)
    n3=g.Rational(2,4)
    n4=g.Rational(2,-4)
    n5=g.Rational(0)
    n6=g.Rational(1)
    n7=g.Rational(3)
    n8=g.Rational(-3)
    assert str(n1*n2) == "1/12"
    assert str(n1*n2) == "1/12"
    assert str(n3) == "1/2"
    assert str(n1*n3) == "1/8"
    assert str(n1+n3) == "3/4"
    assert str(n1+n2) == "7/12"
    assert str(n1+n4) == "-1/4"
    assert str(n4*n4) == "1/4"
    assert str(n4+n2) == "-1/6"
    assert str(n4+n5) == "-1/2"
    assert str(n4*n5) == "0"
    assert str(n3+n4) == "0"
    assert str(n1**n7) == "1/64"
    assert str(n2**n7) == "1/27"
    assert str(n2**n8) == "27"
    assert str(n7**n8) == "1/27"

def test_Rational_cmp():
    n1=g.Rational(1,4)
    n2=g.Rational(1,3)
    n3=g.Rational(2,4)
    n4=g.Rational(2,-4)
    n5=g.Rational(0)
    n6=g.Rational(1)
    n7=g.Rational(3)
    n8=g.Rational(-3)

    assert n8<n5
    assert n5<n6
    assert n6<n7
    assert n8<n7
    assert n7>n8
    assert (n1+1)**n2 < 2
    assert ((n1+n6)/n7) < 1

    assert n4<n3
    assert n2<n3
    assert n1<n2
    assert n3>n1
    assert not n3<n1
    assert not (Rational(-1) > 0)
    assert Rational(-1) < 0
    
def test_Real():
    a = g.Real(2) ** g.Real(3)
    assert a.evalf() == 8.0
    assert abs((g.pi ** -1).evalf() - 0.318309886184) < 0.0000001
    a = g.Real(2) ** g.Real(4)
    assert a.evalf() == g.Real(16.0)


def test_Real_eval():
    a = Real(3.2)
    assert isinstance(a**2, Real)

def test_Infinity():
    assert g.oo == g.oo
    assert g.oo != 1
    assert 1*g.oo == g.oo
    assert 1 != g.oo
    assert g.oo != -g.oo
    assert g.oo != g.Symbol("x")**3
    assert g.oo + 1 == g.oo + 1
    py.test.raises(ArithmeticError, lambda x: x-x, g.oo)
    

def test_powers():
    assert 64**(g.Rational(1)/3)==4
    assert 64**(g.Rational(2)/3)==16
    assert 24*64**(-g.Rational(1)/2)==3

def test_accept_int():
    assert g.Real(4) == 4

def test_accept_str():
    assert g.Real("0.2") == 0.2

def test_complex():
    x=Symbol("x", is_real=True)
    y=Symbol("y", is_real=True)
    a = g.Symbol("a")
    b = g.Symbol("b")
    e=(a+g.I*b)*(a-g.I*b)
    assert e.expand() == a**2+b**2
    assert e.expand() != a**2-b**2

    assert (a+g.I*b).conjugate() !=  a+g.I*b
    assert (a+g.I*b).conjugate() ==  a-g.I*b

    assert str(abs(a))=="abs(a)"

def test_abs1():
    a=Symbol("a", is_real=True)
    b=Symbol("b", is_real=True)
    assert abs(a) == a
    assert abs(-a) == a
    assert abs(-a) != -a
    assert abs(a+g.I*b) == sqrt(a*a+b*b)

def test_abs2():
    a=Symbol("a", is_real=False)
    b=Symbol("b", is_real=False)
    assert abs(a) != a
    assert abs(-a) != a
    assert abs(a+g.I*b) != sqrt(a*a+b*b)

def test_int():
    a=Rational(5)
    assert int(a)==5

def test_real_bug():
    x=g.Symbol("x")
    assert str(2.0*x*x) in ["(2.000*x)*x","2.000*x**2"]
    assert str(2.1*x*x)!="(2.0*x)*x"

def test_bug_sqrt():
    assert ((sqrt(Rational(2))+1)*(sqrt(Rational(2))-1)).expand() == 1

