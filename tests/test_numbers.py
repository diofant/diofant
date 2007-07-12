import sys
sys.path.append(".")

import sympy as g
from sympy import Rational, Symbol, Real, sqrt, oo
from sympy.core.power import integer_nthroot
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
    assert str(Rational("-25")) == "-25"
    assert str(Rational("25/7")) == "25/7"
    assert str(Rational("-123/569")) == "-123/569"
    assert str(Rational("1.25")) == "5/4"
    assert str(Rational("-2.6e-2")) == "-13/500"
    assert str(Rational("0.1[23]")) == "61/495"
    assert str(Rational("5.1[666]")) == "31/6"
    assert str(Rational("-5.1[666]")) == "-31/6"
    assert str(Rational("0.[9]")) == "1"
    assert str(Rational("-0.[9]")) == "-1"

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

def test_Rational_assumptions():
    a = Rational(0)
    b = Rational(1)
    c = Rational(-29)
    d = Rational(41, 37)
    e = Rational(79)

    assert a.is_zero == True
    assert a.is_one == False
    assert a.is_integer == True
    assert a.is_nonnegative == True
    assert a.is_nonpositive == True
    assert a.is_negative == False
    assert a.is_positive == False
    assert a.is_nonzero == False

    assert b.is_zero == False
    assert b.is_one == True
    assert b.is_integer == True
    assert b.is_nonnegative == True
    assert b.is_nonpositive == False
    assert b.is_negative == False
    assert b.is_positive == True
    assert b.is_nonzero == True

    assert c.is_negative_integer == True
    assert c.is_nonpositive_integer == True
    assert c.is_positive_integer == False
    assert c.is_nonnegative_integer == False

    assert d.is_integer == False
    assert d.is_nonzero == True

    assert a.is_even == True
    assert b.is_odd == True

    assert d.is_even == False
    assert d.is_odd == False

    assert e.is_prime == True
    assert d.is_prime == False

    assert a.is_bounded == True

def test_Real():
    def eq(a, b):
        t = Real("1.0E-15")
        return (-t < a-b < t)

    a = g.Real(2) ** g.Real(3)
    assert a.evalf() == g.Real(8)
    assert eq((g.pi ** -1).evalf(), Real("0.31830988618379067"))
    a = g.Real(2) ** g.Real(4)
    assert a.evalf() == g.Real(16)


def test_Real_eval():
    a = Real(3.2)
    assert isinstance(a**2, Real)

def test_Infinity():
    assert oo == oo
    assert oo != 1
    assert 1*oo == oo
    assert 1 != oo
    assert oo != -oo
    assert oo != Symbol("x")**3
    assert oo + 1 == oo + 1
    assert oo + 1 == oo
    assert 2 + oo == oo
    assert 3*oo + 2 == oo
    assert -oo*3 == -oo
    assert oo + oo == oo
    assert -oo + oo*(-5) == -oo
    assert 1/oo  == 0
    assert 1/(-oo)  == 0
    assert 8/oo  == 0
    py.test.raises(ArithmeticError, "oo-oo")
    py.test.raises(ArithmeticError, "oo*0")
    py.test.raises(ArithmeticError, "oo/oo")

def test_powers():
    assert integer_nthroot(1, 2) == (1, True)
    assert integer_nthroot(1, 5) == (1, True)
    assert integer_nthroot(2, 1) == (2, True)
    assert integer_nthroot(2, 2) == (1, False)
    assert integer_nthroot(2, 5) == (1, False)
    assert integer_nthroot(4, 2) == (2, True)
    assert integer_nthroot(123**25, 25) == (123, True)
    assert integer_nthroot(123**25+1, 25) == (123, False)
    assert integer_nthroot(123**25-1, 25) == (122, False)

    assert str(Rational(1,4) ** Rational(1,2)) == "1/2"
    assert str(Rational(1,36) ** Rational(1,2)) == "1/6"
    assert str(Rational(5**3, 8**3) ** Rational(4,3)) == \
        str(Rational(5**4, 8**4))
    assert str(Rational(-4,7) ** Rational(1,2)) == \
        str(g.I*Rational(4,7)**Rational(1,2))

    assert str((123**25) ** Rational(1,25)) == "123"
    assert str((123**25+1)**Rational(1,25)) != "123"
    assert str((123**25-1)**Rational(1,25)) != "123"
    assert str((123**25-1)**Rational(1,25)) != "122"
    assert str(Rational(3,5)**(-Rational(1,2))) == "(3/5)**(-1/2)"
    assert str(Rational(81,36)**(Rational(3,2))) == "27/8"
    assert str(Rational(81,36)**(-Rational(3,2))) == "8/27"

    assert 64**(g.Rational(1)/3)==4
    assert 64**(g.Rational(2)/3)==16
    assert 24*64**(-g.Rational(1)/2)==3

    assert str((-4)**Rational(1,2)) == str(2*g.I)

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
    assert str(abs(Rational(1,6))) == "1/6"
    assert str(abs(Rational(-1,6))) == "1/6"

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




