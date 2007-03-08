import sys
sys.path.append(".")

import sympy as g
from sympy import log,exp

def testeval():
    a=g.Symbol("a")
    b=g.Symbol("b")
    e=a+b+a+b
    s1=str(e)
    e.eval()
    s2=str(e)
    assert s1 == s2

    c=g.Rational(1)
    p=g.Rational(5)
    e=a*b+c+p
    assert e == a*b+6
    e=c+a+p
    assert e == a+6
    e=c+a-p
    assert e == a+(-4)
    e=c+a+b*c+a-p
    assert e == 2*a+b+(-4)
    e=a+a
    assert e == 2*a
    e=a+p+a
    assert e == 2*a+5
    e=a*g.Rational(2)+p+a
    assert e == a*2+5+a
    e=a*g.Rational(2)+p+a
    assert e == 3*a+5
    e=c+p
    assert e == g.Rational(6)
    e=a*g.Rational(2)+a
    assert e == 3*a
    e=b+a-b
    assert e == a

def test_evalpow_bug():
    x=g.Symbol("x")
    e=1/(1/x)
    assert e==x
    e=1/(-1/x)
    assert e==-x


def test_expbug():
    assert exp(-log(3))**(-1) == 3
