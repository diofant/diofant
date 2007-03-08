import sys
sys.path.append(".")

import sympy as s
import py
import decimal

x = s.Symbol('x')

def eq(a,b):
    return abs(a-b)<decimal.Decimal("0.0001")

def test_evalf():
    e = s.log(3)/s.log(2)-1
    assert eq(e.evalf(),0.58496)
    f = 2*x+2
    py.test.raises(ValueError,f.evalf)
    e = (s.Rational(2).sqrt()+1)/3
    assert s.isnumber(e)     

def test_bug1():
    x=s.Symbol('x')
    y=x*x
    assert eq(y.subs(x,s.Real(3.0)).evalf(),9)

def test_bug2():
    a = s.Real(4.)
    x = a + a
    x.eval()

