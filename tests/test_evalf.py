import sys
sys.path.append(".")

from sympy import *
import py
import decimal

x = Symbol('x')

def eq(a,b):
    return abs(a-b)<decimal.Decimal("0.00000000000000001")

def test_evalf():
    e = log(3)/log(2)-1
    assert eq(e.evalf(),decimal.Decimal("0.5849625007250221337163"))
    f = 2*x+2
    py.test.raises(ValueError,f.evalf)
    e = sqrt(Rational(2)+1)/3
    assert e.is_number     

def test_bug1():
    x = Symbol('x')
    y=x*x
    assert eq(y.subs(x,Real(3.0)).evalf(),9)

def test_bug2():
    a = Real(4.)
    x = a + a

