import sys
sys.path.append(".")

import sympy as s
import py

def test_Symbol():
    a = s.Symbol("a")
    x1 = s.Symbol("x")
    x2=s.Symbol("x")
    xdummy1=s.Symbol("x", is_dummy=True)
    xdummy2=s.Symbol("x", is_dummy=True)
    assert a != x1
    assert a != x2
    assert x1 == x2
    assert x1 != xdummy1
    assert xdummy1 != xdummy2

    assert s.Symbol("x")==s.Symbol("x")
    assert s.Symbol("x", is_dummy=True)!=s.Symbol("x", is_dummy=True)

def test_lt_gt():
    x = s.Symbol('x')
    y = s.Symbol('y')
    py.test.raises(NotImplementedError, "x<y")
    py.test.raises(NotImplementedError, "x>y")
    py.test.raises(NotImplementedError, "x>0")
    py.test.raises(NotImplementedError, "x<0")

    # let's check this also on classes Add, Mul, Pow
    py.test.raises(NotImplementedError, "x+1>0")
    py.test.raises(NotImplementedError, "2*x>0")
    py.test.raises(NotImplementedError, "x**2>0")
