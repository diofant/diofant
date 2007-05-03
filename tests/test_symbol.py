import sys
sys.path.append(".")

from sympy import Symbol
import py

def test_Symbol():
    a = Symbol("a")
    x1 = Symbol("x")
    x2 = Symbol("x")
    xdummy1 = Symbol("x", dummy=True)
    xdummy2 = Symbol("x", dummy=True)
    assert a != x1
    assert a != x2
    assert x1 == x2
    assert x1 != xdummy1
    assert xdummy1 != xdummy2

    assert Symbol("x") == Symbol("x")
    assert Symbol("x", dummy=True) != Symbol("x", dummy=True)
    
def test_Symbol_assumptions():
    x = Symbol('x')
    assert x.is_number == False
    assert (1+x).is_number == False
    assert (2*x).is_number == False
    
    assert x.is_dummy != True

def test_lt_gt():
    x = Symbol('x')
    y = Symbol('y')
    py.test.raises(NotImplementedError, "x<y")
    py.test.raises(NotImplementedError, "x>y")
    py.test.raises(NotImplementedError, "x>0")
    py.test.raises(NotImplementedError, "x<0")

    # let's check this also on classes Add, Mul, Pow
    py.test.raises(NotImplementedError, "x+1>0")
    py.test.raises(NotImplementedError, "2*x>0")
    py.test.raises(NotImplementedError, "x**2>0")
