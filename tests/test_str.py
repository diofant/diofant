import sys
sys.path.append(".")

import sym as g

x = g.Symbol('x')
y = g.Symbol('y')
z = g.Symbol('z')
w = g.Symbol('w')

def test_poly_str():
    assert str((2*x-(7*x**2 - 2) + 3*y)) == "2*x-(7*x^2-2)+3*y"
    assert str(x-y) == "x-y" 
    assert str(2+-x) == "2-x"
    assert str(x-2) == "x-2"
    assert str((x-y-z-w)) == "x-y-z-w"
    assert str((x-y-z-w).eval()) == "x-w-y-z"
    assert str((x-z*y**2*z*w).eval()) == "-z^2*y^2*w+x"

def test_bug1():
    e=(x-1*y*x*y)
    #this raised an exception (both lines are needed)
    a=str(e)
    a=str(e.eval())

def test_bug2():
    e=x-y
    a=str(e)
    b=str(e)
    assert a==b
