import sys
sys.path.append(".")

import sym as g

x = g.symbol('x')
y = g.symbol('y')
z = g.symbol('z')
w = g.symbol('w')

def test_poly_str():
    assert str((2*x-(7*x**2 - 2) + 3*y)) == "2*x-(7*x^2-2)+3*y"
    assert str(x-y) == "x-y" 
    assert str(2+-x) == "2-x"
    assert str(x-2) == "x-2"
    assert str((x-y-z-w).eval()) == "x-z-y-w"
    assert str((x-z*y**2*z*w).eval()) == "x-y^2*w*z^2"

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
