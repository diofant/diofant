import sys
sys.path.append(".")

import sym as g

x = g.symbol('x')
y = g.symbol('y')

def test_poly_str():
    assert (2*x-(7*x**2 - 2) + 3*y).__str__() == "2*x-(7*x^2-2)+3*y"
    assert (x-y).__str__() == "x-y" 
    assert (2+-x).__str__() == "2-x"
    assert (x-2).__str__() == "x-2"

def test_bug1():
    e=(x-1*y*x*y)
    #this raised an exception (both lines are needed)
    a=str(e)
    a=str(e.eval())
