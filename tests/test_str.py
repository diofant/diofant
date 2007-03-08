import sys
sys.path.append(".")

import sympy as g

x = g.Symbol('x')
y = g.Symbol('y')
z = g.Symbol('z')
w = g.Symbol('w')

def test_poly_str():
    #if any of these tests fails, it can still be correct, just the terms can
    #be in a different order. That happens for example when we change the 
    #hash algorithm. If it is correct, just add another item in the list [] of
    #correct results.
    #assert str((2*x-(7*x**2 - 2) + 3*y)) == "2*x-(7*x^2-2)+3*y"
    assert str(x-y) in ["x-y", "-y+x"]
    assert str(2+-x) == "2-x"
    #assert str(x-2) == "x-2"
    assert str((x-y-z-w)) in ["x-y-z-w","-w-y-z+x","x-w-y-z"]
    assert str((x-y-z-w).eval()) in ["-w-y-z+x","x-w-y-z","-w+x-z-y",
            "-y-w-z+x","-y+x-z-w","-y+x-w-z"]
    assert str((x-z*y**2*z*w).eval()) in ["-z^2*y^2*w+x", "x-w*y^2*z^2",
            "-y^2*z^2*w+x","x-w*z^2*y^2","x-y^2*z^2*w","x-y^2*w*z^2",
            "x-z^2*y^2*w"]

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
