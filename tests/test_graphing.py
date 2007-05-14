"""
# uncommenting this causes issue 112 found at:
# http://code.google.com/p/sympy/issues/detail?id=112&can=2&q=&sort=-id

import sys
sys.path.append(".")

from sympy import Symbol, Rational, sqrt, log, sin

import sympy.modules.graphing as g

x = Symbol('x')
y = Symbol('y')

def compare_arrays(a, b):
    try:
        if len(a) != len(b):
            return false
        for x in range(len(a)):
            if not compare_arrays(a[x], b[x]):
                return False
        return True
    except: # base case, a isn't enumerable
        return a == b

def test_compare_arrays():
    try:
        from numpy import array
    except:
        print "Numpy required to run graphing tests. Skipping..."
        return

    assert compare_arrays(1,1)
    assert not compare_arrays(1,2)
    assert compare_arrays([1,2],[1,2])
    assert not compare_arrays([1,2,3],[1,2])
    assert compare_arrays([[1,2],[1,2]], [array([1,2]),[1,2]])
    assert not compare_arrays([[1,2],[1,2]], [array([1,2]),[2,1]])
    assert not compare_arrays([[1,2],[1,2]], [array([1,2]),[2,[1]]])

def test_sample2d():
    try:
        from numpy import array
    except:
        return

    s1 = g.sample(x, (x, 0, 2, 2))
    a1 = [[0,1,2],[0,1,2]]
    assert compare_arrays(s1, a1)

    s2 = g.sample(x**2, (x, 1, 3, 2))
    a2 = [[1,2,3], [1,4,9]]
    assert compare_arrays(s2, a2)

def test_sample3d():
    try:
        from numpy import array
    except:
        return

    s1 = g.sample(x+y, (x, 0, 1, 1), (y, 0, 1, 1))
    a1 = ( array([[ 0.,  1.],[ 0.,  1.]]) ,
           array([[ 0.,  0.],[ 1.,  1.]]) ,
           array([[ 0.,  1.],[ 1.,  2.]]) )
    assert compare_arrays(s1, a1)
"""
