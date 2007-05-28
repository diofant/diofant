
import sys
sys.path.append(".")

import py

from sympy import *
from sympy.modules.polynomials import *

def test_ispoly():
    x = Symbol("x")
    y = Symbol("y")
    assert not ispoly( sqrt(x), x )
    assert ispoly( Rational(2), x)
    assert ispoly(x**2+1, x)
    assert ispoly( x, x)
    assert ispoly( x**2, x)
    assert ispoly( x**2 + 3*x - 8, x)
    assert ispoly( x**2 + 3*x*sqrt(y) - 8, x)
    assert not ispoly( x**2 + 3*x*sqrt(y) - 8 , y)
    assert ispoly((x**2)*(y**2) + x*(y**2) + y*x + x + exp(2), (x,y) )

    #assert Rational(1).ispoly(sin(x))
    #assert not exp(x).ispoly(sin(x))

def test_coeff():
    x = Symbol("x")
    assert coeff(x**2, x, 1) == 0
    assert coeff(x**2, x, 2) == 1
    assert coeff(x**2, x, 2) != 0

    assert coeff(2*x+18*x**8, x, 1) == 2
    assert coeff(2*x+18*x**8, x, 4) == 0
    assert coeff(2*x+18*x**8, x, 8) == 18

def test_poly():
    x = Symbol("x")
    y = Symbol("y")
    assert 3*x**2 == poly([(3,2)],x)
    assert 2*x+3*x**2 - 5 == poly([(-5, 0), (2, 1), (3,2)],x)
    assert 2*x**100+3*x**2 - 5 == poly([(-5, 0), (3,2), (2, 100)],x)
    assert 2*x**100+3*x**2 - 6 != poly([(-5, 0), (3,2), (2, 100)],x)

    assert sqrt(y)*x == poly([(sqrt(y),1)],x)
    assert x**2 + 3*x*sqrt(y) - 8 == poly([(-8, 0), (3*sqrt(y), 1),
        (1, 2)],x)

def test_gcd():
    x = Symbol("x")
    assert gcd(x**2, x, x) == x
    assert gcd(3*x**2, x, x) == x
    assert gcd(3*x**2, 3*x, x) == 3*x
    assert gcd(3*x**2, 6*x, x) == 3*x
    assert gcd(x**2+2*x+1, x+1, x) == x+1
    assert gcd(x**2+2*x+2, x+1, x) == 1

    assert gcd(x**2+2*x+1, 2+2*x, x) == 1+x
    assert gcd(x**2+2*x+2, 2+2*x, x) == 1

def test_rep():
    assert rep(101,100) == (1,1)
    assert rep(300,100) == (0,3)
    assert rep(100,100) == (0,1)

    assert rep(100,10) == (0,0,1)

def test_sqf():
    x = Symbol("x")
    assert sqf(3*x**2, x) == 3*x**2
    assert sqf(x**2+2*x+1, x) == (x+1)**2

def test_div():
    x = Symbol("x")
    assert div(x**3-12*x**2-42, x-3, x) == (x**2-9*x-27, -123)
    assert div(x**3-12*x**2-42, x**2+x-3, x) == (x-13, 16*x-81)

    assert div(2+2*x+x**2, 1, x) == (2+2*x+x**2, 0)
    assert div(2+2*x+x**2, 2, x) == (0, 2+2*x+x**2)

    assert div(3*x**3, x**2, x) == (3*x, 0)

def test_resultant():
    x, a, b, c, = [Symbol(y) for y in ['x', 'a', 'b', 'c']]

    s_res = resultant(x**2-1, x**3-x**2+2, x, method='sylvester').expand()
    b_res = resultant(x**2-1, x**3-x**2+2, x, method='bezout').expand()

    assert b_res == s_res == 0

    s_res = resultant(3*x**3-x, 5*x**2+1, x, method='sylvester').expand()
    b_res = resultant(3*x**3-x, 5*x**2+1, x, method='bezout').expand()

    assert b_res == s_res == 64

    s_res = resultant(x**2-2*x+7, x**3-x+5, x, method='sylvester').expand()
    b_res = resultant(x**2-2*x+7, x**3-x+5, x, method='bezout').expand()

    assert b_res == s_res == 265

    s_res = resultant((x-a)**2-2, a**2-3, a, method='sylvester').expand()
    b_res = resultant((x-a)**2-2, a**2-3, a, method='bezout').expand()

    assert b_res == s_res == 1 - 10*x**2 + x**4

    s_res = resultant((x-1)*(x-2)*(x-3), (x-4)*(x-5)*(x-6), x, method='sylvester').expand()
    b_res = resultant((x-1)*(x-2)*(x-3), (x-4)*(x-5)*(x-6), x, method='bezout').expand()

    assert b_res == s_res == -8640

    s_res = resultant((x-1)*(x-2)*(x-3), (x-4)*(x-5)*(x-1), x, method='sylvester').expand()
    b_res = resultant((x-1)*(x-2)*(x-3), (x-4)*(x-5)*(x-1), x, method='bezout').expand()

    assert b_res == s_res == 0

    s_res = resultant(x**3-1, x**3+2*x**2+2*x-1, x, method='sylvester').expand()
    b_res = resultant(x**3-1, x**3+2*x**2+2*x-1, x, method='bezout').expand()

    assert b_res == s_res == 16

    s_res = resultant(3*x**2+2*a*x+3*a**2-2, 3*x**2-2*a*x+3*a**2-2, x, method='sylvester').expand()
    b_res = resultant(3*x**2+2*a*x+3*a**2-2, 3*x**2-2*a*x+3*a**2-2, x, method='bezout').expand()

    assert b_res == s_res == 144*a**4 - 96*a**2

    s_res = resultant((x-a)*(x-b), x-c, x, method='sylvester').expand()
    b_res = resultant((x-a)*(x-b), x-c, x, method='bezout').expand()

    assert b_res == s_res == ((a-c)*(b-c)).expand()

def test_collect():
    x = Symbol('x')
    y = Symbol('y')
    z = Symbol('z')
    a = Symbol('a')

    from sympy.modules.trigonometric import sin, cos

    assert collect(x, [x, y, z]) == ({x: 1}, 0)
    assert collect(x-1, [x, y, z]) == ({x: 1}, -1)
    assert collect(x+y+z, [x, y, z]) == ({x: 1, y: 1, z: 1}, 0)
    assert collect(sin(a)*x-2*cos(a)*y+1024*z-a, [x, y]) \
            == ({x: sin(a), y: -2*cos(a)}, 1024*z-a)
    assert collect(2*x + sin(z)*x + cos(a)*y + z + cos(a) + cos(a)*x + 1, [x, y]) \
            == ({x: 2+sin(z)+cos(a), y: cos(a)}, z+cos(a)+1)
    assert collect(x*y, [x, y]) == None
    assert collect(x*y+2*y+z, [x, y, z]) == None
    assert collect(sin(x)*x+y+z, [x, y, z]) == None
    assert collect(sin(y)*x+y+z, [x, y, z]) == None

def test_coeff_list():
    x = Symbol('x')
    y = Symbol('y')
    z = Symbol('z')

    from sympy.modules.trigonometric import sin

    assert coeff_list(1) == [[1]]
    assert coeff_list(x) == [[1,1]]
    assert coeff_list(x**2+y**3) == [[1,2,0], [1,0,3]]
    assert coeff_list(x**2+y**3, [y,x]) == [[1,3,0], [1,0,2]]
    assert coeff_list(x*y) == [[1,1,1]]
    assert coeff_list(x**2*y**4 + sin(z)*x**3 + x*y**5, [x,y]) \
           == [[sin(z), 3, 0], [1, 2, 4], [1, 1, 5]]
    assert coeff_list(x**2*y**4 + sin(z)*x**3 + x*y**5, [x,y], order='grlex') \
           == [[1, 2, 4], [1, 1, 5], [sin(z), 3, 0]]
    assert coeff_list(x**2*y**4 + sin(z)*x**3 + x*y**5, [x,y],
               order='grevlex') == [[1, 1, 5], [1, 2, 4], [sin(z), 3, 0]]


    py.test.raises(PolynomialException, "coeff_list(sqrt(x),x)")
    py.test.raises(PolynomialException, "coeff_list(sin(x),x)")

def test_div_mv():
    x = Symbol('x')
    y = Symbol('y')

    assert div_mv(1,1) == [1, 0]
    assert div_mv(1,x,[x]) == [0, 1]
    assert div_mv(x*y+2*x+y,x,[x]) == [2+y, y]
    assert div_mv(x*y+2*x+y,x,[y]) == [2+(1+x)*y/x, 0]

    assert div_mv(x*y**2 + 1, [x*y+1, y+1], [x,y]) == [y, -1, 2]
    assert div_mv(x**2*y+x*y**2+y**2, [x*y-1, y**2-1], [x,y]) \
           == [x+y, 1, 1+x+y]
    assert div_mv(x**2*y+x*y**2+y**2, [y**2-1, x*y-1], [x,y]) \
           == [1+x, x, 1+2*x]

def test_groebner():
    x = Symbol('x')
    y = Symbol('y')
    z = Symbol('z')

    assert groebner(y*x, [x]) == [x]
    assert groebner(y*x, [x], reduced=False) == [x*y]
    assert groebner(x*y, [z]) == [1]
    
    # This one already is a Groebner base.
    assert groebner([y-x**2, z-x**3], [y,z,x], 'lex', False) \
           == [-x**2+y, z-x**3]

    assert groebner([x**3-2*x*y, x**2*y-2*y**2+x], [x,y], 'grlex', False) \
           == [x**3-2*x*y, x+x**2*y-2*y**2, -x**2, 2*x*y, 2*y**2-x]
    assert groebner([x**3-2*x*y, x**2*y-2*y**2+x], [x,y], 'grlex', True) \
           == [x**2, x*y, Rational(-1,2)*x+y**2]

def test_lcm_mv():
    x = Symbol('x')
    y = Symbol('y')

    assert lcm_mv(3, 4) == Rational(1)
    assert lcm_mv(4, y) == y
    assert lcm_mv(x, y) == x*y
    assert lcm_mv(y*(x+1), x, [x]) == x+x**2

def test_gcd_mv():
    x = Symbol('x')
    y = Symbol('y')

    assert gcd_mv(3, 4) == Rational(12)
    assert gcd_mv(3, 4, monic=True) == Rational(1)
    assert gcd_mv(x, y) == Rational(1)
    assert gcd_mv(x+y, x**2+2*x*y+y**2) == x+y

