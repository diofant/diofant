
import sys
sys.path.append(".")

import py

from sympy import *
from sympy.modules.polynomials import *

## sympy/modules/polynomials/base.py

def test_Polynomial():
    x = Symbol("x")
    y = Symbol("y")
    f = Polynomial(x+2)
    g = Polynomial(y**2-1)
    h = f + g
    assert f.var == [x]
    assert f.cl == [[1,1], [2,0]]
    assert str(f) == "2+x"
    assert repr(f) == "Polynomial(2+x, [x], 'grevlex', 'sym')"
    assert f.basic == x+2
    assert f == Polynomial(f.cl, f.var, f.var)
    assert h.var == [x, y]
    assert h.cl == [[1, 0, 2], [1, 1, 0], [1, 0, 0]]
    h = f*Polynomial(y,[x])
    assert h.var == [x]
    assert h.cl == [[y, 1], [2*y, 0]]
    h = f*y
    assert h.var == [x, y]
    assert h.cl == [[1, 1, 1], [2, 0, 1]]
    h.var=[y]
    assert h.var == [y]
    assert h.cl == [[2+x, 1]]

def test_coeff_list():
    x = Symbol('x')
    y = Symbol('y')
    z = Symbol('z')

    from sympy.modules.trigonometric import sin

    assert coeff_list(1) == [[1]]
    assert coeff_list(x) == [[1,1]]
    assert coeff_list(x**2+y**3, order='lex') == [[1,2,0], [1,0,3]]
    assert coeff_list(x**2+y**3, [y,x]) == [[1,3,0], [1,0,2]]
    assert coeff_list(x*y) == [[1,1,1]]
    assert coeff_list(x**2*y**4 + sin(z)*x**3 + x*y**5, [x,y], order='lex') \
           == [[sin(z), 3, 0], [1, 2, 4], [1, 1, 5]]
    assert coeff_list(x**2*y**4 + sin(z)*x**3 + x*y**5, [x,y], order='grlex') \
           == [[1, 2, 4], [1, 1, 5], [sin(z), 3, 0]]
    assert coeff_list(x**2*y**4 + sin(z)*x**3 + x**5*y, [x,y],
               order='grevlex') == [[1, 5, 1], [1, 2, 4], [sin(z), 3, 0]]
    assert coeff_list(z*x + x**2*y**2 + x**3*y, [z,x,y], order='1-el') \
           == [[1,1,1,0], [1,0,3,1], [1,0,2,2]]

    #TODO better test that differs between all orders ?

    py.test.raises(PolynomialException, "coeff_list(sqrt(x),x)")
    py.test.raises(PolynomialException, "coeff_list(sin(x),x)")
    
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

## sympy/modules/polynomials/wrapper.py

def test_coeff():
    x = Symbol("x")
    assert coeff(x**2, x, 1) == 0
    assert coeff(x**2, x, 2) == 1
    assert coeff(x**2, x, 2) != 0

    assert coeff(2*x+18*x**8, x, 1) == 2
    assert coeff(2*x+18*x**8, x, 4) == 0
    assert coeff(2*x+18*x**8, x, 8) == 18

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

def test_div():
    x = Symbol("x")
    y = Symbol('y')

    assert div(x**3-12*x**2-42, x-3, x) == (x**2-9*x-27, -123)
    assert div(x**3-12*x**2-42, x**2+x-3, x) == (x-13, 16*x-81)

    assert div(2+2*x+x**2, 1, x) == (2+2*x+x**2, 0)
    assert div(2+2*x+x**2, 2, x, coeff='int') == (1+x, x**2)

    assert div(3*x**3, x**2, x) == (3*x, 0)

    assert div(1,1) == (1, 0)
    assert div(1,x,[x]) == (0, 1)
    assert div(x*y+2*x+y,x,[x]) == (2+y, y)
    assert div(x*y+2*x+y,x,[y]) == (2+(1+1/x)*y, 0)

    assert div(x*y**2 + 1, [x*y+1, y+1], [x,y]) == ([y, -1], 2)
    assert div(x**2*y+x*y**2+y**2, [x*y-1, y**2-1], [x,y]) \
           == ([x+y, 1], 1+x+y)
    assert div(x**2*y+x*y**2+y**2, [y**2-1, x*y-1], [x,y]) \
           == ([1+x, x], 1+2*x)

def test_factor():
    x = Symbol("x")
    assert factor(x**2-1) == (x+1)*(x-1)
    assert factor(x**3-1) == (x-1)*(x**2+x+1)
    assert factor(x**2+2*x+1) == (x+1)**2
    assert factor(x**3-3*x**2+3*x-1) == (x-1)**3
    assert factor(x**3-3*x**2+3*x-1) == (x-1)**3
    assert factor(x**2+x-2) == (x-1)*(x+2)
    assert factor(x**3-x) == x*(x-1)*(x+1)

def test_gcd():
    x = Symbol("x")
    y = Symbol("y")
    z = Symbol("z")
    
    assert gcd(x**2, x, x) == x
    assert gcd(3*x**2, x, x) == x
    assert gcd(3*x**2, 6*x, x, coeff='int') == 3*x
    assert gcd(3*x**2, 6*x, x) == x
    assert gcd(x**2+2*x+1, x+1, x) == x+1
    assert gcd(x**2+2*x+2, x+1, x) == 1

    assert gcd(x**2+2*x+1, 2+2*x, x) == 1+x
    assert gcd(x**2+2*x+2, 2+2*x, x) == 1

    assert gcd(4, 6) == Rational(1)
    assert gcd(6, 4, coeff='int') == Rational(2)
    assert gcd(x, y) == Rational(1)
    assert gcd(sin(z)*(x+y), x**2+2*x*y+y**2, [x, y]) == x+y

def test_groebner():
    x = Symbol('x')
    y = Symbol('y')
    z = Symbol('z')

    assert groebner(y*x, [x]) == [x]
    assert groebner(y*x, [x], reduced=False) == [x*y]
    assert groebner(x*y, [z]) == [1]
    
    # This one already is a Groebner base.
    assert groebner([y-x**2, z-x**3], [y,z,x], 'lex', reduced=False) \
           == [-x**2+y, z-x**3]

    assert groebner([x**3-2*x*y, x**2*y-2*y**2+x], [x,y], 'grlex',
                    reduced=False) \
           == [x**3-2*x*y, x+x**2*y-2*y**2, -x**2, -2*x*y, -2*y**2+x]
    assert groebner([x**3-2*x*y, x**2*y-2*y**2+x], [x,y], 'grlex') \
           == [x**2, x*y, Rational(-1,2)*x+y**2]

def test_lcm():
    x = Symbol('x')
    y = Symbol('y')

    assert lcm(6, 4) == Rational(1)
    assert lcm(6, 4, coeff='int') == Rational(12)
    assert lcm(4, y) == y
    assert lcm(x, y) == x*y
    assert lcm(y*(x+1), x, [x]) ==x+x**2
    assert lcm(2*x, x**2, coeff='int') == 2*x**2 

def test_real_roots():
    x = Symbol('x')

    f = x-1
    assert real_roots(f) == 1
    assert real_roots(f, None, Rational(0)) == 0
    assert real_roots(f, Rational(0), Rational(1)) == 1
    assert real_roots(f, Rational(1), None) == 0
    f = x**2 - 4
    assert real_roots(f) == 2
    assert real_roots(f, None, Rational(0)) == 1
    assert real_roots(f, Rational(-1), Rational(1)) == 0

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

def test_roots():
    x = Symbol("x")
    assert roots(x**2-3*x+2) == [1,2]
    assert roots(x**2-3*x/2+Rational(1,2)) == [1,Rational(1)/2]
    assert roots(2*x**2-3*x+1) == [1,Rational(1)/2]
    assert roots(x**2-1) == [1,-1]
    assert roots(x**2+1) == []
    assert roots(x**3-1) == [1]

    assert roots(x**3) == [0]
    assert roots(x**3-x) == [0,1,-1]
    assert roots(Rational(2),x) == []

def test_sqf():
    x = Symbol("x")
    assert sqf(3*x**2, x) == 3*x**2
    assert sqf(x**2+2*x+1, x) == (x+1)**2

## sympy/modules/polynomials/ideals.py
    
def test_Ideal():
    x = Symbol('x')
    y = Symbol('y')
    z = Symbol('z')

    # TODO: more complete tests?
    assert len(Ideal()) == 1
    assert not x in Ideal()
    I = Ideal([x,y], [x,y])
    assert x*y**2 in I
    assert z*x in I
    assert not z in I
    assert I == I + Ideal()
    assert Ideal() == Ideal()*I
    assert I + I == I
    assert I == I % Ideal()
    assert Ideal() == Ideal(x*y, [x,y]) % I
    assert Ideal(z, [x,y,z]) == Ideal([x,z], [x,y,z]) % Ideal([x,y], [x,y,z])

## sympy/modules/polynomials/roots_.py

def test_sturm():
    from sympy.modules.polynomials import roots_
    x = Symbol('x')

    f = Polynomial(Rational(5), x)
    assert roots_.sturm(f) == [f]
    f = Polynomial(2*x)
    assert roots_.sturm(f) == [f, Polynomial(2, x)]
    f = Polynomial(x**3 - 2*x**2 + 3*x -5)
    assert roots_.sturm(f) == \
           [Polynomial(-5-2*x**2+x**3+3*x), Polynomial(3+3*x**2-4*x),
            Polynomial(Rational(13,3)-Rational(10,9)*x),
            Polynomial(Rational(-3303,100), [x])]
