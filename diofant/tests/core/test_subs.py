import pytest

from diofant import (Add, Basic, Derivative, Dict, E, Eq, Float, Function, I,
                     Integer, Lambda, Min, Piecewise, Rational, RootOf, Subs,
                     Symbol, Tuple, Wild, abc, atan2, cbrt, cos, cot, cse, exp,
                     factor, false, log, nsimplify, oo, pi, sin, sqrt, symbols,
                     tan, zoo)
from diofant.abc import a, b, c, d, e, t, x, y, z
from diofant.core.basic import _aresame
from diofant.core.cache import clear_cache


__all__ = ()


def test_subs():
    n3 = Integer(3)
    e = x
    e = e.subs({x: n3})
    assert e == 3

    e = 2*x
    assert e == 2*x
    e = e.subs({x: n3})
    assert e == 6


def test_trigonometric():
    n3 = Integer(3)
    e = (sin(x)**2).diff(x)
    assert e == 2*sin(x)*cos(x)
    e = e.subs({x: n3})
    assert e == 2*cos(n3)*sin(n3)

    e = (sin(x)**2).diff(x)
    assert e == 2*sin(x)*cos(x)
    e = e.subs({sin(x): cos(x)})
    assert e == 2*cos(x)**2

    i = Symbol('i', integer=True)
    assert tan(x).subs({x: pi/2}) is zoo
    assert cot(x).subs({x: pi}) is zoo
    assert cot(i*x).subs({x: pi}) is zoo
    assert tan(i*x).subs({x: pi/2}) == tan(i*pi/2)
    assert tan(i*x).subs({x: pi/2}).subs({i: 1}) is zoo
    o = Symbol('o', odd=True)
    assert tan(o*x).subs({x: pi/2}) == tan(o*pi/2)


def test_powers():
    assert sqrt(1 - sqrt(x)).subs({x: 4}) == I
    assert (sqrt(1 - x**2)**3).subs({x: 2}) == - 3*I*sqrt(3)
    assert cbrt(x).subs({x: 27}) == 3
    assert cbrt(x).subs({x: -27}) == 3*cbrt(-1)
    assert cbrt(-x).subs({x: 27}) == 3*cbrt(-1)
    n = Symbol('n', negative=True)
    assert (x**n).subs({x: 0}) is zoo
    assert exp(-1).subs({E: 0}) is zoo
    assert (x**(4.0*y)).subs({x**(2.0*y): n}) == n**2.0
    assert (2**(x + 2)).subs({2: 3}) == 3**(x + 3)

    # issue sympy/sympy#10829
    assert (4**x).subs({2**x: y}) == y**2
    assert (9**x).subs({3**x: y}) == y**2

    # issue sympy/sympy#6923
    assert (-2*x*sqrt(2)).subs({2*x: y}) == -sqrt(2)*y


def test_logexppow():   # no eval()
    x = Symbol('x', extended_real=True)
    w = Symbol('w')
    e = (3**(1 + x) + 2**(1 + x))/(3**x + 2**x)
    assert e.subs({2**x: w}) != e
    assert e.subs({exp(x*log(2)): w}) != e


def test_bug():
    x1 = Symbol('x1')
    x2 = Symbol('x2')
    y = x1*x2
    assert y.subs({x1: Float(3.0)}) == Float(3.0)*x2


def test_subbug1():
    # see that they don't fail
    (x**x).subs({x: 1})
    (x**x).subs({x: 1.0})


def test_subbug2():
    # Ensure this does not cause infinite recursion
    assert Float(7.7).epsilon_eq(abs(x).subs({x: -7.7}))


def test_dict_set():
    a, b = map(Wild, 'ab')

    f = 3*cos(4*x)
    r = f.match(a*cos(b*x))
    assert r == {a: 3, b: 4}
    e = a/b*sin(b*x)
    assert e.subs(r) == r[a]/r[b]*sin(r[b]*x)
    assert e.subs(r) == 3*sin(4*x) / 4
    s = set(r.items())
    assert e.subs(s) == r[a]/r[b]*sin(r[b]*x)
    assert e.subs(s) == 3*sin(4*x) / 4

    assert e.subs(r) == r[a]/r[b]*sin(r[b]*x)
    assert e.subs(r) == 3*sin(4*x) / 4
    assert x.subs(Dict((x, 1))) == 1


def test_dict_ambigous():   # see issue sympy/sympy#3566
    y = Symbol('y')
    z = Symbol('z')

    f = x*exp(x)
    g = z*exp(z)

    df = {x: y, exp(x): y}
    dg = {z: y, exp(z): y}

    assert f.subs(df) == y**2
    assert g.subs(dg) == y**2

    # and this is how order can affect the result
    assert f.subs({x: y}).subs({exp(x): y}) == y*exp(y)
    assert f.subs({exp(x): y}).subs({x: y}) == y**2

    # length of args and count_ops are the same so
    # default_sort_key resolves ordering...if one
    # doesn't want this result then an unordered
    # sequence should not be used.
    e = 1 + x*y
    assert e.subs({x: y, y: 2}) == 5
    # here, there are no obviously clashing keys or values
    # but the results depend on the order
    assert exp(x/2 + y).subs({exp(y + 1): 2, x: 2}) == exp(y + 1)


def test_deriv_sub_bug3():
    y = Symbol('y')
    f = Function('f')
    pat = Derivative(f(x), x, x)
    assert pat.subs({y: y**2}) == Derivative(f(x), x, x)
    assert pat.subs({y: y**2}) != Derivative(f(x), x)


def test_equality_subs1():
    f = Function('f')
    x = abc.x
    eq = Eq(f(x)**2, x)
    res = Eq(Integer(16), x)
    assert eq.subs({f(x): 4}) == res


def test_equality_subs2():
    f = Function('f')
    x = abc.x
    eq = Eq(f(x)**2, 16)
    assert bool(eq.subs({f(x): 3})) is False
    assert bool(eq.subs({f(x): 4})) is True


def test_sympyissue_3742():
    y = Symbol('y')

    e = sqrt(x)*exp(y)
    assert e.subs({sqrt(x): 1}) == exp(y)


def test_subs_dict1():
    assert (1 + x*y).subs({x: pi}) == 1 + pi*y
    assert (1 + x*y).subs({x: pi, y: 2}) == 1 + 2*pi

    c2, c3, q1p, q2p, c1, s1, s2, s3 = symbols('c2 c3 q1p q2p c1 s1 s2 s3')
    test = (c2**2*q2p*c3 + c1**2*s2**2*q2p*c3 + s1**2*s2**2*q2p*c3
            - c1**2*q1p*c2*s3 - s1**2*q1p*c2*s3)
    assert (test.subs({c1**2: 1 - s1**2, c2**2: 1 - s2**2, c3**3: 1 - s3**2})
            == c3*q2p*(1 - s2**2) + c3*q2p*s2**2*(1 - s1**2)
            - c2*q1p*s3*(1 - s1**2) + c3*q2p*s1**2*s2**2 - c2*q1p*s3*s1**2)


def test_mul():
    A, B, C = symbols('A B C', commutative=False)
    assert (x*y*z).subs({z*x: y}) == y**2
    assert (z*x).subs({1/x: z}) == z*x
    assert (x*y/z).subs({1/z: a}) == a*x*y
    assert (x*y/z).subs({x/z: a}) == a*y
    assert (x*y/z).subs({y/z: a}) == a*x
    assert (x*y/z).subs({x/z: 1/a}) == y/a
    assert (x*y/z).subs({x: 1/a}) == y/(z*a)
    assert (2*x*y).subs({5*x*y: z}) != 2*z/5
    assert (x*y*A).subs({x*y: a}) == a*A
    assert Subs(x*y*A, (x*y, a)).is_commutative is False
    assert (x**2*y**(3*x/2)).subs({x*y**(x/2): 2}) == 4*y**(x/2)
    assert (x*exp(x*2)).subs({x*exp(x): 2}) == 2*exp(x)
    assert ((x**(2*y))**3).subs({x**y: 2}) == 64
    assert (x*A*B).subs({x*A: y}) == y*B
    assert (x*y*(1 + x)*(1 + x*y)).subs({x*y: 2}) == 6*(1 + x)
    assert ((1 + A*B)*A*B).subs({A*B: x*A*B})
    assert (x*a/z).subs({x/z: A}) == a*A
    assert Subs(x*a/z, (x/z, A)).is_commutative is False
    assert (x**3*A).subs({x**2*A: a}) == a*x
    assert (x**2*A*B).subs({x**2*B: a}) == a*A
    assert (x**2*A*B).subs({x**2*A: a}) == a*B
    assert (b*A**3/(a**3*c**3)).subs({a**4*c**3*A**3/b**4: z}) == \
        b*A**3/(a**3*c**3)
    assert (6*x).subs({2*x: y}) == 3*y
    assert (y*exp(3*x/2)).subs({y*exp(x): 2}) == 2*exp(x/2)
    assert (y*exp(3*x/2)).subs({y*exp(x): 2}) == 2*exp(x/2)
    assert (A**2*B*A**2*B*A**2).subs({A*B*A: C}) == A*C**2*A
    assert (x*A**3).subs({x*A: y}) == y*A**2
    assert (x**2*A**3).subs({x*A: y}) == y**2*A
    assert (x*A**3).subs({x*A: B}) == B*A**2
    assert (x*A*B*A*exp(x*A*B)).subs({x*A: B}) == B**2*A*exp(B*B)
    assert (x**2*A*B*A*exp(x*A*B)).subs({x*A: B}) == B**3*exp(B**2)
    assert (x**3*A*exp(x*A*B)*A*exp(x*A*B)).subs({x*A: B}) == \
        x*B*exp(B**2)*B*exp(B**2)
    assert (x*A*B*C*A*B).subs({x*A*B: C}) == C**2*A*B
    assert (-I*a*b).subs({a*b: 2}) == -2*I

    # issue sympy/sympy#6361
    assert (-8*I*a).subs({-2*a: 1}) == 4*I
    assert (-I*a).subs({-a: 1}) == I

    # issue sympy/sympy#6441
    assert (4*x**2).subs({2*x: y}) == y**2
    assert (2*4*x**2).subs({2*x: y}) == 2*y**2
    assert (-x**3/9).subs({-x/3: z}) == -z**2*x
    assert (-x**3/9).subs({x/3: z}) == -z**2*x
    assert (-2*x**3/9).subs({x/3: z}) == -2*x*z**2
    assert (-2*x**3/9).subs({-x/3: z}) == -2*x*z**2
    assert (-2*x**3/9).subs({-2*x: z}) == z*x**2/9
    assert (-2*x**3/9).subs({2*x: z}) == -z*x**2/9
    assert (2*(3*x/5/7)**2).subs({3*x/5: z}) == 2*Rational(1, 7)**2*z**2
    assert (4*x).subs({-2*x: z}) == 4*x  # try keep subs literal

    assert (A**2*B**2).subs({A*B**3: C}) == A**2*B**2
    assert (A**Rational(5, 3)*B**3).subs({sqrt(A)*B: C}) == A**Rational(5, 3)*B**3
    assert (A**2*B**2*A).subs({A**2*B*A: C}) == A**2*B**2*A


def test_subs_simple():
    a = symbols('a', commutative=True)
    x = symbols('x', commutative=False)

    assert (2*a).subs({1: 3}) == 2*a
    assert (2*a).subs({2: 3}) == 3*a
    assert (2*a).subs({a: 3}) == 6
    assert sin(2).subs({1: 3}) == sin(2)
    assert sin(2).subs({2: 3}) == sin(3)
    assert sin(a).subs({a: 3}) == sin(3)

    assert (2*x).subs({1: 3}) == 2*x
    assert (2*x).subs({2: 3}) == 3*x
    assert (2*x).subs({x: 3}) == 6
    assert sin(x).subs({x: 3}) == sin(3)


def test_subs_constants():
    a, b = symbols('a b', commutative=True)
    x, y = symbols('x y', commutative=False)

    assert (a*b).subs({2*a: 1}) == a*b
    assert (1.5*a*b).subs({a: 1}) == 1.5*b
    assert (2*a*b).subs({2*a: 1}) == b
    assert (2*a*b).subs({4*a: 1}) == 2*a*b

    assert (x*y).subs({2*x: 1}) == x*y
    assert (1.5*x*y).subs({x: 1}) == 1.5*y
    assert (2*x*y).subs({2*x: 1}) == y
    assert (2*x*y).subs({4*x: 1}) == 2*x*y


def test_subs_commutative():
    a, b, c, d, K = symbols('a b c d K', commutative=True)

    assert (a*b).subs({a*b: K}) == K
    assert (a*b*a*b).subs({a*b: K}) == K**2
    assert (a*a*b*b).subs({a*b: K}) == K**2
    assert (a*b*c*d).subs({a*b*c: K}) == d*K
    assert (a*b**c).subs({a: K}) == K*b**c
    assert (a*b**c).subs({b: K}) == a*K**c
    assert (a*b**c).subs({c: K}) == a*b**K
    assert (a*b*c*b*a).subs({a*b: K}) == c*K**2
    assert (a**3*b**2*a).subs({a*b: K}) == a**2*K**2


def test_subs_noncommutative():
    w, x, y, z, L = symbols('w x y z L', commutative=False)

    assert (x*y).subs({x*y: L}) == L
    assert (w*y*x).subs({x*y: L}) == w*y*x
    assert (w*x*y*z).subs({x*y: L}) == w*L*z
    assert (x*y*x*y).subs({x*y: L}) == L**2
    assert (x*x*y).subs({x*y: L}) == x*L
    assert (x*x*y*y).subs({x*y: L}) == x*L*y
    assert (w*x*y).subs({x*y*z: L}) == w*x*y
    assert (x*y**z).subs({x: L}) == L*y**z
    assert (x*y**z).subs({y: L}) == x*L**z
    assert (x*y**z).subs({z: L}) == x*y**L
    assert (w*x*y*z*x*y).subs({x*y*z: L}) == w*L*x*y
    assert (w*x*y*y*w*x*x*y*x*y*y*x*y).subs({x*y: L}) == w*L*y*w*x*L**2*y*L

    # issue sympy/sympy#5284
    A, B = symbols('A B', commutative=False)
    x = symbols('x')
    assert (x*A).subs({x**2*A: B}) == x*A
    assert (A**2).subs({A**3: B}) == A**2
    assert (A**6).subs({A**3: B}) == B**2


def test_subs_basic_funcs():
    a, b, c, d, K = symbols('a b c d K', commutative=True)
    w, x, y, z, L = symbols('w x y z L', commutative=False)

    assert (x + y).subs({x + y: L}) == L
    assert (x - y).subs({x - y: L}) == L
    assert (x/y).subs({x: L}) == L/y
    assert (x**y).subs({x: L}) == L**y
    assert (x**y).subs({y: L}) == x**L
    assert ((a - c)/b).subs({b: K}) == (a - c)/K
    assert (exp(x*y - z)).subs({x*y: L}) == exp(L - z)
    assert (a*exp(x*y - w*z) + b*exp(x*y + w*z)).subs({z: 0}) == \
        a*exp(x*y) + b*exp(x*y)
    assert ((a - b)/(c*d - a*b)).subs({c*d - a*b: K}) == (a - b)/K
    assert (w*exp(a*b - c)*x*y/4).subs({x*y: L}) == w*exp(a*b - c)*L/4


def test_subs_wild():
    R, S, T, U = symbols('R S T U', cls=Wild)

    assert (R*S).subs({R*S: T}) == T
    assert (S*R).subs({R*S: T}) == T
    assert (R + S).subs({R + S: T}) == T
    assert (R**S).subs({R: T}) == T**S
    assert (R**S).subs({S: T}) == R**T
    assert (R*S**T).subs({R: U}) == U*S**T
    assert (R*S**T).subs({S: U}) == R*U**T
    assert (R*S**T).subs({T: U}) == R*S**U


def test_subs_mixed():
    a, b, c, K = symbols('a b c K', commutative=True)
    x, y, z, L = symbols('x y z L', commutative=False)
    R, S, T, U = symbols('R S T U', cls=Wild)

    assert (a*x*y).subs({x*y: L}) == a*L
    assert (a*b*x*y*x).subs({x*y: L}) == a*b*L*x
    assert (R*x*y*exp(x*y)).subs({x*y: L}) == R*L*exp(L)
    assert (a*x*y*y*x - x*y*z*exp(a*b)).subs({x*y: L}) == a*L*y*x - L*z*exp(a*b)
    e = c*y*x*y*x**(R*S - a*b) - T*(a*R*b*S)
    assert e.subs({x*y: L}).subs({a*b: K}).subs({R*S: U}) == \
        c*y*L*x**(U - K) - T*(U*K)


def test_division():
    a, c = symbols('a c', commutative=True)
    x, z = symbols('x z', commutative=True)

    assert (1/a).subs({a: c}) == 1/c
    assert (1/a**2).subs({a: c}) == 1/c**2
    assert (1/a**2).subs({a: -2}) == Rational(1, 4)
    assert (-(1/a**2)).subs({a: -2}) == -Rational(1, 4)

    assert (1/x).subs({x: z}) == 1/z
    assert (1/x**2).subs({x: z}) == 1/z**2
    assert (1/x**2).subs({x: -2}) == Rational(1, 4)
    assert (-(1/x**2)).subs({x: -2}) == -Rational(1, 4)

    # issue sympy/sympy#5360
    assert (1/x).subs({x: 0}) == 1/Integer(0)


def test_add():
    assert (a**2 - b - c).subs({a**2 - b: d}) in [d - c, a**2 - b - c]
    assert (a**2 - c).subs({a**2 - c: d}) == d
    assert (a**2 - b - c).subs({a**2 - c: d}) in [d - b, a**2 - b - c]
    assert (a**2 - x - c).subs({a**2 - c: d}) in [d - x, a**2 - x - c]
    assert (a**2 - b - sqrt(a)).subs({a**2 - sqrt(a): c}) == c - b
    assert (a + b + exp(a + b)).subs({a + b: c}) == c + exp(c)
    assert (c + b + exp(c + b)).subs({c + b: a}) == a + exp(a)
    assert (a + b + c + d).subs({b + c: x}) == a + d + x
    assert (a + b + c + d).subs({-b - c: x}) == a + d - x
    assert ((x + 1)*y).subs({x + 1: t}) == t*y
    assert ((-x - 1)*y).subs({x + 1: t}) == -t*y
    assert ((x - 1)*y).subs({x + 1: t}) == y*(t - 2)
    assert ((-x + 1)*y).subs({x + 1: t}) == y*(-t + 2)

    # this should work everytime:
    e = a**2 - b - c
    assert e.subs({Add(*e.args[:2]): d}) == d + e.args[2]
    assert e.subs({a**2 - c: d}) == d - b

    # the fallback should recognize when a change has
    # been made; while .1 == Rational(1, 10) they are not the same
    # and the change should be made
    assert (0.1 + a).subs({0.1: Rational(1, 10)}) == Rational(1, 10) + a

    e = (-x*(-y + 1) - y*(y - 1))
    ans = (-x*x - y*(-x)).expand()
    assert e.subs({-y + 1: x}) == ans


def test_sympyissue_4009():
    assert (I*Symbol('a')).subs({1: 2}) == I*Symbol('a')


def test_functions_subs():
    f, g = symbols('f g', cls=Function)
    l = Lambda((x, y), sin(x) + y)
    assert (g(y, x) + cos(x)).subs({g: l}) == sin(y) + x + cos(x)
    assert (f(x)**2).subs({f: sin}) == sin(x)**2
    assert (f(x, y)).subs({f: log}) == log(x, y)
    assert (f(x, y)).subs({f: sin}) == f(x, y)
    assert (sin(x) + atan2(x, y)).subs({atan2: f, sin: g}) == \
        f(x, y) + g(x)
    assert (g(f(x + y, x))).subs({f: l, g: Lambda(x, exp(x))}) == exp(x + sin(x + y))


def test_derivative_subs():
    y = Symbol('y')
    f = Function('f')
    assert Derivative(f(x), x).subs({f(x): y}) != 0
    assert Derivative(f(x), x).subs({f(x): y}).subs({y: f(x)}) == \
        Derivative(f(x), x)
    # issues sympy/sympy#5085, sympy/sympy#5037
    assert cse(Derivative(f(x), x) + f(x))[1][0].has(Derivative)
    assert cse(Derivative(f(x, y), x) +
               Derivative(f(x, y), y))[1][0].has(Derivative)


def test_derivative_subs2():
    f_func, g_func = symbols('f g', cls=Function)
    f, g = f_func(x, y, z), g_func(x, y, z)
    assert Derivative(f, x, y).subs({Derivative(f, x, y): g}) == g
    assert Derivative(f, y, x).subs({Derivative(f, x, y): g}) == g
    assert Derivative(f, x, y).subs({Derivative(f, x): g}) == Derivative(g, y)
    assert Derivative(f, x, y).subs({Derivative(f, y): g}) == Derivative(g, x)
    assert (Derivative(f, x, y, z).subs({Derivative(f, x, z): g}) == Derivative(g, y))
    assert (Derivative(f, x, y, z).subs({Derivative(f, z, y): g}) == Derivative(g, x))
    assert Derivative(f, x, y, z).subs({Derivative(f, z, y, x): g}) == g
    assert (Derivative(sin(x), (x, 2)).subs({Derivative(sin(x), f_func(x)): g_func}) ==
            Derivative(sin(x), (x, 2)))

    # issue sympy/sympy#9135
    assert (Derivative(f, x, x, y).subs({Derivative(f, y, y): g}) == Derivative(f, x, x, y))
    assert (Derivative(f, x, y, y, z).subs({Derivative(f, x, y, y, y): g}) == Derivative(f, x, y, y, z))

    assert Derivative(f, x, y).subs({Derivative(f_func(x), x, y): g}) == Derivative(f, x, y)


def test_derivative_subs3():
    x = Symbol('x')
    dex = Derivative(exp(x), x)
    assert Derivative(dex, x).subs({dex: exp(x)}) == dex
    assert dex.subs({exp(x): dex}) == Derivative(exp(x), x, x)


def test_subs_iter():
    assert x.subs(reversed([[x, y]])) == y
    it = iter([[x, y]])
    assert x.subs(it) == y
    assert x.subs(Tuple((x, y))) == y


def test_subs_dict():
    assert (2*x + y + z).subs({x: 1, y: 2}) == 4 + z

    l = [(sin(x), 2), (x, 1)]
    assert (sin(x)).subs(l) == \
           (sin(x)).subs(dict(l)) == 2
    assert sin(x).subs(reversed(l)) == sin(1)

    expr = sin(2*x) + sqrt(sin(2*x))*cos(2*x)*sin(exp(x)*x)
    reps = {sin(2*x): c, sqrt(sin(2*x)): a, cos(2*x): b,
            exp(x): e, x: d}
    assert expr.subs(reps) == c + a*b*sin(d*e)

    l = [(x, 3), (y, x**2)]
    assert (x + y).subs(l) == 3 + x**2
    assert (x + y).subs(reversed(l)) == 12

    # If changes are made to convert lists into dictionaries and do
    # a dictionary-lookup replacement, these tests will help to catch
    # some logical errors that might occur
    l = [(y, z + 2), (1 + z, 5), (z, 2)]
    assert (y - 1 + 3*x).subs(l) == 5 + 3*x
    l = [(y, z + 2), (z, 3)]
    assert (y - 2).subs(l) == 3


def test_no_arith_subs_on_floats():
    assert (x + 3).subs({x + 3: a}) == a
    assert (x + 3).subs({x + 2: a}) == a + 1

    assert (x + y + 3).subs({x + 3: a}) == a + y
    assert (x + y + 3).subs({x + 2: a}) == a + y + 1

    assert (x + 3.0).subs({x + 3.0: a}) == a
    assert (x + 3.0).subs({x + 2.0: a}) == x + 3.0

    assert (x + y + 3.0).subs({x + 3.0: a}) == a + y
    assert (x + y + 3.0).subs({x + 2.0: a}) == x + y + 3.0


def test_sympyissue_5651():
    a, b, c, K = symbols('a b c K', commutative=True)
    assert (a/(b*c)).subs({b*c: K}) == a/K
    assert (a/(b**2*c**3)).subs({b*c: K}) == a/(c*K**2)
    assert (1/(x*y)).subs({x*y: 2}) == Rational(1, 2)
    assert ((1 + x*y)/(x*y)).subs({x*y: 1}) == 2
    assert (x*y*z).subs({x*y: 2}) == 2*z
    assert ((1 + x*y)/(x*y)/z).subs({x*y: 1}) == 2/z


def test_sympyissue_6075():
    assert Tuple(1, True).subs({1: 2}) == Tuple(2, True)


def test_sympyissue_6079():
    # since x + 2.0 == x + 2 we can't do a simple equality test
    assert _aresame((x + 2.0).subs({2: 3}), x + 2.0)
    assert _aresame((x + 2.0).subs({2.0: 3}), x + 3)
    assert not _aresame(x + 2, x + 2.0)
    assert not _aresame(Basic(cos, 1), Basic(cos, 1.))
    assert _aresame(cos, cos)
    assert not _aresame(1, Integer(1))
    assert not _aresame(x, symbols('x', positive=True))


def test_sympyissue_4680():
    N = Symbol('N')
    assert N.subs({N: 3}) == 3


def test_sympyissue_6158():
    assert (x - 1).subs({1: y}) == x - y
    assert (x - 1).subs({-1: y}) == x + y
    assert (x - oo).subs({oo: y}) == x - y
    assert (x - oo).subs({-oo: y}) == x + y


def test_Function_subs():
    f, g, h, i = symbols('f g h i', cls=Function)
    p = Piecewise((g(f(x, y)), x < -1), (g(x), x <= 1))
    assert p.subs({g: h}) == Piecewise((h(f(x, y)), x < -1), (h(x), x <= 1))
    assert (f(y) + g(x)).subs({f: h, g: i}) == i(x) + h(y)


def test_simultaneous_subs():
    reps = {x: 0, y: 0}
    assert (x/y).subs(reps) != (y/x).subs(reps)
    assert (x/y).subs(reps, simultaneous=True) == \
        (y/x).subs(reps, simultaneous=True)
    reps = reps.items()
    assert (x/y).subs(reps) != (y/x).subs(reps)
    assert (x/y).subs(reps, simultaneous=True) == \
        (y/x).subs(reps, simultaneous=True)
    assert Derivative(x, y, z).subs(reps, simultaneous=True) == \
        Subs(Derivative(0, y, z), (y, 0))


def test_sympyissue_6419_6421():
    assert (1/(1 + x/y)).subs({x/y: x}) == 1/(1 + x)
    assert (-2*I).subs({2*I: x}) == -x
    assert (-I*x).subs({I*x: x}) == -x
    assert (-3*I*y**4).subs({3*I*y**2: x}) == -x*y**2


def test_sympyissue_6559():
    assert (-12*x + y).subs({-x: 1}) == 12 + y
    # though this involves cse it generated a failure in Mul._eval_subs
    x0 = Symbol('x0')
    e = -log(-12*sqrt(2) + 17)/24 - log(-2*sqrt(2) + 3)/12 + sqrt(2)/3
    # XXX modify cse so x1 is eliminated and x0 = -sqrt(2)?
    assert cse(e) == (
        [(x0, sqrt(2))], [x0/3 - log(-12*x0 + 17)/24 - log(-2*x0 + 3)/12])


def test_sympyissue_5261():
    x = symbols('x', extended_real=True)
    e = I*x
    assert exp(e).subs({exp(x): y}) == y**I
    assert (2**e).subs({2**x: y}) == y**I
    eq = (-2)**e
    assert eq.subs({(-2)**x: y}) == eq


@pytest.mark.xfail
def test_mul2():
    # When this fails, remove things labelled "2-arg hack"
    # 1) remove special handling in the fallback of subs that
    # was added in the same commit as this test
    # 2) remove the special handling in Mul.flatten
    assert (2*(x + 1)).is_Mul


def test_noncommutative_subs():
    x, y = symbols('x, y', commutative=False)
    assert (x*y*x).subs({x: x*y, y: x}, simultaneous=True) == (x*y*x**2*y)


def test_sympyissue_2877():
    f = Float(2.0)
    assert (x + f).subs({f: 2}) == x + 2

    def r(a, b, c):
        return factor(a*x**2 + b*x + c)
    e = r(5/6, 10, 5)
    assert nsimplify(e) == 5*x**2/6 + 10*x + 5


def test_sympyissue_5910():
    t = Symbol('t')
    assert (1/(1 - t)).subs({t: 1}) == zoo
    n = t
    d = t - 1
    assert (n/d).subs({t: 1}) == zoo
    assert (-n/-d).subs({t: 1}) == zoo


def test_sympyissue_5217():
    s = Symbol('s')
    z = (1 - 2*x*x)
    w = (1 + 2*x*x)
    q = 2*x*x*2*y*y
    sub = {2*x*x: s}
    assert w.subs(sub) == 1 + s
    assert z.subs(sub) == 1 - s
    assert q == 4*x**2*y**2
    assert q.subs(sub) == 2*y**2*s


def test_pow_eval_subs_no_cache():
    s = 1/sqrt(x**2)
    # This bug only appeared when the cache was turned off.
    clear_cache()

    # This used to fail with a wrong result.
    # It incorrectly returned 1/sqrt(x**2) before.
    result = s.subs({sqrt(x**2): y})
    assert result == 1/y


def test_issue_124():
    n = Symbol('n', integer=True)
    assert exp(n*x).subs({exp(x): x}) == x**n


def test_sympyissue_11159():
    exp1 = E
    exp0 = exp1*exp1
    assert exp0.subs({exp1: exp0}) == E**4


def test_sympyissue_10092():
    r = RootOf(x**3 - 17*x**2 + 81*x - 118, 0)
    assert (x < r).subs({x: r}) is false


def test_sympyissue_11746():
    x = Symbol('x', real=True)
    assert (1/x).subs({x**2: 1}) != 1


def test_diff_subs():
    # issue diofant/diofant#376
    f = symbols('f', cls=Function)
    e1 = Subs(Derivative(f(x), x), (x, y*z))
    e2 = Subs(Derivative(f(t), t), (t, y*z))
    assert (x*e1).subs({e2: y}) == x*y
    e3 = Subs(Derivative(f(t), t), (t, y*z**2))
    assert e3.subs({e2: y}) == e3
    e4 = Subs(Derivative(f(t), t, t), (t, y*z))
    assert e4.subs({e2: y}) == e4


def test_aresame_ordering():
    e = Min(3, z)
    s = {Min(z, 3): 3}
    assert e.subs(s) == 3


def test_underscores():
    _0, _1 = symbols('_0 _1')
    e = Subs(_0 + _1, (_0, 1), (_1, 0))
    assert e._expr == Symbol('__0') + Symbol('__1')
