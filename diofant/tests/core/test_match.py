import pytest

from diofant import (Add, Derivative, FiniteSet, Float, Function, I, Integer,
                     Mul, Rational, Symbol, Wild, WildFunction, cos, diff, exp,
                     log, meijerg, oo, pi, root, sin, sqrt, symbols)
from diofant.abc import X, Y, Z, a, b, c, gamma, mu, x, y


__all__ = ()


def test_symbol():
    a, c = map(Wild, 'ac')

    e = x
    assert e.match(x) == {}
    assert e.match(a) == {a: x}

    e = Integer(5)
    assert e.match(c) == {c: 5}
    assert e.match(e) == {}
    assert e.match(e + 1) is None


def test_add():
    p = Wild('p')

    e = a + b
    assert e.match(p + b) == {p: a}
    assert e.match(p + a) == {p: b}

    e = 1 + b
    assert e.match(p + b) == {p: 1}

    e = a + b + c
    assert e.match(a + p + c) == {p: b}
    assert e.match(b + p + c) == {p: a}

    e = a + b + c + x
    assert e.match(a + p + x + c) == {p: b}
    assert e.match(b + p + c + x) == {p: a}
    assert e.match(b) is None
    assert e.match(b + p) == {p: a + c + x}
    assert e.match(a + p + c) == {p: b + x}
    assert e.match(b + p + c) == {p: a + x}

    e = 4*x + 5
    assert e.match(4*x + p) == {p: 5}
    assert e.match(3*x + p) == {p: x + 5}
    assert e.match(p*x + 5) == {p: 4}


def test_power():
    p, q, r = map(Wild, 'pqr')

    e = (x + y)**a
    assert e.match(p**q) == {p: x + y, q: a}
    assert e.match(p**p) is None

    e = (x + y)**(x + y)
    assert e.match(p**p) == {p: x + y}
    assert e.match(p**q) == {p: x + y, q: x + y}

    e = (2*x)**2
    assert e.match(p*q**r) == {p: 4, q: x, r: 2}

    e = Integer(1)
    assert e.match(x**p) == {p: 0}


def test_match_exclude():
    p = Wild('p')
    q = Wild('q')
    r = Wild('r')

    e = Integer(6)
    assert e.match(2*p) == {p: 3}

    e = 3/(4*x + 5)
    assert e.match(3/(p*x + q)) == {p: 4, q: 5}

    e = 3/(4*x + 5)
    assert e.match(p/(q*x + r)) == {p: 3, q: 4, r: 5}

    e = 2/(x + 1)
    assert e.match(p/(q*x + r)) == {p: 2, q: 1, r: 1}

    e = 1/(x + 1)
    assert e.match(p/(q*x + r)) == {p: 1, q: 1, r: 1}

    e = 4*x + 5
    assert e.match(p*x + q) == {p: 4, q: 5}

    e = 4*x + 5*y + 6
    assert e.match(p*x + q*y + r) == {p: 4, q: 5, r: 6}

    a = Wild('a', exclude=[x])

    e = 3*x
    assert e.match(p*x) == {p: 3}
    assert e.match(a*x) == {a: 3}

    e = 3*x**2
    assert e.match(p*x) == {p: 3*x}
    assert e.match(a*x) is None

    e = 3*x + 3 + 6/x
    assert e.match(p*x**2 + p*x + 2*p) == {p: 3/x}
    assert e.match(a*x**2 + a*x + 2*a) is None


def test_mul():
    p, q = map(Wild, 'pq')

    e = 4*x
    assert e.match(p*x) == {p: 4}
    assert e.match(p*y) is None
    assert e.match(e + p*y) == {p: 0}

    e = a*x*b*c
    assert e.match(p*x) == {p: a*b*c}
    assert e.match(c*p*x) == {p: a*b}

    e = (a + b)*(a + c)
    assert e.match((p + b)*(p + c)) == {p: a}

    e = x
    assert e.match(p*x) == {p: 1}

    e = exp(x)
    assert e.match(x**p*exp(x*q)) == {p: 0, q: 1}

    e = I*x
    assert e.match(I*p) == {p: x}


def test_mul_noncommutative():
    A, B = symbols('A B', commutative=False)
    u, v = symbols('u v', cls=Wild)
    w, z = symbols('u v', cls=Wild, commutative=False)

    assert x.match(u*v) in ({v: x, u: 1}, {u: x, v: 1})
    assert (x*y).match(u*v) in ({v: y, u: x}, {u: y, v: x})
    assert A.match(u*v) is None
    assert (A*B).match(u*v) is None
    assert (x*A).match(u*v) is None
    assert (x*y*A).match(u*v) is None
    assert (x*A*B).match(u*v) is None
    assert (x*y*A*B).match(u*v) is None

    assert x.match(v*w) is None
    assert (x*y).match(v*w) is None
    assert A.match(v*w) == {w: A, v: 1}
    assert (A*B).match(v*w) == {w: A*B, v: 1}
    assert (x*A).match(v*w) == {w: A, v: x}
    assert (x*y*A).match(v*w) == {w: A, v: x*y}
    assert (x*A*B).match(v*w) == {w: A*B, v: x}
    assert (x*y*A*B).match(v*w) == {w: A*B, v: x*y}

    assert (-x).match(v*w) is None
    assert (-x*y).match(v*w) is None
    assert (-A).match(v*w) == {w: A, v: -1}
    assert (-A*B).match(v*w) == {w: A*B, v: -1}
    assert (-x*A).match(v*w) == {w: A, v: -x}
    assert (-x*y*A).match(v*w) == {w: A, v: -x*y}
    assert (-x*A*B).match(v*w) == {w: A*B, v: -x}
    assert (-x*y*A*B).match(v*w) == {w: A*B, v: -x*y}

    assert (x*y*A).match(u*v*w) == {u: x, v: y, w: A}
    assert (x*A).match(w*z) is None


def test_complex():
    x, y = map(Wild, 'xy')

    assert (1 + I).match(x + I) == {x: 1}
    assert (a + I).match(x + I) == {x: a}
    assert (2*I).match(x*I) == {x: 2}
    assert (a*I).match(x*I) == {x: a}
    assert (a*I).match(x*y) == {x: I, y: a}
    assert (2*I).match(x*y) == {x: 2, y: I}
    assert (a + b*I).match(x + y*I) == {x: a, y: b}


def test_functions():
    g = WildFunction('g')
    p = Wild('p')
    q = Wild('q')

    f = cos(5*x)
    notf = x
    assert f.match(p*cos(q*x)) == {p: 1, q: 5}
    assert f.match(p*g) == {p: 1, g: cos(5*x)}
    assert notf.match(g) is None

    F = WildFunction('F', nargs=2)
    assert F.nargs == FiniteSet(2)
    f = Function('f')
    assert f(x).match(F) is None

    F = WildFunction('F', nargs=(1, 2))
    assert F.nargs == FiniteSet(1, 2)

    # issue sympy/sympy#2711
    f = meijerg(((), ()), ((0,), ()), x)
    a = Wild('a')
    b = Wild('b')

    assert f.find(a) == {0: 1, x: 1, meijerg(((), ()), ((0,), ()), x): 1,
                         (): 3, (0,): 1, ((), ()): 1, ((0,), ()): 1}
    assert f.find(a + b) == {0: 1, x: 1, meijerg(((), ()), ((0,), ()), x): 1}
    assert f.find(a**2) == {x: 1, meijerg(((), ()), ((0,), ()), x): 1}


@pytest.mark.xfail
def test_functions_X1():
    g = WildFunction('g')
    p = Wild('p')
    q = Wild('q')

    f = cos(5*x)
    assert f.match(p*g(q*x)) == {p: 1, g: cos, q: 5}


def test_interface():
    p, q = map(Wild, 'pq')

    assert (x + 1).match(p + 1) == {p: x}
    assert (x*3).match(p*3) == {p: x}
    assert (x**3).match(p**3) == {p: x}
    assert (x*cos(y)).match(p*cos(q)) == {p: x, q: y}

    assert (x*y).match(p*q) in [{p: x, q: y}, {p: y, q: x}]
    assert (x + y).match(p + q) in [{p: x, q: y}, {p: y, q: x}]
    assert (x*y + 1).match(p*q) in [{p: 1, q: 1 + x*y}, {p: 1 + x*y, q: 1}]


def test_derivative1():
    p, q = map(Wild, 'pq')

    f = Function('f', nargs=1)
    fd = Derivative(f(x), x)

    assert fd.match(p) == {p: fd}
    assert (fd + 1).match(p + 1) == {p: fd}
    assert fd.match(fd) == {}
    assert (3*fd).match(p*fd) is not None
    assert (3*fd - 1).match(p*fd + q) == {p: 3, q: -1}


def test_derivative_bug1():
    f = Function('f')
    a = Wild('a', exclude=[f, x])
    b = Wild('b', exclude=[f])
    pattern = a * Derivative(f(x), x, x) + b
    expr = Derivative(f(x), x) + x**2
    d1 = {b: x**2}
    d2 = expr.match(pattern.xreplace(d1))
    assert d2 is None


def test_derivative2():
    f = Function('f')
    a = Wild('a', exclude=[f, x])
    b = Wild('b', exclude=[f])
    e = Derivative(f(x), x)
    assert e.match(Derivative(f(x), x)) == {}
    assert e.match(Derivative(f(x), x, x)) is None
    e = Derivative(f(x), x, x)
    assert e.match(Derivative(f(x), x)) is None
    assert e.match(Derivative(f(x), x, x)) == {}
    e = Derivative(f(x), x) + x**2
    assert e.match(a*Derivative(f(x), x) + b) == {a: 1, b: x**2}
    assert e.match(a*Derivative(f(x), x, x) + b) is None
    e = Derivative(f(x), x, x) + x**2
    assert e.match(a*Derivative(f(x), x) + b) is None
    assert e.match(a*Derivative(f(x), x, x) + b) == {a: 1, b: x**2}


def test_match_deriv_bug1():
    n = Function('n')
    l = Function('l')

    p = Wild('p')

    e = diff(l(x), x)/x - diff(diff(n(x), x), x)/2 - \
        diff(n(x), x)**2/4 + diff(n(x), x)*diff(l(x), x)/4
    e = e.subs({n(x): -l(x)}).doit()
    t = x*exp(-l(x))
    t2 = t.diff(x, x)/t
    assert e.match((p*t2).expand()) == {p: -Rational(1, 2)}


def test_match_bug2():
    p, q, r = map(Wild, 'pqr')
    res = (x + y).match(p + q + r)
    assert (p + q + r).subs(res) == x + y


def test_match_bug3():
    p = Wild('p')
    assert (b*x*exp(a*x)).match(x*exp(p*x)) is None


def test_match_bug4():
    p = Wild('p')
    e = x
    assert e.match(-p*x) == {p: -1}


def test_match_bug5():
    p = Wild('p')
    e = -x
    assert e.match(-p*x) == {p: 1}


def test_match_bug6():
    p = Wild('p')
    e = x
    assert e.match(3*p*x) == {p: Rational(1, 3)}


def test_match_polynomial():
    a = Wild('a', exclude=[x])
    b = Wild('b', exclude=[x])
    c = Wild('c', exclude=[x])
    d = Wild('d', exclude=[x])

    eq = 4*x**3 + 3*x**2 + 2*x + 1
    pattern = a*x**3 + b*x**2 + c*x + d
    assert eq.match(pattern) == {a: 4, b: 3, c: 2, d: 1}
    assert (eq - 3*x**2).match(pattern) == {a: 4, b: 0, c: 2, d: 1}
    assert (x + sqrt(2) + 3).match(a + b*x + c*x**2) == \
        {b: 1, a: sqrt(2) + 3, c: 0}


def test_exclude():
    p = Wild('p', exclude=[1, x])
    q = Wild('q')
    r = Wild('r', exclude=[sin, y])

    assert sin(x).match(r) is None
    assert cos(y).match(r) is None

    e = 3*x**2 + y*x + a
    assert e.match(p*x**2 + q*x + r) == {p: 3, q: y, r: a}

    e = x + 1
    assert e.match(x + p) is None
    assert e.match(p + 1) is None
    assert e.match(x + 1 + p) == {p: 0}

    e = cos(x) + 5*sin(y)
    assert e.match(r) is None
    assert e.match(cos(y) + r) is None
    assert e.match(r + p*sin(q)) == {r: cos(x), p: 5, q: y}


def test_floats():
    a, b = map(Wild, 'ab')

    e = cos(0.12345, evaluate=False)**2
    r = e.match(a*cos(b)**2)
    assert r == {a: 1, b: Float(0.12345)}


def test_Derivative_bug1():
    f = Function('f')
    a = Wild('a', exclude=[f(x)])
    b = Wild('b', exclude=[f(x)])
    eq = f(x).diff(x)
    assert eq.match(a*Derivative(f(x), x) + b) == {a: 1, b: 0}


def test_match_wild_wild():
    p = Wild('p')
    q = Wild('q')
    r = Wild('r')

    assert p.match(q + r) in [{q: p, r: 0}, {q: 0, r: p}]
    assert p.match(q * r) in [{q: p, r: 1}, {q: 1, r: p}]

    p = Wild('p')
    q = Wild('q', exclude=[p])
    r = Wild('r')

    assert p.match(q + r) == {q: 0, r: p}
    assert p.match(q*r) == {q: 1, r: p}

    p = Wild('p')
    q = Wild('q', exclude=[p])
    r = Wild('r', exclude=[p])

    assert p.match(q + r) is None
    assert p.match(q*r) is None


def test_combine_inverse():
    assert Mul._combine_inverse(x*I*y, x*I) == y
    assert Mul._combine_inverse(x*I*y, y*I) == x
    assert Mul._combine_inverse(oo*I*y, y*I) == oo
    assert Mul._combine_inverse(oo*I*y, oo*I) == y
    assert Add._combine_inverse(oo, oo) == Integer(0)
    assert Add._combine_inverse(oo*I, oo*I) == Integer(0)


def test_sympyissue_3773():
    z, phi, r = symbols('z phi r')
    c, A, B, N = symbols('c A B N', cls=Wild)
    l = Wild('l', exclude=(0,))

    eq = z * sin(2*phi) * r**7
    matcher = c * sin(phi*N)**l * r**A * log(r)**B

    assert eq.match(matcher) == {c: z, l: 1, N: 2, A: 7, B: 0}
    assert (-eq).match(matcher) == {c: -z, l: 1, N: 2, A: 7, B: 0}
    assert (x*eq).match(matcher) == {c: x*z, l: 1, N: 2, A: 7, B: 0}
    assert (-7*x*eq).match(matcher) == {c: -7*x*z, l: 1, N: 2, A: 7, B: 0}

    matcher = c*sin(phi*N)**l * r**A

    assert eq.match(matcher) == {c: z, l: 1, N: 2, A: 7}
    assert (-eq).match(matcher) == {c: -z, l: 1, N: 2, A: 7}
    assert (x*eq).match(matcher) == {c: x*z, l: 1, N: 2, A: 7}
    assert (-7*x*eq).match(matcher) == {c: -7*x*z, l: 1, N: 2, A: 7}


def test_sympyissue_3883():
    f = (-gamma * (x - mu)**2 - log(gamma) + log(2*pi))/2
    a, b, c = symbols('a b c', cls=Wild, exclude=(gamma,))

    assert f.match(a * log(gamma) + b * gamma + c) == \
        {a: -Rational(1, 2), b: -(mu - x)**2/2, c: log(2*pi)/2}
    assert f.expand().collect(gamma).match(a * log(gamma) + b * gamma + c) == \
        {a: -Rational(1, 2), b: (-(x - mu)**2/2).expand(), c: (log(2*pi)/2).expand()}
    g1 = Wild('g1', exclude=[gamma])
    g2 = Wild('g2', exclude=[gamma])
    g3 = Wild('g3', exclude=[gamma])
    assert f.expand().match(g1 * log(gamma) + g2 * gamma + g3) == \
        {g3: log(2)/2 + log(pi)/2, g1: -Rational(1, 2), g2: -mu**2/2 + mu*x - x**2/2}


def test_sympyissue_4418():
    a, b, c = symbols('a b c', cls=Wild, exclude=(x,))
    f, g = symbols('f g', cls=Function)

    eq = diff(g(x)*f(x).diff(x), x)

    assert eq.match(
        g(x).diff(x)*f(x).diff(x) + g(x)*f(x).diff(x, x) + c) == {c: 0}
    assert eq.match(a*g(x).diff(
        x)*f(x).diff(x) + b*g(x)*f(x).diff(x, x) + c) == {a: 1, b: 1, c: 0}


def test_sympyissue_4700():
    f = Function('f')
    a, b = symbols('a b', cls=Wild, exclude=(f(x),))

    p = a*f(x) + b
    eq1 = sin(x)
    eq2 = f(x) + sin(x)
    eq3 = f(x) + x + sin(x)
    eq4 = x + sin(x)

    assert eq1.match(p) == {a: 0, b: sin(x)}
    assert eq2.match(p) == {a: 1, b: sin(x)}
    assert eq3.match(p) == {a: 1, b: x + sin(x)}
    assert eq4.match(p) == {a: 0, b: x + sin(x)}


def test_sympyissue_5168():
    a, b, c = symbols('a b c', cls=Wild)
    f = Function('f')

    assert x.match(a) == {a: x}
    assert x.match(a*f(x)**c) == {a: x, c: 0}
    assert x.match(a*b) == {a: 1, b: x}
    assert x.match(a*b*f(x)**c) == {a: 1, b: x, c: 0}

    assert (-x).match(a) == {a: -x}
    assert (-x).match(a*f(x)**c) == {a: -x, c: 0}
    assert (-x).match(a*b) == {a: -1, b: x}
    assert (-x).match(a*b*f(x)**c) == {a: -1, b: x, c: 0}

    assert (2*x).match(a) == {a: 2*x}
    assert (2*x).match(a*f(x)**c) == {a: 2*x, c: 0}
    assert (2*x).match(a*b) == {a: 2, b: x}
    assert (2*x).match(a*b*f(x)**c) == {a: 2, b: x, c: 0}

    assert (-2*x).match(a) == {a: -2*x}
    assert (-2*x).match(a*f(x)**c) == {a: -2*x, c: 0}
    assert (-2*x).match(a*b) == {a: -2, b: x}
    assert (-2*x).match(a*b*f(x)**c) == {a: -2, b: x, c: 0}


def test_sympyissue_4559():
    e = Symbol('e')
    w = Wild('w', exclude=[x])
    y = Wild('y')

    # this is as it should be

    assert (3/x).match(w/y) == {w: 3, y: x}
    assert (3*x).match(w*y) == {w: 3, y: x}
    assert (x/3).match(y/w) == {w: 3, y: x}
    assert (3*x).match(y/w) == {w: Rational(1, 3), y: x}

    # these could be allowed to fail

    assert (x/3).match(w/y) == {w: Rational(1, 3), y: 1/x}
    assert (3*x).match(w/y) == {w: 3, y: 1/x}
    assert (3/x).match(w*y) == {w: 3, y: 1/x}

    # Note that solve will give
    # multiple roots but match only gives one:
    #
    # >>> solve(x**r - y**2, y)
    # [{y: -x**(r/2)}, {y: x**(r/2)}]

    r = Symbol('r', rational=True)
    assert (x**r).match(y**2) == {y: x**(r/2)}
    assert (x**e).match(y**2) == {y: sqrt(x**e)}

    # since (x**i = y) -> x = y**(1/i) where i is an integer
    # the following should also be valid as long as y is not
    # zero when i is negative.

    a = Wild('a')

    e = Integer(0)
    assert e.match(a) == {a: e}
    assert e.match(1/a) is None
    assert e.match(a**.3) is None

    e = Integer(3)
    assert e.match(1/a) == {a: 1/e}
    assert e.match(1/a**2) == {a: 1/sqrt(e)}
    e = pi
    assert e.match(1/a) == {a: 1/e}
    assert e.match(1/a**2) == {a: 1/sqrt(e)}
    assert (-e).match(sqrt(a)) is None
    assert (-e).match(a**2) == {a: I*sqrt(pi)}


def test_sympyissue_4883():
    a = Wild('a')

    e = [i**2 for i in (x - 2, 2 - x)]
    p = [i**2 for i in (x - a, a - x)]
    for eq in e:
        for pat in p:
            assert eq.match(pat) == {a: 2}


def test_sympyissue_4319():
    p = -x*(Rational(1, 8) - y)
    ans = {0, y - Rational(1, 8)}

    def ok(pat):
        assert set(p.match(pat).values()) == ans

    ok(Wild('coeff', exclude=[x])*x + Wild('rest'))
    ok(Wild('w', exclude=[x])*x + Wild('rest'))
    ok(Wild('coeff', exclude=[x])*x + Wild('rest'))
    ok(Wild('w', exclude=[x])*x + Wild('rest'))
    ok(Wild('e', exclude=[x])*x + Wild('rest'))
    ok(Wild('ress', exclude=[x])*x + Wild('rest'))
    ok(Wild('resu', exclude=[x])*x + Wild('rest'))


def test_sympyissue_3778():
    p, c, q = symbols('p c q', cls=Wild)

    assert (sin(x)**2).match(sin(p)*sin(q)*c) == {q: x, c: 1, p: x}
    assert (2*sin(x)).match(sin(p) + sin(q) + c) == {q: x, c: 0, p: x}


def test_sympyissue_6103():
    a = Wild('a')
    assert (-I*x*oo).match(I*a*oo) == {a: -x}


def test_sympyissue_3539():
    a = Wild('a')
    assert (x - 2).match(a - x) is None
    assert (6/x).match(a*x) is None
    assert (6/x**2).match(a/x) == {a: 6/x}


def test_issue_423():
    a1, b1, c1, d1, a2, b2, c2, d2 = symbols('a1 b1 c1 d1 a2 b2 c2 d2',
                                             cls=Wild, exclude=[x])
    pat = (a1*x + b1)/(c1*x + d1) + (a2*x + b2)/(c2*x + d2)
    expr = (2*x + 9)/(x + 4) + (4*x + 3)/(5*x + 12)
    ans = {a1: 2, a2: 4, b1: 9, b2: 3, c1: 1, c2: 5, d1: 4, d2: 12}
    assert expr.match(pat) == ans


def test_sympyissue_8694():
    theta1, theta2, rho = symbols('theta1, theta2, rho')
    S1, C1 = sin(theta1), cos(theta1)
    X1 = Wild('X1', exclude=[rho, theta1, theta2])
    Y1 = Wild('Y1', exclude=[rho, theta1, theta2])
    Z1 = Wild('Z1', exclude=[rho, theta1, theta2])
    eq = -Y + (-X + Z)*cos(theta1) + (X + Y)*sin(theta1)
    assert eq.match(X1*C1 + Y1*S1 + Z1) == {X1: Z - X, Y1: X + Y, Z1: -Y}
    eq = -Y + Z*cos(theta1) + (X + Y)*sin(theta1)
    assert eq.match(X1*C1 + Y1*S1 + Z1) == {X1: Z, Y1: X + Y, Z1: -Y}


def test_issue_462():
    w = Wild('w', exclude=[x])
    assert (-x).match(2*w*x) == {w: Rational(-1, 2)}
    assert (-x**2).match(2*w*x**2) == {w: Rational(-1, 2)}
    # see also sympy/sympy#12238
    a11, a12, a22, a13, a23, a33 = symbols('a11 a12 a22 a13 a23 a33',
                                           exclude=[x, y], cls=Wild)
    eq = -x**2 + 12*x*y - 8*x - 36*y**2 - y - 4
    tmpl = a11*x**2 + 2*a12*x*y + a22*y**2 + 2*a13*x + 2*a23*y + a33
    assert eq.match(tmpl) == {a11: -1, a12: 6, a22: -36, a13: -4,
                              a23: Rational(-1, 2), a33: -4}


def test_sympyissue_16774():
    a1, b1, c1, d1, a2, b2, c2, d2 = symbols('a1 b1 c1 d1 a2 b2 c2 d2',
                                             cls=Wild, exclude=[a])
    e = (4*a + 3)/(5*a + 12) - (2*a + 9)/(a + 4)
    t = (a*a1 + b1)/(a*c1 + d1) + (a*a2 + b2)/(a*c2 + d2)
    assert e.match(t) == {a1: -2, a2: 4, b1: -9, b2: 3, c1: 1,
                          c2: 5, d1: 4, d2: 12}


def test_sympyissue_21466():
    a, b, k, m, n = symbols('a b k m n', cls=Wild)
    pattern = x**n * (a + b * x**k)**m
    expr = root(x, 3)*root(3 - x**2, 3)
    assert expr.match(pattern) == {k: 2, b: -1, a: 3, m: Rational(1, 3),
                                   n: Rational(1, 3)}
