"""Tests for the implementation of RootOf class and related tools."""

import pytest

from diofant import (Eq, Float, Function, GeneratorsNeededError, I, Integer,
                     Lambda, MultivariatePolynomialError, PolynomialError, Pow,
                     PurePoly, Rational, RootOf, RootSum, Symbol, conjugate,
                     exp, expand_func, false, legendre_poly, log, oo, root,
                     solve, sqrt, tan, true)
from diofant.abc import a, b, r, x, y, z


__all__ = ()


def test_RootOf___new__():
    assert RootOf(x, 0) == 0
    assert RootOf(x, -1) == 0

    assert RootOf(x - 1, 0) == 1
    assert RootOf(x - 1, -1) == 1

    assert RootOf(x + 1, 0) == -1
    assert RootOf(x + 1, -1) == -1

    assert RootOf(x**2 + 2*x + 3, 0) == -1 - I*sqrt(2)
    assert RootOf(x**2 + 2*x + 3, 1) == -1 + I*sqrt(2)
    assert RootOf(x**2 + 2*x + 3, -1) == -1 + I*sqrt(2)
    assert RootOf(x**2 + 2*x + 3, -2) == -1 - I*sqrt(2)

    r = RootOf(x**2 + 2*x + 3, 0, radicals=False)
    assert isinstance(r, RootOf) is True

    r = RootOf(x**2 + 2*x + 3, 1, radicals=False)
    assert isinstance(r, RootOf) is True

    r = RootOf(x**2 + 2*x + 3, -1, radicals=False)
    assert isinstance(r, RootOf) is True

    r = RootOf(x**2 + 2*x + 3, -2, radicals=False)
    assert isinstance(r, RootOf) is True

    assert RootOf((x - 1)*(x + 1), 0, radicals=False) == -1
    assert RootOf((x - 1)*(x + 1), 1, radicals=False) == 1
    assert RootOf((x - 1)*(x + 1), -1, radicals=False) == 1
    assert RootOf((x - 1)*(x + 1), -2, radicals=False) == -1

    assert RootOf((x - 1)*(x + 1), 0, radicals=True) == -1
    assert RootOf((x - 1)*(x + 1), 1, radicals=True) == 1
    assert RootOf((x - 1)*(x + 1), -1, radicals=True) == 1
    assert RootOf((x - 1)*(x + 1), -2, radicals=True) == -1

    assert RootOf((x - 1)*(x**3 + x + 3), 0) == RootOf(x**3 + x + 3, 0)
    assert RootOf((x - 1)*(x**3 + x + 3), 1) == 1
    assert RootOf((x - 1)*(x**3 + x + 3), 2) == RootOf(x**3 + x + 3, 1)
    assert RootOf((x - 1)*(x**3 + x + 3), 3) == RootOf(x**3 + x + 3, 2)
    assert RootOf((x - 1)*(x**3 + x + 3), -1) == RootOf(x**3 + x + 3, 2)
    assert RootOf((x - 1)*(x**3 + x + 3), -2) == RootOf(x**3 + x + 3, 1)
    assert RootOf((x - 1)*(x**3 + x + 3), -3) == 1
    assert RootOf((x - 1)*(x**3 + x + 3), -4) == RootOf(x**3 + x + 3, 0)

    assert RootOf(x**4 + 3*x**3, 0) == -3
    assert RootOf(x**4 + 3*x**3, 1) == 0
    assert RootOf(x**4 + 3*x**3, 2) == 0
    assert RootOf(x**4 + 3*x**3, 3) == 0

    pytest.raises(GeneratorsNeededError, lambda: RootOf(0, 0))
    pytest.raises(GeneratorsNeededError, lambda: RootOf(1, 0))

    pytest.raises(PolynomialError, lambda: RootOf(Integer(0).as_poly(x), 0))
    pytest.raises(PolynomialError, lambda: RootOf(Integer(1).as_poly(x), 0))

    pytest.raises(PolynomialError, lambda: RootOf(x - y, 0))

    pytest.raises(IndexError, lambda: RootOf(x**2 - 1, -4))
    pytest.raises(IndexError, lambda: RootOf(x**2 - 1, -3))
    pytest.raises(IndexError, lambda: RootOf(x**2 - 1, 2))
    pytest.raises(IndexError, lambda: RootOf(x**2 - 1, 3))
    pytest.raises(ValueError, lambda: RootOf(x**2 - 1, x))
    pytest.raises(NotImplementedError,
                  lambda: RootOf(Symbol('a', nonzero=False)*x**5 +
                                 2*x - 1, x, 0))
    pytest.raises(NotImplementedError,
                  lambda: (Symbol('a', nonzero=False)*x**5 +
                           2*x - 1).as_poly(x).all_roots())

    assert RootOf((x - y).as_poly(x), 0) == y

    assert RootOf((x**2 - y).as_poly(x), 0) == -sqrt(y)
    assert RootOf((x**2 - y).as_poly(x), 1) == sqrt(y)

    assert isinstance(RootOf(x**3 - y, x, 0), RootOf)
    p = Symbol('p', positive=True)
    assert RootOf(x**3 - p, x, 0) == root(p, 3)*RootOf(x**3 - 1, 0)

    assert RootOf(y*x**3 + y*x + 2*y, x, 0) == -1

    assert RootOf(x**3 + x + 1, 0).is_commutative is True

    e = RootOf(x**2 - 4, x, 1, evaluate=False)
    assert isinstance(e, RootOf)
    assert e.doit() == 2
    assert e.args == (x**2 - 4, x, 1)
    assert e.poly == PurePoly(x**2 - 4, x)
    assert e.index == 1

    assert RootOf(x**7 - 0.1*x + 1, 0) == RootOf(10*x**7 - x + 10, 0)

    e = x**7 - x
    p = (x**7 - x).as_poly(modulus=7)
    F7 = p.domain
    assert (RootOf(p, 1) == RootOf(e, 1, modulus=7) ==
            RootOf(e, x, 1, modulus=7) == RootOf(p, 1, evaluate=False) ==
            RootOf(e, 1, domain=F7, evaluate=False))


def test_RootOf_attributes():
    r = RootOf(x**3 + x + 3, 0)
    assert r.is_number
    assert r.free_symbols == set()

    r = RootOf(x**3 + y*x + 1, x, 0)
    assert isinstance(r, RootOf)
    assert r.expr == x**3 + y*x + 1
    assert r.free_symbols == {y}
    assert r.is_number is False


def test_RootOf___eq__():
    assert (RootOf(x**3 + x + 3, 0) == RootOf(x**3 + x + 3, 0)) is True
    assert (RootOf(x**3 + x + 3, 0) == RootOf(x**3 + x + 3, 1)) is False
    assert (RootOf(x**3 + x + 3, 1) == RootOf(x**3 + x + 3, 1)) is True
    assert (RootOf(x**3 + x + 3, 1) == RootOf(x**3 + x + 3, 2)) is False
    assert (RootOf(x**3 + x + 3, 2) == RootOf(x**3 + x + 3, 2)) is True

    assert (RootOf(x**3 + x + 3, 0) == RootOf(y**3 + y + 3, 0)) is True
    assert (RootOf(x**3 + x + 3, 0) == RootOf(y**3 + y + 3, 1)) is False
    assert (RootOf(x**3 + x + 3, 1) == RootOf(y**3 + y + 3, 1)) is True
    assert (RootOf(x**3 + x + 3, 1) == RootOf(y**3 + y + 3, 2)) is False
    assert (RootOf(x**3 + x + 3, 2) == RootOf(y**3 + y + 3, 2)) is True


def test_RootOf___eval_Eq__():
    f = Function('f')
    r = RootOf(x**3 + x + 3, 2)
    r1 = RootOf(x**3 + x + 3, 1)
    assert Eq(r, r1) is false
    assert Eq(r, r) is true
    assert Eq(r, x) is false
    assert Eq(r, 0) is false
    assert Eq(r, oo) is false
    assert Eq(r, I) is false
    assert Eq(r, f(0)) is false
    assert Eq(r, f(0)) is false
    sol = solve(r.expr, x)
    for s in sol:
        if s[x].is_real:
            assert Eq(r, s[x]) is false
    r = RootOf(r.expr, 0)
    for s in sol:
        if s[x].is_real:
            assert Eq(r, s[x]) is true
    eq = x**3 + x + 1
    assert ([Eq(RootOf(eq, i), j[x])
             for i in range(3) for j in solve(eq)] ==
            [False, False, True, False, True, False, True, False, False])
    assert Eq(RootOf(eq, 0), 1 + I) is false


def test_RootOf_is_real():
    assert RootOf(x**3 + x + 3, 0).is_real is True
    assert RootOf(x**3 + x + 3, 1).is_real is False
    assert RootOf(x**3 + x + 3, 2).is_real is False

    r = RootOf(x**3 + y*x + 1, x, 0)
    assert r.is_real is None

    assert RootOf(x**3 + I*x + 2, 0).is_real is False


def test_RootOf_is_imaginary():
    assert RootOf(x**3 + x + 3, 0).is_imaginary is False
    assert RootOf(x**3 + x + 3, 1).is_imaginary is False
    assert RootOf(x**3 + y*x + 1, x, 0).is_imaginary is None

    assert RootOf(x**3 + I*x + 2, 0).is_real is False

    assert RootOf(x**4 + 10*x**2 + 1, 2).is_imaginary is True


def test_RootOf_is_complex():
    assert RootOf(x**3 + x + 3, 0).is_complex is True
    assert RootOf(x**3 + y*x + 3, x, 0).is_complex is None

    assert RootOf(x**3 + y*x + 3, x, 0).is_commutative

    assert RootOf(x**3 + I*x + 2, 0).is_complex is True


def test_RootOf_is_algebraic():
    assert RootOf(x**3 + x + 3, 0).is_algebraic is True
    assert RootOf(x**3 + y*x + 3, x, 0).is_algebraic is None


def test_RootOf_power():
    e = RootOf(y**3 - x, y, 0)
    assert e**3 == x
    assert e**2 == Pow(e, 2, evaluate=False)
    e2 = RootOf(y**3 - x*y, y, 0)
    assert e2**3 == Pow(e2, 3, evaluate=False)
    e3 = RootOf(3*x**5 + 2*x - 1, 0)
    assert e3**5 == -2*e3/3 + Rational(1, 3)  # issue sympy/sympy#8543
    assert e3**4 == Pow(e3, 4, evaluate=False)
    assert e3**-1 == 3*e3**4 + 2


def test_RootOf_conjugate():
    p = x**7 + x + 1
    assert RootOf(p, 0).conjugate() == RootOf(p, 0)
    assert RootOf(p, 1).conjugate() == RootOf(p, 2)
    assert RootOf(p, 2).conjugate() == RootOf(p, 1)
    assert RootOf(p, 6).conjugate() == RootOf(p, 5)

    p2 = p*(x - 123)
    assert RootOf(p2, 0).conjugate() == RootOf(p2, 0)
    assert RootOf(p2, 1).conjugate() == RootOf(p2, 1)
    assert RootOf(p2, 2).conjugate() == RootOf(p2, 3)
    assert RootOf(p2, 3).conjugate() == RootOf(p2, 2)
    assert RootOf(p2, 7).conjugate() == RootOf(p2, 6)

    p3 = (x**7 + x*y + 1).as_poly(x)
    assert RootOf(p3, x, 0).conjugate() == conjugate(RootOf(p3, x, 0),
                                                     evaluate=False)

    p4 = x**12 - 4*x**8 + 2*x**6 + 4*x**4 + 4*x**2 + 1
    r4 = RootOf(p4, 4)
    r5 = RootOf(p4, 5)
    assert r4.conjugate() == r5
    assert r4.evalf() == -r5.evalf()


def test_RootOf_subs():
    assert RootOf(x**3 + x + 1, 0).subs({x: y}) == RootOf(y**3 + y + 1, 0)
    eq = -x + RootOf(y**3 - x**3 + 3*x**2, y, 0) + 1
    assert eq.subs({x: Rational(1, 3)}) == 0


def test_RootOf_diff():
    assert RootOf(x**3 + x + 1, 0).diff(x) == 0
    assert RootOf(x**3 + x + 1, 0).diff(y) == 0

    r = RootOf(x**7 + x*y + 1, x, 0)
    assert r.diff(y) == -r/(y + 7*r**6)
    assert r.diff(x) == 0


def test_RootOf_evalf():
    real = RootOf(x**3 + x + 3, 0).evalf(20)

    assert real.epsilon_eq(Float('-1.2134116627622296341'))

    re, im = RootOf(x**3 + x + 3, 1).evalf(20).as_real_imag()

    assert re.epsilon_eq(+Float('0.60670583138111481707'))
    assert im.epsilon_eq(-Float('1.45061224918844152650'))

    re, im = RootOf(x**3 + x + 3, 2).evalf(20).as_real_imag()

    assert re.epsilon_eq(Float('0.60670583138111481707'))
    assert im.epsilon_eq(Float('1.45061224918844152650'))

    p = legendre_poly(4, x, polys=True)
    roots = [str(r.evalf(17)) for r in p.real_roots()]
    assert roots == [
        '-0.86113631159405258',
        '-0.33998104358485626',
        '0.33998104358485626',
        '0.86113631159405258',
    ]

    re = RootOf(x**5 - 5*x + 12, 0).evalf(20)
    assert re.epsilon_eq(Float('-1.84208596619025438271'))

    re, im = RootOf(x**5 - 5*x + 12, 1).evalf(20).as_real_imag()
    assert re.epsilon_eq(Float('-0.351854240827371999559'))
    assert im.epsilon_eq(Float('-1.709561043370328882010'))

    re, im = RootOf(x**5 - 5*x + 12, 2).evalf(20).as_real_imag()
    assert re.epsilon_eq(Float('-0.351854240827371999559'))
    assert im.epsilon_eq(Float('+1.709561043370328882010'))

    re, im = RootOf(x**5 - 5*x + 12, 3).evalf(20).as_real_imag()
    assert re.epsilon_eq(Float('+1.272897223922499190910'))
    assert im.epsilon_eq(Float('-0.719798681483861386681'))

    re, im = RootOf(x**5 - 5*x + 12, 4).evalf(20).as_real_imag()
    assert re.epsilon_eq(Float('+1.272897223922499190910'))
    assert im.epsilon_eq(Float('+0.719798681483861386681'))

    # issue sympy/sympy#6393
    assert str(RootOf(x**5 + 2*x**4 + x**3 - 68719476736, 0).evalf(3)) == '147.'
    eq = (531441*x**11 + 3857868*x**10 + 13730229*x**9 + 32597882*x**8 +
          55077472*x**7 + 60452000*x**6 + 32172064*x**5 - 4383808*x**4 -
          11942912*x**3 - 1506304*x**2 + 1453312*x + 512)
    a, b = RootOf(eq, 1).evalf(2).as_real_imag()
    c, d = RootOf(eq, 2).evalf(2).as_real_imag()
    assert a == c
    assert b < d
    assert b == -d
    # issue sympy/sympy#6451
    r = RootOf(legendre_poly(64, x), 7)
    assert r.evalf(2) == r.evalf(100).evalf(2)
    # issue sympy/sympy#8617
    ans = [w[x].evalf(2) for w in solve(x**3 - x - 4)]
    assert RootOf(exp(x)**3 - exp(x) - 4, 0).evalf(2) in ans
    # issue sympy/sympy#9019
    r0 = RootOf(x**2 + 1, 0, radicals=False)
    r1 = RootOf(x**2 + 1, 1, radicals=False)
    assert r0.evalf(4, chop=True) == -1.0*I
    assert r1.evalf(4, chop=True) == +1.0*I

    # make sure verification is used in case a max/min traps the "root"
    assert str(RootOf(4*x**5 + 16*x**3 + 12*x**2 + 7, 0).evalf(3)) == '-0.976'

    assert isinstance(RootOf(x**3 + y*x + 1, x, 0).evalf(2), RootOf)

    assert RootOf(x**3 + I*x + 2, 0).evalf(7) == (Float('-1.260785326', dps=7) +
                                                  I*Float('0.2684419416', dps=7))

    r = RootOf(x**2 - 4456178*x + 60372201703370, 0, radicals=False)
    assert r.evalf(2) == Float('2.2282e+6', dps=2) - I*Float('7.4465e+6', dps=2)


def test_RootOf_evalf_caching_bug():
    r = RootOf(x**5 - 5*x + 12, 1)
    r.evalf()
    a = r.interval
    r = RootOf(x**5 - 5*x + 12, 1)
    r.evalf()
    b = r.interval
    assert a == b


def test_RootOf_real_roots():
    assert (x**5 + x + 1).as_poly().real_roots() == [RootOf(x**3 - x**2 + 1, 0)]
    assert (x**5 + x + 1).as_poly().real_roots(radicals=False) == [RootOf(x**3 - x**2 + 1, 0)]
    assert (x**7 - 0.1*x + 1).as_poly(x).real_roots() == [RootOf(10*x**7 - x + 10, 0)]
    assert ((x - 2)*(x - 1)*(x**2 - 2)).as_poly().real_roots() == [-sqrt(2), 1, sqrt(2), 2]
    assert ((x + 1)**3*(3*x + 1)).as_poly().real_roots(multiple=False) == [(-1, 3), (Rational(-1, 3), 1)]


def test_RootOf_all_roots():
    assert (x**5 + x + 1).as_poly().all_roots() == [
        RootOf(x**3 - x**2 + 1, 0),
        -Rational(1, 2) - sqrt(3)*I/2,
        -Rational(1, 2) + sqrt(3)*I/2,
        RootOf(x**3 - x**2 + 1, 1),
        RootOf(x**3 - x**2 + 1, 2),
    ]

    assert (x**5 + x + 1).as_poly().all_roots(radicals=False) == [
        RootOf(x**3 - x**2 + 1, 0),
        RootOf(x**2 + x + 1, 0, radicals=False),
        RootOf(x**2 + x + 1, 1, radicals=False),
        RootOf(x**3 - x**2 + 1, 1),
        RootOf(x**3 - x**2 + 1, 2),
    ]

    r = ((x**3 + x + 20)*(x**3 + x + 21)).as_poly().all_roots()

    assert r[0].is_real
    assert r[1].is_real
    assert all(not _.is_real for _ in r[2:])

    assert r == [RootOf(x**3 + x + 21, 0), RootOf(x**3 + x + 20, 0),
                 RootOf(x**3 + x + 20, 1), RootOf(x**3 + x + 20, 2),
                 RootOf(x**3 + x + 21, 1), RootOf(x**3 + x + 21, 2)]


def test_RootOf_eval_rational():
    p = legendre_poly(4, x, polys=True)
    roots = [r.eval_rational(Rational(1, 10)**20) for r in p.real_roots()]
    for r in roots:
        assert isinstance(r, Rational)
    # All we know is that the Rational instance will be at most 1/10^20 from
    # the exact root. So if we evaluate to 17 digits, it must be exactly equal
    # to:
    roots = [str(r.evalf(17)) for r in roots]
    assert roots == [
        '-0.86113631159405258',
        '-0.33998104358485626',
        '0.33998104358485626',
        '0.86113631159405258',
    ]

    pytest.raises(NotImplementedError,
                  lambda: RootOf(x**3 + x + 3, 1).eval_rational(1e-3))


def test_RootSum___new__():
    f = x**3 + x + 3

    g = Lambda(r, log(r*x))
    s = RootSum(f, g)

    assert isinstance(s, RootSum) is True

    assert RootSum(f**2, g) == 2*RootSum(f, g)
    assert RootSum((x - 7)*f**3, g) == log(7*x) + 3*RootSum(f, g)

    # issue sympy/sympy#5571
    assert hash(RootSum((x - 7)*f**3, g)) == hash(log(7*x) + 3*RootSum(f, g))

    pytest.raises(MultivariatePolynomialError, lambda: RootSum(x**3 + x + y))
    pytest.raises(ValueError, lambda: RootSum(x**2 + 3, lambda x: x))

    assert RootSum(f, log) == RootSum(f, Lambda(x, log(x)))

    assert isinstance(RootSum(f, auto=False), RootSum) is True

    assert RootSum(f) == 0
    assert RootSum(f, Lambda(x, x)) == 0
    assert RootSum(f, Lambda(x, x**2)) == -2

    assert RootSum(f, Lambda(x, 1)) == 3
    assert RootSum(f, Lambda(x, 2)) == 6

    assert RootSum(f, auto=False).is_commutative is True

    assert RootSum(f, Lambda(x, 1/(x + x**2))) == Rational(11, 3)
    assert RootSum(f, Lambda(x, y/(x + x**2))) == Rational(11, 3)*y

    assert RootSum(x**2 - 1, Lambda(x, 3*x**2), x) == 6
    assert RootSum(x**2 - y, Lambda(x, 3*x**2), x) == 6*y

    assert RootSum(x**2 - 1, Lambda(x, z*x**2), x) == 2*z
    assert RootSum(x**2 - y, Lambda(x, z*x**2), x) == 2*z*y

    assert RootSum(
        x**2 - 1, Lambda(x, exp(x)), quadratic=True) == exp(-1) + exp(1)

    assert RootSum(x**3 + a*x + a**3, tan, x) == \
        RootSum(x**3 + x + 1, Lambda(x, tan(a*x)))
    assert RootSum(a**3*x**3 + a*x + 1, tan, x) == \
        RootSum(x**3 + x + 1, Lambda(x, tan(x/a)))

    assert isinstance(RootSum(x**7 + 2*x + 1,
                              Lambda(x, log(x))).doit(),
                      RootSum)


def test_RootSum_free_symbols():
    assert RootSum(x**3 + x + 3, Lambda(r, exp(r))).free_symbols == set()
    assert RootSum(x**3 + x + 3, Lambda(r, exp(a*r))).free_symbols == {a}
    assert RootSum(
        x**3 + x + y, Lambda(r, exp(a*r)), x).free_symbols == {a, y}


def test_RootSum___eq__():
    f = Lambda(x, exp(x))

    assert (RootSum(x**3 + x + 1, f) == RootSum(x**3 + x + 1, f)) is True
    assert (RootSum(x**3 + x + 1, f) == RootSum(y**3 + y + 1, f)) is True

    assert (RootSum(x**3 + x + 1, f) == RootSum(x**3 + x + 2, f)) is False
    assert (RootSum(x**3 + x + 1, f) == RootSum(y**3 + y + 2, f)) is False


def test_RootSum_doit():
    rs = RootSum(x**2 + 1, Lambda(x, exp(x)))

    assert isinstance(rs, RootSum) is True
    assert rs.doit() == exp(-I) + exp(I)

    rs = RootSum(x**2 + a, Lambda(x, exp(x)), x)

    assert isinstance(rs, RootSum) is True
    assert rs.doit() == exp(-sqrt(-a)) + exp(sqrt(-a))


def test_RootSum_evalf():
    rs = RootSum(x**2 + 1, Lambda(x, exp(x)))

    assert rs.evalf(20, chop=True).epsilon_eq(
        Float('1.0806046117362794348', 20), Float('1e-20')) is true
    assert rs.evalf(15, chop=True).epsilon_eq(
        Float('1.08060461173628', 15), Float('1e-15')) is true

    rs = RootSum(x**2 + a, Lambda(x, exp(x)), x)

    assert rs.evalf() == rs


def test_RootSum_diff():
    f = x**3 + x + 3

    g = Lambda(r, exp(r*x))
    h = Lambda(r, r*exp(r*x))

    assert RootSum(f, g).diff(x) == RootSum(f, h)


def test_RootSum_subs():
    f = x**3 + x + 3
    g = Lambda(r, exp(r*x))

    F = y**3 + y + 3
    G = Lambda(r, exp(r*y))

    assert RootSum(f, g).subs({y: 1}) == RootSum(f, g)
    assert RootSum(f, g).subs({x: y}) == RootSum(F, G)


def test_RootSum_rational():
    assert RootSum(
        z**5 - z + 1, Lambda(z, z/(x - z))) == (4*x - 5)/(x**5 - x + 1)

    f = 161*z**3 + 115*z**2 + 19*z + 1
    g = Lambda(z, z*log(
        -3381*z**4/4 - 3381*z**3/4 - 625*z**2/2 - 125*z/2 - 5 + exp(x)))

    assert RootSum(f, g).diff(x) == -(
        (5*exp(2*x) - 6*exp(x) + 4)*exp(x)/(exp(3*x) - exp(2*x) + 1))/7


def test_RootSum_independent():
    f = (x**3 - a)**2*(x**4 - b)**3

    g = Lambda(x, 5*tan(x) + 7)
    h = Lambda(x, tan(x))

    r0 = RootSum(x**3 - a, h, x)
    r1 = RootSum(x**4 - b, h, x)

    assert RootSum(f, g, x).as_ordered_terms() == [10*r0, 15*r1, 126]


def test_sympyissue_7876():
    l1 = (x**6 - x + 1).as_poly().all_roots()
    l2 = [RootOf(x**6 - x + 1, i) for i in range(6)]
    assert frozenset(l1) == frozenset(l2)


def test_sympyissue_8316():
    f = (7*x**8 - 9).as_poly()
    assert len(f.all_roots()) == 8
    f = (7*x**8 - 10).as_poly()
    assert len(f.all_roots()) == 8


def test_rewrite():
    r3 = RootOf(x**3 + x - 1, 0)
    assert r3.evalf() == r3.rewrite(Pow).evalf()
    assert r3.rewrite(Pow) == (-1/(3*root(Rational(1, 2) + sqrt(93)/18, 3)) +
                               root(Rational(1, 2) + sqrt(93)/18, 3))
    r4 = RootOf(x**4 - x + 5, 0)
    assert r4.evalf() == r4.rewrite(Pow).evalf()
    r11 = RootOf(x**11 + x - 3, 0)
    assert r11.rewrite(Pow) == r11


def test_RootOf_expand_func1():
    r0 = RootOf(x**3 + x + 1, 0)
    assert expand_func(r0) == r0
    r1 = RootOf(x**3 - sqrt(2)*x + I, 1)
    assert expand_func(r1) == RootOf(x**12 - 4*x**8 + 2*x**6 +
                                     4*x**4 + 4*x**2 + 1, 7)


@pytest.mark.slow
def test_RootOf_expand_func2():
    r0 = RootOf(x**3 + I*x + 2, 0)
    assert expand_func(r0) == RootOf(x**6 + 4*x**3 + x**2 + 4, 1)
    r1 = RootOf(x**3 + I*x + 2, 1)
    assert expand_func(r1) == RootOf(x**6 + 4*x**3 + x**2 + 4, 3)
    r2 = RootOf(x**4 + sqrt(2)*x**3 - I*x + 1, 0)
    assert expand_func(r2) == RootOf(x**16 - 4*x**14 + 8*x**12 - 6*x**10 +
                                     10*x**8 + 5*x**4 + 2*x**2 + 1, 1)
    r3 = RootOf(x**3 - I*sqrt(2)*x + 5, 1)
    assert expand_func(r3) == RootOf(x**6 + 10*x**3 + 2*x**2 + 25, 2)


@pytest.mark.slow
def test_RootOf_algebraic():
    e = RootOf(sqrt(2)*x**4 + sqrt(2)*x**3 - I*x + sqrt(2), x, 0)
    assert e.interval.as_tuple() == ((Rational(-201, 100), 0),
                                     (Rational(-201, 200), Rational(201, 200)))
    assert e.evalf(7) == Float('-1.22731258', dps=7) + I*Float('0.6094138324', dps=7)

    t = RootOf(x**5 + 4*x + 2, 0)
    e = RootOf(x**4 + t*x + 1, 0)
    assert e.interval.as_tuple() == ((Rational(-201, 200), Rational(-201, 200)),
                                     (Rational(-201, 400), Rational(-201, 400)))
    assert e.evalf(7) == Float('-0.7123350278', dps=7) - I*Float('0.8248345032', dps=7)


@pytest.mark.timeout(10)
def test_issue_730():
    e = RootOf(x**3 + 10*x**2 + 1, 2)
    assert e.is_real is False
    assert e.is_imaginary is False
    assert e.evalf(3) == Float('0.00498962', dps=3) + I*Float('0.31604', dps=3)
    assert e.conjugate().conjugate() == e


@pytest.mark.timeout(150)
@pytest.mark.slow
def test_issue_723():
    p = x**5 + sqrt(3)*x - 2
    for _ in range(20):
        for j in (1, 2):
            RootOf(p, j)


def test_sympyissue_15413():
    assert (sqrt(2)*x**3 + x).as_poly(x).all_roots() == [0, -I*root(2, -4),
                                                         I*root(2, -4)]
