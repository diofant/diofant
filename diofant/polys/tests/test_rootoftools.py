"""Tests for the implementation of RootOf class and related tools. """

import pytest

from diofant import (Eq, Float, Function, I, Lambda, Pow, Rational, Symbol,
                     exp, false, legendre_poly, log, oo, root, solve, sqrt,
                     tan, true)
from diofant.abc import a, b, r, x, y, z
from diofant.polys.polyerrors import (GeneratorsNeeded,
                                      MultivariatePolynomialError,
                                      PolynomialError)
from diofant.polys.polytools import Poly, PurePoly
from diofant.polys.rootoftools import RootOf, RootSum


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

    pytest.raises(GeneratorsNeeded, lambda: RootOf(0, 0))
    pytest.raises(GeneratorsNeeded, lambda: RootOf(1, 0))

    pytest.raises(PolynomialError, lambda: RootOf(Poly(0, x), 0))
    pytest.raises(PolynomialError, lambda: RootOf(Poly(1, x), 0))

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
                  lambda: Poly(Symbol('a', nonzero=False)*x**5 +
                               2*x - 1, x).all_roots())

    assert RootOf(Poly(x - y, x), 0) == y

    assert RootOf(Poly(x**2 - y, x), 0) == -sqrt(y)
    assert RootOf(Poly(x**2 - y, x), 1) == sqrt(y)

    assert isinstance(RootOf(x**3 - y, x, 0), RootOf)
    p = Symbol('p', positive=True)
    assert RootOf(x**3 - p, x, 0) == root(p, 3)*RootOf(x**3 - 1, x, 0)

    assert RootOf(y*x**3 + y*x + 2*y, x, 0) == -1

    assert RootOf(x**3 + x + 1, 0).is_commutative is True

    e = RootOf(x**2 - 4, x, 1, evaluate=False)
    assert isinstance(e, RootOf)
    assert e.doit() == 2
    assert e.args == (x**2 - 4, x, 1)
    assert e.poly == PurePoly(x**2 - 4, x)
    assert e.index == 1

    assert RootOf(x**7 - 0.1*x + 1, x, 0) == RootOf(10*x**7 - x + 10, x, 0)


def test_RootOf_attributes():
    r = RootOf(x**3 + x + 3, 0)
    assert r.is_number
    assert r.free_symbols == set()

    r = RootOf(x**3 + y*x + 1, x, 0)
    assert isinstance(r, RootOf) and r.expr == x**3 + y*x + 1
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


def test_RootOf_is_complex():
    assert RootOf(x**3 + x + 3, 0).is_complex is True
    assert RootOf(x**3 + y*x + 3, x, 0).is_complex is None

    assert RootOf(x**3 + y*x + 3, x, 0).is_commutative


def test_RootOf_is_algebraic():
    assert RootOf(x**3 + x + 3, 0).is_algebraic is True
    assert RootOf(x**3 + y*x + 3, x, 0).is_algebraic is None


def test_RootOf_power():
    e = RootOf(y**3 - x, y, 0)
    assert e**3 == x
    assert e**2 == Pow(e, 2, evaluate=False)
    e2 = RootOf(y**3 - x*y, y, 0)
    assert e2**3 == Pow(e2, 3, evaluate=False)


def test_RootOf_subs():
    assert RootOf(x**3 + x + 1, 0).subs(x, y) == RootOf(y**3 + y + 1, 0)
    eq = -x + RootOf(y**3 - x**3 + 3*x**2, y, 0) + 1
    assert eq.subs(x, Rational(1, 3)) == 0


def test_RootOf_diff():
    assert RootOf(x**3 + x + 1, 0).diff(x) == 0
    assert RootOf(x**3 + x + 1, 0).diff(y) == 0


def test_RootOf_evalf():
    real = RootOf(x**3 + x + 3, 0).evalf(n=20)

    assert real.epsilon_eq(Float("-1.2134116627622296341"))

    re, im = RootOf(x**3 + x + 3, 1).evalf(n=20).as_real_imag()

    assert re.epsilon_eq( Float("0.60670583138111481707"))
    assert im.epsilon_eq(-Float("1.45061224918844152650"))

    re, im = RootOf(x**3 + x + 3, 2).evalf(n=20).as_real_imag()

    assert re.epsilon_eq(Float("0.60670583138111481707"))
    assert im.epsilon_eq(Float("1.45061224918844152650"))

    p = legendre_poly(4, x, polys=True)
    roots = [str(r.n(17)) for r in p.real_roots()]
    assert roots == [
        "-0.86113631159405258",
        "-0.33998104358485626",
        "0.33998104358485626",
        "0.86113631159405258",
    ]

    re = RootOf(x**5 - 5*x + 12, 0).evalf(n=20)
    assert re.epsilon_eq(Float("-1.84208596619025438271"))

    re, im = RootOf(x**5 - 5*x + 12, 1).evalf(n=20).as_real_imag()
    assert re.epsilon_eq(Float("-0.351854240827371999559"))
    assert im.epsilon_eq(Float("-1.709561043370328882010"))

    re, im = RootOf(x**5 - 5*x + 12, 2).evalf(n=20).as_real_imag()
    assert re.epsilon_eq(Float("-0.351854240827371999559"))
    assert im.epsilon_eq(Float("+1.709561043370328882010"))

    re, im = RootOf(x**5 - 5*x + 12, 3).evalf(n=20).as_real_imag()
    assert re.epsilon_eq(Float("+1.272897223922499190910"))
    assert im.epsilon_eq(Float("-0.719798681483861386681"))

    re, im = RootOf(x**5 - 5*x + 12, 4).evalf(n=20).as_real_imag()
    assert re.epsilon_eq(Float("+1.272897223922499190910"))
    assert im.epsilon_eq(Float("+0.719798681483861386681"))

    # issue sympy/sympy#6393
    assert str(RootOf(x**5 + 2*x**4 + x**3 - 68719476736, 0).n(3)) == '147.'
    eq = (531441*x**11 + 3857868*x**10 + 13730229*x**9 + 32597882*x**8 +
          55077472*x**7 + 60452000*x**6 + 32172064*x**5 - 4383808*x**4 -
          11942912*x**3 - 1506304*x**2 + 1453312*x + 512)
    a, b = RootOf(eq, 1).n(2).as_real_imag()
    c, d = RootOf(eq, 2).n(2).as_real_imag()
    assert a == c
    assert b < d
    assert b == -d
    # issue sympy/sympy#6451
    r = RootOf(legendre_poly(64, x), 7)
    assert r.n(2) == r.n(100).n(2)
    # issue sympy/sympy#8617
    ans = [w[x].n(2) for w in solve(x**3 - x - 4)]
    assert RootOf(exp(x)**3 - exp(x) - 4, 0).n(2) in ans
    # issue sympy/sympy#9019
    r0 = RootOf(x**2 + 1, 0, radicals=False)
    r1 = RootOf(x**2 + 1, 1, radicals=False)
    assert r0.n(4) == -1.0*I
    assert r1.n(4) == 1.0*I

    # make sure verification is used in case a max/min traps the "root"
    assert str(RootOf(4*x**5 + 16*x**3 + 12*x**2 + 7, 0).n(3)) == '-0.976'

    assert isinstance(RootOf(x**3 + y*x + 1, x, 0).n(2), RootOf)


def test_RootOf_evalf_caching_bug():
    r = RootOf(x**5 - 5*x + 12, 1)
    r.n()
    a = r.interval
    r = RootOf(x**5 - 5*x + 12, 1)
    r.n()
    b = r.interval
    assert a == b


def test_RootOf_real_roots():
    assert Poly(x**5 + x + 1).real_roots() == [RootOf(x**3 - x**2 + 1, 0)]
    assert Poly(x**5 + x + 1).real_roots(radicals=False) == [RootOf(
        x**3 - x**2 + 1, 0)]
    assert Poly(x**7 - 0.1*x + 1, x).real_roots() == [RootOf(10*x**7 - x + 10, x, 0)]


def test_RootOf_all_roots():
    assert Poly(x**5 + x + 1).all_roots() == [
        RootOf(x**3 - x**2 + 1, 0),
        -Rational(1, 2) - sqrt(3)*I/2,
        -Rational(1, 2) + sqrt(3)*I/2,
        RootOf(x**3 - x**2 + 1, 1),
        RootOf(x**3 - x**2 + 1, 2),
    ]

    assert Poly(x**5 + x + 1).all_roots(radicals=False) == [
        RootOf(x**3 - x**2 + 1, 0),
        RootOf(x**2 + x + 1, 0, radicals=False),
        RootOf(x**2 + x + 1, 1, radicals=False),
        RootOf(x**3 - x**2 + 1, 1),
        RootOf(x**3 - x**2 + 1, 2),
    ]


def test_RootOf_eval_rational():
    p = legendre_poly(4, x, polys=True)
    roots = [r.eval_rational(Rational(1, 10)**20) for r in p.real_roots()]
    for r in roots:
        assert isinstance(r, Rational)
    # All we know is that the Rational instance will be at most 1/10^20 from
    # the exact root. So if we evaluate to 17 digits, it must be exactly equal
    # to:
    roots = [str(r.n(17)) for r in roots]
    assert roots == [
        "-0.86113631159405258",
        "-0.33998104358485626",
        "0.33998104358485626",
        "0.86113631159405258",
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

    assert rs.evalf(n=20, chop=True).epsilon_eq(
        Float("1.0806046117362794348", 20), Float("1e-20")) is true
    assert rs.evalf(n=15, chop=True).epsilon_eq(
        Float("1.08060461173628", 15), Float("1e-15")) is true

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

    assert RootSum(f, g).subs(y, 1) == RootSum(f, g)
    assert RootSum(f, g).subs(x, y) == RootSum(F, G)


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
    l1 = Poly(x**6 - x + 1, x).all_roots()
    l2 = [RootOf(x**6 - x + 1, i) for i in range(6)]
    assert frozenset(l1) == frozenset(l2)


def test_sympyissue_8316():
    f = Poly(7*x**8 - 9)
    assert len(f.all_roots()) == 8
    f = Poly(7*x**8 - 10)
    assert len(f.all_roots()) == 8
