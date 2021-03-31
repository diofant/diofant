"""Limit computation tests."""

import itertools

import pytest

from diofant import (E, Float, Function, I, Integral, Limit, Piecewise,
                     PoleError, Rational, Sum, Symbol, acos, asin, atan,
                     besselk, binomial, cbrt, ceiling, cos, cosh, cot, diff,
                     digamma, erf, erfi, exp, factorial, floor, gamma,
                     integrate, limit, log, nan, oo, pi, polygamma, root, sign,
                     simplify, sin, sinh, sqrt, subfactorial, symbols, tan)
from diofant.abc import a, b, c, n, x, y, z
from diofant.series.limits import heuristics
from diofant.series.order import O


__all__ = ()


def test_basic1():
    assert limit(x, x, oo) == oo
    assert limit(x, x, -oo) == -oo
    assert limit(-x, x, oo) == -oo
    assert limit(x**2, x, -oo) == oo
    assert limit(-x**2, x, oo) == -oo
    assert limit(x*log(x), x, 0, dir='+') == 0
    assert limit(1/x, x, oo) == 0  # issue sympy/sympy#11667
    assert limit(exp(x), x, oo) == oo
    assert limit(-exp(x), x, oo) == -oo
    assert limit(exp(x)/x, x, oo) == oo
    assert limit(1/x - exp(-x), x, oo) == 0
    assert limit(x + 1/x, x, oo) == oo
    assert limit(x - x**2, x, oo) == -oo
    assert limit((1 + x)**(1 + sqrt(2)), x, 0) == 1
    assert limit((1 + x)**oo, x, 0) == oo
    assert limit((1 + x)**oo, x, 0, dir='-') == 0
    assert limit((1 + x + y)**oo, x, 0, dir='-') == (1 + y)**oo
    assert limit(y/x/log(x), x, 0) == -oo*sign(y)
    assert limit(cos(x + y)/x, x, 0) == sign(cos(y))*oo
    limit(Sum(1/x, (x, 1, y)) - log(y), y, oo)
    limit(Sum(1/x, (x, 1, y)) - 1/y, y, oo)
    assert limit(gamma(1/x + 3), x, oo) == 2
    assert limit(nan, x, -oo) == nan
    assert limit(O(2)*x, x, nan) == nan
    assert limit(sin(O(x)), x, 0) == 0
    assert limit(1/(x - 1), x, 1, dir='+') == oo
    assert limit(1/(x - 1), x, 1, dir='-') == -oo
    assert limit(1/(5 - x)**3, x, 5, dir='+') == -oo
    assert limit(1/(5 - x)**3, x, 5, dir='-') == oo
    assert limit(1/sin(x), x, pi, dir='+') == -oo
    assert limit(1/sin(x), x, pi, dir='-') == oo
    assert limit(1/cos(x), x, pi/2, dir='+') == -oo
    assert limit(1/cos(x), x, pi/2, dir='-') == oo
    assert limit(1/tan(x**3), x, cbrt(2*pi), dir='+') == oo
    assert limit(1/tan(x**3), x, cbrt(2*pi), dir='-') == -oo
    assert limit(1/cot(x)**3, x, 3*pi/2, dir='+') == -oo
    assert limit(1/cot(x)**3, x, 3*pi/2, dir='-') == oo

    # approaching 0
    # from dir='+'
    assert limit(1 + 1/x, x, 0) == oo
    # from dir='-'
    # Add
    assert limit(1 + 1/x, x, 0, dir='-') == -oo
    # Pow
    assert limit(x**(-2), x, 0, dir='-') == oo
    assert limit(x**(-3), x, 0, dir='-') == -oo
    assert limit(1/sqrt(x), x, 0, dir='-') == (-oo)*I
    assert limit(x**2, x, 0, dir='-') == 0
    assert limit(sqrt(x), x, 0, dir='-') == 0
    assert limit(x**-pi, x, 0, dir='-') == oo*sign((-1)**(-pi))
    assert limit((1 + cos(x))**oo, x, 0) == oo

    assert limit(x**2, x, 0, dir='real') == 0
    assert limit(exp(x), x, 0, dir='real') == 1

    pytest.raises(PoleError, lambda: limit(1/x, x, 0, dir='real'))

    # issue diofant/diofant#74
    assert limit(sign(log(1 - 1/x)), x, oo) == -1


def test_basic2():
    assert limit(x**x, x, 0, dir='+') == 1
    assert limit((exp(x) - 1)/x, x, 0) == 1
    assert limit(1 + 1/x, x, oo) == 1
    assert limit(-exp(1/x), x, oo) == -1
    assert limit(x + exp(-x), x, oo) == oo
    assert limit(x + exp(-x**2), x, oo) == oo
    assert limit(x + exp(-exp(x)), x, oo) == oo
    assert limit(13 + 1/x - exp(-x), x, oo) == 13


def test_basic3():
    assert limit(1/x, x, 0, dir='+') == oo
    assert limit(1/x, x, 0, dir='-') == -oo


def test_basic4():
    assert limit(2*x + y*x, x, 0) == 0
    assert limit(2*x + y*x, x, 1) == 2 + y
    assert limit(2*x**8 + y*x**(-3), x, -2) == 512 - y/8
    assert limit(sqrt(x + 1) - sqrt(x), x, oo) == 0
    assert integrate(1/(x**3 + 1), (x, 0, oo)) == 2*pi*sqrt(3)/9

    # coverage test
    l = Limit(Piecewise((x, x > 1), (0, True)), x, -1)
    assert l.doit() == l

    # issue sympy/sympy#16714
    e = ((n**(n + 1) + (n + 1)**n)/n**(n + 1))**n
    assert limit(e, n, oo) == E**E


def test_basic5():
    class MyFunction(Function):
        @classmethod
        def eval(cls, arg):
            if arg is oo:
                return nan
    assert limit(MyFunction(x), x, oo) == Limit(MyFunction(x), x, oo)

    assert limit(4/x > 8, x, 0)  # relational test
    assert limit(MyFunction(x) > 0, x, oo) == Limit(MyFunction(x) > 0, x, oo)

    # from https://groups.google.com/forum/#!topic/sympy/LkTMQKC_BOw
    # fix bisected to ade6d20 and c459d18
    a = Symbol('a', positive=True)
    f = exp(x*(-a - 1)) * sinh(x)
    assert limit(f, x, oo) == 0

    assert limit(O(x), x, x**2) == Limit(O(x), x, x**2)


def test_sympyissue_3885():
    assert limit(x*y + x*z, z, 2) == x*y + 2*x


def test_Limit():
    assert Limit(sin(x)/x, x, 0) != 1
    assert Limit(sin(x)/x, x, 0).doit() == 1


def test_floor():
    assert limit(floor(x), x, -2, '+') == -2
    assert limit(floor(x), x, -2, '-') == -3
    assert limit(floor(x), x, -1, '+') == -1
    assert limit(floor(x), x, -1, '-') == -2
    assert limit(floor(x), x, 0, '+') == 0
    assert limit(floor(x), x, 0, '-') == -1
    assert limit(floor(x), x, 1, '+') == 1
    assert limit(floor(x), x, 1, '-') == 0
    assert limit(floor(x), x, 2, '+') == 2
    assert limit(floor(x), x, 2, '-') == 1
    assert limit(floor(x), x, 248, '+') == 248
    assert limit(floor(x), x, 248, '-') == 247


def test_floor_requires_robust_assumptions():
    assert limit(floor(sin(x)), x, 0, '+') == 0
    assert limit(floor(sin(x)), x, 0, '-') == -1
    assert limit(floor(cos(x)), x, 0, '+') == 0
    assert limit(floor(cos(x)), x, 0, '-') == 0
    assert limit(floor(5 + sin(x)), x, 0, '+') == 5
    assert limit(floor(5 + sin(x)), x, 0, '-') == 4
    assert limit(floor(5 + cos(x)), x, 0, '+') == 5
    assert limit(floor(5 + cos(x)), x, 0, '-') == 5


def test_ceiling():
    assert limit(ceiling(x), x, -2, '+') == -1
    assert limit(ceiling(x), x, -2, '-') == -2
    assert limit(ceiling(x), x, -1, '+') == 0
    assert limit(ceiling(x), x, -1, '-') == -1
    assert limit(ceiling(x), x, 0, '+') == 1
    assert limit(ceiling(x), x, 0, '-') == 0
    assert limit(ceiling(x), x, 1, '+') == 2
    assert limit(ceiling(x), x, 1, '-') == 1
    assert limit(ceiling(x), x, 2, '+') == 3
    assert limit(ceiling(x), x, 2, '-') == 2
    assert limit(ceiling(x), x, 248, '+') == 249
    assert limit(ceiling(x), x, 248, '-') == 248


def test_ceiling_requires_robust_assumptions():
    assert limit(ceiling(sin(x)), x, 0, '+') == 1
    assert limit(ceiling(sin(x)), x, 0, '-') == 0
    assert limit(ceiling(cos(x)), x, 0, '+') == 1
    assert limit(ceiling(cos(x)), x, 0, '-') == 1
    assert limit(ceiling(5 + sin(x)), x, 0, '+') == 6
    assert limit(ceiling(5 + sin(x)), x, 0, '-') == 5
    assert limit(ceiling(5 + cos(x)), x, 0, '+') == 6
    assert limit(ceiling(5 + cos(x)), x, 0, '-') == 6


def test_atan():
    assert limit(atan(x)*sin(1/x), x, 0) == 0
    assert limit(atan(x) + sqrt(x + 1) - sqrt(x), x, oo) == pi/2


def test_abs():
    assert limit(abs(x), x, 0) == 0
    assert limit(abs(sin(x)), x, 0) == 0
    assert limit(abs(cos(x)), x, 0) == 1
    assert limit(abs(sin(x + 1)), x, 0) == sin(1)

    # sympy/sympy#12398
    assert limit(abs(log(n)/n**3), n, oo) == 0
    expr = abs(log(n)/n**3)
    expr2 = expr.subs({n: n + 1})
    assert limit(n*(expr/expr2 - 1), n, oo) == 3


def test_heuristic():
    x = Symbol('x', extended_real=True)
    assert heuristics(sin(1/x) + atan(x), x, 0, '+') == sin(oo)
    assert heuristics(log(2 + sqrt(atan(x))*sin(1/x)), x, 0, '+') == log(2)
    assert heuristics(tan(tan(1/x)), x, 0, '+') is None


def test_sympyissue_3871():
    z = Symbol('z', positive=True)
    f = -1/z*exp(-z*x)
    assert limit(f, x, oo) == 0
    assert f.limit(x, oo) == 0


def test_exponential():
    x = Symbol('x', extended_real=True)
    assert limit((1 + x/n)**n, n, oo) == exp(x)
    assert limit((1 + x/(2*n))**n, n, oo) == exp(x/2)
    assert limit((1 + x/(2*n + 1))**n, n, oo) == exp(x/2)
    assert limit(((x - 1)/(x + 1))**x, x, oo) == exp(-2)
    assert limit(1 + (1 + 1/x)**x, x, oo) == 1 + E


@pytest.mark.xfail
def test_exponential2():
    assert limit((1 + x/(n + sin(n)))**n, n, oo) == exp(x)


def test_doit():
    f = Integral(2 * x, x)
    l = Limit(f, x, oo)
    assert l.doit() == oo


def test_doit2():
    f = Integral(2 * x, x)
    l = Limit(f, x, oo)
    # limit() breaks on the contained Integral.
    assert l.doit(deep=False) == l


def test_sympyissue_3792():
    assert limit((1 - cos(x))/x**2, x, Rational(1, 2)) == 4 - 4*cos(Rational(1, 2))
    assert limit(sin(sin(x + 1) + 1), x, 0) == sin(1 + sin(1))
    assert limit(abs(sin(x + 1) + 1), x, 0) == 1 + sin(1)


def test_sympyissue_4090():
    assert limit(1/(x + 3), x, 2) == Rational(1, 5)
    assert limit(1/(x + pi), x, 2) == 1/(2 + pi)
    assert limit(log(x)/(x**2 + 3), x, 2) == log(2)/7
    assert limit(log(x)/(x**2 + pi), x, 2) == log(2)/(4 + pi)


def test_sympyissue_4547():
    assert limit(cot(x), x, 0, dir='+') == oo
    assert limit(cot(x), x, pi/2, dir='+') == 0


def test_sympyissue_5164():
    assert limit(x**0.5, x, oo) == oo**0.5 == oo
    assert limit(x**0.5, x, 16) == 4.0
    assert limit(x**0.5, x, 0) == 0
    assert limit(x**(-0.5), x, oo) == 0
    assert limit(x**(-0.5), x, 4) == 0.5


def test_sympyissue_5183():
    # using list(...) so pytest can recalculate values
    tests = list(itertools.product([x, -x],
                                   [-1, 1],
                                   [2, 3, Rational(1, 2), Rational(2, 3)],
                                   ['-', '+']))
    results = (oo, oo, -oo, oo, -oo*I, oo, -oo*sign(cbrt(-1)), oo,
               0, 0, 0, 0, 0, 0, 0, 0,
               oo, oo, oo, -oo, oo, -oo*I, oo, -oo*sign(cbrt(-1)),
               0, 0, 0, 0, 0, 0, 0, 0)
    assert len(tests) == len(results)
    for i, (args, res) in enumerate(zip(tests, results)):
        y, s, e, d = args
        eq = y**(s*e)
        assert limit(eq, x, 0, dir=d) == res


def test_sympyissue_5184():
    assert limit(sin(x)/x, x, oo) == 0
    assert limit(atan(x), x, oo) == pi/2
    assert limit(gamma(x), x, oo) == oo
    assert limit(cos(x)/x, x, oo) == 0
    assert limit(gamma(x), x, Rational(1, 2)) == sqrt(pi)

    r = Symbol('r', real=True)
    assert limit(r*sin(1/r), r, 0) == 0


def test_sympyissue_5229():
    assert limit((1 + y)**(1/y) - E, y, 0) == 0


def test_sympyissue_4546():
    # using list(...) so pytest can recalculate values
    tests = list(itertools.product([cot, tan],
                                   [-pi/2, 0, pi/2, pi, 3*pi/2],
                                   ['-', '+']))
    results = (0, 0, -oo, oo, 0, 0, -oo, oo, 0, 0,
               oo, -oo, 0, 0, oo, -oo, 0, 0, oo, -oo)
    assert len(tests) == len(results)
    for i, (args, res) in enumerate(zip(tests, results)):
        f, l, d = args
        eq = f(x)
        assert limit(eq, x, l, dir=d) == res


def test_sympyissue_3934():
    assert limit((1 + x**log(3))**(1/x), x, 0) == 1
    assert limit((5**(1/x) + 3**(1/x))**x, x, 0) == 5


def test_compute_leading_term():
    assert limit(x**Rational(77, 3)/(1 + x**Rational(77, 3)), x, oo) == 1
    assert limit(x**Rational('101.1')/(1 + x**Rational('101.1')), x, oo) == 1


def test_sympyissue_5955():
    assert limit((x**16)/(1 + x**16), x, oo) == 1
    assert limit((x**100)/(1 + x**100), x, oo) == 1
    assert limit((x**1885)/(1 + x**1885), x, oo) == 1
    assert limit((x**100/((x + 1)**100 + exp(-x))), x, oo) == 1


def test_newissue():
    assert limit(exp(1/sin(x))/exp(cot(x)), x, 0) == 1


def test_extended_real_line():
    assert limit(x - oo, x, oo) == -oo
    assert limit(oo - x, x, -oo) == oo
    assert limit(x**2/(x - 5) - oo, x, oo) == -oo
    assert limit(1/(x + sin(x)) - oo, x, 0) == -oo
    assert limit(oo/x, x, oo) == oo
    assert limit(x - oo + 1/x, x, oo) == -oo
    assert limit(x - oo + 1/x, x, 0) == -oo


@pytest.mark.xfail
def test_order_oo():
    x = Symbol('x', positive=True, finite=True)
    assert O(x)*oo != O(1, x)
    assert limit(oo/(x**2 - 4), x, oo) == oo


def test_sympyissue_5436():
    # also issue sympy/sympy#13312 (but see diofant/diofant#425!)
    assert limit(exp(x*y), x, oo) == exp(oo*sign(y))
    assert limit(exp(-x*y), x, oo) == exp(-oo*sign(y))


def test_Limit_dir():
    pytest.raises(TypeError, lambda: Limit(x, x, 0, dir=0))
    pytest.raises(ValueError, lambda: Limit(x, x, 0, dir='0'))


def test_polynomial():
    assert limit((x + 1)**1000/((x + 1)**1000 + 1), x, oo) == 1
    assert limit((x + 1)**1000/((x + 1)**1000 + 1), x, -oo) == 1


def test_rational():
    assert limit(1/y - (1/(y + x) + x/(y + x)/y)/z, x, oo) == (z - 1)/(y*z)
    assert limit(1/y - (1/(y + x) + x/(y + x)/y)/z, x, -oo) == (z - 1)/(y*z)


def test_sympyissue_5740():
    assert limit(log(x)*z - log(2*x)*y, x, 0) == oo*sign(y - z)


def test_sympyissue_6366():
    n = Symbol('n', integer=True, positive=True)
    r = (n + 1)*x**(n + 1)/(x**(n + 1) - 1) - x/(x - 1)
    assert limit(r, x, 1).simplify() == n/2


def test_factorial():
    f = factorial(x)
    assert limit(f, x, oo) == oo
    assert limit(x/f, x, oo) == 0
    # see Stirling's approximation:
    # https://en.wikipedia.org/wiki/Stirling's_approximation
    assert limit(f/(sqrt(2*pi*x)*(x/E)**x), x, oo) == 1
    assert limit(f, x, -oo) == factorial(-oo)
    assert (limit(f, x, x**2) - factorial(x**2)).simplify() == 0
    assert (limit(f, x, -x**2) - factorial(-x**2)).simplify() == 0


def test_sympyissue_6560():
    e = 5*x**3/4 - 3*x/4 + (y*(3*x**2/2 - Rational(1, 2)) +
                            35*x**4/8 - 15*x**2/4 + Rational(3, 8))/(2*(y + 1))
    assert limit(e, y, oo) == (5*x**3 + 3*x**2 - 3*x - 1)/4


def test_sympyissue_5172():
    r = Symbol('r', positive=True)
    p = Symbol('p', positive=True)
    m = Symbol('m', negative=True)
    expr = ((2*n*(n - r + 1)/(n + r*(n - r + 1)))**c +
            (r - 1)*(n*(n - r + 2)/(n + r*(n - r + 1)))**c - n)/(n**c - n)
    expr = expr.subs({c: c + 1})
    assert limit(expr.subs({c: m}), n, oo) == 1
    assert limit(expr.subs({c: p}), n, oo).simplify() == \
        (2**(p + 1) + r - 1)/(r + 1)**(p + 1)


def test_sympyissue_7088():
    assert limit(sqrt(x/(x + a)), x, oo) == 1


def test_sympyissue_6364():
    e = z/(1 - sqrt(1 + z)*sin(a)**2 - sqrt(1 - z)*cos(a)**2)
    assert (limit(e, z, 0) - 2/cos(2*a)).simplify() == 0


def test_sympyissue_4099():
    assert limit(a/x, x, 0) == oo*sign(a)
    assert limit(-a/x, x, 0) == -oo*sign(a)
    assert limit(-a*x, x, oo) == -oo*sign(a)
    assert limit(a*x, x, oo) == oo*sign(a)


def test_sympyissue_4503():
    dx = Symbol('dx')
    assert limit((sqrt(1 + exp(x + dx)) - sqrt(1 + exp(x)))/dx, dx, 0) == \
        exp(x)/(2*sqrt(exp(x) + 1))


def test_sympyissue_8730():
    assert limit(subfactorial(x), x, oo) == oo


def test_diofantissue_55():
    assert limit((x + exp(x))/(x - 1), x, -oo) == 1
    assert limit((x*exp(x))/(exp(x) - 1), x, -oo) == 0  # issue sympy/sympy#2929


def test_sympyissue_8061():
    assert limit(4**(acos(1/(1 + x**2))**2)/log(1 + x, 4), x, 0) == oo


def test_sympyissue_8229():
    assert limit((root(x, 4) - 2)/(sqrt(x) - 4)**Rational(2, 3),
                 x, 16) == 0


def test_sympyissue_9205():
    assert Limit(x, x, a).free_symbols == {a}
    assert Limit(x, x, a, '-').free_symbols == {a}
    assert Limit(x + y, x + y, a).free_symbols == {a}
    assert Limit(-x**2 + y, x**2, a).free_symbols == {y, a}


def test_sympyissue_10610():
    assert limit(3**x*3**(-x - 1)*(x + 1)**2/x**2, x, oo) == Rational(1, 3)
    assert limit(2**x*2**(-x - 1)*(x + 1)*(y - 1)**(-x) *
                 (y - 1)**(x + 1)/(x + 2), x, oo) == y/2 - Rational(1, 2)


def test_sympyissue_9075():
    assert limit((6**(x + 1) + x + 1)/(6**x + x), x, oo) == 6


def test_sympyissue_8634():
    p = Symbol('p', positive=True)
    assert limit(x**p, x, -oo) == oo*sign((-1)**p)


def test_sympyissue_9558():
    assert limit(sin(x)**15, x, 0, '-') == 0  # should be fast


def test_diofantissue_296():
    e = log(exp(1/x)/Float(2) + exp(-1/x)/2)*x**2
    assert e.limit(x, oo) == 0.5


def test_sympyissue_5383():
    e = (1.0 + 1.0*x)**(1.0/x)
    assert e.limit(x, 0) == E.evalf()


def test_sympyissue_6171():
    e = Piecewise((0, x < 0), (1, True))
    assert e.limit(x, 0) == 1
    assert e.limit(x, 0, '-') == 0


def test_sympyissue_11526():
    df = diff(1/(a*log((x - b)/(x - c))), x)
    res = -1/(-a*c + a*b)
    assert limit(df, x, oo) == res
    assert limit(simplify(df), x, oo) == res

    e = log((1/x - b)/(1/x - c))
    assert e.as_leading_term(x) == x*(c - b)


def test_sympyissue_11672():
    assert limit(Rational(-1, 2)**x, x, oo) == 0
    assert limit(1/(-2)**x, x, oo) == 0


def test_sympyissue_8635():
    n = Symbol('n', integer=True, positive=True)

    k = 0
    assert limit(x**n - x**(n - k), x, oo) == 0
    k = 1
    assert limit(x**n - x**(n - k), x, oo) == oo
    k = 2
    assert limit(x**n - x**(n - k), x, oo) == oo
    k = 3
    assert limit(x**n - x**(n - k), x, oo) == oo


def test_sympyissue_8157():
    n = Symbol('n', integer=True)
    limit(cos(pi*n), n, oo)  # not raises


def test_sympyissue_5415():
    assert limit(polygamma(2 + 1/x, 3 + exp(-x)), x, oo) == polygamma(2, 3)


def test_sympyissue_2865():
    l1 = limit(O(1/x, (x, oo)), x, 0)
    assert l1 != 0 and isinstance(l1, Limit)
    l2 = limit(O(x, (x, oo)), x, 0)
    assert l2 != 0 and isinstance(l2, Limit)


def test_sympyissue_11879():
    assert limit(((x + y)**n - x**n)/y, y, 0).powsimp() == n*x**(n-1)


def test_sympyissue_12555():
    assert limit((3**x + 2*x**10)/(x**10 + E**x), x, -oo) == 2


def test_sympyissue_12769():
    r, z, x = symbols('r,z,x', real=True)
    a, b, s0, K, F0, s, T = symbols('a,b,s0,K,F0,s,T',
                                    positive=True, real=True)
    fx = (F0**b*K**b*r*s0 -
          sqrt((F0**2*K**(2*b)*a**2*(b - 1) +
                F0**(2*b)*K**2*a**2*(b - 1) +
                F0**(2*b)*K**(2*b)*s0**2*(b - 1)*(b**2 - 2*b + 1) -
                2*F0**(2*b)*K**(b + 1)*a*r*s0*(b**2 - 2*b + 1) +
                2*F0**(b + 1)*K**(2*b)*a*r*s0*(b**2 - 2*b + 1) -
                2*F0**(b + 1)*K**(b + 1)*a**2*(b - 1))/((b - 1)*(b**2 - 2*b + 1))))*(b*r - b - r + 1)
    assert limit(fx, K, F0) == (F0**(2*b)*b*r**2*s0 - 2*F0**(2*b)*b*r*s0 +
                                F0**(2*b)*b*s0 - F0**(2*b)*r**2*s0 +
                                2*F0**(2*b)*r*s0 - F0**(2*b)*s0)


def test_sympyissue_13332():
    assert limit(sqrt(30)*5**(-5*n - 1)*(46656*n)**n *
                 (5*n + 2)**(5*n + Rational(5, 2)) *
                 (6*n + 2)**(-6*n - Rational(5, 2)), n, oo) == Rational(25, 36)


def test_sympyissue_13382():
    assert limit(n*(((n + 1)**2 + 1)/(n**2 + 1) - 1), n, oo) == 2


def test_sympyissue_13403():
    assert limit(n*(-1 + (n + log(n + 1) + 1)/(n + log(n))), n, oo) == 1


def test_sympyissue_13416():
    assert limit((-n**3*log(n)**3 +
                  (n - 1)*(n + 1)**2*log(n + 1)**3)/(n**2*log(n)**3),
                 n, oo) == 1


def test_sympyissue_13462():
    assert limit(n**2*(2*n*(-(1 - 1/(2*n))**x + 1) -
                 x - (-x**2/4 + x/4)/n), n, oo) == x/12 - x**2/8 + x**3/24


def test_sympyissue_13575():
    assert limit(acos(erfi(x)), x, 1) == pi/2 + I*log(sqrt(erf(I)**2 + 1) +
                                                      erf(I))


def test_diofantissue_558():
    n = Symbol('n')
    r = Symbol('r', positive=True)
    c = Symbol('c')
    expr = ((2*n*(n - r + 1)/(n + r*(n - r + 1)))**c +
            (r - 1)*(n*(n - r + 2)/(n + r*(n - r + 1)))**c - n)/(n**c - n)
    expr = expr.subs({c: c + 1})
    assert limit(expr, n, oo) == Limit(expr, n, oo)


def test_sympyissue_14393():
    assert limit((x**b - y**b)/(x**a - y**a), x, y) == b*y**b/y**a/a


def test_sympyissue_14590():
    assert limit((n**3*((n + 1)/n)**n)/((n + 1)*(n + 2)*(n + 3)), n, oo) == E


def test_sympyissue_14793():
    e = ((x + Rational(1, 2))*log(x) - x +
         log(2*pi)/2 - log(factorial(x)) + 1/(12*x))*x**3
    assert limit(e, x, oo) == Rational(1, 360)


def test_sympyissue_14811():
    assert limit(((1 + Rational(2, 3)**(x + 1))**2**x)/(2**Rational(4, 3)**(x - 1)), x, oo) == oo


def test_sympyissue_15055():
    assert limit(n**3*((-n - 1)*sin(1/n) + (n + 2)*sin(1/(n + 1)))/(-n + 1), n, oo) == 1


def test_sympyissue_15146():
    assert limit((n/2)*(-2*n**3 - 2*(n**3 - 1)*n**2*digamma(n**3 + 1) +
                        2*(n**3 - 1)*n**2*digamma(n**3 + n + 1) +
                        n + 3), n, oo) == Rational(1, 3)


def test_sympyissue_15323():
    assert limit(((1 - 1/x)**x).diff(x), x, 1) == 1


def test_sympyissue_15984():
    assert limit((-x + log(exp(x) + 1))/x, x, oo, dir='-') == 0


def test_sympyissue_16222():
    assert limit(exp(x), x, 10000000) == exp(10000000)
    assert limit(exp(x), x, 100000000) == exp(100000000)
    assert limit(exp(x), x, 1000000000) == exp(1000000000)


@pytest.mark.timeout(20)
def test_sympyissue_15282():
    assert limit((x**2000 - (x + 1)**2000)/x**1999, x, oo) == -2000


def test_sympyissue_16722():
    n, z = symbols('n z')
    assert isinstance(limit(binomial(n + z, n)*n**-z, n, oo), Limit)

    z = symbols('z', positive=True)
    assert limit(binomial(n + z, n)*n**-z, n, oo) == 1/gamma(z + 1)

    n = symbols('n', positive=True, integer=True)
    z = symbols('z', positive=True)
    assert limit(binomial(n + z, n)*n**-z, n, oo) == 1/gamma(z + 1)

    n, z = symbols('n z', positive=True, integer=True)
    assert limit(binomial(n + z, n)*n**-z, n, oo) == 1/gamma(z + 1)


def test_sympyissue_15673():
    p = symbols('p')
    alpha = symbols('α', positive=True)

    e = Limit(4*pi*p**(-alpha)*(p**3 - p**alpha)/(alpha - 3), p, 0)
    assert isinstance(e.doit(), Limit)  # but see diofant/diofant#425


def test_sympyissue_17380():
    assert limit(x*(((x + 1)**2 + 1)/(x**2 + 1) - 1), x, oo) == 2


def test_sympyissue_17431():
    assert limit(((n + 1) + 1)/(((n + 1) + 2)*factorial(n + 1)) *
                 (n + 2)*factorial(n)/(n + 1), n, oo) == 0
    assert limit((n + 2)**2*factorial(n)/((n + 1)*(n + 3)*factorial(n + 1)),
                 n, oo) == 0

    # test from sympy/sympy#17434 (see also diofant/diofant#425):
    y = symbols('y', integer=True, positive=True)
    assert isinstance(limit(x*factorial(x)/factorial(x + y), x, oo), Limit)


def test_sympyissue_17792():
    n = Symbol('n', positive=True, integer=True)
    assert limit(factorial(n)/sqrt(n)*(E/n)**n, n, oo) == sqrt(2*pi)


def test_sympyissue_18118():
    assert limit(sign(x), x, 0, '+') == +1
    assert limit(sign(x), x, 0, '-') == -1

    assert limit(sign(sin(x)), x, 0, '+') == +1
    assert limit(sign(sin(x)), x, 0, '-') == -1


def test_sympyissue_6599():
    assert limit((x + cos(x))/x, x, oo) == 1


def test_sympyissue_18176():
    x = Symbol('x', real=True, positive=True)
    n = Symbol('n', integer=True, positive=True)
    k = Symbol('k')
    e = x**n - x**(n - k)
    assert limit(e.subs({k: 0}), x, oo) == 0
    assert limit(e.subs({k: 1}), x, oo) == oo


def test_sympyissue_18306():
    assert limit(sin(sqrt(x))/sqrt(sin(x)), x, 0) == 1


def test_sympyissue_18378():
    assert limit(log(exp(3*x) + x)/log(exp(x) + x**100), x, oo) == 3


def test_sympyissue_18399():
    assert limit((1 - x/2)**(3*x), x, oo) == oo
    assert limit((-x)**x, x, oo) == oo


def test_sympyissue_18452():
    assert limit(abs(log(x))**x, x, 0, '+') == 1
    assert limit(abs(log(x))**x, x, 0, '-') == 1


def test_sympyissue_18482():
    assert limit(sqrt(x**2 + 6*x) + (x**3 + x**2)/(x**2 + 1), x, -oo) == -2
    assert limit((x**3 + x**2 + sqrt(x*(x + 6))*(x**2 + 1))/(x**2 + 1),
                 x, -oo) == -2
    assert limit((2*exp(3*x)/(exp(2*x) + 1))**(1/x), x, oo) == E


def test_sympyissue_18501():
    assert limit(abs(log(x - 1)**3 - 1), x, 1, '+') == oo


def test_sympyissue_18508():
    assert limit(sin(x)/sqrt(1 - cos(x)), x, 0) == sqrt(2)


def test_sympyissue_18707():
    p = Symbol('p', positive=True, real=True)

    assert limit(1/p**n, n, oo) == p**-oo


def test_sympyissue_18997():
    assert limit(abs(log(x)), x, 0) == oo
    assert limit(abs(log(abs(x))), x, 0) == oo


def test_sympyissue_18992():
    assert limit(x/(factorial(x)**(1/x)), x, oo) == E

    n = Symbol('n', positive=True, integer=True)

    assert limit(n/(factorial(n)**(1/n)), n, oo) == E


def test_sympyissue_19026():
    assert limit(abs(log(x) + 1)/log(x), x, oo) == 1


def test_sympyissue_19770():
    m = Symbol('m', extended_real=True)

    assert limit(cos(x*m)/x, x, oo) == 0


def test_sympyissue_19766():
    assert limit(2**(-x)*sqrt(4**(x + 1) + 1), x, oo) == 2


def test_sympyissue_14874():
    assert limit(besselk(0, x), x, oo) == 0


def test_sympyissue_20365():
    assert limit(((x + 1)**(1/x) - E)/x, x, 0) == -E/2


def test_sympyissue_20704():
    assert limit(x*(abs(1/x + y) - abs(y - 1/x))/2, x, 0) == 0


def test_sympyissue_21031():
    assert limit(((1 + x)**(1/x) -
                  (1 + 2*x)**(1/(2*x)))/asin(x), x, 0) == E/2


def test_sympyissue_21038():
    assert limit(sin(pi*x)/(3*x - 12), x, 4) == pi/3


def test_sympyissue_21029():
    e = (sinh(x) + cosh(x) - 1)/x/4
    ans0 = Rational(1, 4)
    ansz = (exp(z) - 1)/z/4

    assert limit(e, x, 0) == ans0
    assert limit(e, x, z) == ansz
    assert limit(ansz, z, 0) == ans0


def test_sympyissue_20578():
    e = abs(x)*sin(1/x)

    assert all(_ == 0 for _ in [limit(e, x, 0, '+'),
                                limit(e, x, 0, '-'),
                                limit(e, x, 0, 'real')])


def test_sympyissue_19453():
    beta = Symbol('beta', real=True, positive=True)
    h = Symbol('h', real=True, positive=True)
    m = Symbol('m', real=True, positive=True)
    w = Symbol('omega', real=True, positive=True)
    g = Symbol('g', real=True, positive=True)

    q = 3*h**2*beta*g*exp(h*beta*w/2)
    p = m**2*w**2
    s = exp(h*beta*w) - 1
    z = (-q/(4*p*s) - q/(2*p*s**2) -
         q*(exp(h*beta*w) + 1)/(2*p*s**3) + exp(h*beta*w/2)/s)
    e = -diff(log(z), beta)

    assert limit(e - h*w/2, beta, oo) == 0
    assert limit(e.simplify() - h*w/2, beta, oo) == 0


def test_sympyissue_19442():
    pytest.raises(PoleError, lambda: limit(1/x, x, 0, 'real'))
