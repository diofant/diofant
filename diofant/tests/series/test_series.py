import pytest

from diofant import (Derivative, E, Function, Integer, Integral, O, Rational,
                     Subs, Symbol, cos, exp, log, oo, pi, series, sin, sqrt,
                     symbols)
from diofant.abc import x, y


__all__ = ()


def test_sin():
    e1 = sin(x).series(x, 0)
    e2 = series(sin(x), x, 0)
    assert e1 == e2

    # issue sympy/sympy#5223:
    assert ((1/sin(x))**oo).series() == oo


def test_cos():
    e1 = cos(x).series(x, 0)
    e2 = series(cos(x), x, 0)
    assert e1 == e2


def test_exp():
    e1 = exp(x).series(x, 0)
    e2 = series(exp(x), x, 0)
    assert e1 == e2


def test_exp2():
    e1 = exp(cos(x)).series(x, 0)
    e2 = series(exp(cos(x)), x, 0)
    assert e1 == e2


def test_simple():
    # issue sympy/sympy#5223
    assert series(1, x) == 1
    pytest.raises(ValueError, lambda: cos(x + y).series())
    pytest.raises(ValueError, lambda: x.series(dir=''))

    D = Derivative
    assert D(x**2 + x**3*y**2, (x, 2), (y, 1)).series(x).doit().removeO() == 12*x*y

    assert (1 + x).getn() is None

    # issue sympy/sympy#8805
    assert series(1, n=8) == 1

    # issue sympy/sympy#5223
    assert (cos(x).series(x, 1) -
            cos(x + 1).series(x).subs({x: x - 1})).removeO() == 0

    assert abs(x).series(x, oo, n=5, dir='+') == x
    assert abs(x).series(x, -oo, n=5, dir='-') == -x
    assert abs(-x).series(x, oo, n=5, dir='+') == x
    assert abs(-x).series(x, -oo, n=5, dir='-') == -x

    # issue sympy/sympy#7203
    assert series(cos(x), x, pi, 3) == -1 + (x - pi)**2/2 + O((x - pi)**3,
                                                              (x, pi))


def test_sympyissue_5223():
    assert next(Integer(0).lseries(x)) == 0
    assert cos(x).series() == cos(x).series(x)

    e = cos(x).series(x, 1, n=None)
    assert [next(e) for i in range(2)] == [cos(1), -((x - 1)*sin(1))]
    e = cos(x).series(x, 1, n=None, dir='-')
    assert [next(e) for i in range(2)] == [cos(1), (1 - x)*sin(1)]
    # the following test is exact so no need for x -> x - 1 replacement
    assert abs(x).series(x, 1, dir='-') == x
    assert exp(x).series(x, 1, dir='-', n=3).removeO() == \
        E - E*(-x + 1) + E*(-x + 1)**2/2

    D = Derivative
    assert next(D(cos(x), x).lseries()) == D(1, x)
    assert D(exp(x), x).series(n=3) == (D(1, x) + D(x, x) + D(x**2/2, x) +
                                        D(x**3/6, x) + O(x**3))

    assert Integral(x, (x, 1, 3), (y, 1, x)).series(x) == -4 + 4*x

    assert (1 + x + O(x**2)).getn() == 2

    assert ((sin(x))**y).nseries(x, n=1) == x**y + O(x**(y + 2), x)

    assert sin(1/x).series(x, oo, n=5) == 1/x - 1/(6*x**3) + O(x**(-5), (x, oo))

    assert exp(x*log(x)).series(n=3) == \
        1 + x*log(x) + x**2*log(x)**2/2 + x**3*log(x)**3/6 + O(x**3)

    p = Symbol('p', positive=True)
    assert exp(sqrt(p)**3*log(p)).series(n=3) == \
        1 + p**3*log(p)**2/2 + p**Rational(3, 2)*log(p) + O(p**3)

    assert exp(sin(x)*log(x)).series(n=2) == \
        1 + x*log(x) + x**2*log(x)**2/2 + O(x**2)


def test_sympyissue_3978():
    f = Function('f')
    assert f(x).series(x, 0, 3, dir='-') == \
        f(0) + x*Subs(Derivative(f(x), x), (x, 0)) + \
        x**2*Subs(Derivative(f(x), x, x), (x, 0))/2 + O(x**3)
    assert f(x).series(x, 0, 3) == \
        f(0) + x*Subs(Derivative(f(x), x), (x, 0)) + \
        x**2*Subs(Derivative(f(x), x, x), (x, 0))/2 + O(x**3)
    assert f(x**2).series(x, 0, 3) == \
        f(0) + x**2*Subs(Derivative(f(x), x), (x, 0)) + O(x**3)
    assert f(x**2+1).series(x, 0, 3) == \
        f(1) + x**2*Subs(Derivative(f(x), x), (x, 1)) + O(x**3)

    class TestF(Function):
        pass

    assert TestF(x).series(x, 0, 3) == TestF(0) + \
        x*Subs(Derivative(TestF(x), x), (x, 0)) + \
        x**2*Subs(Derivative(TestF(x), x, x), (x, 0))/2 + O(x**3)


def test_sympyissue_5852():
    assert series(1/cos(x/log(x)), x, 0) == 1 + x**2/(2*log(x)**2) + \
        5*x**4/(24*log(x)**4) + O(x**6)


def test_sympyissue_4583():
    assert cos(1 + x + x**2).series(x, 0, 5) == cos(1) - x*sin(1) + \
        x**2*(-sin(1) - cos(1)/2) + x**3*(-cos(1) + sin(1)/6) + \
        x**4*(-11*cos(1)/24 + sin(1)/2) + O(x**5)


def test_sympyissue_6318():
    eq = (1/x)**Rational(2, 3)
    assert (eq + 1).as_leading_term(x) == eq


def test_x_is_base_detection():
    eq = (x**2)**Rational(2, 3)
    assert eq.series() == x**Rational(4, 3)


def test_sin_power():
    e = sin(x)**1.2
    assert e.compute_leading_term(x) == x**1.2


@pytest.mark.xfail(reason='https://github.com/diofant/diofant/pull/158')
def test_exp_product_positive_factors():
    a, b = symbols('a, b', positive=True)
    x = a * b
    series(exp(x), x, n=8)
    # (1 + a*b + a**2*b**2/2 +
    #  a**3*b**3/6 + a**4*b**4/24 + a**5*b**5/120 + a**6*b**6/720 +
    #  a**7*b**7/5040 + O(a**8*b**8))


def test_sympyissue_9173():
    p_0, p_1, p_2, p_3, b_0, b_1, b_2 = symbols('p_0:4, b_0:3')
    Q = (p_0 + (p_1 + (p_2 + p_3/y)/y)/y)/(1 + ((p_3/(b_0*y) +
                                                 (b_0*p_2 - b_1*p_3)/b_0**2)/y + (b_0**2*p_1 - b_0*b_1*p_2 -
                                                                                  p_3*(b_0*b_2 - b_1**2))/b_0**3)/y)
    assert Q.series(y, n=3) == b_2*y**2 + b_1*y + b_0 + O(y**3)


@pytest.mark.slow
def test_sympyissue_9549():
    e = (x**2 + x + 1)/(x**3 + x**2)
    r = e.series(x, oo)
    assert r == x**(-5) - 1/x**4 + x**(-3) + 1/x + O(x**(-6), (x, oo))
    assert e.series(x, oo, n=8) + O(1/x**6, (x, oo)) == r


def test_sympyissue_10761():
    e = 1/(x**-2 + x**-3)
    assert e.series(x) == x**3 - x**4 + x**5 + O(x**6)
    # more tests from https://github.com/sympy/sympy/pull/10762
    assert e.series(x, n=10) == (x**3 - x**4 + x**5 - x**6 + x**7
                                 - x**8 + x**9 + O(x**10))
    assert e.series(x, n=20) == (x**3 - x**4 + x**5 - x**6 + x**7
                                 - x**8 + x**9 - x**10 + x**11 - x**12
                                 + x**13 - x**14 + x**15 - x**16
                                 + x**17 - x**18 + x**19 + O(x**20))


def test_sympyissue_11407():
    a, b, c = symbols('a, b, c')
    assert sqrt(a + b + c*x).series(x, 0, 1) == sqrt(a + b) + O(x)
    assert sqrt(a + b + c + c*x).series(x, 0, 1) == sqrt(a + b + c) + O(x)


def test_sympyissue_6179():
    assert (sin(x)*log(x)).series(x, 0, 4) == (x*log(x) -
                                               x**3*log(x)/6 + O(x**4))
    assert ((x**2*(x**3 + x**2 + 1)*log(x)).series(x, 0, 4) ==
            x**2*log(x) + x**4*log(x) + O(x**4))


@pytest.mark.slow
def test_sympyissue_11722():
    t, g = symbols('t g')
    good = -g**4*t**4/4 + 7*g**3*t**4/3 + g**3*t**3/3 - 27*g**2*t**4/4 - 2*g**2*t**3 - g**2*t**2/2 + 15*g*t**4/2 + 19*g*t**3/6 + 3*g*t**2/2 + g*t + g - 2009*t**4/720 - 13*t**3/9 - 5*t**2/6 - t/2 - (g + log(-t + 1) - 1 + (g + log(-t + 1))/(-1 + 1/t) - 1/(2*(-1 + 1/t)) - (g + log(-t + 1))**2/(2*(-1 + 1/t)**2) + 3*(g + log(-t + 1))/(2*(-1 + 1/t)**2) - 5/(6*(-1 + 1/t)**2) + (g + log(-t + 1))**3/(3*(-1 + 1/t)**3) - 2*(g + log(-t + 1))**2/(-1 + 1/t)**3 + 19*(g + log(-t + 1))/(6*(-1 + 1/t)**3) - 13/(9*(-1 + 1/t)**3) - (g + log(-t + 1))**4/(4*(-1 + 1/t)**4) + 7*(g + log(-t + 1))**3/(3*(-1 + 1/t)**4) - 27*(g + log(-t + 1))**2/(4*(-1 + 1/t)**4) + 15*(g + log(-t + 1))/(2*(-1 + 1/t)**4) - 2009/(720*(-1 + 1/t)**4) + 1/t)/(1 - 1/(g + log(-t + 1) - 1 + (g + log(-t + 1))/(-1 + 1/t) - 1/(2*(-1 + 1/t)) - (g + log(-t + 1))**2/(2*(-1 + 1/t)**2) + 3*(g + log(-t + 1))/(2*(-1 + 1/t)**2) - 5/(6*(-1 + 1/t)**2) + (g + log(-t + 1))**3/(3*(-1 + 1/t)**3) - 2*(g + log(-t + 1))**2/(-1 + 1/t)**3 + 19*(g + log(-t + 1))/(6*(-1 + 1/t)**3) - 13/(9*(-1 + 1/t)**3) - (g + log(-t + 1))**4/(4*(-1 + 1/t)**4) + 7*(g + log(-t + 1))**3/(3*(-1 + 1/t)**4) - 27*(g + log(-t + 1))**2/(4*(-1 + 1/t)**4) + 15*(g + log(-t + 1))/(2*(-1 + 1/t)**4) - 2009/(720*(-1 + 1/t)**4) + 1/t)) + 1/t
    bad = good.subs({g: log(1/t)})
    assert bad.series(t, x0=0, n=5) == O(t**5)


def test_sympyissue_11884():
    assert O(x).subs({x: x - 1}) + 1 == 1 + O(x - 1, (x, 1))
    assert cos(x).series(x, x0=1, n=1) == cos(1) + O(x - 1, (x, 1))


def test_sympyissue_12375():
    s = (x + 1).series(x, 2, 1)
    assert s == 3 + O(x - 2, (x, 2))
    assert s.removeO() == 3


def test_sympyissue_12747():
    assert exp(x).series(x, y, n=1) == exp(y) + O(x - y, (x, y))


def test_sympyissue_14384():
    assert (x**y).series(x) == x**y


def test_sympyissue_14885():
    assert ((x**Rational(-3, 2)*exp(x)).series(x) ==
            (x**Rational(-3, 2) + 1/sqrt(x) + sqrt(x)/2 + x**Rational(3, 2)/6 +
             x**Rational(5, 2)/24 + x**Rational(7, 2)/120 +
             x**Rational(9, 2)/720 + x**Rational(11, 2)/5040 + O(x**6)))


def test_sympyissue_15539():
    assert series(exp(x), x, x0=-oo) == exp(-x)


def test_sympyissue_18008():
    e = x*(x*(-x + 1) + 1)/(x*(-x + 1) - (-x + 1)**2 + 1)
    es = e.simplify()
    s = e.series(x, x0=oo, n=4)
    ss = es.series(x, x0=oo, n=4)
    assert s == ss


def test_sympyissue_20697():
    p0, p1, p2, p3 = symbols('p:4')
    b0, b1, b2 = symbols('b:3')

    e = ((p0 + (p1 + (p2 + p3/y)/y)/y) /
         (1 + ((p3/(b0*y) + (b0*p2 - b1*p3)/b0**2)/y +
               (b0**2*p1 - b0*b1*p2 - p3*(b0*b2 - b1**2))/b0**3)/y))

    assert e.series(y, n=3) == b2*y**2 + b1*y + b0 + O(y**3)
