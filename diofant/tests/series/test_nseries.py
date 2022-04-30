import pytest

from diofant import (Derivative, E, I, O, PoleError, Rational, Symbol, acosh,
                     acoth, asin, asinh, atanh, besselk, cbrt, ceiling, cos,
                     cosh, cot, coth, exp, floor, limit, ln, log, oo, pi, sign,
                     sin, sinh, sqrt, tan, tanh)
from diofant.abc import a, b, l, w, x, y, z


__all__ = ()


def test_simple_1():
    assert x.series(x, n=5) == x
    assert y.series(x, n=5) == y
    assert (1/(x*y)).series(y, n=5) == 1/(x*y)
    assert Rational(3, 4).series(x, n=5) == Rational(3, 4)
    assert x.series(x) == x

    # issue sympy/sympy#5183
    assert ((1 + x)**2).series(x) == 1 + 2*x + x**2
    assert (1 + 1/x).series() == 1 + 1/x
    assert (Derivative(exp(x).series(), x).doit() ==
            1 + x + x**2/2 + x**3/6 + x**4/24 + Derivative(O(x**6), x))


def test_mul_0():
    assert (x*ln(x)).series(x, n=5) == x*ln(x)


def test_mul_1():
    assert (x*ln(2 + x)).series(x, n=5) == x*log(2) + x**2/2 - x**3/8 + \
        x**4/24 + O(x**5)
    assert (x*ln(1 + x)).series(x, n=5) == x**2 - x**3/2 + x**4/3 + O(x**5)


def test_pow_0():
    assert (x**2).series(x, n=5) == x**2
    assert (1/x).series(x, n=5) == 1/x
    assert (1/x**2).series(x, n=5) == 1/x**2
    assert (x**Rational(2, 3)).series(x, n=5) == (x**Rational(2, 3))
    assert (sqrt(x)**3).series(x, n=5) == (sqrt(x)**3)


def test_pow_1():
    assert ((1 + x)**2).series(x, n=5) == 1 + 2*x + x**2


def test_geometric_1():
    assert (1/(1 - x)).series(x, n=5) == 1 + x + x**2 + x**3 + x**4 + O(x**5)
    assert (x/(1 - x)).series(x, n=7) == x + x**2 + x**3 + x**4 + x**5 + \
        x**6 + O(x**7)
    assert (x**3/(1 - x)).series(x, n=8) == x**3 + x**4 + x**5 + x**6 + \
        x**7 + O(x**8)


def test_sqrt_1():
    assert sqrt(1 + x).series(x, n=5) == 1 + x/2 - x**2/8 + x**3/16 - 5*x**4/128 + O(x**5)


def test_exp_1():
    assert exp(x).series(x, n=5) == 1 + x + x**2/2 + x**3/6 + x**4/24 + O(x**5)
    assert exp(x).series(x, n=12) == 1 + x + x**2/2 + x**3/6 + x**4/24 + x**5/120 +  \
        x**6/720 + x**7/5040 + x**8/40320 + x**9/362880 + x**10/3628800 +  \
        x**11/39916800 + O(x**12)
    assert exp(1/x).series(x, n=5) == exp(1/x)
    assert exp(1/(1 + x)).series(x, n=4) ==  \
        (E*(1 - x - 13*x**3/6 + 3*x**2/2)).expand() + O(x**4)
    assert exp(2 + x).series(x, n=5) ==  \
        (exp(2)*(1 + x + x**2/2 + x**3/6 + x**4/24)).expand() + O(x**5)


def test_exp_sqrt_1():
    assert exp(1 + sqrt(x)).series(x, n=2) == \
        E + E*x/2 + E*sqrt(x) + E*x**Rational(3, 2)/6 + O(x**2)


def test_power_x_x1():
    assert (exp(x*ln(x))).series(x, n=3) == \
        1 + x*log(x) + x**2*log(x)**2/2 + x**3*log(x)**3/6 + O(x**3)


def test_power_x_x2():
    assert (x**x).series(x, n=3) == \
        1 + x*log(x) + x**2*log(x)**2/2 + x**3*log(x)**3/6 + O(x**3)


def test_log_singular1():
    assert log(1 + 1/x).series(x, n=5) == x - log(x) - x**2/2 + x**3/3 - \
        x**4/4 + O(x**5)


def test_log_power1():
    e = 1 / (1/x + x ** (log(3)/log(2)))
    assert e.series(x, n=3) == x - x**(2 + log(3)/log(2)) + O(x**3)


def test_log_series():
    e = 1/(1 - log(x))
    assert e.series(x, n=5, logx=l) == 1/(1 - l)


def test_log2():
    e = log(-1/x)
    assert e.series(x, n=5) == -log(x) + log(-1)


def test_log3():
    e = 1/log(-1/x)
    assert e.series(x, n=4, logx=l) == 1/(-l + log(-1))


def test_log4():
    assert (log(1 + x).series(x, x0=I*oo) ==
            1/(5*x**5) - 1/(4*x**4) + 1/(3*x**3) - 1/(2*x**2) + 1/x +
            I*pi/2 + log(-I*x) + O(x**(-6), (x, oo*I)))


def test_x0():
    # issue sympy/sympy#5654
    assert ((1/(x**2 + a**2)**2).series(x, x0=I*a, n=0) ==
            -I/(4*a**3*(-I*a + x)) - 1/(4*a**2*(-I*a + x)**2) + O(1, (x, I*a)))
    assert ((1/(x**2 + a**2)**2).series(x, x0=I*a, n=1) ==
            3/(16*a**4) - I/(4*a**3*(-I*a + x)) - 1/(4*a**2*(-I*a + x)**2) +
            O(-I*a + x, (x, I*a)))


def test_series1():
    e = sin(x)
    assert e.series(x, n=0) != 0
    assert e.series(x, n=1) == O(x)
    assert e.series(x, n=2) == x + O(x**2)
    assert e.series(x, n=3) == x + O(x**3)
    assert e.series(x, n=5) == x - x**3/6 + O(x**5)

    e = (exp(x) - 1)/x
    assert e.series(x, n=3) == 1 + x/2 + x**2/6 + O(x**3)

    assert x.series(x, n=2) == x


def test_series1_failing():
    assert x.series(x, n=0) == O(1, x)
    assert x.series(x, n=1) == O(x)


def test_seriesbug1():
    assert (1/x).series(x, n=3) == 1/x
    assert (x + 1/x).series(x, n=3) == x + 1/x


def test_series2x():
    assert ((x + 1)**(-2)).series(x, n=4) == 1 - 2*x + 3*x**2 - 4*x**3 + O(x**4)
    assert ((x + 1)**(-1)).series(x, n=4) == 1 - x + x**2 - x**3 + O(x**4)
    assert ((x + 1)**0).series(x, n=3) == 1
    assert ((x + 1)**1).series(x, n=3) == 1 + x
    assert ((x + 1)**2).series(x, n=3) == 1 + 2*x + x**2
    assert ((x + 1)**3).series(x, n=3) == 1 + 3*x + 3*x**2 + O(x**3)

    assert (1/(1 + x)).series(x, n=4) == 1 - x + x**2 - x**3 + O(x**4)
    assert (x + 3/(1 + 2*x)).series(x, n=4) == 3 - 5*x + 12*x**2 - 24*x**3 + O(x**4)

    assert ((1/x + 1)**3).series(x, n=4) == 1 + x**(-3) + 3*x**(-2) + 3/x
    assert (1/(1 + 1/x)).series(x, n=4) == x - x**2 + x**3 - O(x**4)
    assert (1/(1 + 1/x**2)).series(x) == x**2 - x**4 + O(x**6)


def test_bug2():  # 1/log(0) * log(0) problem
    e = (w**(-1) + w**(
        -log(3)*log(2)**(-1)))**(-1)*(3*w**(-log(3)*log(2)**(-1)) + 2*w**(-1))
    e = e.expand()
    assert e.series(w, n=4).subs({w: 0}) == 3


def test_exp():
    e = (1 + x)**(1/x)
    assert e.series(x, n=2) == exp(1) - x*exp(1)/2 + O(x**2)


def test_exp2():
    e = w**(1 - log(x)/(log(2) + log(x)))
    assert e.series(w, n=1) == e


def test_bug3():
    e = (2/x + 3/x**2)/(1/x + 1/x**2)
    assert e.series(x, n=3) == 3 - x + x**2 + O(x**3)


def test_generalexponent():
    p = 2
    e = (2/x + 3/x**p)/(1/x + 1/x**p)
    assert e.series(x, n=3) == 3 - x + x**2 + O(x**3)
    p = Rational(1, 2)
    e = (2/x + 3/x**p)/(1/x + 1/x**p)
    assert e.series(x, n=2) == 2 - x + sqrt(x) + x**Rational(3, 2) + O(x**2)

    e = 1 + sqrt(x)
    assert e.series(x, n=4) == 1 + sqrt(x)

# more complicated example


def test_genexp_x():
    e = 1/(1 + sqrt(x))
    assert e.series(x, n=2) == 1 + x - sqrt(x) - sqrt(x)**3 + O(x**2)

# more complicated example


def test_genexp_x2():
    p = Rational(3, 2)
    e = (2/x + 3/x**p)/(1/x + 1/x**p)
    assert e.series(x, n=2) == 3 + x - sqrt(x) - x**p + O(x**2)


def test_seriesbug2():
    # simple case (1):
    e = ((2*w)/w)**(1 + w)
    assert e.series(w, n=1) == 2 + O(w)
    assert e.series(w, n=1).subs({w: 0}) == 2


def test_seriesbug2b():
    # test sin
    e = sin(2*w)/w
    assert e.series(w, n=2) == 2 + O(w**2)


def test_seriesbug2d():
    w = Symbol('w', extended_real=True)
    e = log(sin(2*w)/w)
    assert e.series(w, n=5) == log(2) - 2*w**2/3 - 4*w**4/45 + O(w**5)


def test_seriesbug2c():
    w = Symbol('w', extended_real=True)
    # more complicated case, but sin(x)~x, so the result is the same as in (1)
    e = (sin(2*w)/w)**(1 + w)
    assert e.series(w, 0, 1) == 2 + O(w)
    assert e.series(w, 0, 3) == 2 + 2*w*log(2) + \
        w**2*(-Rational(4, 3) + log(2)**2) + O(w**3)
    assert e.series(w, 0, 2).subs({w: 0}) == 2


def test_expbug4():
    x = Symbol('x', extended_real=True)
    assert (log(
        sin(2*x)/x)*(1 + x)).series(x, 0, 2) == log(2) + x*log(2) + O(x**2)
    assert exp(
        log(sin(2*x)/x)*(1 + x)).series(x, 0, 2) == 2 + 2*x*log(2) + O(x**2)

    assert exp(log(2) + O(x)).series(x, n=1) == 2 + O(x)
    assert ((2 + O(x))**(1 + x)).series(x, n=1) == 2 + O(x)


def test_logbug4():
    assert log(2 + O(x)).series(x, n=1) == log(2) + O(x)


def test_expbug5():
    assert exp(log(1 + x)/x).series(x, n=2) == exp(1) + -exp(1)*x/2 + O(x**2)

    assert exp(O(x)).series(x, n=1) == 1 + O(x)


def test_sinsinbug():
    assert sin(sin(x)).series(x, n=9) == x - x**3/3 + x**5/10 - 8*x**7/315 + O(x**9)


def test_sympyissue_3258():
    a = x/(exp(x) - 1)
    assert a.series(x, n=5) == 1 - x/2 - x**4/720 + x**2/12 + O(x**5)


def test_sympyissue_3204():
    x = Symbol('x', nonnegative=True)
    f = cbrt(sin(x**3))
    assert f.series(x, n=19) == x - x**7/18 - x**13/3240 + O(x**19)


def test_sympyissue_3224():
    f = sqrt(1 - sqrt(y))
    assert f.series(y, n=2) == 1 - sqrt(y)/2 - y/8 - sqrt(y)**3/16 + O(y**2)


def test_sympyissue_3463():
    r = log(5)/log(3)
    p = w**(-1 + r)
    e = 1/x*(-log(w**(1 + r)) + log(w + w**r))
    e_ser = -r*log(w)/x + p/x - p**2/(2*x) + O(p**3)
    assert (e.series(w, n=3) - e_ser).removeO().simplify() == 0


def test_sin():
    assert sin(8*x).series(x, n=4) == 8*x - 256*x**3/3 + O(x**4)
    assert sin(x + y).series(x, n=1) == sin(y) + O(x)
    assert sin(x + y).series(x, n=2) == sin(y) + cos(y)*x + O(x**2)
    assert sin(x + y).series(x, n=5) == sin(y) + cos(y)*x - sin(y)*x**2/2 - \
        cos(y)*x**3/6 + sin(y)*x**4/24 + O(x**5)


def test_sympyissue_3515():
    e = sin(8*x)/x
    assert e.series(x) == 8 - 256*x**2/3 + 4096*x**4/15 + O(x**6)


def test_sympyissue_3505():
    e = sin(x)**(-4)*(sqrt(cos(x))*sin(x)**2 - cbrt(cos(x))*sin(x)**2)
    assert e.series(x) == -Rational(1, 12) - 7*x**2/288 - \
        43*x**4/10368 + O(x**6)


def test_sympyissue_3501():
    e = x**(-2)*(x*sin(a + x) - x*sin(a))
    assert e.series(x, n=4) == cos(a) - sin(a)*x/2 - cos(a)*x**2/6 + \
        sin(a)*x**3/24 + O(x**4)
    e = x**(-2)*(x*cos(a + x) - x*cos(a))
    assert e.series(x, n=5) == -sin(a) - cos(a)*x/2 + sin(a)*x**2/6 + \
        cos(a)*x**3/24 - x**4*sin(a)/120 + O(x**5)


def test_sympyissue_3502():
    e = sin(5*x)/sin(2*x)
    assert e.series(x, n=2) == Rational(5, 2) + O(x**2)
    assert e.series(x) == \
        Rational(5, 2) - 35*x**2/4 + 329*x**4/48 + O(x**6)


def test_sympyissue_3503():
    e = sin(2 + x)/(2 + x)
    assert e.series(x, n=2) == sin(2)/2 + x*(-sin(2)/4 + cos(2)/2) + O(x**2)


def test_sympyissue_3506():
    e = (x + sin(3*x))**(-2)*(x*(x + sin(3*x)) - (x + sin(3*x))*sin(2*x))
    assert e.series(x) == -Rational(1, 4) + 5*x**2/96 + 91*x**4/768 + O(x**6)


def test_sympyissue_3508():
    x = Symbol('x', extended_real=True)
    assert log(sin(x)).series(x, n=5) == log(x) - x**2/6 - x**4/180 + O(x**5)
    e = -log(x) + x*(-log(x) + log(sin(2*x))) + log(sin(2*x))
    assert e.series(x, n=5) == \
        log(2) + log(2)*x - 2*x**2/3 - 2*x**3/3 - 4*x**4/45 + O(x**5)


def test_sympyissue_3507():
    e = x**(-4)*(x**2 - x**2*sqrt(cos(x)))
    assert e.series(x) == \
        Rational(1, 4) + x**2/96 + 19*x**4/5760 + O(x**6)


def test_sympyissue_3639():
    assert sin(cos(x)).series(x, n=5) == \
        sin(1) - x**2*cos(1)/2 + x**4*(-sin(1)/8 + cos(1)/24) + O(x**5)


def test_hyperbolic():
    assert sinh(x).series(x, n=7) == x + x**3/6 + x**5/120 + O(x**7)
    assert cosh(x).series(x) == 1 + x**2/2 + x**4/24 + O(x**6)
    assert tanh(x).series(x, n=7) == x - x**3/3 + 2*x**5/15 + O(x**7)
    assert coth(x).series(x, n=7) == \
        1/x - x**3/45 + x/3 + 2*x**5/945 + O(x**7)
    assert asinh(x).series(x, n=7) == x - x**3/6 + 3*x**5/40 + O(x**7)
    assert acosh(x).series(x, n=7) == \
        pi*I/2 - I*x - 3*I*x**5/40 - I*x**3/6 + O(x**7)
    assert atanh(x).series(x, n=7) == x + x**3/3 + x**5/5 + O(x**7)
    assert acoth(x).series(x, n=7) == -I*pi/2 + x + x**3/3 + x**5/5 + O(x**7)


def test_series2():
    w = Symbol('w', extended_real=True)
    x = Symbol('x', extended_real=True)
    e = w**(-2)*(w*exp(1/x - w) - w*exp(1/x))
    assert e.series(w, n=2) == -exp(1/x) + w*exp(1/x)/2 + O(w**2)


def test_series3():
    w = Symbol('w', extended_real=True)
    e = w**(-6)*(w**3*tan(w) - w**3*sin(w))
    assert e.series(w, n=2) == Rational(1, 2) + O(w**2)


def test_bug4():
    e = x/(w**4 + x**2*w**4 + 2*x*w**4)*w**4
    assert e.series(w, n=2).simplify() in [x/(1 + 2*x + x**2),
                                           1/(1 + x/2 + 1/x/2)/2, 1/x/(1 + 2/x + x**(-2))]


def test_bug5():
    e = (-log(w) + log(1 + w*log(x)))**(-2)*w**(-2)*((-log(w) +
                                                      log(1 + x*w))*(-log(w) + log(1 + w*log(x)))*w - x*(-log(w) +
                                                                                                         log(1 + w*log(x)))*w)
    assert e.series(w, n=0, logx=l) == (x/l + 1)/w + O(1, w)
    assert e.series(w, n=1, logx=l) == x*log(x)/l**2 + log(x)/l - \
        x/l + (1 + x/l)/w + O(w)


def test_sympyissue_4115():
    assert (sin(x)/(1 - cos(x))).series(x, n=1) == 2/x + O(x)
    assert (sin(x)**2/(1 - cos(x))).series(x, n=2) == 2 + O(x**2)


def test_pole():
    pytest.raises(PoleError, lambda: sin(1/x).series(x, 0, 5))
    pytest.raises(PoleError, lambda: sin(1 + 1/x).series(x, 0, 5))
    pytest.raises(PoleError, lambda: (x*sin(1/x)).series(x, 0, 5))
    pytest.raises(PoleError, lambda: besselk(0, x).series(x, 0, 2))


def test_expsinbug():
    assert exp(sin(x)).series(x, 0, 0) == O(1, x)
    assert exp(sin(x)).series(x, 0, 1) == 1 + O(x)
    assert exp(sin(x)).series(x, 0, 2) == 1 + x + O(x**2)
    assert exp(sin(x)).series(x, 0, 3) == 1 + x + x**2/2 + O(x**3)
    assert exp(sin(x)).series(x, 0, 4) == 1 + x + x**2/2 + O(x**4)
    assert exp(sin(x)).series(x, 0, 5) == 1 + x + x**2/2 - x**4/8 + O(x**5)


def test_floor():
    x = Symbol('x')
    assert floor(x).series(x) == 0
    assert floor(-x).series(x) == -1
    assert floor(sin(x)).series(x) == 0
    assert floor(sin(-x)).series(x) == -1
    assert floor(x**3).series(x) == 0
    assert floor(-x**3).series(x) == -1
    assert floor(cos(x)).series(x) == 0
    assert floor(cos(-x)).series(x) == 0
    assert floor(5 + sin(x)).series(x) == 5
    assert floor(5 + sin(-x)).series(x) == 4

    assert floor(x).series(x, 2) == 2
    assert floor(-x).series(x, 2) == -3

    x = Symbol('x', negative=True)
    assert floor(x + 1.5).series(x) == 1


def test_ceiling():
    assert ceiling(x).series(x) == 1
    assert ceiling(-x).series(x) == 0
    assert ceiling(sin(x)).series(x) == 1
    assert ceiling(sin(-x)).series(x) == 0
    assert ceiling(1 - cos(x)).series(x) == 1
    assert ceiling(1 - cos(-x)).series(x) == 1
    assert ceiling(x).series(x, 2) == 3
    assert ceiling(-x).series(x, 2) == -2


def test_abs():
    assert abs(x).series(x, n=4) == x
    assert abs(-x).series(x, n=4) == x
    assert abs(x + 1).series(x, n=4) == x + 1
    assert abs(sin(x)).series(x, n=4) == x - x**3/6 + O(x**4)
    assert abs(sin(-x)).series(x, n=4) == x - x**3/6 + O(x**4)
    assert abs(x - a).series(x, 1) == (x - a)/sign(1 - a)

    # issue sympy/sympy#5183
    assert abs(x + x**2).series(n=1) == O(x)
    assert abs(x + x**2).series(n=2) == x + O(x**2)


def test_dir():
    assert abs(x).series(x, 0, dir=-1) == x
    assert abs(x).series(x, 0, dir=+1) == -x
    assert floor(x + 2).series(x, 0, dir=-1) == 2
    assert floor(x + 2).series(x, 0, dir=+1) == 1
    assert floor(x + 2.2).series(x, 0, dir=+1) == 2
    assert ceiling(x + 2.2).series(x, 0, dir=+1) == 3
    assert sin(x + y).series(x, 0, dir=+1) == sin(x + y).series(x, 0, dir=-1)


def test_sympyissue_3504():
    e = asin(a*x)/x
    assert e.series(x, 4, n=2).removeO() == \
        (x - 4)*(a/(4*sqrt(-16*a**2 + 1)) - asin(4*a)/16) + asin(4*a)/4


def test_sympyissue_4441():
    f = 1/(1 + a*x)
    assert f.series(x, 0, 5) == 1 - a*x + a**2*x**2 - a**3*x**3 + \
        a**4*x**4 + O(x**5)
    f = 1/(1 + (a + b)*x)
    assert f.series(x, 0, 3) == 1 + x*(-a - b) + x**2*(a**2 + 2*a*b + b**2) + O(x**3)


def test_sympyissue_4329():
    assert tan(x).series(x, pi/2, n=3).removeO() == \
        -pi/6 + x/3 - 1/(x - pi/2)
    assert cot(x).series(x, pi, n=3).removeO() == \
        -x/3 + pi/3 + 1/(x - pi)
    assert limit(tan(x)**tan(2*x), x, pi/4) == exp(-1)


def test_sympyissue_5925():
    sx = sqrt(x + z).series(z, 0, 1)
    sxy = sqrt(x + y + z).series(z, 0, 1)
    s1, s2 = sx.subs({x: x + y}), sxy
    assert (s1 - s2).expand().removeO().simplify() == 0

    sx = sqrt(x + z).series(z, 0, 1)
    sxy = sqrt(x + y + z).series(z, 0, 1)
    assert sxy.subs({x: 1, y: 2}) == sx.subs({x: 3})


def test_sympyissue_6235():
    q = Symbol('q', positive=True)
    assert (((x - 1)**q + 1)/(x**q - 1)).nseries(x, n=2).removeO() == \
        (-1 - x**q + (-1)**(q + 1) + (-1)**(q + 1)*x**q +
         (-1)**q*q*x**(q + 1) + (-1)**q*q*x)


def test_sympyissue_6236():
    q = Symbol('q', positive=True)
    assert (((x - 1)**q)/(x**q - 1)).nseries(x, n=2).removeO() == \
        (-1)**(q + 1) + (-1)**(q + 1)*x**q + (-1)**q*q*x**(q + 1) + (-1)**q*q*x


def test_issue_210():
    assert cos(x**6).series(x, n=12) == 1 + O(x**12)
    assert cos(x**6).series(x, n=23) == 1 - x**12/2 + O(x**23)
    assert cos(x**6).series(x, n=24) == 1 - x**12/2 + O(x**24)
    assert cos(x**6).series(x, n=36) == 1 - x**12/2 + x**24/24 + O(x**36)

    # issue sympy/sympy#10503
    f = exp(x**3)*cos(x**6)
    assert f.series(x, n=14) == (1 + x**3 + x**6/2 +
                                 x**9/6 - 11*x**12/24 + O(x**14))
    assert f.series(x, n=15) == (1 + x**3 + x**6/2 +
                                 x**9/6 - 11*x**12/24 + O(x**15))
    assert f.series(x, n=16) == (1 + x**3 + x**6/2 + x**9/6 - 11*x**12/24 -
                                 59*x**15/120 + O(x**16))


def test_sympyissue_21075():
    e = (sqrt(x) + cbrt(x))**2
    assert e.series(x) == e.expand()


def test_sympyissue_21227():
    f = log(x)

    assert f.series(x, logx=y) == y
    assert f.series(x, logx=-x) == -x

    f = log(-log(x))

    assert f.series(x, logx=y) == log(-y)
    assert f.series(x, logx=-x) == log(x)

    f = log(log(x))

    assert f.series(x, logx=y) == log(y)
    assert f.series(x, logx=-x) == log(-x)
    assert f.series(x, logx=x) == log(x)
