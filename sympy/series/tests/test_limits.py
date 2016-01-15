from itertools import product as cartes

import pytest

from sympy import (limit, exp, oo, log, sqrt, Limit, sin, floor, cos,
                   acos, ceiling, atan, gamma, Symbol, S, pi, Integral,
                   cot, Rational, I, zoo, tan, cot, integrate, Sum, sign,
                   Function, subfactorial, PoleError, Integer)
from sympy.series.limits import heuristics
from sympy.series.order import Order

from sympy.abc import x, y, z, a


def test_basic1():
    assert limit(x, x, oo) == oo
    assert limit(x, x, -oo) == -oo
    assert limit(-x, x, oo) == -oo
    assert limit(x**2, x, -oo) == oo
    assert limit(-x**2, x, oo) == -oo
    assert limit(x*log(x), x, 0, dir="+") == 0
    assert limit(1/x, x, oo) == 0
    assert limit(exp(x), x, oo) == oo
    assert limit(-exp(x), x, oo) == -oo
    assert limit(exp(x)/x, x, oo) == oo
    assert limit(1/x - exp(-x), x, oo) == 0
    assert limit(x + 1/x, x, oo) == oo
    assert limit(x - x**2, x, oo) == -oo
    assert limit((1 + x)**(1 + sqrt(2)), x, 0) == 1
    assert limit((1 + x)**oo, x, 0) == oo
    assert limit((1 + x)**oo, x, 0, dir='-') == 0
    assert limit((1 + x + y)**oo, x, 0, dir='-') == (1 + y)**(oo)
    assert limit(y/x/log(x), x, 0) == -oo*sign(y)
    assert limit(cos(x + y)/x, x, 0) == sign(cos(y))*oo
    limit(Sum(1/x, (x, 1, y)) - log(y), y, oo)
    limit(Sum(1/x, (x, 1, y)) - 1/y, y, oo)
    assert limit(gamma(1/x + 3), x, oo) == 2
    assert limit(S.NaN, x, -oo) == S.NaN
    assert limit(Order(2)*x, x, S.NaN) == S.NaN
    assert limit(1/(x - 1), x, 1, dir="+") == oo
    assert limit(1/(x - 1), x, 1, dir="-") == -oo
    assert limit(1/(5 - x)**3, x, 5, dir="+") == -oo
    assert limit(1/(5 - x)**3, x, 5, dir="-") == oo
    assert limit(1/sin(x), x, pi, dir="+") == -oo
    assert limit(1/sin(x), x, pi, dir="-") == oo
    assert limit(1/cos(x), x, pi/2, dir="+") == -oo
    assert limit(1/cos(x), x, pi/2, dir="-") == oo
    assert limit(1/tan(x**3), x, (2*pi)**Rational(1, 3), dir="+") == oo
    assert limit(1/tan(x**3), x, (2*pi)**Rational(1, 3), dir="-") == -oo
    assert limit(1/cot(x)**3, x, (3*pi/2), dir="+") == -oo
    assert limit(1/cot(x)**3, x, (3*pi/2), dir="-") == oo

    # approaching 0
    # from dir="+"
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


def test_basic2():
    assert limit(x**x, x, 0, dir="+") == 1
    assert limit((exp(x) - 1)/x, x, 0) == 1
    assert limit(1 + 1/x, x, oo) == 1
    assert limit(-exp(1/x), x, oo) == -1
    assert limit(x + exp(-x), x, oo) == oo
    assert limit(x + exp(-x**2), x, oo) == oo
    assert limit(x + exp(-exp(x)), x, oo) == oo
    assert limit(13 + 1/x - exp(-x), x, oo) == 13


def test_basic3():
    assert limit(1/x, x, 0, dir="+") == oo
    assert limit(1/x, x, 0, dir="-") == -oo


def test_basic4():
    assert limit(2*x + y*x, x, 0) == 0
    assert limit(2*x + y*x, x, 1) == 2 + y
    assert limit(2*x**8 + y*x**(-3), x, -2) == 512 - y/8
    assert limit(sqrt(x + 1) - sqrt(x), x, oo) == 0
    assert integrate(1/(x**3 + 1), (x, 0, oo)) == 2*pi*sqrt(3)/9


def test_basic5():
    class my(Function):
        @classmethod
        def eval(cls, arg):
            if arg is S.Infinity:
                return S.NaN
    assert limit(my(x), x, oo) == Limit(my(x), x, oo)


def test_issue_3885():
    assert limit(x*y + x*z, z, 2) == x*y + 2*x


def test_Limit():
    assert Limit(sin(x)/x, x, 0) != 1
    assert Limit(sin(x)/x, x, 0).doit() == 1


def test_floor():
    assert limit(floor(x), x, -2, "+") == -2
    assert limit(floor(x), x, -2, "-") == -3
    assert limit(floor(x), x, -1, "+") == -1
    assert limit(floor(x), x, -1, "-") == -2
    assert limit(floor(x), x, 0, "+") == 0
    assert limit(floor(x), x, 0, "-") == -1
    assert limit(floor(x), x, 1, "+") == 1
    assert limit(floor(x), x, 1, "-") == 0
    assert limit(floor(x), x, 2, "+") == 2
    assert limit(floor(x), x, 2, "-") == 1
    assert limit(floor(x), x, 248, "+") == 248
    assert limit(floor(x), x, 248, "-") == 247


def test_floor_requires_robust_assumptions():
    assert limit(floor(sin(x)), x, 0, "+") == 0
    assert limit(floor(sin(x)), x, 0, "-") == -1
    assert limit(floor(cos(x)), x, 0, "+") == 0
    assert limit(floor(cos(x)), x, 0, "-") == 0
    assert limit(floor(5 + sin(x)), x, 0, "+") == 5
    assert limit(floor(5 + sin(x)), x, 0, "-") == 4
    assert limit(floor(5 + cos(x)), x, 0, "+") == 5
    assert limit(floor(5 + cos(x)), x, 0, "-") == 5


def test_ceiling():
    assert limit(ceiling(x), x, -2, "+") == -1
    assert limit(ceiling(x), x, -2, "-") == -2
    assert limit(ceiling(x), x, -1, "+") == 0
    assert limit(ceiling(x), x, -1, "-") == -1
    assert limit(ceiling(x), x, 0, "+") == 1
    assert limit(ceiling(x), x, 0, "-") == 0
    assert limit(ceiling(x), x, 1, "+") == 2
    assert limit(ceiling(x), x, 1, "-") == 1
    assert limit(ceiling(x), x, 2, "+") == 3
    assert limit(ceiling(x), x, 2, "-") == 2
    assert limit(ceiling(x), x, 248, "+") == 249
    assert limit(ceiling(x), x, 248, "-") == 248


def test_ceiling_requires_robust_assumptions():
    assert limit(ceiling(sin(x)), x, 0, "+") == 1
    assert limit(ceiling(sin(x)), x, 0, "-") == 0
    assert limit(ceiling(cos(x)), x, 0, "+") == 1
    assert limit(ceiling(cos(x)), x, 0, "-") == 1
    assert limit(ceiling(5 + sin(x)), x, 0, "+") == 6
    assert limit(ceiling(5 + sin(x)), x, 0, "-") == 5
    assert limit(ceiling(5 + cos(x)), x, 0, "+") == 6
    assert limit(ceiling(5 + cos(x)), x, 0, "-") == 6


def test_atan():
    assert limit(atan(x)*sin(1/x), x, 0) == 0
    assert limit(atan(x) + sqrt(x + 1) - sqrt(x), x, oo) == pi/2


def test_abs():
    assert limit(abs(x), x, 0) == 0
    assert limit(abs(sin(x)), x, 0) == 0
    assert limit(abs(cos(x)), x, 0) == 1
    assert limit(abs(sin(x + 1)), x, 0) == sin(1)


def test_heuristic():
    x = Symbol("x", extended_real=True)
    assert heuristics(sin(1/x) + atan(x), x, 0, '+') == sin(oo)
    assert heuristics(log(2 + sqrt(atan(x))*sin(1/x)), x, 0, '+') == log(2)


def test_issue_3871():
    z = Symbol("z", positive=True)
    f = -1/z*exp(-z*x)
    assert limit(f, x, oo) == 0
    assert f.limit(x, oo) == 0


def test_exponential():
    n = Symbol('n')
    x = Symbol('x', extended_real=True)
    assert limit((1 + x/n)**n, n, oo) == exp(x)
    assert limit((1 + x/(2*n))**n, n, oo) == exp(x/2)
    assert limit((1 + x/(2*n + 1))**n, n, oo) == exp(x/2)
    assert limit(((x - 1)/(x + 1))**x, x, oo) == exp(-2)
    assert limit(1 + (1 + 1/x)**x, x, oo) == 1 + S.Exp1


@pytest.mark.xfail
def test_exponential2():
    n = Symbol('n')
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


def test_issue_3792():
    assert limit( (1 - cos(x))/x**2, x, Rational(1, 2)) == 4 - 4*cos(Rational(1, 2))
    assert limit(sin(sin(x + 1) + 1), x, 0) == sin(1 + sin(1))
    assert limit(abs(sin(x + 1) + 1), x, 0) == 1 + sin(1)


def test_issue_4090():
    assert limit(1/(x + 3), x, 2) == Rational(1, 5)
    assert limit(1/(x + pi), x, 2) == Integer(1)/(2 + pi)
    assert limit(log(x)/(x**2 + 3), x, 2) == log(2)/7
    assert limit(log(x)/(x**2 + pi), x, 2) == log(2)/(4 + pi)


def test_issue_4547():
    assert limit(cot(x), x, 0, dir='+') == oo
    assert limit(cot(x), x, pi/2, dir='+') == 0


def test_issue_5164():
    assert limit(x**0.5, x, oo) == oo**0.5 == oo
    assert limit(x**0.5, x, 16) == Integer(2)**2.0
    assert limit(x**0.5, x, 0) == 0
    assert limit(x**(-0.5), x, oo) == 0
    assert limit(x**(-0.5), x, 4) == Integer(2)**(-1.0)


def test_issue_5183():
    # using list(...) so py.test can recalculate values
    tests = list(cartes([x, -x],
                        [-1, 1],
                        [2, 3, Rational(1, 2), Rational(2, 3)],
                        ['-', '+']))
    results = (oo, oo, -oo, oo, -oo*I, oo, -oo*sign((-1)**Rational(1, 3)), oo,
               0, 0, 0, 0, 0, 0, 0, 0,
               oo, oo, oo, -oo, oo, -oo*I, oo, -oo*sign((-1)**Rational(1, 3)),
               0, 0, 0, 0, 0, 0, 0, 0)
    assert len(tests) == len(results)
    for i, (args, res) in enumerate(zip(tests, results)):
        y, s, e, d = args
        eq = y**(s*e)
        assert limit(eq, x, 0, dir=d) == res


def test_issue_5184():
    assert limit(sin(x)/x, x, oo) == 0
    assert limit(atan(x), x, oo) == pi/2
    assert limit(gamma(x), x, oo) == oo
    assert limit(cos(x)/x, x, oo) == 0
    assert limit(gamma(x), x, Rational(1, 2)) == sqrt(pi)

    r = Symbol('r', extended_real=True, finite=True)
    assert limit(r*sin(1/r), r, 0) == 0


def test_issue_5229():
    assert limit((1 + y)**(1/y) - S.Exp1, y, 0) == 0


def test_issue_4546():
    # using list(...) so py.test can recalculate values
    tests = list(cartes([cot, tan],
                        [-pi/2, 0, pi/2, pi, 3*pi/2],
                        ['-', '+']))
    results = (0, 0, -oo, oo, 0, 0, -oo, oo, 0, 0,
               oo, -oo, 0, 0, oo, -oo, 0, 0, oo, -oo)
    assert len(tests) == len(results)
    for i, (args, res) in enumerate(zip(tests, results)):
        f, l, d = args
        eq = f(x)
        assert limit(eq, x, l, dir=d) == res


def test_issue_3934():
    assert limit((1 + x**log(3))**(1/x), x, 0) == 1
    assert limit((5**(1/x) + 3**(1/x))**x, x, 0) == 5


def test_compute_leading_term():
    # needs series to go to n = 32
    assert limit(x**Rational(77, 3)/(1 + x**Rational(77, 3)), x, oo) == 1
    # needs series to go to n = 128
    assert limit(x**101.1/(1 + x**101.1), x, oo) == 1


def test_issue_5955():
    assert limit((x**16)/(1 + x**16), x, oo) == 1
    assert limit((x**100)/(1 + x**100), x, oo) == 1
    assert limit((x**1885)/(1 + x**1885), x, oo) == 1
    assert limit((x**1000/((x + 1)**1000 + exp(-x))), x, oo) == 1


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
    assert Order(x)*oo != Order(1, x)
    assert limit(oo/(x**2 - 4), x, oo) == oo


def test_issue_5436():
    limit(exp(x*y), x, oo)
    limit(exp(-x*y), x, oo)


def test_Limit_dir():
    pytest.raises(TypeError, lambda: Limit(x, x, 0, dir=0))
    pytest.raises(ValueError, lambda: Limit(x, x, 0, dir='0'))


def test_polynomial():
    assert limit((x + 1)**1000/((x + 1)**1000 + 1), x, oo) == 1
    assert limit((x + 1)**1000/((x + 1)**1000 + 1), x, -oo) == 1


def test_rational():
    assert limit(1/y - (1/(y + x) + x/(y + x)/y)/z, x, oo) == (z - 1)/(y*z)
    assert limit(1/y - (1/(y + x) + x/(y + x)/y)/z, x, -oo) == (z - 1)/(y*z)


def test_issue_5740():
    assert limit(log(x)*z - log(2*x)*y, x, 0) == oo*sign(y - z)


def test_issue_6366():
    n = Symbol('n', integer=True, positive=True)
    r = (n + 1)*x**(n + 1)/(x**(n + 1) - 1) - x/(x - 1)
    assert limit(r, x, 1).simplify() == n/2


def test_factorial():
    from sympy import factorial, E
    f = factorial(x)
    assert limit(f, x, oo) == oo
    assert limit(x/f, x, oo) == 0
    # see Stirling's approximation:
    # http://en.wikipedia.org/wiki/Stirling's_approximation
    assert limit(f/(sqrt(2*pi*x)*(x/E)**x), x, oo) == 1
    assert limit(f, x, -oo) == factorial(-oo)
    assert (limit(f, x, x**2) - factorial(x**2)).simplify() == 0
    assert (limit(f, x, -x**2) - factorial(-x**2)).simplify() == 0


def test_issue_6560():
    e = 5*x**3/4 - 3*x/4 + (y*(3*x**2/2 - Rational(1, 2)) +
        35*x**4/8 - 15*x**2/4 + Rational(3, 8))/(2*(y + 1))
    assert limit(e, y, oo) == (5*x**3 + 3*x**2 - 3*x - 1)/4


def test_issue_5172():
    n = Symbol('n')
    r = Symbol('r', positive=True)
    c = Symbol('c')
    p = Symbol('p', positive=True)
    m = Symbol('m', negative=True)
    expr = ((2*n*(n - r + 1)/(n + r*(n - r + 1)))**c +
        (r - 1)*(n*(n - r + 2)/(n + r*(n - r + 1)))**c - n)/(n**c - n)
    expr = expr.subs(c, c + 1)
    assert limit(expr.subs(c, m), n, oo) == 1
    assert limit(expr.subs(c, p), n, oo).simplify() == \
        (2**(p + 1) + r - 1)/(r + 1)**(p + 1)


def test_issue_7088():
    a = Symbol('a')
    assert limit(sqrt(x/(x + a)), x, oo) == 1


def test_issue_6364():
    a = Symbol('a')
    e = z/(1 - sqrt(1 + z)*sin(a)**2 - sqrt(1 - z)*cos(a)**2)
    assert (limit(e, z, 0) - 2/cos(2*a)).simplify() == 0


def test_issue_4099():
    a = Symbol('a')
    assert limit(a/x, x, 0) == oo*sign(a)
    assert limit(-a/x, x, 0) == -oo*sign(a)
    assert limit(-a*x, x, oo) == -oo*sign(a)
    assert limit(a*x, x, oo) == oo*sign(a)


def test_issue_4503():
    dx = Symbol('dx')
    assert limit((sqrt(1 + exp(x + dx)) - sqrt(1 + exp(x)))/dx, dx, 0) == \
        exp(x)/(2*sqrt(exp(x) + 1))


def test_issue_8730():
    assert limit(subfactorial(x), x, oo) == oo


def test_omgissue_55():
    assert limit((x + exp(x))/(x - 1), x, -oo) == 1
    assert limit((x*exp(x))/(exp(x) - 1), x, -oo) == 0  # issue 2929


def test_issue_8061():
    assert limit(4**(acos(1/(1 + x**2))**2)/log(1 + x, 4), x, 0) == oo


def test_issue_8229():
    assert limit((x**Rational(1, 4) - 2)/(sqrt(x) - 4)**Rational(2, 3),
                 x, 16) == 0


def test_issue_9205():
    assert Limit(x, x, a).free_symbols == {a}
    assert Limit(x, x, a, '-').free_symbols == {a}
    assert Limit(x + y, x + y, a).free_symbols == {a}
    assert Limit(-x**2 + y, x**2, a).free_symbols == {y, a}
