from sympy import (S, symbols, I, atan, log, Poly, sqrt,
                   simplify, integrate, Integer, Rational)
from sympy.integrals.rationaltools import ratint, ratint_logpart, log_to_atan

from sympy.abc import a, b, x, t

half = Rational(1, 2)


def test_ratint():
    assert ratint(Integer(0), x) == 0
    assert ratint(Integer(7), x) == 7*x

    assert ratint(x, x) == x**2/2
    assert ratint(2*x, x) == x**2
    assert ratint(-2*x, x) == -x**2

    assert ratint(8*x**7 + 2*x + 1, x) == x**8 + x**2 + x

    f = Integer(1)
    g = x + 1

    assert ratint(f / g, x) == log(x + 1)
    assert ratint((f, g), x) == log(x + 1)

    f = x**3 - x
    g = x - 1

    assert ratint(f/g, x) == x**3/3 + x**2/2

    f = x
    g = (x - a)*(x + a)

    assert ratint(f/g, x) == log(x**2 - a**2)/2

    f = Integer(1)
    g = x**2 + 1

    assert ratint(f/g, x, extended_real=None) == atan(x)
    assert ratint(f/g, x, extended_real=True) == atan(x)

    assert ratint(f/g, x, extended_real=False) == I*log(x + I)/2 - I*log(x - I)/2

    f = Integer(36)
    g = x**5 - 2*x**4 - 2*x**3 + 4*x**2 + x - 2

    assert ratint(f/g, x) == \
        -4*log(x + 1) + 4*log(x - 2) + (12*x + 6)/(x**2 - 1)

    f = x**4 - 3*x**2 + 6
    g = x**6 - 5*x**4 + 5*x**2 + 4

    assert ratint(f/g, x) == \
        atan(x) + atan(x**3) + atan(x/2 - 3*x**3/2 + Rational(1, 2)*x**5)

    f = x**7 - 24*x**4 - 4*x**2 + 8*x - 8
    g = x**8 + 6*x**6 + 12*x**4 + 8*x**2

    assert ratint(f/g, x) == \
        (4 + 6*x + 8*x**2 + 3*x**3)/(4*x + 4*x**3 + x**5) + log(x)

    assert ratint((x**3*f)/(x*g), x) == \
        -(12 - 16*x + 6*x**2 - 14*x**3)/(4 + 4*x**2 + x**4) - \
        5*sqrt(2)*atan(x*sqrt(2)/2) + Rational(1, 2)*x**2 - 3*log(2 + x**2)

    f = x**5 - x**4 + 4*x**3 + x**2 - x + 5
    g = x**4 - 2*x**3 + 5*x**2 - 4*x + 4

    assert ratint(f/g, x) == \
        x + Rational(1, 2)*x**2 + Rational(1, 2)*log(2 - x + x**2) - (4*x - 9)/(14 - 7*x + 7*x**2) + \
        13*sqrt(7)*atan(-Rational(1, 7)*sqrt(7) + 2*x*sqrt(7)/7)/49

    assert ratint(1/(x**2 + x + 1), x) == \
        2*sqrt(3)*atan(sqrt(3)/3 + 2*x*sqrt(3)/3)/3

    assert ratint(1/(x**3 + 1), x) == \
        -log(1 - x + x**2)/6 + log(1 + x)/3 + sqrt(3)*atan(-sqrt(3)
             / 3 + 2*x*sqrt(3)/3)/3

    assert ratint(1/(x**2 + x + 1), x, extended_real=False) == \
        -I*3**half*log(half + x - half*I*3**half)/3 + \
        I*3**half*log(half + x + half*I*3**half)/3

    assert ratint(1/(x**3 + 1), x, extended_real=False) == log(1 + x)/3 + \
        (-Rational(1, 6) + I*3**half/6)*log(-half + x + I*3**half/2) + \
        (-Rational(1, 6) - I*3**half/6)*log(-half + x - I*3**half/2)

    # issue 4991
    assert ratint(1/(x*(a + b*x)**3), x) == \
        (3*a + 2*b*x)/(2*a**4 + 4*a**3*b*x + 2*a**2*b**2*x**2) + (
            log(x) - log(a/b + x))/a**3

    assert ratint(x/(1 - x**2), x) == -log(x**2 - 1)/2
    assert ratint(-x/(1 - x**2), x) == log(x**2 - 1)/2

    assert ratint((x/4 - 4/(1 - x)).diff(x), x) == x/4 + 4/(x - 1)

    ans = atan(x)
    assert ratint(1/(x**2 + 1), x, symbol=x) == ans
    assert ratint(1/(x**2 + 1), x, symbol='x') == ans
    assert ratint(1/(x**2 + 1), x, symbol=a) == ans


def test_ratint_logpart():
    assert ratint_logpart(x, x**2 - 9, x, t) == \
        [(Poly(x**2 - 9, x), Poly(-2*t + 1, t))]
    assert ratint_logpart(x**2, x**3 - 5, x, t) == \
        [(Poly(x**3 - 5, x), Poly(-3*t + 1, t))]


def test_issue_5414():
    assert ratint(1/(x**2 + 16), x) == atan(x/4)/4


def test_issue_5249():
    assert ratint(
        1/(x**2 + a**2), x) == (-I*log(-I*a + x)/2 + I*log(I*a + x)/2)/a


def test_issue_5817():
    a, b, c = symbols('a,b,c', positive=True)

    assert simplify(ratint(a/(b*c*x**2 + a**2 + b*a), x)) == \
        sqrt(a)*atan(sqrt(
            b)*sqrt(c)*x/(sqrt(a)*sqrt(a + b)))/(sqrt(b)*sqrt(c)*sqrt(a + b))


def test_issue_5981():
    u = symbols('u')
    assert integrate(1/(u**2 + 1)) == atan(u)


def test_log_to_atan():
    f, g = (Poly(x + Rational(1, 2), x, domain='QQ'), Poly(sqrt(3)/2, x, domain='EX'))
    fg_ans = 2*atan(2*sqrt(3)*x/3 + sqrt(3)/3)
    assert log_to_atan(f, g) == fg_ans
    assert log_to_atan(g, f) == -fg_ans
