import pytest

from diofant import FF, Rational, erf, log, pi, ratsimp, ratsimpmodprime, sqrt
from diofant.abc import a, b, c, d, e, t, x, y, z


__all__ = ()


def test_ratsimp():
    f, g = 1/x + 1/y, (x + y)/(x*y)

    assert f != g
    assert ratsimp(f) == g

    f, g = 1/(1 + 1/x), 1 - 1/(x + 1)

    assert f != g
    assert ratsimp(f) == g

    f, g = x/(x + y) + y/(x + y), 1

    assert f != g
    assert ratsimp(f) == g

    f, g = -x - y - y**2/(x + y) + x**2/(x + y), -2*y

    assert f != g
    assert ratsimp(f) == g

    f = (a*c*x*y + a*c*z - b*d*x*y - b*d*z - b*t*x*y - b*t*x - b*t*z +
         e*x)/(x*y + z)
    G = [a*c - b*d - b*t + (-b*t*x + e*x)/(x*y + z),
         a*c - b*d - b*t - ( b*t*x - e*x)/(x*y + z)]

    assert f != g
    assert ratsimp(f) in G

    A = sqrt(pi)

    B = log(erf(x) - 1)
    C = log(erf(x) + 1)

    D = 8 - 8*erf(x)

    f = A*B/D - A*C/D + A*C*erf(x)/D - A*B*erf(x)/D + 2*A/D

    assert ratsimp(f) == A*B/8 - A*C/8 - A/(4*erf(x) - 4)


def test_ratsimpmodprime():
    a = y**5 + x + y
    b = x - y
    F = [x*y**5 - x - y]
    assert ratsimpmodprime(a/b, F, x, y, order='lex') == \
        (x**2 + x*y + x + y)/(x**2 - x*y)

    a = x + y**2 - 2
    b = x + y**2 - y - 1
    F = [x*y - 1]
    assert ratsimpmodprime(a/b, F, x, y, order='lex') == \
        (x - y - 1)/(x - y)

    a = 5*x**3 + 21*x**2 + 4*x*y + 23*x + 12*y + 15
    b = 7*x**3 - y*x**2 + 31*x**2 + 2*x*y + 15*y + 37*x + 21
    F = [x**2 + y**2 - 1]
    assert ratsimpmodprime(a/b, F, x, y, order='lex') == \
        (5*x - 5*y - 1)/(6*x - 8*y)

    a = x*y - x - 2*y + 4
    b = x + y**2 - 2*y
    F = [x - 2, y - 3]
    assert ratsimpmodprime(a/b, F, x, y, order='lex') == \
        Rational(2, 5)

    # Test a bug where denominators would be dropped
    assert ratsimpmodprime(x, [y - 2*x], order='lex') == \
        y/2

    # coverage tests
    assert ratsimpmodprime(x, [x - y/2], x, y, order='lex') == 2*y
    pytest.raises(ValueError, lambda: ratsimpmodprime((x + y)/(x - y),
                                                      [y**2 - 1, x**2 - 1]))

    a = x**5 + 2*x**4 + 2*x**3 + 2*x**2 + x + 2/x + x**(-2)
    assert ratsimpmodprime(a, [x + 1], domain=FF(2)) == 1
    assert ratsimpmodprime(a, [x + 1], domain=FF(3)) == 2
