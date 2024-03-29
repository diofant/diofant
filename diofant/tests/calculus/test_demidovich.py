# Numbers listed with the tests refer to problem numbers in books
# :cite:`demidovich1970problems` and :cite:`anti-demidovich2001` (with '*').

from diofant import (Rational, asin, cos, exp, limit, log, oo, pi, root, sign,
                     sin, sqrt, tan)
from diofant.abc import a, h, k, m, n, x


__all__ = ()


def test_Limits_simple_0():
    assert limit((2**(x + 1) + 3**(x + 1))/(2**x + 3**x), x, oo) == 3  # 175


def test_Limits_simple_1():
    assert limit((n + 1)*(n + 2)*(n + 3)/n**3, n, oo) == 1  # 172
    assert limit(sqrt(n + 1) - sqrt(n), n, oo) == 0  # 179
    assert limit((2*x - 3)*(3*x + 5)*(4*x - 6)/(3*x**3 + x - 1),
                 x, oo) == 8  # Example 1
    assert limit(x/root(x**3 + 10, 3), x, oo) == 1  # Example 2
    assert limit((x + 1)**2/(x**2 + 1), x, oo) == 1  # 181


def test_Limits_simple_2():
    assert limit(1000*x/(x**2 - 1), x, oo) == 0  # 182
    assert limit((x**2 - 5*x + 1)/(3*x + 7), x, oo) == oo  # 183
    assert limit((2*x**2 - x + 3)/(x**3 - 8*x + 5), x, oo) == 0  # 184
    assert limit((2*x**2 - 3*x - 4)/sqrt(x**4 + 1), x, oo) == 2  # 186
    assert limit((2*x + 3)/(x + root(x, 3)), x, oo) == 2  # 187
    assert limit(x**2/(10 + x*sqrt(x)), x, oo) == oo  # 188
    assert limit(root(x**2 + 1, 3)/(x + 1), x, oo) == 0  # 189
    assert limit(sqrt(x)/sqrt(x + sqrt(x + sqrt(x))), x, oo) == 1  # 190


def test_Limits_simple_3a():
    # issue sympy/sympy#3513
    assert limit((x**2 - (a + 1)*x + a)/(x**3 - a**3),
                 x, a) == (a - 1)/(3*a**2)  # 196


def test_Limits_simple_3b():
    assert limit(((x + h)**3 - x**3)/h, h, 0) == 3*x**2  # 197
    assert limit(1/(1 - x) - 3/(1 - x**3), x, 1) == -1  # 198
    assert limit((sqrt(1 + x) - 1)/(root(1 + x, 3) - 1),
                 x, 0) == Rational(3, 2)  # Example 4
    assert limit((sqrt(x) - 1)/(x - 1), x, 1) == Rational(1, 2)  # 199
    assert limit((sqrt(x) - 8)/(root(x, 3) - 4), x, 64) == 3  # 200
    assert limit((root(x, 3) - 1)/(root(x, 4) - 1),
                 x, 1) == Rational(4, 3)  # 201
    assert limit((root(x**2, 3) - 2*root(x, 3) + 1)/(x - 1)**2,
                 x, 1) == Rational(1, 9)  # 202


def test_Limits_simple_4a():
    assert limit((sqrt(x) - sqrt(a))/(x - a),
                 x, a) == 1/(2*sqrt(a))  # Example 5
    assert limit((sqrt(x) - 1)/(root(x, 3) - 1), x, 1) == Rational(3, 2)  # 205
    assert limit((sqrt(1 + x) - sqrt(1 - x))/x, x, 0) == 1  # 207
    assert limit(sqrt(x**2 - 5*x + 6) - x, x, oo) == -Rational(5, 2)  # 213


def test_limits_simple_4aa():
    assert limit(x*(sqrt(x**2 + 1) - x), x, oo) == Rational(1, 2)  # 214


def test_Limits_simple_4b():
    # issue sympy/sympy#3511
    assert limit(x - root(x**3 - 1, 3), x, oo) == 0
    assert limit(x + root(1 - x**3, 3),
                 x, oo) == oo*sign(1 + root(-1, 3))  # 215


def test_Limits_simple_4c():
    assert limit(log(1 + exp(x))/x, x, -oo) == 0  # 267a
    assert limit(log(1 + exp(x))/x, x, oo) == 1  # 267b


def test_bounded():
    assert limit(sin(x)/x, x, oo) == 0  # 216b
    assert limit(x*sin(1/x), x, 0) == 0  # 227a


def test_f1a():
    # issue sympy/sympy#3508, Example 7
    assert limit((sin(2*x)/x)**(1 + x), x, 0) == 2


def test_f1a2():
    # issue sympy/sympy#3509, Example 9
    assert limit(((x - 1)/(x + 1))**x, x, oo) == exp(-2)
    assert limit((1 + k/x)**x, x, oo) == exp(k)


def test_f1b():
    assert limit(sin(x)/x, x, 2) == sin(2)/2  # 216a
    assert limit(sin(3*x)/x, x, 0) == 3  # 217
    assert limit(sin(5*x)/sin(2*x), x, 0) == Rational(5, 2)  # 218
    assert limit(sin(pi*x)/sin(3*pi*x), x, 0) == Rational(1, 3)  # 219
    assert limit(n*sin(pi/n), n, oo) == pi  # 220
    assert limit((1 - cos(x))/x**2, x, 0) == Rational(1, 2)  # 221
    assert limit(x*sin(1/x), x, oo) == 1  # 227b
    assert limit((cos(m*x) - cos(n*x))/x**2, x, 0) == ((n**2 - m**2)/2)  # 232
    assert limit((tan(x) - sin(x))/x**3, x, 0) == Rational(1, 2)  # 233
    assert limit((x - sin(2*x))/(x + sin(3*x)), x, 0) == -Rational(1, 4)  # 237
    assert limit((1 - sqrt(cos(x)))/x**2, x, 0) == Rational(1, 4)  # 239
    assert limit((sqrt(1 + sin(x)) - sqrt(1 - sin(x)))/x, x, 0) == 1  # 240

    assert limit((sin(x) - sin(a))/(x - a), x, a) == cos(a)  # 222, *176
    assert limit((cos(x) - cos(a))/(x - a), x, a) == -sin(a)  # 223
    assert limit((sin(x + h) - sin(x))/h, h, 0) == cos(x)  # 225


def test_f2a():
    assert limit(((x + 1)/(2*x + 1))**x**2, x, oo) == 0  # Primer 8


def test_f2():
    assert limit((sqrt(cos(x)) - root(cos(x), 3))/sin(x)**2,
                 x, 0) == -Rational(1, 12)  # *184


def test_f3():
    # issue sympy/sympy#3504
    assert limit(asin(a*x)/x, x, 0) == a  # *222
