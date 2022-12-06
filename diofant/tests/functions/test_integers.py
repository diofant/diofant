import pytest

from diofant import (E, Float, I, Rational, Symbol, ceiling, exp, factorial,
                     false, floor, log, nan, oo, pi, sin, sqrt, symbols, true)
from diofant.abc import x, y


__all__ = ()


def test_floor():
    i = Symbol('i', imaginary=True)
    y = Symbol('y', extended_real=True)
    r = Symbol('r', real=True)
    k, n = symbols('k,n', integer=True)

    assert floor(y).is_extended_real
    assert floor(x).is_extended_real is None
    assert floor(r).is_finite
    assert floor(y).is_finite is None
    assert floor(r).is_integer
    assert floor(y).is_integer is None

    assert floor(nan) == nan

    assert floor(oo) == oo
    assert floor(-oo) == -oo

    assert floor(0) == 0

    assert floor(1) == 1
    assert floor(-1) == -1

    assert floor(E) == 2
    assert floor(-E) == -3

    assert floor(2*E) == 5
    assert floor(-2*E) == -6

    assert floor(pi) == 3
    assert floor(-pi) == -4

    assert floor(Rational(1, 2)) == 0
    assert floor(-Rational(1, 2)) == -1

    assert floor(Rational(7, 3)) == 2
    assert floor(-Rational(7, 3)) == -3

    assert floor(Float(17.0)) == 17
    assert floor(-Float(17.0)) == -17

    assert floor(Float(7.69)) == 7
    assert floor(-Float(7.69)) == -8

    assert floor(I) == I
    assert floor(-I) == -I
    e = floor(i)
    assert isinstance(e, floor)
    assert e.args[0] == i

    assert floor(oo*I) == oo*I
    assert floor(-oo*I) == -oo*I

    assert floor(2*I) == 2*I
    assert floor(-2*I) == -2*I

    assert floor(I/2) == 0
    assert floor(-I/2) == -I

    assert floor(E + 17) == 19
    assert floor(pi + 2) == 5

    assert floor(E + pi) == floor(E + pi)
    assert floor(I + pi) == floor(I + pi)

    assert floor(floor(pi)) == 3
    assert floor(floor(y)) == floor(y)
    assert floor(floor(x)) == floor(floor(x))

    assert floor(x) == floor(x)
    assert floor(2*x) == floor(2*x)
    assert floor(k*x) == floor(k*x)

    assert floor(k) == k
    assert floor(2*k) == 2*k
    assert floor(k*n) == k*n

    assert floor(k/2) == floor(k/2)

    assert floor(x + y) == floor(x + y)

    assert floor(x + 3) == floor(x + 3)
    assert floor(x + k) == floor(x + k)

    assert floor(y + 3) == floor(y) + 3
    assert floor(y + k) == floor(y) + k

    assert floor(3 + I*y + pi) == 6 + floor(y)*I

    assert floor(k + n) == k + n

    assert floor(x*I) == floor(x*I)
    assert floor(k*I) == k*I

    assert floor(Rational(23, 10) - E*I) == 2 - 3*I

    assert floor(sin(1)) == 0
    assert floor(sin(-1)) == -1

    assert floor(exp(2)) == 7

    assert floor(log(8)/log(2)) != 2
    assert int(floor(log(8)/log(2)).evalf(chop=True)) == 3

    assert floor(factorial(50)/exp(1)) == \
        11188719610782480504630258070757734324011354208865721592720336800

    assert (floor(y) <= y) is true
    assert (floor(y) > y) is false
    assert (floor(x) <= x).is_Relational  # x could be non-real
    assert (floor(x) > x).is_Relational
    assert (floor(x) <= y).is_Relational  # arg is not same as rhs
    assert (floor(x) > y).is_Relational

    assert floor(x).as_leading_term(x) == floor(x)

    # issue sympy/sympy#11207
    assert floor(floor(x)) == floor(x)
    assert floor(ceiling(x)) == ceiling(x)


def test_ceiling():
    i = Symbol('i', imaginary=True)
    y = Symbol('y', extended_real=True)
    k, n = symbols('k,n', integer=True)

    assert ceiling(nan) == nan

    assert ceiling(oo) == oo
    assert ceiling(-oo) == -oo

    assert ceiling(0) == 0

    assert ceiling(1) == 1
    assert ceiling(-1) == -1

    assert ceiling(E) == 3
    assert ceiling(-E) == -2

    assert ceiling(2*E) == 6
    assert ceiling(-2*E) == -5

    assert ceiling(pi) == 4
    assert ceiling(-pi) == -3

    assert ceiling(Rational(1, 2)) == 1
    assert ceiling(-Rational(1, 2)) == 0

    assert ceiling(Rational(7, 3)) == 3
    assert ceiling(-Rational(7, 3)) == -2

    assert ceiling(Float(17.0)) == 17
    assert ceiling(-Float(17.0)) == -17

    assert ceiling(Float(7.69)) == 8
    assert ceiling(-Float(7.69)) == -7

    assert ceiling(I) == I
    assert ceiling(-I) == -I
    e = ceiling(i)
    assert isinstance(e, ceiling)
    assert e.args[0] == i

    assert ceiling(oo*I) == oo*I
    assert ceiling(-oo*I) == -oo*I

    assert ceiling(2*I) == 2*I
    assert ceiling(-2*I) == -2*I

    assert ceiling(I/2) == I
    assert ceiling(-I/2) == 0

    assert ceiling(E + 17) == 20
    assert ceiling(pi + 2) == 6

    assert ceiling(E + pi) == ceiling(E + pi)
    assert ceiling(I + pi) == ceiling(I + pi)

    assert ceiling(ceiling(pi)) == 4
    assert ceiling(ceiling(y)) == ceiling(y)
    assert ceiling(ceiling(x)) == ceiling(ceiling(x))

    assert ceiling(x) == ceiling(x)
    assert ceiling(2*x) == ceiling(2*x)
    assert ceiling(k*x) == ceiling(k*x)

    assert ceiling(k) == k
    assert ceiling(2*k) == 2*k
    assert ceiling(k*n) == k*n

    assert ceiling(k/2) == ceiling(k/2)

    assert ceiling(x + y) == ceiling(x + y)

    assert ceiling(x + 3) == ceiling(x + 3)
    assert ceiling(x + k) == ceiling(x + k)

    assert ceiling(y + 3) == ceiling(y) + 3
    assert ceiling(y + k) == ceiling(y) + k

    assert ceiling(3 + pi + y*I) == 7 + ceiling(y)*I

    assert ceiling(k + n) == k + n

    assert ceiling(x*I) == ceiling(x*I)
    assert ceiling(k*I) == k*I

    assert ceiling(Rational(23, 10) - E*I) == 3 - 2*I

    assert ceiling(sin(1)) == 1
    assert ceiling(sin(-1)) == 0

    assert ceiling(exp(2)) == 8

    assert ceiling(-log(8)/log(2)) != -2
    assert int(ceiling(-log(8)/log(2)).evalf(chop=True)) == -3

    assert ceiling(factorial(50)/exp(1)) == \
        11188719610782480504630258070757734324011354208865721592720336801

    assert (ceiling(y) >= y) is true
    assert (ceiling(y) < y) is false
    assert (ceiling(x) >= x).is_Relational  # x could be non-real
    assert (ceiling(x) < x).is_Relational
    assert (ceiling(x) >= y).is_Relational  # arg is not same as rhs
    assert (ceiling(x) < y).is_Relational

    # issue sympy/sympy#11207
    assert ceiling(floor(x)) == floor(x)
    assert ceiling(ceiling(x)) == ceiling(x)


def test_series():
    assert floor(x).series(x, y, 100) == floor(y)
    assert ceiling(x).series(x, y, 100) == ceiling(y)
    assert floor(x).series(x, pi, 100) == 3
    assert ceiling(x).series(x, pi, 100) == 4
    assert floor(x).series(x, n=100) == 0
    assert ceiling(x).series(x, n=100) == 1
    assert floor(-x).series(x, n=100) == -1
    assert ceiling(-x).series(x, n=100) == 0


@pytest.mark.xfail
def test_sympyissue_4149():
    y = Symbol('y', real=True)
    assert floor(3 + pi*I + y*I) == 3 + floor(pi + y)*I
    assert floor(3*I + pi*I + y*I) == floor(3 + pi + y)*I
    assert floor(3 + E + pi*I + y*I) == 5 + floor(pi + y)*I


def test_issue_1055():
    e = 1/(sqrt(2) - 1) - sqrt(2)
    e1 = e - Rational(1, 10**15)
    e2 = e - 1
    e3 = e - Rational(1, 10**1000)

    assert floor(e) == 1
    assert floor(e1) == 0
    assert floor(e2) == e2
    assert floor(e3) == floor(e3, evaluate=False)
