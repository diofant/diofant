"""Tests for algorithms for partial fraction decomposition of rational
functions.
"""

import pytest

from diofant import (ZZ, Dummy, E, Eq, Expr, I, Integer, Lambda, Mul, Rational,
                     RootSum, Symbol, apart, apart_list,
                     assemble_partfrac_list, factor, numbered_symbols, pi,
                     sqrt, symbols, together)
from diofant.abc import a, b, c, x, y
from diofant.polys.partfrac import apart_undetermined_coeffs


__all__ = ()


def test_apart():
    assert apart(1) == 1
    assert apart(1, x) == 1

    f, g = (x**2 + 1)/(x + 1), 2/(x + 1) + x - 1

    assert apart(f, full=False) == g
    assert apart(f, full=True) == g

    f, g = 1/(x + 2)/(x + 1), 1/(1 + x) - 1/(2 + x)

    assert apart(f, full=False) == g
    assert apart(f, full=True) == g

    f, g = 1/(x + 1)/(x + 5), -1/(5 + x)/4 + 1/(1 + x)/4

    assert apart(f, full=False) == g
    assert apart(f, full=True) == g

    assert apart((E*x + 2)/(x - pi)*(x - 1), x) == \
        2 - E + E*pi + E*x + (E*pi + 2)*(pi - 1)/(x - pi)

    assert apart(Eq((x**2 + 1)/(x + 1), x), x) == Eq(x - 1 + 2/(x + 1), x)

    assert apart(x/2, y) == x/2  # issue sympy/sympy#9123

    # issue sympy/sympy#12177

    f, g = (x + y)/(2*x - y), 3*y/(2*x - y)/2 + Rational(1, 2)

    assert apart(f, x, full=False) == g
    assert apart(f, x, full=True) == g

    f, g = (x + y)/(2*x - y), 3*x/(2*x - y) - 1

    assert apart(f, y, full=False) == g
    assert apart(f, y, full=True) == g


def test_apart_symbolic():
    f = a*x**4 + (2*b + 2*a*c)*x**3 + (4*b*c - a**2 + a*c**2)*x**2 + \
        (-2*a*b + 2*b*c**2)*x - b**2
    g = a**2*x**4 + (2*a*b + 2*c*a**2)*x**3 + (4*a*b*c + b**2 +
                                               a**2*c**2)*x**2 + (2*c*b**2 + 2*a*b*c**2)*x + b**2*c**2

    assert apart(f/g, x) == 1/a - 1/(x + c)**2 - b**2/(a*(a*x + b)**2)

    assert apart(1/((x + a)*(x + b)*(x + c)), x) == \
        1/((a - c)*(b - c)*(c + x)) - 1/((a - b)*(b - c)*(b + x)) + \
        1/((a - b)*(a - c)*(a + x))


def test_apart_extension():
    f = 2/(x**2 + 1)
    g = I/(x + I) - I/(x - I)

    assert apart(f, extension=I) == g
    assert apart(f, gaussian=True) == g

    f = x/((x - 2)*(x + I))

    assert factor(together(apart(f))) == f


def test_apart_full():
    f = 1/(x**2 + 1)

    assert apart(f, full=False) == f
    assert apart(f, full=True) == \
        -RootSum(x**2 + 1, Lambda(a, a/(x - a)), auto=False)/2

    f = 1/(x**3 + x + 1)

    assert apart(f, full=False) == f
    assert apart(f, full=True) == \
        RootSum(x**3 + x + 1,
                Lambda(a, (6*a**2/31 - 9*a/31 + Rational(4, 31))/(x - a)), auto=False)

    f = 1/(x**5 + 1)

    assert apart(f, full=False) == \
        (-Rational(1, 5))*((x**3 - 2*x**2 + 3*x - 4)/(x**4 - x**3 + x**2 -
                                                      x + 1)) + Rational(1, 5)/(x + 1)
    assert apart(f, full=True) == \
        -RootSum(x**4 - x**3 + x**2 - x + 1,
                 Lambda(a, a/(x - a)), auto=False)/5 + Rational(1, 5)/(x + 1)


def test_apart_undetermined_coeffs():
    p = (2*x - 3).as_poly()
    q = (x**9 - x**8 - x**6 + x**5 - 2*x**2 + 3*x - 1).as_poly()
    r = (-x**7 - x**6 - x**5 + 4)/(x**8 - x**5 - 2*x + 1) + 1/(x - 1)

    assert apart_undetermined_coeffs(p, q) == r

    dom = ZZ.inject(a, b)
    p = Integer(1).as_poly(x, domain=dom)
    q = ((x + a)*(x + b)).as_poly(x, domain=dom)
    r = 1/((a - b)*(b + x)) - 1/((a - b)*(a + x))

    assert apart_undetermined_coeffs(p, q) == r


def test_apart_list():
    assert apart_list(1) == 1

    w0, w1, w2 = Symbol('w0'), Symbol('w1'), Symbol('w2')
    _a = Dummy('a')

    f = (-2*x - 2*x**2) / (3*x**2 - 6*x)
    assert (apart_list(f, x, dummies=numbered_symbols('w')) ==
            (-1, Rational(2, 3).as_poly(x),
             [((w0 - 2).as_poly(), Lambda(_a, 2), Lambda(_a, -_a + x), 1)]))

    assert (apart_list(2/(x**2-2), x, dummies=numbered_symbols('w')) ==
            (1, Integer(0).as_poly(x),
             [((w0**2 - 2).as_poly(), Lambda(_a, _a/2), Lambda(_a, -_a + x), 1)]))

    f = 36 / (x**5 - 2*x**4 - 2*x**3 + 4*x**2 + x - 2)
    assert (apart_list(f, x, dummies=numbered_symbols('w')) ==
            (1, Integer(0).as_poly(x),
             [((w0 - 2).as_poly(), Lambda(_a, 4), Lambda(_a, -_a + x), 1),
              ((w1**2 - 1).as_poly(), Lambda(_a, -3*_a - 6), Lambda(_a, -_a + x), 2),
              ((w2 + 1).as_poly(), Lambda(_a, -4), Lambda(_a, -_a + x), 1)]))

    f = 1/(2*(x - 1)**2)
    assert (apart_list(f, x, dummies=numbered_symbols('w')) ==
            (1, Integer(0).as_poly(x),
             [(Integer(2).as_poly(w0), Lambda(_a, 0), Lambda(_a, x - _a), 1),
              ((w1 - 1).as_poly(), Lambda(_a, Rational(1, 2)), Lambda(_a, x - _a), 2),
              (Integer(1).as_poly(w2), Lambda(_a, 0), Lambda(_a, x - _a), 1)]))


def test_assemble_partfrac_list():
    f = 36 / (x**5 - 2*x**4 - 2*x**3 + 4*x**2 + x - 2)
    pfd = apart_list(f)
    assert assemble_partfrac_list(pfd) == -4/(x + 1) - 3/(x + 1)**2 - 9/(x - 1)**2 + 4/(x - 2)

    a = Dummy('a')
    pfd = (1, Integer(0).as_poly(x), [([sqrt(2), -sqrt(2)], Lambda(a, a/2), Lambda(a, -a + x), 1)])
    assert assemble_partfrac_list(pfd) == -1/(sqrt(2)*(x + sqrt(2))) + 1/(sqrt(2)*(x - sqrt(2)))


def test_noncommutative_pseudomultivariate():
    class Foo(Expr):
        is_commutative = False
    e = x/(x + x*y)
    c = 1/(1 + y)
    fe, fc = map(Foo, [e, c])
    assert apart(e + fe) == c + fc
    assert apart(e*fe) == c*fc

    A, B = symbols('A, B', commutative=False)
    assert apart(A*B) == A*B


def test_sympyissue_5798():
    assert apart(
        2*x/(x**2 + 1) - (x - 1)/(2*(x**2 + 1)) + 1/(2*(x + 1)) - 2/x) == \
        (3*x + 1)/(x**2 + 1)/2 + 1/(x + 1)/2 - 2/x


@pytest.mark.timeout(100)
def test_sympyissue_18531():
    e = (x**2 + 1)**3/((x - 1)**2*(x + 1)**2*(-x**2 + 2*x + 1)*(x**2 + 2*x - 1))

    def mul2(x):
        return Mul(2, x, evaluate=False)

    r = (1/mul2(x - sqrt(2) + 1) - 1/mul2(x - sqrt(2) - 1) +
         1/mul2(x + 1 + sqrt(2)) - 1/mul2(x - 1 + sqrt(2)) +
         1/mul2((x + 1)**2) + 1/mul2((x - 1)**2))

    assert apart(e, x, extension=sqrt(2)) == r


@pytest.mark.slow
def test_sympyissue_20163():
    e = 6/(x**6 + 1)
    r = ((sqrt(3) + I)/(2*x + sqrt(3) + I) -
         (-sqrt(3) + I)/(2*x + sqrt(3) - I) +
         (-sqrt(3) + I)/(2*x - sqrt(3) + I) -
         (sqrt(3) + I)/(2*x - sqrt(3) - I) + I/(x + I) - I/(x - I))

    assert apart(e, extension=sqrt(3) + I) == r
