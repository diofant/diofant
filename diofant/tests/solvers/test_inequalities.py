"""Tests for tools for solving inequalities and systems of inequalities."""

import pytest

from diofant import (E, Eq, FiniteSet, Float, Integer, Integral, Interval, Max,
                     Min, Ne, Piecewise, PurePoly, Rational, RootOf, S, Symbol,
                     Union, false, log, oo, pi, reduce_inequalities, root, sin,
                     solve, sqrt, true)
from diofant.abc import x, y
from diofant.solvers.inequalities import (reduce_piecewise_inequality,
                                          reduce_rational_inequalities,
                                          solve_poly_inequalities)
from diofant.solvers.inequalities import solve_poly_inequality as psolve
from diofant.solvers.inequalities import solve_univariate_inequality as isolve


__all__ = ()

inf = oo.evalf()


def test_solve_linear_inequalities():
    eqs = [x >= 0, 2*x + 4*y <= 14, x - 2*y <= 1]
    ans = ((x >= Integer(0)) & (x <= Integer(4)) &
           (y >= x/2 - Rational(1, 2)) & (y <= -x/2 + Rational(7, 2)))

    assert reduce_inequalities(eqs) == ans

    eqs = [x + y >= 4, x <= 1, y <= 1]

    assert reduce_inequalities(eqs) == false

    eqs = [x + 2*y <= 3, 2*x + y <= 5]

    assert reduce_inequalities(eqs) == (y <= Min(-2*x + 5,
                                                 -x/2 + Rational(3, 2)))

    eqs = [x + 2*y < 3, 2*x + y < 5]

    assert reduce_inequalities(eqs) == (y < Min(-2*x + 5,
                                                -x/2 + Rational(3, 2)))

    eqs = [x + 2*y <= 3, 2*x + y < 5]
    ans = (((y <= -x/2 + Rational(3, 2)) & (-x/2 + Rational(3, 2) < -2*x + 5)) |
           ((y < -2*x + 5) & (-2*x + 5 <= -x/2 + Rational(3, 2))))

    assert reduce_inequalities(eqs) == ans

    eqs = [x + 2*y >= 3, 2*x + y >= 5]

    assert reduce_inequalities(eqs) == (y >= Max(-2*x + 5,
                                                 -x/2 + Rational(3, 2)))

    eqs = [x + 2*y > 3, 2*x + y > 5]

    assert reduce_inequalities(eqs) == (y > Max(-2*x + 5,
                                                -x/2 + Rational(3, 2)))

    eqs = [x + 2*y >= 3, 2*x + y > 5]
    ans = (((y >= -x/2 + Rational(3, 2)) & (-x/2 + Rational(3, 2) > -2*x + 5)) |
           ((y > -2*x + 5) & (-2*x + 5 >= -x/2 + Rational(3, 2))))

    assert reduce_inequalities(eqs) == ans


def test_solve_poly_inequality():
    assert psolve(Integer(0).as_poly(x), '==') == [S.ExtendedReals]
    assert psolve(Integer(1).as_poly(x), '==') == [S.EmptySet]
    assert psolve(PurePoly(x + 1, x), '>') == [Interval(-1, oo, True, False)]
    pytest.raises(ValueError, lambda: psolve(x, '=='))
    pytest.raises(ValueError, lambda: psolve(x.as_poly(), '??'))

    assert (solve_poly_inequalities((((x**2 - 3).as_poly(), '>'),
                                     ((-x**2 + 1).as_poly(), '>'))) ==
            Union(Interval(-oo, -sqrt(3), False, True),
                  Interval(-1, 1, True, True),
                  Interval(sqrt(3), oo, True, False)))


def test_reduce_poly_inequalities_real_interval():
    assert reduce_rational_inequalities(
        [[Eq(x**2, 0)]], x, relational=False) == FiniteSet(0)
    assert reduce_rational_inequalities(
        [[x**2 <= 0]], x, relational=False) == FiniteSet(0)
    assert reduce_rational_inequalities(
        [[x**2 < 0]], x, relational=False) == S.EmptySet
    assert reduce_rational_inequalities(
        [[x**2 >= 0]], x, relational=False) == S.ExtendedReals
    assert reduce_rational_inequalities(
        [[x**2 > 0]], x, relational=False) == \
        FiniteSet(0).complement(S.ExtendedReals)
    assert reduce_rational_inequalities(
        [[Ne(x**2, 0)]], x, relational=False) == \
        FiniteSet(0).complement(S.ExtendedReals)

    assert reduce_rational_inequalities(
        [[Eq(x**2, 1)]], x, relational=False) == FiniteSet(-1, 1)
    assert reduce_rational_inequalities(
        [[x**2 <= 1]], x, relational=False) == Interval(-1, 1)
    assert reduce_rational_inequalities(
        [[x**2 < 1]], x, relational=False) == Interval(-1, 1, True, True)
    assert reduce_rational_inequalities(
        [[x**2 >= 1]], x, relational=False) == \
        Union(Interval(-oo, -1, False), Interval(1, oo, False, False))
    assert reduce_rational_inequalities(
        [[x**2 > 1]], x, relational=False) == \
        Interval(-1, 1).complement(S.ExtendedReals)
    assert reduce_rational_inequalities(
        [[Ne(x**2, 1)]], x, relational=False) == \
        FiniteSet(-1, 1).complement(S.ExtendedReals)
    assert reduce_rational_inequalities([[Eq(
        x**2, 1.0)]], x, relational=False) == FiniteSet(-1.0, 1.0).evalf()
    assert reduce_rational_inequalities(
        [[x**2 <= 1.0]], x, relational=False) == Interval(-1.0, 1.0)
    assert reduce_rational_inequalities([[
        x**2 < 1.0]], x, relational=False) == Interval(-1.0, 1.0, True, True)
    assert reduce_rational_inequalities(
        [[x**2 >= 1.0]], x, relational=False) == \
        Union(Interval(-inf, -1.0), Interval(1.0, inf))
    assert reduce_rational_inequalities(
        [[x**2 > 1.0]], x, relational=False) == \
        Union(Interval(-inf, -1.0, False, True),
              Interval(1.0, inf, True))
    assert reduce_rational_inequalities([[Ne(
        x**2, 1.0)]], x, relational=False) == \
        FiniteSet(-1.0, 1.0).complement(S.ExtendedReals)

    s = sqrt(2)

    assert reduce_rational_inequalities([[
        x**2 - 1 < 0, x**2 - 1 > 0]], x, relational=False) == S.EmptySet
    assert reduce_rational_inequalities([[x**2 - 1 <= 0,
                                          x**2 - 1 >= 0]], x,
                                        relational=False) == FiniteSet(-1, 1)
    assert reduce_rational_inequalities(
        [[x**2 - 2 <= 0, x**2 - 1 >= 0]], x, relational=False
    ) == Union(Interval(-s, -1, False, False), Interval(1, s, False, False))
    assert reduce_rational_inequalities(
        [[x**2 - 2 <= 0, x**2 - 1 > 0]], x, relational=False
    ) == Union(Interval(-s, -1, False, True), Interval(1, s, True, False))
    assert reduce_rational_inequalities(
        [[x**2 - 2 < 0, x**2 - 1 >= 0]], x, relational=False
    ) == Union(Interval(-s, -1, True, False), Interval(1, s, False, True))
    assert reduce_rational_inequalities(
        [[x**2 - 2 < 0, x**2 - 1 > 0]], x, relational=False
    ) == Union(Interval(-s, -1, True, True), Interval(1, s, True, True))
    assert reduce_rational_inequalities(
        [[x**2 - 2 < 0, Ne(x**2 - 1, 0)]], x, relational=False
    ) == Union(Interval(-s, -1, True, True), Interval(-1, 1, True, True),
               Interval(1, s, True, True))

    # issue sympy/sympy#10237
    assert reduce_rational_inequalities(
        [[x < oo, x >= 0, -oo < x]], x, relational=False) == Interval(0, oo, False, True)

    assert reduce_rational_inequalities([[Eq((x + 1)/(x**2 - 1),
                                             0)]], x) is false


def test_reduce_poly_inequalities_complex_relational():
    assert reduce_rational_inequalities(
        [[Eq(x**2, 0)]], x, relational=True) == Eq(x, 0)
    assert reduce_rational_inequalities(
        [[x**2 <= 0]], x, relational=True) == Eq(x, 0)
    assert reduce_rational_inequalities(
        [[x**2 < 0]], x, relational=True) is false
    assert reduce_rational_inequalities(
        [[x**2 >= 0]], x, relational=True) is true
    assert reduce_rational_inequalities(
        [[x**2 > 0]], x, relational=True) == \
        (x < 0) | (Integer(0) < x)
    assert reduce_rational_inequalities(
        [[Ne(x**2, 0)]], x, relational=True) == \
        (x < 0) | (Integer(0) < x)

    for one in (Integer(1), Float(1.0)):
        assert reduce_rational_inequalities(
            [[Eq(x**2, one)]], x, relational=True) == \
            Eq(x, -one) | Eq(x, one)
        assert reduce_rational_inequalities(
            [[x**2 <= one]], x, relational=True) == \
            (-one <= x) & (x <= one)
        assert reduce_rational_inequalities(
            [[x**2 < one]], x, relational=True) == \
            (-one < x) & (x < one)
        assert reduce_rational_inequalities(
            [[x**2 >= one]], x, relational=True) == \
            (one <= x) | (x <= -one)
        assert reduce_rational_inequalities(
            [[x**2 > one]], x, relational=True) == \
            (one < x) | (x < -one)
        assert reduce_rational_inequalities(
            [[Ne(x**2, one)]], x, relational=True) == \
            (x < -one) | ((-one < x) & (x < one)) | (one < x)


def test_reduce_rational_inequalities_real_relational():
    assert reduce_rational_inequalities([], x) is false
    assert reduce_rational_inequalities(
        [[(x**2 + 3*x + 2)/(x**2 - 16) >= 0]], x, relational=False) == \
        Union(Interval.Ropen(-oo, -4), Interval(-2, -1), Interval.Lopen(4, oo))

    assert reduce_rational_inequalities(
        [[((-2*x - 10)*(3 - x))/((x**2 + 5)*(x - 2)**2) < 0]], x,
        relational=False) == \
        Union(Interval.open(-5, 2), Interval.open(2, 3))

    assert reduce_rational_inequalities([[(x + 1)/(x - 5) <= 0]], x,
                                        relational=False) == \
        Interval.Ropen(-1, 5)

    assert reduce_rational_inequalities([[(x**2 + 4*x + 3)/(x - 1) > 0]], x,
                                        relational=False) == \
        Union(Interval.open(-3, -1), Interval.Lopen(1, oo))

    assert reduce_rational_inequalities([[(x**2 - 16)/(x - 1)**2 < 0]], x,
                                        relational=False) == \
        Union(Interval.open(-4, 1), Interval.open(1, 4))

    assert reduce_rational_inequalities([[(3*x + 1)/(x + 4) >= 1]], x,
                                        relational=False) == \
        Union(Interval.Ropen(-oo, -4), Interval(Rational(3, 2), oo))

    assert reduce_rational_inequalities([[(x - 8)/x <= 3 - x]], x,
                                        relational=False) == \
        Union(Interval(-oo, -2), Interval.Lopen(0, 4))


def test_reduce_piecewise_inequalities():
    e = abs(x - 5) < 3
    ans = (Integer(2) < x) & (x < 8)
    assert reduce_inequalities(e) == ans
    assert reduce_inequalities(e, x) == ans
    assert reduce_inequalities(abs(x - 5)) == Eq(x, 5)
    assert reduce_inequalities(
        abs(2*x + 3) >= 8) == ((Rational(5, 2) <= x) |
                               (x <= -Rational(11, 2)))
    assert reduce_inequalities(abs(x - 4) + abs(
        3*x - 5) < 7) == (Rational(1, 2) < x) & (x < 4)
    assert reduce_inequalities(abs(x - 4) + abs(3*abs(x) - 5) < 7) == \
        ((Integer(-2) < x) & (x < -1)) | ((Rational(1, 2) < x) & (x < 4))

    nr = Symbol('nr', extended_real=False)
    pytest.raises(TypeError, lambda: reduce_inequalities(abs(nr - 5) < 3))

    # sympy/sympy#10198
    assert reduce_inequalities(-1 + 1/abs(1/x - 1) < 0) == \
        ((Integer(0) < x) & (x < Rational(1, 2))) | (x < 0)
    assert reduce_inequalities(-1 + 1/abs(-1/x - 1) < 0) == \
        ((-Rational(1, 2) < x) & (x < 0)) | (Integer(0) < x)

    # sympy/sympy#10255
    assert reduce_inequalities(Piecewise((1, x < 1), (3, True)) > 1) == \
        (Integer(1) <= x)
    assert reduce_inequalities(Piecewise((x**2, x < 0), (2*x, True)) < 1) == \
        (Integer(-1) < x) & (x < Rational(1, 2))


def test_solve_inequalities():
    eqs = [x**2 - 2 < 0, x**2 - 1 > 0]
    assert reduce_inequalities(eqs) == (((-sqrt(2) < x) & (x < -1)) |
                                        ((Integer(1) < x) & (x < sqrt(2))))

    # issue sympy/sympy#6627, sympy/sympy#6547
    assert reduce_inequalities((x - 3)/(x - 2) < 0) == (Integer(2) < x) & (x < 3)
    assert reduce_inequalities(x/(x + 1) > 1, x) == (x < -1)

    assert reduce_inequalities(sin(x) > Rational(1, 2)) == (pi/6 < x) & (x < 5*pi/6)


def test_reduce_inequalities_general():
    assert reduce_inequalities(sqrt(2)*x >= 1) == (sqrt(2)/2 <= x)
    assert reduce_inequalities(PurePoly(x + 1, x) > 0) == (Integer(-1) < x)

    # issue sympy/sympy#10196
    assert reduce_inequalities(x**2 >= 0)
    assert reduce_inequalities(x**2 < 0) is false


def test_reduce_inequalities_boolean():
    assert reduce_inequalities(
        [Eq(x**2, 0), True]) == Eq(x, 0)
    assert reduce_inequalities([Eq(x**2, 0), False]) is false


def test_reduce_inequalities_multivariate():
    assert (reduce_inequalities([x**2 >= 1, y**2 >= 1]) ==
            ((Integer(1) <= x) | (x <= -1)) & ((Integer(1) <= y) | (y <= -1)))


def test_reduce_inequalities_errors():
    pytest.raises(NotImplementedError, lambda: reduce_inequalities(sin(x) + x >= 1))
    pytest.raises(NotImplementedError, lambda: reduce_inequalities(x**2*y + y >= 1))


def test_hacky_inequalities():
    y = Symbol('y', real=True)
    assert reduce_inequalities(x + y < 1, symbols=[x]) == (x < -y + 1)
    assert reduce_inequalities(x + y >= 1, symbols=[x]) == (-y + 1 <= x)


def test_sympyissue_10203():
    y = Symbol('y', extended_real=True)
    assert reduce_inequalities(Eq(0, x - y), symbols=[x]) == Eq(x, y)
    assert reduce_inequalities(Ne(0, x - y), symbols=[x]) == \
        (y < x) | (x < y)


def test_sympyissue_6343():
    eq = -3*x**2/2 - 45*x/4 + Rational(33, 2) > 0
    assert reduce_inequalities(eq) == \
        (x < -Rational(15, 4) + sqrt(401)/4) & (-sqrt(401)/4 - Rational(15, 4) < x)


def test_sympyissue_8235():
    assert reduce_inequalities(x**2 - 1 < 0) == \
        (Integer(-1) < x) & (x < Integer(1))
    assert reduce_inequalities(x**2 - 1 <= 0) == \
        (Integer(-1) <= x) & (x <= 1)
    assert reduce_inequalities(x**2 - 1 > 0) == \
        (x < -1) | (Integer(1) < x)
    assert reduce_inequalities(x**2 - 1 >= 0) == \
        (x <= Integer(-1)) | (Integer(1) <= x)

    eq = x**8 + x - 9  # we want RootOf solns here
    sol = reduce_inequalities(eq >= 0)
    tru = (RootOf(eq, 1) <= x) | (x <= RootOf(eq, 0))
    assert sol == tru

    # recast vanilla as real
    assert reduce_inequalities(sqrt((-x + 1)**2) < 1) == (Integer(0) < x) & (x < 2)


def test_sympyissue_5526():
    assert reduce_inequalities(Integer(0) <=
                               x + Integral(y**2, (y, 1, 3)) - 1, [x]) == \
        (-Integral(y**2, (y, 1, 3)) + 1 <= x)


def test_solve_univariate_inequality():
    assert isolve(x**2 >= 4, x, relational=False) == Union(Interval(-oo, -2),
                                                           Interval(2, oo))
    assert isolve(x**2 >= 4, x) == (Integer(2) <= x) | (x <= -2)
    assert isolve((x - 1)*(x - 2)*(x - 3) >= 0, x, relational=False) == \
        Union(Interval(1, 2), Interval(3, oo))
    assert isolve((x - 1)*(x - 2)*(x - 3) >= 0, x) == \
        ((Integer(1) <= x) & (x <= 2)) | (Integer(3) <= x)
    # issue sympy/sympy#2785:
    assert isolve(x**3 - 2*x - 1 > 0, x, relational=False) == \
        Union(Interval(-1, -sqrt(5)/2 + Rational(1, 2), True, True),
              Interval(Rational(1, 2) + sqrt(5)/2, oo, True))
    # issue sympy/sympy#2794:
    assert isolve(x**3 - x**2 + x - 1 > 0, x, relational=False) == \
        Interval(1, oo, True)

    # XXX should be limited in domain, e.g. between 0 and 2*pi
    assert isolve(sin(x) < Rational(1, 2), x) == (x < pi/6) | (5*pi/6 < x)
    assert isolve(sin(x) > Rational(1, 2), x) == (pi/6 < x) & (x < 5*pi/6)

    # numerical testing in valid() is needed
    assert isolve(x**7 - x - 2 > 0, x) == (RootOf(x**7 - x - 2, 0) < x)

    # handle numerator and denominator; although these would be handled as
    # rational inequalities, these test confirm that the right thing is done
    # when the domain is EX (e.g. when 2 is replaced with sqrt(2))
    assert isolve(1/(x - 2) > 0, x) == (Integer(2) < x)
    den = ((x - 1)*(x - 2)).expand()
    assert isolve((x - 1)/den <= 0, x) == (x < 1) | ((Integer(1) < x) & (x < 2))

    assert isolve(x > oo, x) is false

    # issue sympy/sympy#10268
    assert reduce_inequalities(log(x) < 300) == (-oo < x) & (x < E**300)


def test_slow_general_univariate():
    assert reduce_inequalities(sqrt(x) + 1/root(x, 3) > 1) == (Integer(0) < x)


def test_sympyissue_8545():
    eq = 1 - x - abs(1 - x)
    ans = Integer(1) < x
    assert reduce_piecewise_inequality(eq, '<', x) == ans
    eq = 1 - x - sqrt((1 - x)**2)
    assert reduce_inequalities(eq < 0) == ans


def test_sympyissue_8974():
    assert isolve(-oo < x, x) == (-oo < x)
    assert isolve(oo > x, x) == (x < oo)


def test_issue_453():
    x = Symbol('x', real=True)
    assert isolve(abs((x - 1)/(x - 5)) <= Rational(1, 3),
                  x) == (Integer(-1) <= x) & (x <= 2)
    assert solve(abs((x - 1)/(x - 5)) - Rational(1, 3), x) == [{x: -1}, {x: 2}]


def test_issue_836():
    assert reduce_inequalities(x + y > 1) == (y > -x + 1)

    pytest.raises(NotImplementedError,
                  lambda: reduce_inequalities([x**2*y**2 <= x**2*y,
                                               x**2*y**2 > x**2*y]))


def test_sympyissue_20861():
    assert reduce_inequalities([3/x < 0, x >= 2, x >= 7], x) is false


def test_sympyissue_20902():
    eq = y/((1 + y)**2)

    assert (reduce_inequalities(eq.subs({y: 3*x + 2}).diff(x) > 0) ==
            (Integer(-1) < x) & (x < Rational(-1, 3)))
    assert (reduce_inequalities(eq.subs({y: 3*x + 3}).diff(x) > 0) ==
            (Rational(-4, 3) < x) & (x < Rational(-2, 3)))
    assert (reduce_inequalities(eq.subs({y: 3*x + 4}).diff(x) > 0) ==
            (Rational(-5, 3) < x) & (x < Integer(-1)))
    assert (reduce_inequalities(eq.subs({y: 3*x + 2}).diff(x) > 0) ==
            (Integer(-1) < x) & (x < Rational(-1, 3)))
