"""Tests for tools for solving inequalities and systems of inequalities. """

import pytest

from diofant import (And, E, Eq, FiniteSet, Float, Ge, Gt, Integer, Integral,
                     Interval, Le, Lt, Ne, Or, Piecewise, Poly, PurePoly,
                     Rational, S, Symbol, Union, false, log, oo, pi, root, sin,
                     sqrt)
from diofant.abc import x, y
from diofant.polys.rootoftools import RootOf
from diofant.solvers.inequalities import (reduce_inequalities,
                                          reduce_piecewise_inequality,
                                          reduce_rational_inequalities)
from diofant.solvers.inequalities import solve_poly_inequality as psolve
from diofant.solvers.inequalities import solve_univariate_inequality as isolve
from diofant.solvers.solvers import solve


__all__ = ()

inf = oo.evalf()


def test_solve_poly_inequality():
    assert psolve(Poly(0, x), '==') == [S.Reals]
    assert psolve(Poly(1, x), '==') == [S.EmptySet]
    assert psolve(PurePoly(x + 1, x), ">") == [Interval(-1, oo, True, True)]
    pytest.raises(ValueError, lambda: psolve(x, '=='))
    pytest.raises(ValueError, lambda: psolve(Poly(x, x), '??'))


def test_reduce_poly_inequalities_real_interval():
    assert reduce_rational_inequalities(
        [[Eq(x**2, 0)]], x, relational=False) == FiniteSet(0)
    assert reduce_rational_inequalities(
        [[Le(x**2, 0)]], x, relational=False) == FiniteSet(0)
    assert reduce_rational_inequalities(
        [[Lt(x**2, 0)]], x, relational=False) == S.EmptySet
    assert reduce_rational_inequalities(
        [[Ge(x**2, 0)]], x, relational=False) == \
        S.Reals if x.is_extended_real else Interval(-oo, oo)
    assert reduce_rational_inequalities(
        [[Gt(x**2, 0)]], x, relational=False) == \
        FiniteSet(0).complement(S.Reals)
    assert reduce_rational_inequalities(
        [[Ne(x**2, 0)]], x, relational=False) == \
        FiniteSet(0).complement(S.Reals)

    assert reduce_rational_inequalities(
        [[Eq(x**2, 1)]], x, relational=False) == FiniteSet(-1, 1)
    assert reduce_rational_inequalities(
        [[Le(x**2, 1)]], x, relational=False) == Interval(-1, 1)
    assert reduce_rational_inequalities(
        [[Lt(x**2, 1)]], x, relational=False) == Interval(-1, 1, True, True)
    assert reduce_rational_inequalities(
        [[Ge(x**2, 1)]], x, relational=False) == \
        Union(Interval(-oo, -1, True), Interval(1, oo, False, True))
    assert reduce_rational_inequalities(
        [[Gt(x**2, 1)]], x, relational=False) == \
        Interval(-1, 1).complement(S.Reals)
    assert reduce_rational_inequalities(
        [[Ne(x**2, 1)]], x, relational=False) == \
        FiniteSet(-1, 1).complement(S.Reals)
    assert reduce_rational_inequalities([[Eq(
        x**2, 1.0)]], x, relational=False) == FiniteSet(-1.0, 1.0).evalf()
    assert reduce_rational_inequalities(
        [[Le(x**2, 1.0)]], x, relational=False) == Interval(-1.0, 1.0)
    assert reduce_rational_inequalities([[Lt(
        x**2, 1.0)]], x, relational=False) == Interval(-1.0, 1.0, True, True)
    assert reduce_rational_inequalities(
        [[Ge(x**2, 1.0)]], x, relational=False) == \
        Union(Interval(-inf, -1.0, True), Interval(1.0, inf, False, True))
    assert reduce_rational_inequalities(
        [[Gt(x**2, 1.0)]], x, relational=False) == \
        Union(Interval(-inf, -1.0, True, True),
              Interval(1.0, inf, True, True))
    assert reduce_rational_inequalities([[Ne(
        x**2, 1.0)]], x, relational=False) == \
        FiniteSet(-1.0, 1.0).complement(S.Reals)

    s = sqrt(2)

    assert reduce_rational_inequalities([[Lt(
        x**2 - 1, 0), Gt(x**2 - 1, 0)]], x, relational=False) == S.EmptySet
    assert reduce_rational_inequalities([[Le(x**2 - 1, 0), Ge(
        x**2 - 1, 0)]], x, relational=False) == FiniteSet(-1, 1)
    assert reduce_rational_inequalities(
        [[Le(x**2 - 2, 0), Ge(x**2 - 1, 0)]], x, relational=False
    ) == Union(Interval(-s, -1, False, False), Interval(1, s, False, False))
    assert reduce_rational_inequalities(
        [[Le(x**2 - 2, 0), Gt(x**2 - 1, 0)]], x, relational=False
    ) == Union(Interval(-s, -1, False, True), Interval(1, s, True, False))
    assert reduce_rational_inequalities(
        [[Lt(x**2 - 2, 0), Ge(x**2 - 1, 0)]], x, relational=False
    ) == Union(Interval(-s, -1, True, False), Interval(1, s, False, True))
    assert reduce_rational_inequalities(
        [[Lt(x**2 - 2, 0), Gt(x**2 - 1, 0)]], x, relational=False
    ) == Union(Interval(-s, -1, True, True), Interval(1, s, True, True))
    assert reduce_rational_inequalities(
        [[Lt(x**2 - 2, 0), Ne(x**2 - 1, 0)]], x, relational=False
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
        [[Le(x**2, 0)]], x, relational=True) == Eq(x, 0)
    assert reduce_rational_inequalities(
        [[Lt(x**2, 0)]], x, relational=True) is false
    assert reduce_rational_inequalities(
        [[Ge(x**2, 0)]], x, relational=True) == And(Lt(-oo, x), Lt(x, oo))
    assert reduce_rational_inequalities(
        [[Gt(x**2, 0)]], x, relational=True) == \
        And(Or(And(Lt(-oo, x), Lt(x, 0)), And(Lt(0, x), Lt(x, oo))))
    assert reduce_rational_inequalities(
        [[Ne(x**2, 0)]], x, relational=True) == \
        And(Or(And(Lt(-oo, x), Lt(x, 0)), And(Lt(0, x), Lt(x, oo))))

    for one in (Integer(1), Float(1.0)):
        inf = one*oo
        assert reduce_rational_inequalities(
            [[Eq(x**2, one)]], x, relational=True) == \
            Or(Eq(x, -one), Eq(x, one))
        assert reduce_rational_inequalities(
            [[Le(x**2, one)]], x, relational=True) == \
            And(And(Le(-one, x), Le(x, one)))
        assert reduce_rational_inequalities(
            [[Lt(x**2, one)]], x, relational=True) == \
            And(And(Lt(-one, x), Lt(x, one)))
        assert reduce_rational_inequalities(
            [[Ge(x**2, one)]], x, relational=True) == \
            And(Or(And(Le(one, x), Lt(x, inf)), And(Le(x, -one), Lt(-inf, x))))
        assert reduce_rational_inequalities(
            [[Gt(x**2, one)]], x, relational=True) == \
            And(Or(And(Lt(-inf, x), Lt(x, -one)), And(Lt(one, x), Lt(x, inf))))
        assert reduce_rational_inequalities(
            [[Ne(x**2, one)]], x, relational=True) == \
            Or(And(Lt(-inf, x), Lt(x, -one)),
               And(Lt(-one, x), Lt(x, one)),
               And(Lt(one, x), Lt(x, inf)))


def test_reduce_rational_inequalities_real_relational():
    assert reduce_rational_inequalities([], x) is false
    assert reduce_rational_inequalities(
        [[(x**2 + 3*x + 2)/(x**2 - 16) >= 0]], x, relational=False) == \
        Union(Interval.open(-oo, -4), Interval(-2, -1), Interval.open(4, oo))

    assert reduce_rational_inequalities(
        [[((-2*x - 10)*(3 - x))/((x**2 + 5)*(x - 2)**2) < 0]], x,
        relational=False) == \
        Union(Interval.open(-5, 2), Interval.open(2, 3))

    assert reduce_rational_inequalities([[(x + 1)/(x - 5) <= 0]], x,
                                        relational=False) == \
        Interval.Ropen(-1, 5)

    assert reduce_rational_inequalities([[(x**2 + 4*x + 3)/(x - 1) > 0]], x,
                                        relational=False) == \
        Union(Interval.open(-3, -1), Interval.open(1, oo))

    assert reduce_rational_inequalities([[(x**2 - 16)/(x - 1)**2 < 0]], x,
                                        relational=False) == \
        Union(Interval.open(-4, 1), Interval.open(1, 4))

    assert reduce_rational_inequalities([[(3*x + 1)/(x + 4) >= 1]], x,
                                        relational=False) == \
        Union(Interval.open(-oo, -4), Interval.Ropen(Rational(3, 2), oo))

    assert reduce_rational_inequalities([[(x - 8)/x <= 3 - x]], x,
                                        relational=False) == \
        Union(Interval.Lopen(-oo, -2), Interval.Lopen(0, 4))


def test_reduce_piecewise_inequalities():
    e = abs(x - 5) < 3
    ans = And(Lt(2, x), Lt(x, 8))
    assert reduce_inequalities(e) == ans
    assert reduce_inequalities(e, x) == ans
    assert reduce_inequalities(abs(x - 5)) == Eq(x, 5)
    assert reduce_inequalities(
        abs(2*x + 3) >= 8) == Or(And(Le(Rational(5, 2), x), Lt(x, oo)),
                                 And(Le(x, -Rational(11, 2)), Lt(-oo, x)))
    assert reduce_inequalities(abs(x - 4) + abs(
        3*x - 5) < 7) == And(Lt(Rational(1, 2), x), Lt(x, 4))
    assert reduce_inequalities(abs(x - 4) + abs(3*abs(x) - 5) < 7) == \
        Or(And(Integer(-2) < x, x < -1), And(Rational(1, 2) < x, x < 4))

    nr = Symbol('nr', extended_real=False)
    pytest.raises(TypeError, lambda: reduce_inequalities(abs(nr - 5) < 3))

    # sympy/sympy#10198
    assert reduce_inequalities(-1 + 1/abs(1/x - 1) < 0) == \
        Or(And(Lt(0, x), x < Rational(1, 2)), And(-oo < x, x < 0))

    # sympy/sympy#10255
    assert reduce_inequalities(Piecewise((1, x < 1), (3, True)) > 1) == \
        And(Le(1, x), x < oo)
    assert reduce_inequalities(Piecewise((x**2, x < 0), (2*x, x >= 0)) < 1) == \
        And(Lt(-1, x), x < Rational(1, 2))


def test_reduce_inequalities_general():
    assert reduce_inequalities(Ge(sqrt(2)*x, 1)) == And(sqrt(2)/2 <= x, x < oo)
    assert reduce_inequalities(PurePoly(x + 1, x) > 0) == And(Integer(-1) < x, x < oo)


def test_reduce_inequalities_boolean():
    assert reduce_inequalities(
        [Eq(x**2, 0), True]) == Eq(x, 0)
    assert reduce_inequalities([Eq(x**2, 0), False]) is false


def test_reduce_inequalities_multivariate():
    assert reduce_inequalities([Ge(x**2, 1), Ge(y**2, 1)]) == And(
        Or(And(Le(1, x), Lt(x, oo)), And(Le(x, -1), Lt(-oo, x))),
        Or(And(Le(1, y), Lt(y, oo)), And(Le(y, -1), Lt(-oo, y))))


def test_reduce_inequalities_errors():
    pytest.raises(NotImplementedError, lambda: reduce_inequalities(Ge(sin(x) + x, 1)))
    pytest.raises(NotImplementedError, lambda: reduce_inequalities(Ge(x**2*y + y, 1)))


def test_hacky_inequalities():
    y = Symbol('y', real=True)
    assert reduce_inequalities(x + y < 1, symbols=[x]) == And(-oo < x, x < -y + 1)
    assert reduce_inequalities(x + y >= 1, symbols=[x]) == And(-y + 1 <= x, x < oo)


def test_sympyissue_10203():
    y = Symbol('y', extended_real=True)
    assert reduce_inequalities(Eq(0, x - y), symbols=[x]) == Eq(x, y)
    assert reduce_inequalities(Ne(0, x - y), symbols=[x]) == \
        Or(And(-oo < x, x < y), And(x < oo, y < x))


def test_sympyissue_6343():
    eq = -3*x**2/2 - 45*x/4 + Rational(33, 2) > 0
    assert reduce_inequalities(eq) == \
        And(x < -Rational(15, 4) + sqrt(401)/4, -sqrt(401)/4 - Rational(15, 4) < x)


def test_sympyissue_8235():
    assert reduce_inequalities(x**2 - 1 < 0) == \
        And(Integer(-1) < x, x < Integer(1))
    assert reduce_inequalities(x**2 - 1 <= 0) == \
        And(Integer(-1) <= x, x <= 1)
    assert reduce_inequalities(x**2 - 1 > 0) == \
        Or(And(-oo < x, x < -1), And(x < oo, Integer(1) < x))
    assert reduce_inequalities(x**2 - 1 >= 0) == \
        Or(And(-oo < x, x <= Integer(-1)), And(Integer(1) <= x, x < oo))

    eq = x**8 + x - 9  # we want RootOf solns here
    sol = reduce_inequalities(eq >= 0)
    tru = Or(And(RootOf(eq, 1) <= x, x < oo), And(-oo < x, x <= RootOf(eq, 0)))
    assert sol == tru

    # recast vanilla as real
    assert reduce_inequalities(sqrt((-x + 1)**2) < 1) == And(Integer(0) < x, x < 2)


def test_sympyissue_5526():
    assert reduce_inequalities(Integer(0) <=
                               x + Integral(y**2, (y, 1, 3)) - 1, [x]) == \
        And(-Integral(y**2, (y, 1, 3)) + 1 <= x, x < oo)


def test_solve_univariate_inequality():
    assert isolve(x**2 >= 4, x, relational=False) == Union(Interval(-oo, -2, True),
                                                           Interval(2, oo, False, True))
    assert isolve(x**2 >= 4, x) == Or(And(Le(2, x), Lt(x, oo)), And(Le(x, -2),
                                                                    Lt(-oo, x)))
    assert isolve((x - 1)*(x - 2)*(x - 3) >= 0, x, relational=False) == \
        Union(Interval(1, 2), Interval(3, oo, False, True))
    assert isolve((x - 1)*(x - 2)*(x - 3) >= 0, x) == \
        Or(And(Le(1, x), Le(x, 2)), And(Le(3, x), Lt(x, oo)))
    # issue sympy/sympy#2785:
    assert isolve(x**3 - 2*x - 1 > 0, x, relational=False) == \
        Union(Interval(-1, -sqrt(5)/2 + Rational(1, 2), True, True),
              Interval(Rational(1, 2) + sqrt(5)/2, oo, True, True))
    # issue sympy/sympy#2794:
    assert isolve(x**3 - x**2 + x - 1 > 0, x, relational=False) == \
        Interval(1, oo, True, True)

    # XXX should be limited in domain, e.g. between 0 and 2*pi
    assert isolve(sin(x) < Rational(1, 2), x) == \
        Or(And(-oo < x, x < pi/6), And(5*pi/6 < x, x < oo))
    assert isolve(sin(x) > Rational(1, 2), x) == And(pi/6 < x, x < 5*pi/6)

    # numerical testing in valid() is needed
    assert isolve(x**7 - x - 2 > 0, x) == \
        And(RootOf(x**7 - x - 2, 0) < x, x < oo)

    # handle numerator and denominator; although these would be handled as
    # rational inequalities, these test confirm that the right thing is done
    # when the domain is EX (e.g. when 2 is replaced with sqrt(2))
    assert isolve(1/(x - 2) > 0, x) == And(Integer(2) < x, x < oo)
    den = ((x - 1)*(x - 2)).expand()
    assert isolve((x - 1)/den <= 0, x) == \
        Or(And(-oo < x, x < 1), And(Integer(1) < x, x < 2))

    assert isolve(x > oo, x) is false


def test_slow_general_univariate():
    r = RootOf(x**5 - x**2 + 1, 0)
    assert reduce_inequalities(sqrt(x) + 1/root(x, 3) > 1) == \
        Or(And(Integer(0) < x, x < r**6), And(r**6 < x, x < oo))


def test_sympyissue_8545():
    eq = 1 - x - abs(1 - x)
    ans = And(Lt(1, x), Lt(x, oo))
    assert reduce_piecewise_inequality(eq, '<', x) == ans
    eq = 1 - x - sqrt((1 - x)**2)
    assert reduce_inequalities(eq < 0) == ans


def test_sympyissue_8974():
    assert isolve(-oo < x, x) == And(-oo < x, x < oo)
    assert isolve(oo > x, x) == And(-oo < x, x < oo)


def test_sympyissue_10196():
    assert reduce_inequalities(x**2 >= 0)
    assert reduce_inequalities(x**2 < 0) is false


def test_sympyissue_10268():
    assert reduce_inequalities(log(x) < 300) == And(-oo < x, x < E**300)


def test_diofantissue_453():
    x = Symbol('x', real=True)
    assert isolve(abs((x - 1)/(x - 5)) <= Rational(1, 3),
                  x) == And(Integer(-1) <= x, x <= 2)
    assert solve(abs((x - 1)/(x - 5)) - Rational(1, 3), x) == [{x: -1}, {x: 2}]
