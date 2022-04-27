"""Special sets tests."""

import itertools

import pytest

from diofant import (Basic, Contains, EmptySet, FiniteSet, I, ImageSet,
                     Integer, Intersection, Interval, Lambda, Range, Rational,
                     S, Set, Symbol, cos, exp, imageset, log, oo, pi, sin,
                     sqrt, symbols, tan)
from diofant.abc import m, n, t, x, y


__all__ = ()


def test_naturals():
    N = S.Naturals
    assert 5 in N
    assert -5 not in N
    assert 5.5 not in N
    assert (1, 2) not in N  # issue sympy/sympy#11732
    ni = iter(N)
    a, b, c, d = next(ni), next(ni), next(ni), next(ni)
    assert (a, b, c, d) == (1, 2, 3, 4)
    assert isinstance(a, Basic)

    assert N.intersection(Interval(-5, 5)) == Range(1, 6)
    assert N.intersection(Interval(-5, 5, True, True)) == Range(1, 5)
    assert N.intersection(Set(x)) == Intersection(N, Set(x),
                                                  evaluate=False)

    assert N.boundary == N

    assert N.inf == 1
    assert N.sup == oo


def test_naturals0():
    N = S.Naturals0
    assert 0 in N
    assert -1 not in N
    assert (1, 2) not in N  # issue sympy/sympy#11732
    assert next(iter(N)) == 0


def test_integers():
    Z = S.Integers
    assert 5 in Z
    assert -5 in Z
    assert 5.5 not in Z
    assert (1, 2) not in Z  # issue sympy/sympy#11732
    zi = iter(Z)
    a, b, c, d = next(zi), next(zi), next(zi), next(zi)
    assert (a, b, c, d) == (0, 1, -1, 2)
    assert isinstance(a, Basic)

    assert Z.intersection(Interval(-5, 5)) == Range(-5, 6)
    assert Z.intersection(Interval(-5, 5, True, True)) == Range(-4, 5)
    assert Z.intersection(Set(x)) == Intersection(Z, Set(x),
                                                  evaluate=False)

    assert Z.inf == -oo
    assert Z.sup == oo

    assert Z.boundary == Z

    assert imageset(Lambda((x, y), x*y), Z) == ImageSet(Lambda((x, y), x*y), Z)


def test_rationals():
    Q = S.Rationals
    assert 5 in Q
    assert -5 in Q
    assert Rational(1, 3) in Q
    assert pi not in Q
    assert Q.contains(x) == Contains(x, Q, evaluate=False)

    assert Q.inf == -oo
    assert Q.sup == oo

    assert Q.boundary == Q

    assert (list(itertools.islice(S.Rationals, 17)) ==
            [Integer(0), Integer(1), Rational(1, 2), Rational(1, 3),
             Integer(-1), Rational(-1, 2), Rational(-1, 3), Rational(1, 4),
             Rational(-1, 4), Integer(2), Rational(2, 3), Rational(1, 5),
             Rational(-1, 5), Rational(2, 5), Integer(-2),
             Rational(-2, 3), Rational(-2, 5)])


def test_ImageSet():
    squares = ImageSet(Lambda(x, x**2), S.Naturals)
    assert 4 in squares
    assert 5 not in squares
    assert FiniteSet(*range(10)).intersection(squares) == FiniteSet(1, 4, 9)

    assert 16 not in squares.intersection(Interval(0, 10))

    si = iter(squares)
    a, b, c, d = next(si), next(si), next(si), next(si)
    assert (a, b, c, d) == (1, 4, 9, 16)

    harmonics = ImageSet(Lambda(x, 1/x), S.Naturals)
    assert Rational(1, 5) in harmonics
    assert Rational(.25) in harmonics
    assert 0.25 not in harmonics
    assert Rational(.3) not in harmonics

    assert harmonics.is_iterable


def test_image_is_ImageSet():
    assert isinstance(imageset(x, sqrt(sin(x)), Range(5)), ImageSet)


@pytest.mark.xfail
def test_halfcircle():
    r, th = symbols('r, theta', extended_real=True)
    L = Lambda((r, th), (r*cos(th), r*sin(th)))
    halfcircle = ImageSet(L, Interval(0, 1)*Interval(0, pi))

    assert (1, 0) in halfcircle
    assert (0, -1) not in halfcircle
    assert (0, 0) in halfcircle

    assert not halfcircle.is_iterable


def test_ImageSet_iterator_not_injetive():
    L = Lambda(x, x - x % 2)  # produces 0, 2, 2, 4, 4, 6, 6, ...
    evens = ImageSet(L, S.Naturals)
    i = iter(evens)
    # No repeats here
    assert (next(i), next(i), next(i), next(i)) == (0, 2, 4, 6)


def test_Range():
    assert Range(5) == Range(0, 5) == Range(0, 5, 1)

    r = Range(10, 20, 2)
    assert 12 in r
    assert 8 not in r
    assert 11 not in r
    assert 30 not in r

    assert list(Range(0, 5)) == list(range(5))
    assert list(Range(5, 0, -1)) == list(range(1, 6))

    assert Range(0, 10, -1) == S.EmptySet

    assert Range(5, 15).sup == 14
    assert Range(5, 15).inf == 5
    assert Range(15, 5, -1).sup == 15
    assert Range(15, 5, -1).inf == 6
    assert Range(10, 67, 10).sup == 60
    assert Range(60, 7, -10).inf == 10

    assert len(Range(10, 38, 10)) == 3
    assert Range(0, 0, 5) == S.EmptySet

    assert Range(1, 1) == S.EmptySet
    pytest.raises(ValueError, lambda: Range(0, oo, oo))
    pytest.raises(ValueError, lambda: Range(0, pi, 1))

    assert 5 in Range(0, oo, 5)
    assert -5 in Range(-oo, 0, 5)

    assert Range(0, oo)
    assert Range(-oo, 0)
    assert Range(0, -oo, -1)
    assert Range(1, oo)
    assert Range(-oo, oo)

    assert Range(0, oo, 2)._last_element is oo
    assert Range(-oo, 1, 1)._last_element is Integer(0)

    it = iter(Range(-oo, 0, 2))
    assert (next(it), next(it)) == (-2, -4)

    assert Range(-1, 10, 1).intersection(S.Integers) == Range(-1, 10, 1)
    assert Range(-1, 10, 1).intersection(S.Naturals) == Range(1, 10, 1)
    assert (Range(-1, 10, 1).intersection(Set(x)) ==
            Intersection(Range(-1, 10, 1), Set(x), evaluate=False))

    assert Range(1, 10, 1)._ith_element(5) == 6  # the index starts from zero
    assert Range(1, 10, 1)._last_element == 9

    assert Range(1, 10, 1).boundary == Range(1, 10, 1)

    assert Range(range(10)) == Range(10)
    assert Range(range(1, 10)) == Range(1, 10)
    assert Range(range(1, 10, 2)) == Range(1, 10, 2)
    assert Range(range(1000000000000)) == Range(1000000000000)


def test_range_interval_intersection():
    # Intersection with intervals
    assert FiniteSet(*Range(0, 10, 1).intersection(Interval(2, 6))) == \
        FiniteSet(2, 3, 4, 5, 6)

    # Open Intervals are removed
    assert (FiniteSet(*Range(0, 10, 1).intersection(Interval(2, 6,
                                                             True, True))) ==
            FiniteSet(3, 4, 5))

    # Try this with large steps
    assert (FiniteSet(*Range(0, 100, 10).intersection(Interval(15, 55))) ==
            FiniteSet(20, 30, 40, 50))

    # Going backwards
    assert (FiniteSet(*Range(10, -9, -3).intersection(Interval(-5, 6))) ==
            FiniteSet(-5, -2, 1, 4))
    assert (FiniteSet(*Range(10, -9, -3).intersection(Interval(-5, 6, True))) ==
            FiniteSet(-2, 1, 4))


def test_fun():
    assert (FiniteSet(*ImageSet(Lambda(x, sin(pi*x/4)), Range(-10, 11))) ==
            FiniteSet(-1, -sqrt(2)/2, 0, sqrt(2)/2, 1))


def test_reals():
    assert 5 in S.Reals
    assert pi in S.Reals
    assert -sqrt(2) in S.Reals
    assert (2, 5) not in S.Reals
    assert sqrt(-1) not in S.Reals
    assert S.Reals == Interval.open(-oo, oo)
    assert hash(S.Reals) == hash(Interval.open(-oo, oo))
    assert S.Reals.is_open is True
    assert S.Reals.is_closed is True
    assert S.ExtendedReals == Interval(-oo, oo)
    assert hash(S.ExtendedReals) == hash(Interval(-oo, oo))
    assert S.ExtendedReals.is_closed is True


def test_intersections():
    assert S.Integers.intersection(S.Reals) == S.Integers
    assert 5 in S.Integers.intersection(S.Reals)
    assert 5 in S.Integers.intersection(S.Reals)
    assert -5 not in S.Naturals.intersection(S.Reals)
    assert 5.5 not in S.Integers.intersection(S.Reals)
    assert 5 in S.Integers.intersection(Interval(3, oo))
    assert -5 in S.Integers.intersection(Interval(-oo, 3))
    assert all(x.is_Integer
               for x in itertools.islice(S.Integers.intersection(Interval(3, oo)), 10))


def test_infinitely_indexed_set_1():
    assert imageset(Lambda(n, n), S.Integers) == imageset(Lambda(m, m), S.Integers)

    assert (imageset(Lambda(n, 2*n),
                     S.Integers).intersection(imageset(Lambda(m, 2*m + 1),
                                                       S.Integers)) ==
            EmptySet())

    assert (imageset(Lambda(n, 2*n),
                     S.Integers).intersection(imageset(Lambda(n, 2*n + 1),
                                                       S.Integers)) ==
            EmptySet())

    assert (imageset(Lambda(m, 2*m),
                     S.Integers).intersection(imageset(Lambda(n, 3*n),
                                                       S.Integers)) ==
            ImageSet(Lambda(t, 6*t), S.Integers))


def test_infinitely_indexed_set_2():
    a = Symbol('a', integer=True)
    assert imageset(Lambda(n, n), S.Integers) == imageset(Lambda(n, n + a), S.Integers)
    assert imageset(Lambda(n, n), S.Integers) == imageset(Lambda(n, -n + a), S.Integers)
    assert imageset(Lambda(n, -6*n), S.Integers) == ImageSet(Lambda(n, 6*n), S.Integers)
    assert imageset(Lambda(n, 2*n + pi), S.Integers) == ImageSet(Lambda(n, 2*n + pi), S.Integers)
    assert imageset(Lambda(n, pi*n + pi), S.Integers) == ImageSet(Lambda(n, pi*n + pi), S.Integers)
    assert imageset(Lambda(n, exp(n)), S.Integers) != imageset(Lambda(n, n), S.Integers)


def test_imageset_intersection_real():
    assert (imageset(Lambda(n, n + (n - 1)*(n + 1)*I),
                     S.Integers).intersection(S.Reals) ==
            FiniteSet(-1, 1))

    s = ImageSet(Lambda(n, -I*(I*(2*pi*n - pi/4) + log(abs(sqrt(-I))))), S.Integers)
    assert s.intersection(S.Reals) == imageset(Lambda(n, 2*n*pi - pi/4), S.Integers)


def test_infinitely_indexed_diophantine():
    assert (imageset(Lambda(m, 2*pi*m),
                     S.Integers).intersection(imageset(Lambda(n, 3*pi*n),
                                                       S.Integers)) ==
            ImageSet(Lambda(t, 6*pi*t), S.Integers))


@pytest.mark.xfail
def test_infinitely_indexed_set_3():
    assert imageset(Lambda(n, 2*n + 1), S.Integers) == imageset(Lambda(n, 2*n - 1), S.Integers)
    assert imageset(Lambda(n, 3*n + 2), S.Integers) == imageset(Lambda(n, 3*n - 1), S.Integers)


def test_ImageSet_simplification():
    assert imageset(Lambda(n, n), S.Integers) == S.Integers
    assert (imageset(Lambda(n, sin(n)),
                     imageset(Lambda(m, tan(m)), S.Integers)) ==
            imageset(Lambda(m, sin(tan(m))), S.Integers))


def test_sympyissue_10497():
    # __iter__ for infinite product set
    assert (list(itertools.islice(iter(S.Integers**2), 7)) ==
            [(0, 0), (0, 1), (1, 0), (1, 1), (0, -1), (1, -1), (-1, 0)])


def test_sympyissue_11732():
    interval12 = Interval(1, 2)
    finiteset1234 = FiniteSet(1, 2, 3, 4)

    assert (interval12 in S.Naturals) is False
    assert (interval12 in S.Naturals0) is False
    assert (interval12 in S.Integers) is False

    assert (finiteset1234 in S.Naturals) is False
    assert (finiteset1234 in S.Naturals0) is False
    assert (finiteset1234 in S.Integers) is False
