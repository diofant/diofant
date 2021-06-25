"""Generic set theory tests."""

import pytest
from mpmath import mpi

from diofant import (And, Complement, Contains, E, EmptySet, Eq, FiniteSet,
                     Float, I, ImageSet, Integer, Intersection, Interval,
                     Lambda, Le, LessThan, Lt, Max, Min, Or, Piecewise, Pow,
                     ProductSet, Rational, RootOf, S, Set, Sum, Symbol,
                     SymmetricDifference, Union, cos, false, imageset, log,
                     nan, oo, pi, sin, sqrt, sympify, true, zoo)
from diofant.abc import a, b, x, y, z


__all__ = ()


def test_interval_arguments():
    assert Interval(0, oo).right_open is false
    assert Interval(-oo, 0).left_open is false
    assert Interval(oo, -oo) == S.EmptySet

    assert isinstance(Interval(1, 1), FiniteSet)
    e = Sum(x, (x, 1, 3))
    assert isinstance(Interval(e, e), FiniteSet)

    assert Interval(1, 0) == S.EmptySet
    assert Interval(1, 1).measure == 0

    assert Interval(1, 1, False, True) == S.EmptySet
    assert Interval(1, 1, True, False) == S.EmptySet
    assert Interval(1, 1, True, True) == S.EmptySet

    assert isinstance(Interval(0, Symbol('a')), Interval)
    assert Interval(Symbol('a', extended_real=True, positive=True), 0) == S.EmptySet
    pytest.raises(ValueError, lambda: Interval(0, I))
    pytest.raises(ValueError, lambda: Interval(0, Symbol('z', extended_real=False)))

    pytest.raises(NotImplementedError, lambda: Interval(0, 1, And(x, y)))
    pytest.raises(NotImplementedError, lambda: Interval(0, 1, False, And(x, y)))
    pytest.raises(NotImplementedError, lambda: Interval(0, 1, z, And(x, y)))


def test_interval_evalf():
    assert (Interval(E, pi).evalf() ==
            Interval(Float('2.7182818284590451', dps=15),
                     Float('3.1415926535897931', dps=15), false, false))


def test_interval_symbolic_end_points():
    a = Symbol('a', extended_real=True)

    assert Union(Interval(0, a), Interval(0, 3)).sup == Max(a, 3)
    assert Union(Interval(a, 0), Interval(-3, 0)).inf == Min(-3, a)

    assert Interval(0, a).contains(1) == LessThan(1, a)


def test_union():
    assert Union(Interval(1, 2), Interval(2, 3)) == Interval(1, 3)
    assert Union(Interval(1, 2), Interval(2, 3, True)) == Interval(1, 3)
    assert Union(Interval(1, 3), Interval(2, 4)) == Interval(1, 4)
    assert Union(Interval(1, 2), Interval(1, 3)) == Interval(1, 3)
    assert Union(Interval(1, 3), Interval(1, 2)) == Interval(1, 3)
    assert Union(Interval(1, 3, False, True), Interval(1, 2)) == \
        Interval(1, 3, False, True)
    assert Union(Interval(1, 3), Interval(1, 2, False, True)) == Interval(1, 3)
    assert Union(Interval(1, 2, True), Interval(1, 3)) == Interval(1, 3)
    assert Union(Interval(1, 2, True), Interval(1, 3, True)) == \
        Interval(1, 3, True)
    assert Union(Interval(1, 2, True), Interval(1, 3, True, True)) == \
        Interval(1, 3, True, True)
    assert Union(Interval(1, 2, True, True), Interval(1, 3, True)) == \
        Interval(1, 3, True)
    assert Union(Interval(1, 3), Interval(2, 3)) == Interval(1, 3)
    assert Union(Interval(1, 3, False, True), Interval(2, 3)) == \
        Interval(1, 3)
    assert Union(Interval(1, 2, False, True), Interval(2, 3, True)) != \
        Interval(1, 3)
    assert Union(Interval(1, 2), S.EmptySet) == Interval(1, 2)
    assert Union(S.EmptySet) == S.EmptySet

    assert Union(Interval(0, 1), [FiniteSet(1.0/n) for n in range(1, 10)]) == \
        Interval(0, 1)

    assert Interval(1, 2).union(Interval(2, 3)) == \
        Interval(1, 2) + Interval(2, 3)

    assert Interval(1, 2).union(Interval(2, 3)) == Interval(1, 3)

    assert Union(Set()) == Set()

    assert FiniteSet(1) + FiniteSet(2) + FiniteSet(3) == FiniteSet(1, 2, 3)
    assert FiniteSet('ham') + FiniteSet('eggs') == FiniteSet('ham', 'eggs')
    assert FiniteSet(1, 2, 3) + S.EmptySet == FiniteSet(1, 2, 3)

    assert FiniteSet(1, 2, 3) & FiniteSet(2, 3, 4) == FiniteSet(2, 3)
    assert FiniteSet(1, 2, 3) | FiniteSet(2, 3, 4) == FiniteSet(1, 2, 3, 4)

    assert S.EmptySet | FiniteSet(x, FiniteSet(y, z)) == \
        FiniteSet(x, FiniteSet(y, z))

    # Test that Intervals and FiniteSets play nicely
    assert Interval(1, 3) + FiniteSet(2) == Interval(1, 3)
    assert Interval(1, 3, True, True) + FiniteSet(3) == \
        Interval(1, 3, True, False)
    X = Interval(1, 3) + FiniteSet(5)
    Y = Interval(1, 2) + FiniteSet(3)
    XandY = X.intersection(Y)
    assert 2 in X and 3 in X and 3 in XandY
    assert XandY.is_subset(X) and XandY.is_subset(Y)

    pytest.raises(TypeError, lambda: Union(1, 2, 3))

    assert X.is_iterable is False
    Z = Union(FiniteSet(1, 2)*FiniteSet(3, 4), FiniteSet(1, 2, 3, 4))
    assert Z.is_iterable

    # issue sympy/sympy#7843
    assert Union(S.EmptySet, FiniteSet(-sqrt(-I), sqrt(-I))) == FiniteSet(-sqrt(-I), sqrt(-I))

    assert Union(ProductSet(FiniteSet(1), FiniteSet(2)), Interval(0, 1)).is_Union
    assert Union(ProductSet(FiniteSet(1), FiniteSet(2)),
                 ProductSet(FiniteSet(1), FiniteSet(2), FiniteSet(3))).is_Union

    assert list(Union(FiniteSet(1, 2), FiniteSet(3, 4), evaluate=False)) == [1, 3, 2, 4]
    pytest.raises(TypeError, lambda: iter(Union(FiniteSet(1, 2), Interval(0, 1))))

    assert (Union(FiniteSet(E), FiniteSet(pi), evaluate=False).evalf() ==
            FiniteSet(Float('2.7182818284590451', dps=15),
                      Float('3.1415926535897931', dps=15)))


def test_difference():
    assert Interval(1, 3) - Interval(1, 2) == Interval(2, 3, True)
    assert Interval(1, 3) - Interval(2, 3) == Interval(1, 2, False, True)
    assert Interval(1, 3, True) - Interval(2, 3) == Interval(1, 2, True, True)
    assert Interval(1, 3, True) - Interval(2, 3, True) == \
        Interval(1, 2, True, False)
    assert Interval(0, 2) - FiniteSet(1) == \
        Union(Interval(0, 1, False, True), Interval(1, 2, True, False))

    assert FiniteSet(1, 2, 3) - FiniteSet(2) == FiniteSet(1, 3)
    assert (FiniteSet('ham', 'eggs') - FiniteSet('eggs') ==
            Complement(FiniteSet('ham'), FiniteSet('eggs')))
    assert FiniteSet(1, 2, 3, 4) - Interval(2, 10, True, False) == \
        FiniteSet(1, 2)
    assert FiniteSet(1, 2, 3, 4) - S.EmptySet == FiniteSet(1, 2, 3, 4)
    assert Union(Interval(0, 2), FiniteSet(2, 3, 4)) - Interval(1, 3) == \
        Union(Interval(0, 1, False, True), FiniteSet(4))

    assert -1 in S.Reals - S.Naturals


def test_Complement():
    assert Complement(Interval(1, 3), Interval(1, 2)) == Interval(2, 3, True)
    assert Complement(FiniteSet(1, 3, 4), FiniteSet(3, 4)) == FiniteSet(1)
    assert Complement(Union(Interval(0, 2),
                            FiniteSet(2, 3, 4)), Interval(1, 3)) == \
        Union(Interval(0, 1, False, True), FiniteSet(4))

    assert 3 not in Complement(Interval(0, 5), Interval(1, 4), evaluate=False)
    assert -1 in Complement(S.Reals, S.Naturals, evaluate=False)
    assert 1 not in Complement(S.Reals, S.Naturals, evaluate=False)

    assert Complement(S.Integers, S.UniversalSet) == EmptySet()
    assert S.UniversalSet.complement(S.Integers) == EmptySet()

    assert (0 not in S.Reals.intersection(S.Integers - FiniteSet(0)))

    assert S.EmptySet - S.Integers == S.EmptySet

    assert (S.Integers - FiniteSet(0)) - FiniteSet(1) == S.Integers - FiniteSet(0, 1)

    assert (S.Reals - Union(S.Naturals, FiniteSet(pi)) ==
            Intersection(S.Reals - S.Naturals, S.Reals - FiniteSet(pi)))


def test_complement():
    assert Interval(0, 1).complement(S.Reals) == \
        Union(Interval(-oo, 0, True, True), Interval(1, oo, True, True))
    assert Interval(0, 1, True, False).complement(S.Reals) == \
        Union(Interval(-oo, 0, True, False), Interval(1, oo, True, True))
    assert Interval(0, 1, False, True).complement(S.Reals) == \
        Union(Interval(-oo, 0, True, True), Interval(1, oo, False, True))
    assert Interval(0, 1, True, True).complement(S.Reals) == \
        Union(Interval(-oo, 0, True, False), Interval(1, oo, False, True))

    assert S.UniversalSet.complement(S.EmptySet) == S.EmptySet
    assert S.UniversalSet.complement(S.Reals) == S.EmptySet
    assert S.UniversalSet.complement(S.UniversalSet) == S.EmptySet

    assert S.EmptySet.complement(S.Reals) == S.Reals

    assert Union(Interval(0, 1), Interval(2, 3)).complement(S.Reals) == \
        Union(Interval(-oo, 0, True, True), Interval(1, 2, True, True),
              Interval(3, oo, True, True))

    assert FiniteSet(0).complement(S.Reals) ==  \
        Union(Interval(-oo, 0, True, True), Interval(0, oo, True, True))

    assert (FiniteSet(5) + Interval(-oo,
                                    0)).complement(S.Reals) == \
        Interval(0, 5, True, True) + Interval(5, oo, True, True)

    assert FiniteSet(1, 2, 3).complement(S.Reals) == \
        Interval(-oo, 1, True, True) + \
        Interval(1, 2, True, True) + Interval(2, 3, True, True) +\
        Interval(3, oo, True, True)

    assert FiniteSet(x).complement(S.Reals) == Complement(S.Reals, FiniteSet(x))

    assert FiniteSet(0, x).complement(S.Reals) == Complement(Interval(-oo, 0, True, True) +
                                                             Interval(0, oo, True, True),
                                                             FiniteSet(x), evaluate=False)

    square = Interval(0, 1) * Interval(0, 1)
    notsquare = square.complement(S.Reals*S.Reals)

    assert all(pt in square for pt in [(0, 0), (.5, .5), (1, 0), (1, 1)])
    assert not any(
        pt in notsquare for pt in [(0, 0), (.5, .5), (1, 0), (1, 1)])
    assert not any(pt in square for pt in [(-1, 0), (1.5, .5), (10, 10)])
    assert all(pt in notsquare for pt in [(-1, 0), (1.5, .5), (10, 10)])


def test_intersection():
    pytest.raises(TypeError, lambda: Intersection(1))

    assert Intersection() == S.UniversalSet

    assert Interval(0, 2).intersection(Interval(1, 2)) == Interval(1, 2)
    assert Interval(0, 2).intersection(Interval(1, 2, True)) == \
        Interval(1, 2, True)
    assert Interval(0, 2, True).intersection(Interval(1, 2)) == \
        Interval(1, 2, False, False)
    assert Interval(0, 2, True, True).intersection(Interval(1, 2)) == \
        Interval(1, 2, False, True)
    assert Interval(0, 2).intersection(Union(Interval(0, 1), Interval(2, 3))) == \
        Union(Interval(0, 1), Interval(2, 2))

    x = Symbol('x')

    assert Intersection(Interval(0, 1), FiniteSet(x)).is_iterable
    assert not Intersection(Interval(0, 1), Interval(0, x)).is_iterable

    assert FiniteSet(1, 2, x).intersection(FiniteSet(x)) == FiniteSet(x)
    assert FiniteSet('ham', 'eggs').intersection(FiniteSet('ham')) == \
        FiniteSet('ham')
    assert FiniteSet(1, 2, 3, 4, 5).intersection(S.EmptySet) == S.EmptySet

    assert Interval(0, 5).intersection(FiniteSet(1, 3)) == FiniteSet(1, 3)
    assert Interval(0, 1, True, True).intersection(FiniteSet(1)) == S.EmptySet

    assert Union(Interval(0, 1), Interval(2, 3)).intersection(Interval(1, 2)) == \
        Union(Interval(1, 1), Interval(2, 2))
    assert Union(Interval(0, 1), Interval(2, 3)).intersection(Interval(0, 2)) == \
        Union(Interval(0, 1), Interval(2, 2))
    assert Union(Interval(0, 1), Interval(2, 3)).intersection(Interval(1, 2, True, True)) == \
        S.EmptySet
    assert Union(Interval(0, 1), Interval(2, 3)).intersection(S.EmptySet) == \
        S.EmptySet
    assert Union(Interval(0, 5), FiniteSet('ham')).intersection(FiniteSet(2, 3, 4, 5, 6)) == \
        Union(FiniteSet(2, 3, 4, 5), Intersection(FiniteSet(6), Union(Interval(0, 5), FiniteSet('ham'))))

    # issue sympy/sympy#8217
    assert Intersection(FiniteSet(x), FiniteSet(y)) == \
        Intersection(FiniteSet(x), FiniteSet(y), evaluate=False)
    assert FiniteSet(x).intersection(S.Reals) == \
        Intersection(S.Reals, FiniteSet(x), evaluate=False)

    # tests for the intersection alias
    assert Interval(0, 5).intersection(FiniteSet(1, 3)) == FiniteSet(1, 3)
    assert Interval(0, 1, True, True).intersection(FiniteSet(1)) == S.EmptySet

    assert Union(Interval(0, 1), Interval(2, 3)).intersection(Interval(1, 2)) == \
        Union(Interval(1, 1), Interval(2, 2))

    assert ProductSet(FiniteSet(1),
                      FiniteSet(2)).intersection(Interval(1, 2)).is_Intersection

    # iterable
    i = Intersection(FiniteSet(1, 2, 3), Interval(2, 5), evaluate=False)
    assert i.is_iterable
    assert set(i) == {Integer(2), Integer(3)}

    # challenging intervals
    x = Symbol('x', real=True)
    i = Intersection(Interval(0, 3), Interval(x, 6))
    assert (5 in i) is False
    pytest.raises(TypeError, lambda: 2 in i)

    # Singleton special cases
    assert Intersection(Interval(0, 1), S.EmptySet) == S.EmptySet
    assert Intersection(S.Reals, Interval(-oo, x, True)) == Interval(-oo, x, True)
    assert Intersection(S.Reals, Interval(-oo, x)) == Interval.Lopen(-oo, x)
    assert Intersection(S.Reals, Interval(x, oo)) == Interval.Ropen(x, oo)
    assert Intersection(S.Reals, Interval(-oo, 1)) == Interval.Lopen(-oo, 1)
    assert Intersection(S.Reals, Interval(1, oo)) == Interval.Ropen(1, oo)
    assert Intersection(S.ExtendedReals, Interval(1, 2)) == Interval(1, 2)

    # Products
    line = Interval(0, 5)
    i = Intersection(line**2, line**3, evaluate=False)
    assert (2, 2) not in i
    assert (2, 2, 2) not in i
    pytest.raises(ValueError, lambda: list(i))

    assert (Intersection(Intersection(S.Integers, S.Naturals, evaluate=False),
                         S.Reals, evaluate=False) ==
            Intersection(S.Integers, S.Naturals, S.Reals, evaluate=False))

    assert (imageset(Lambda(x, x**2),
                     Intersection(FiniteSet(1, 2), FiniteSet(2, 3),
                                  evaluate=False)) == FiniteSet(4))


def test_is_disjoint():
    assert Interval(0, 2).is_disjoint(Interval(1, 2)) is False
    assert Interval(0, 2).is_disjoint(Interval(3, 4)) is True
    assert Interval(0, 2).isdisjoint(Interval(3, 4)) is True


def test_ProductSet_of_single_arg_is_arg():
    assert ProductSet(Interval(0, 1)) == Interval(0, 1)


def test_pow_args():
    pytest.raises(ValueError, lambda: Interval(0, 1)**-1)
    pytest.raises(TypeError, lambda: ProductSet(2))


def test_ProductSet_iter():
    pytest.raises(TypeError, lambda: iter(ProductSet(Set(1), Set(2))))


def test_interval_subs():
    a = Symbol('a', extended_real=True)

    assert Interval(0, a).subs({a: 2}) == Interval(0, 2)
    assert Interval(a, 0).subs({a: 2}) == S.EmptySet


def test_interval_to_mpi():
    assert Interval(0, 1).to_mpi() == mpi(0, 1)
    assert Interval(0, 1, True, False).to_mpi() == mpi(0, 1)
    assert isinstance(Interval(0, 1).to_mpi(), type(mpi(0, 1)))


def test_measure():
    a = Symbol('a', extended_real=True)

    assert Interval(1, 3).measure == 2
    assert Interval(0, a).measure == a
    assert Interval(1, a).measure == a - 1

    assert Union(Interval(1, 2), Interval(3, 4)).measure == 2
    assert Union(Interval(1, 2), Interval(3, 4), FiniteSet(5, 6, 7)).measure \
        == 2

    assert FiniteSet(1, 2, oo, a, -oo, -5).measure == 0

    assert S.EmptySet.measure == 0

    square = Interval(0, 10) * Interval(0, 10)
    offsetsquare = Interval(5, 15) * Interval(5, 15)
    band = Interval(-oo, oo) * Interval(2, 4)

    assert square.measure == offsetsquare.measure == 100
    assert (square + offsetsquare).measure == 175  # there is some overlap
    assert (square - offsetsquare).measure == 75
    assert (square * FiniteSet(1, 2, 3)).measure == 0
    assert (square.intersection(band)).measure == 20
    assert (square + band).measure == oo
    assert (band * FiniteSet(1, 2, 3)).measure == nan


def test_is_subset():
    assert Interval(0, 1).is_subset(Interval(0, 2)) is True
    assert Interval(0, 3).is_subset(Interval(0, 2)) is False

    assert FiniteSet(1, 2).is_subset(FiniteSet(1, 2, 3, 4))
    assert FiniteSet(4, 5).is_subset(FiniteSet(1, 2, 3, 4)) is False
    assert FiniteSet(1).is_subset(Interval(0, 2))
    assert FiniteSet(1, 2).is_subset(Interval(0, 2, True, True)) is False
    assert (Interval(1, 2) + FiniteSet(3)).is_subset(
        (Interval(0, 2, False, True) + FiniteSet(2, 3)))

    assert Interval(3, 4).is_subset(Union(Interval(0, 1), Interval(2, 5))) is True
    assert Interval(3, 6).is_subset(Union(Interval(0, 1), Interval(2, 5))) is False

    assert FiniteSet(1, 2, 3, 4).is_subset(Interval(0, 5)) is True
    assert S.EmptySet.is_subset(FiniteSet(1, 2, 3)) is True

    assert Interval(0, 1).is_subset(S.EmptySet) is False
    assert S.EmptySet.is_subset(S.EmptySet) is True

    pytest.raises(ValueError, lambda: S.EmptySet.is_subset(1))

    # tests for the issubset alias
    assert FiniteSet(1, 2, 3, 4).issubset(Interval(0, 5)) is True
    assert S.EmptySet.issubset(FiniteSet(1, 2, 3)) is True


def test_is_proper_subset():
    assert Interval(0, 1).is_proper_subset(Interval(0, 2)) is True
    assert Interval(0, 3).is_proper_subset(Interval(0, 2)) is False
    assert S.EmptySet.is_proper_subset(FiniteSet(1, 2, 3)) is True

    pytest.raises(ValueError, lambda: Interval(0, 1).is_proper_subset(0))


def test_is_superset():
    assert Interval(0, 1).is_superset(Interval(0, 2)) is False
    assert Interval(0, 3).is_superset(Interval(0, 2))

    assert FiniteSet(1, 2).is_superset(FiniteSet(1, 2, 3, 4)) is False
    assert FiniteSet(4, 5).is_superset(FiniteSet(1, 2, 3, 4)) is False
    assert FiniteSet(1).is_superset(Interval(0, 2)) is False
    assert FiniteSet(1, 2).is_superset(Interval(0, 2, True, True)) is False
    assert (Interval(1, 2) + FiniteSet(3)).is_superset(
        (Interval(0, 2, False, True) + FiniteSet(2, 3))) is False

    assert Interval(3, 4).is_superset(Union(Interval(0, 1), Interval(2, 5))) is False

    assert FiniteSet(1, 2, 3, 4).is_superset(Interval(0, 5)) is False
    assert S.EmptySet.is_superset(FiniteSet(1, 2, 3)) is False

    assert Interval(0, 1).is_superset(S.EmptySet) is True
    assert S.EmptySet.is_superset(S.EmptySet) is True

    pytest.raises(ValueError, lambda: S.EmptySet.is_superset(1))

    # tests for the issuperset alias
    assert Interval(0, 1).issuperset(S.EmptySet) is True
    assert S.EmptySet.issuperset(S.EmptySet) is True


def test_is_proper_superset():
    assert Interval(0, 1).is_proper_superset(Interval(0, 2)) is False
    assert Interval(0, 3).is_proper_superset(Interval(0, 2)) is True
    assert FiniteSet(1, 2, 3).is_proper_superset(S.EmptySet) is True

    pytest.raises(ValueError, lambda: Interval(0, 1).is_proper_superset(0))


def test_contains():
    assert Interval(0, 2).contains(1) is true
    assert Interval(0, 2).contains(3) is false
    assert Interval(0, 2, True, False).contains(0) is false
    assert Interval(0, 2, True, False).contains(2) is true
    assert Interval(0, 2, False, True).contains(0) is true
    assert Interval(0, 2, False, True).contains(2) is false
    assert Interval(0, 2, True, True).contains(0) is false
    assert Interval(0, 2, True, True).contains(2) is false

    assert Interval(0, 2) not in Interval(0, 2)

    # issue sympy/sympy#10326
    assert S.Reals.contains(oo) is false
    assert S.Reals.contains(-oo) is false
    assert Interval(-oo, oo, True).contains(oo) is true
    assert Interval(-oo, oo).contains(-oo) is true
    bad = [EmptySet(), FiniteSet(1), Interval(1, 2), zoo,
           I, oo, nan, -oo]
    assert all(i not in Interval(0, 5) for i in bad)

    assert FiniteSet(1, 2, 3).contains(2) is true
    assert FiniteSet(1, 2, x).contains(x) is true

    # issue sympy/sympy#8197
    assert isinstance(FiniteSet(b).contains(-a), Contains)
    assert isinstance(FiniteSet(b).contains(a), Contains)
    assert isinstance(FiniteSet(a).contains(1), Contains)
    pytest.raises(TypeError, lambda: 1 in FiniteSet(a))

    # issue sympy/sympy#8209
    rad1 = Integer(4)
    rad2 = Pow(2, 2, evaluate=False)
    s1 = FiniteSet(rad1)
    s2 = FiniteSet(rad2)
    assert s1 - s2 == S.EmptySet

    items = [1, 2, oo, Symbol('ham'), -1.1]
    fset = FiniteSet(*items)
    assert all(item in fset for item in items)
    assert all(fset.contains(item) is true for item in items)

    assert Union(Interval(0, 1), Interval(2, 5)).contains(3) is true
    assert Union(Interval(0, 1), Interval(2, 5)).contains(6) is false
    assert Union(Interval(0, 1), FiniteSet(2, 5)).contains(3) is false

    assert S.EmptySet.contains(1) is false
    assert FiniteSet(RootOf(x**3 + x - 1, 0)).contains(oo) is false

    assert RootOf(x**5 + x**3 + 1, 0) in S.Reals
    assert not RootOf(x**5 + x**3 + 1, 1) in S.Reals


def test_interval_symbolic():
    e = Interval(0, 1)
    assert e.contains(x) == And(0 <= x, x <= 1)
    pytest.raises(TypeError, lambda: x in e)
    e = Interval(0, 1, True, True)
    assert e.contains(x) == And(0 < x, x < 1)


def test_union_contains():
    i1 = Interval(0, 1)
    i2 = Interval(2, 3)
    i3 = Union(i1, i2)
    pytest.raises(TypeError, lambda: x in i3)
    e = i3.contains(x)
    assert e == Or(And(0 <= x, x <= 1), And(2 <= x, x <= 3))
    assert e.subs({x: -0.5}) is false
    assert e.subs({x: 0.5}) is true
    assert e.subs({x: 1.5}) is false
    assert e.subs({x: 2.5}) is true
    assert e.subs({x: 3.5}) is false

    U = Interval(0, 2, True, True) + Interval(10, oo) + FiniteSet(-1, 2, 5, 6)
    assert all(el not in U for el in [0, 4, -oo])
    assert all(el in U for el in [2, 5, 10])


def test_is_number():
    assert Interval(0, 1).is_number is False
    assert Set().is_number is False


def test_Interval_is_left_unbounded():
    assert Interval(3, 4).is_left_unbounded is False
    assert Interval(-oo, 3).is_left_unbounded is True
    assert Interval(Float('-inf'), 3).is_left_unbounded is True


def test_Interval_is_right_unbounded():
    assert Interval(3, 4).is_right_unbounded is False
    assert Interval(3, oo).is_right_unbounded is True
    assert Interval(3, Float('+inf')).is_right_unbounded is True


def test_Interval_as_relational():
    x = Symbol('x')

    assert Interval(-1, 2, False, False).as_relational(x) == \
        And(Le(-1, x), Le(x, 2))
    assert Interval(-1, 2, True, False).as_relational(x) == \
        And(Lt(-1, x), Le(x, 2))
    assert Interval(-1, 2, False, True).as_relational(x) == \
        And(Le(-1, x), Lt(x, 2))
    assert Interval(-1, 2, True, True).as_relational(x) == \
        And(Lt(-1, x), Lt(x, 2))

    assert Interval(-oo, 2, True, False).as_relational(x) == And(Lt(-oo, x), Le(x, 2))
    assert Interval(-oo, 2, True, True).as_relational(x) == And(Lt(-oo, x), Lt(x, 2))

    assert Interval(-oo, 2, False, False).as_relational(x) == Le(x, 2)
    assert Interval(-oo, 2, False, True).as_relational(x) == Lt(x, 2)

    assert Interval(-2, oo, False, True).as_relational(x) == And(Le(-2, x), Lt(x, oo))
    assert Interval(-2, oo, True, True).as_relational(x) == And(Lt(-2, x), Lt(x, oo))

    assert Interval(-2, oo, False, False).as_relational(x) == Le(-2, x)
    assert Interval(-2, oo, True, False).as_relational(x) == Lt(-2, x)

    assert Interval(-oo, oo, True, True).as_relational(x) == And(Lt(-oo, x), Lt(x, oo))
    assert Interval(-oo, oo).as_relational(x) == true

    x = Symbol('x', extended_real=True)
    y = Symbol('y', extended_real=True)
    assert Interval(x, y).as_relational(x) == (x <= y)
    assert Interval(y, x).as_relational(x) == (y <= x)


def test_Finite_as_relational():
    assert FiniteSet(1, 2).as_relational(x) == Or(Eq(x, 1), Eq(x, 2))
    assert FiniteSet(y, -5).as_relational(x) == Or(Eq(x, y), Eq(x, -5))


def test_Union_as_relational():
    assert (Interval(0, 1) + FiniteSet(2)).as_relational(x) == \
        Or(And(Le(0, x), Le(x, 1)), Eq(x, 2))
    assert (Interval(0, 1, True, True) + FiniteSet(1)).as_relational(x) == \
        And(Lt(0, x), Le(x, 1))


def test_Intersection_as_relational():
    assert (Intersection(Interval(0, 1), FiniteSet(2),
                         evaluate=False).as_relational(x) ==
            And(And(Le(0, x), Le(x, 1)), Eq(x, 2)))


def test_EmptySet():
    assert S.EmptySet.as_relational(x) is false
    assert S.EmptySet.intersection(S.UniversalSet) == S.EmptySet
    assert S.EmptySet.boundary == S.EmptySet
    assert Interval(0, 1).symmetric_difference(S.EmptySet) == Interval(0, 1)
    assert Interval(1, 2).intersection(S.EmptySet) == S.EmptySet


def test_finite_basic():
    assert isinstance(FiniteSet(evaluate=False), FiniteSet)

    A = FiniteSet(1, 2, 3)
    B = FiniteSet(3, 4, 5)
    AorB = Union(A, B)
    AandB = A.intersection(B)
    assert A.is_subset(AorB) and B.is_subset(AorB)
    assert AandB.is_subset(A)
    assert AandB == FiniteSet(3)

    assert A.inf == 1 and A.sup == 3
    assert AorB.inf == 1 and AorB.sup == 5
    assert FiniteSet(x, 1, 5).sup == Max(x, 5)
    assert FiniteSet(x, 1, 5).inf == Min(x, 1)

    # issue sympy/sympy#7335
    assert FiniteSet(S.EmptySet) != S.EmptySet
    assert FiniteSet(FiniteSet(1, 2, 3)) != FiniteSet(1, 2, 3)
    assert FiniteSet((1, 2, 3)) != FiniteSet(1, 2, 3)

    # Ensure a variety of types can exist in a FiniteSet
    assert FiniteSet((1, 2), Float, A, -5, x, 'eggs', x**2, Interval)

    assert (A > B) is False
    assert (A >= B) is False
    assert (A < B) is False
    assert (A <= B) is False
    assert AorB > A and AorB > B
    assert AorB >= A and AorB >= B
    assert A >= A and A <= A
    assert A >= AandB and B >= AandB
    assert A > AandB and B > AandB

    assert (FiniteSet(pi, E).evalf() ==
            FiniteSet(Float('2.7182818284590451', dps=15),
                      Float('3.1415926535897931', dps=15)))

    # issue sympy/sympy#10337
    assert (FiniteSet(2) == 3) is False
    assert (FiniteSet(2) != 3) is True

    pytest.raises(TypeError, lambda: FiniteSet(2) < 3)
    pytest.raises(TypeError, lambda: FiniteSet(2) <= 3)
    pytest.raises(TypeError, lambda: FiniteSet(2) > 3)
    pytest.raises(TypeError, lambda: FiniteSet(2) >= 3)


def test_powerset():
    # EmptySet
    A = FiniteSet()
    pset = A.powerset()
    assert len(pset) == 1
    assert pset == FiniteSet(S.EmptySet)

    # FiniteSets
    A = FiniteSet(1, 2)
    pset = A.powerset()
    assert len(pset) == 2**len(A)
    assert pset == FiniteSet(FiniteSet(), FiniteSet(1),
                             FiniteSet(2), A)
    # Not finite sets
    I = Interval(0, 1)
    pytest.raises(NotImplementedError, I.powerset)


def test_product_basic():
    H, T = 'H', 'T'
    unit_line = Interval(0, 1)
    d6 = FiniteSet(1, 2, 3, 4, 5, 6)
    d4 = FiniteSet(1, 2, 3, 4)
    coin = FiniteSet(H, T)

    square = unit_line * unit_line

    assert (0, 0) in square
    assert 0 not in square
    assert (H, T) in coin ** 2
    assert (.5, .5, .5) in square * unit_line
    assert (H, 3, 3) in coin*d6*d6
    HH, TT = sympify(H), sympify(T)
    assert set(coin**2) == {(HH, HH), (HH, TT), (TT, HH), (TT, TT)}

    assert (d4*d4).is_subset(d6*d6)

    assert (square.complement(Interval(-oo, oo)*Interval(-oo, oo)) ==
            Union((Interval(-oo, 0, False, True) +
                   Interval(1, oo, True))*Interval(-oo, oo),
                  Interval(-oo, oo)*(Interval(-oo, 0, False, True) +
                                     Interval(1, oo, True))))

    assert (Interval(-5, 5)**3).is_subset(Interval(-10, 10)**3)
    assert not (Interval(-10, 10)**3).is_subset(Interval(-5, 5)**3)
    assert not (Interval(-5, 5)**2).is_subset(Interval(-10, 10)**3)

    assert (Interval(.2, .5)*FiniteSet(.5)).is_subset(square)  # segment in square

    assert len(coin*coin*coin) == 8
    assert len(S.EmptySet*S.EmptySet) == 0
    assert len(S.EmptySet*coin) == 0
    pytest.raises(TypeError, lambda: len(coin*Interval(0, 2)))


def test_real():
    x = Symbol('x', real=True)

    I = Interval(0, 5)
    J = Interval(10, 20)
    A = FiniteSet(1, 2, 30, x, pi)
    B = FiniteSet(-4, 0)
    C = FiniteSet(100)
    D = FiniteSet('Ham', 'Eggs')

    assert all(s.is_subset(S.Reals) for s in [I, J, A, B, C])
    assert not D.is_subset(S.Reals)
    assert all((a + b).is_subset(S.Reals) for a in [I, J, A, B, C] for b in [I, J, A, B, C])
    assert not any((a + D).is_subset(S.Reals) for a in [I, J, A, B, C, D])

    assert not (I + A + D).is_subset(S.Reals)


def test_supinf():
    x = Symbol('x', extended_real=True)
    y = Symbol('y', extended_real=True)

    assert (Interval(0, 1) + FiniteSet(2)).sup == 2
    assert (Interval(0, 1) + FiniteSet(2)).inf == 0
    assert (Interval(0, 1) + FiniteSet(x)).sup == Max(1, x)
    assert (Interval(0, 1) + FiniteSet(x)).inf == Min(0, x)
    assert FiniteSet(5, 1, x).sup == Max(5, x)
    assert FiniteSet(5, 1, x).inf == Min(1, x)
    assert FiniteSet(5, 1, x, y).sup == Max(5, x, y)
    assert FiniteSet(5, 1, x, y).inf == Min(1, x, y)
    assert FiniteSet(5, 1, x, y, oo, -oo).sup == +oo
    assert FiniteSet(5, 1, x, y, oo, -oo).inf == -oo
    assert FiniteSet('Ham', 'Eggs').sup == Max('Ham', 'Eggs')


def test_universalset():
    U = S.UniversalSet
    assert U.as_relational(x) is true
    assert U.union(Interval(2, 4)) == U

    assert U.intersection(Interval(2, 4)) == Interval(2, 4)
    assert U.measure == oo
    assert U.boundary == S.EmptySet
    assert U.contains(0) is true
    assert Interval(0, 1).symmetric_difference(U) == Interval(0, 1)
    assert U.complement(U) == S.EmptySet


def test_Union_of_ProductSets_shares():
    line = Interval(0, 2)
    points = FiniteSet(0, 1, 2)
    assert Union(line * line, line * points) == line * line


def test_Interval_free_symbols():
    # issue sympy/sympy#6211
    assert Interval(0, 1).free_symbols == set()
    x = Symbol('x', extended_real=True)
    assert Interval(0, x).free_symbols == {x}


def test_image_interval():
    x = Symbol('x', extended_real=True)
    a = Symbol('a', extended_real=True)
    assert imageset(x, 2*x, Interval(-2, 1)) == Interval(-4, 2)
    assert imageset(x, 2*x, Interval(-2, 1, True, False)) == \
        Interval(-4, 2, True, False)
    assert imageset(x, x**2, Interval(-2, 1, True, False)) == \
        Interval(0, 4, False, True)
    assert imageset(x, x**2, Interval(-2, 1)) == Interval(0, 4)
    assert imageset(x, x**2, Interval(-2, 1, True, False)) == \
        Interval(0, 4, False, True)
    assert imageset(x, x**2, Interval(-2, 1, True, True)) == \
        Interval(0, 4, False, True)
    assert imageset(x, (x - 2)**2, Interval(1, 3)) == Interval(0, 1)
    assert imageset(x, 3*x**4 - 26*x**3 + 78*x**2 - 90*x, Interval(0, 4)) == \
        Interval(-35, 0)  # Multiple Maxima
    assert imageset(x, x + 1/x, Interval(-oo, oo)) == Interval(-oo, -2) \
        + Interval(2, oo)  # Single Infinite discontinuity
    assert imageset(x, 1/x + 1/(x - 1)**2, Interval(0, 2, True, False)) == \
        Interval(Rational(3, 2), oo, False, True)  # Multiple Infinite discontinuities

    # issue sympy/sympy#10113
    assert imageset(x, x**2/(x**2 - 4),
                    Interval(-2, 2)) == Interval(-oo, 0, True)

    # Test for Python lambda
    assert imageset(lambda x: 2*x, Interval(-2, 1)) == Interval(-4, 2)

    assert (imageset(Lambda(x, a*x), Interval(0, 1)) ==
            ImageSet(Lambda(x, a*x), Interval(0, 1)))

    assert (imageset(Lambda(x, sin(cos(x))), Interval(0, 1)) ==
            ImageSet(Lambda(x, sin(cos(x))), Interval(0, 1)))

    assert imageset(x, -(x - 2)*(x + 2),
                    Interval(-3, 4)) == Interval(-12, 4)

    assert (imageset(z, 2*z, ImageSet(Lambda((x, y), x*y),
                                      Interval(0, 2))) ==
            ImageSet(Lambda(z, 2*z), ImageSet(Lambda((x, y), x*y),
                                              Interval(0, 2))))


def test_image_piecewise():
    f = Piecewise((x, x <= -1), (1/x**2, x <= 5), (x**3, True))
    f1 = Piecewise((0, x <= 1), (1, x <= 2), (2, True))
    f2 = Piecewise((x, x <= -1), (x**3, True))
    assert imageset(x, f, Interval(-5, 5)) == Union(Interval(-5, -1),
                                                    Interval(Rational(1, 25),
                                                             oo, false, true))
    assert imageset(x, f1, Interval(1, 2)) == FiniteSet(0, 1)
    assert imageset(x, f2, Interval(-2, 2)) == Interval(-2, 8)


@pytest.mark.xfail  # See: https://github.com/sympy/sympy/pull/2723#discussion_r8659826
def test_image_Intersection():
    x = Symbol('x', extended_real=True)
    y = Symbol('y', extended_real=True)
    assert (imageset(x, x**2, Interval(-2, 0).intersection(Interval(x, y))) ==
            Interval(0, 4).intersection(Interval(Min(x**2, y**2),
                                                 Max(x**2, y**2))))


def test_image_FiniteSet():
    x = Symbol('x', extended_real=True)
    assert imageset(x, 2*x, FiniteSet(1, 2, 3)) == FiniteSet(2, 4, 6)


def test_image_Union():
    x = Symbol('x', extended_real=True)
    assert (imageset(x, x**2, Interval(-2, 0) + FiniteSet(1, 2, 3)) ==
            Interval(0, 4) + FiniteSet(9))


def test_image_EmptySet():
    x = Symbol('x', extended_real=True)
    assert imageset(x, 2*x, S.EmptySet) == S.EmptySet


def test_sympyissue_5724_7680():
    assert I not in S.Reals  # issue sympy/sympy#7680
    assert Interval(-oo, oo).contains(I) is false


def test_boundary():
    assert FiniteSet(1).boundary == FiniteSet(1)
    assert all(Interval(0, 1, left_open, right_open).boundary == FiniteSet(0, 1)
               for left_open in (true, false) for right_open in (true, false))


def test_boundary_Union():
    assert (Interval(0, 1) + Interval(2, 3)).boundary == FiniteSet(0, 1, 2, 3)
    assert ((Interval(0, 1, False, True) +
             Interval(1, 2, True, False)).boundary == FiniteSet(0, 1, 2))

    assert (Interval(0, 1) + FiniteSet(2)).boundary == FiniteSet(0, 1, 2)
    assert (Union(Interval(0, 10), Interval(5, 15), evaluate=False).boundary ==
            FiniteSet(0, 15))

    assert (Union(Interval(0, 10), Interval(0, 1), evaluate=False).boundary ==
            FiniteSet(0, 10))
    assert (Union(Interval(0, 10, True, True),
                  Interval(10, 15, True, True), evaluate=False).boundary ==
            FiniteSet(0, 10, 15))


@pytest.mark.xfail
def test_union_boundary_of_joining_sets():
    assert (Union(Interval(0, 10), Interval(10, 15), evaluate=False).boundary ==
            FiniteSet(0, 15))


def test_boundary_ProductSet():
    open_square = Interval(0, 1, True, True) ** 2
    assert open_square.boundary == (FiniteSet(0, 1) * Interval(0, 1) +
                                    Interval(0, 1) * FiniteSet(0, 1))

    second_square = Interval(1, 2, True, True) * Interval(0, 1, True, True)
    assert ((open_square + second_square).boundary ==
            (FiniteSet(0, 1)*Interval(0, 1) + FiniteSet(1, 2)*Interval(0, 1) +
             Interval(0, 1)*FiniteSet(0, 1) + Interval(1, 2)*FiniteSet(0, 1)))


def test_boundary_ProductSet_line():
    line_in_r2 = Interval(0, 1) * FiniteSet(0)
    assert line_in_r2.boundary == line_in_r2


def test_is_open():
    assert not Interval(0, 1, False, False).is_open
    assert not Interval(0, 1, True, False).is_open
    assert Interval(0, 1, True, True).is_open
    assert not FiniteSet(1, 2, 3).is_open


def test_is_closed():
    assert Interval(0, 1, False, False).is_closed
    assert not Interval(0, 1, True, False).is_closed
    assert FiniteSet(1, 2, 3).is_closed


def test_closure():
    assert Interval(0, 1, False, True).closure == Interval(0, 1, False, False)


def test_interior():
    assert Interval(0, 1, False, True).interior == Interval(0, 1, True, True)


def test_sympyissue_7841():
    pytest.raises(TypeError, lambda: x in S.Reals)


def test_Eq():
    assert Eq(Interval(0, 1), Interval(0, 1))
    assert Eq(Interval(0, 1), Union(Interval(2, 3),
                                    Interval(4, 5))).is_Equality
    assert Eq(Interval(0, 1), Interval(0, 2)) is false

    s1 = FiniteSet(0, 1)
    s2 = FiniteSet(1, 2)

    assert Eq(s1, s1)
    assert Eq(s1, s2) is false
    assert Eq(FiniteSet(1, 2), FiniteSet(3, 4, 5)) is false

    assert Eq(s1*s2, s1*s2)
    assert Eq(s1*s2, s2*s1) is false

    p1 = ProductSet(FiniteSet(1), FiniteSet(2))
    assert Eq(p1, s1).is_Equality
    p2 = ProductSet(FiniteSet(1), FiniteSet(2), FiniteSet(3))
    assert Eq(p1, p2) is false


def test_SymmetricDifference():
    assert (SymmetricDifference(FiniteSet(0, 1, 2, 3, 4, 5),
                                FiniteSet(2, 4, 6, 8, 10)) ==
            FiniteSet(0, 1, 3, 5, 6, 8, 10))
    assert (SymmetricDifference(FiniteSet(2, 3, 4),
                                FiniteSet(2, 3, 4, 5)) ==
            FiniteSet(5))
    assert FiniteSet(2).symmetric_difference(FiniteSet(2, 5)) == FiniteSet(5)
    assert (FiniteSet(1, 2, 3, 4, 5) ^ FiniteSet(1, 2, 5, 6) ==
            FiniteSet(3, 4, 6))
    assert (Set(Integer(1), Integer(2),
                Integer(3)) ^ Set(Integer(2), Integer(3), Integer(4)) ==
            Union(Set(Integer(1), Integer(2),
                      Integer(3)) - Set(Integer(2), Integer(3), Integer(4)),
                  Set(Integer(2), Integer(3),
                      Integer(4)) - Set(Integer(1), Integer(2), Integer(3))))
    assert (Interval(0, 4) ^ Interval(2, 5) ==
            Union(Interval(0, 4) - Interval(2, 5),
                  Interval(2, 5) - Interval(0, 4)))

    class TestSet(Set):
        def _symmetric_difference(self, other):
            return

    t1, t2 = TestSet(1), TestSet(2)
    assert t1.symmetric_difference(t2) == SymmetricDifference(t1, t2,
                                                              evaluate=False)


def test_sympyissue_9956():
    assert Union(Interval(-oo, oo), FiniteSet(1)) == Interval(-oo, oo)
    assert Interval(-oo, oo).contains(1) is true


def test_sympyissue_9536():
    a = Symbol('a', real=True)
    assert FiniteSet(log(a)).intersection(S.Reals) == Intersection(S.Reals, FiniteSet(log(a)))


def test_sympyissue_9637():
    n = Symbol('n')
    a = FiniteSet(n)
    b = FiniteSet(2, n)
    assert Complement(S.Reals, a) == Complement(S.Reals, a, evaluate=False)
    assert Complement(Interval(1, 3), a) == Complement(Interval(1, 3), a, evaluate=False)
    assert (Complement(Interval(1, 3), b) ==
            Complement(Union(Interval(1, 2, right_open=True),
                             Interval(2, 3, left_open=True)), a, evaluate=False))
    assert Complement(a, S.Reals) == Complement(a, S.Reals, evaluate=False)
    assert Complement(a, Interval(1, 3)) == Complement(a, Interval(1, 3), evaluate=False)


def test_sympyissue_10113():
    f = x**2/(x**2 - 4)
    assert imageset(x, f, S.Reals) == Union(Interval(-oo, 0, True),
                                            Interval(1, oo, True, True))
    assert imageset(x, f, Interval(-2, 3)) == Union(Interval(-oo, 0, True),
                                                    Interval(Rational(9, 5),
                                                             oo, False, True))


def test_sympyissue_9808():
    assert (Complement(FiniteSet(y), FiniteSet(1)) ==
            Complement(FiniteSet(y), FiniteSet(1), evaluate=False))
    assert (Complement(FiniteSet(1, 2, x), FiniteSet(x, y, 2, 3)) ==
            Complement(FiniteSet(1), FiniteSet(y), evaluate=False))


def test_sympyissue_9447():
    a = Interval(0, 1) + Interval(2, 3)
    assert (Complement(S.UniversalSet, a) ==
            Complement(S.UniversalSet,
                       Union(Interval(0, 1), Interval(2, 3)), evaluate=False))
    # issue sympy/sympy#10305:
    assert (Complement(S.Naturals, a) ==
            Complement(S.Naturals,
                       Union(Interval(0, 1), Interval(2, 3)), evaluate=False))


def test_sympyissue_2799():
    U = S.UniversalSet
    a = Symbol('a', real=True)
    inf_interval = Interval(a, oo)
    R = S.Reals

    assert U + inf_interval == inf_interval + U
    assert U + R == R + U
    assert R + inf_interval == inf_interval + R


def test_sympyissue_9706():
    assert Interval(-oo, 0).closure == Interval(-oo, 0)
    assert Interval(0, oo).closure == Interval(0, oo)
    assert Interval(-oo, oo).closure == Interval(-oo, oo)


def test_sympyissue_8257():
    assert Interval(-oo, oo) + FiniteSet(+oo) == Interval(-oo, oo)
    assert FiniteSet(+oo) + Interval(-oo, oo) == Interval(-oo, oo)
    assert Interval(-oo, oo) + FiniteSet(-oo) == Interval(-oo, oo)
    assert FiniteSet(-oo) + Interval(-oo, oo) == Interval(-oo, oo)


def test_sympyissue_10931():
    assert S.Integers - S.Integers == S.EmptySet
    assert S.Integers - S.Reals == S.EmptySet
