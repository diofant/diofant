from mpmath import mpi
import pytest

from sympy import (Symbol, Set, Union, Interval, oo, S, sympify, nan,
                   LessThan, Max, Min, And, Or, Eq, Le,
                   Lt, Float, FiniteSet, Intersection, imageset, I, true, false,
                   ProductSet, sqrt, Complement, EmptySet, sin, cos, Lambda,
                   ImageSet, pi, Eq, Pow, Contains, Sum, RootOf,
                   SymmetricDifference, Integer, Rational)

from sympy.abc import a, b, x, y, z


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
    pytest.raises(ValueError, lambda: Interval(0, S.ImaginaryUnit))
    pytest.raises(ValueError, lambda: Interval(0, Symbol('z', extended_real=False)))

    pytest.raises(NotImplementedError, lambda: Interval(0, 1, And(x, y)))
    pytest.raises(NotImplementedError, lambda: Interval(0, 1, False, And(x, y)))
    pytest.raises(NotImplementedError, lambda: Interval(0, 1, z, And(x, y)))


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

    x = Symbol("x")
    y = Symbol("y")
    z = Symbol("z")
    assert S.EmptySet | FiniteSet(x, FiniteSet(y, z)) == \
        FiniteSet(x, FiniteSet(y, z))

    # Test that Intervals and FiniteSets play nicely
    assert Interval(1, 3) + FiniteSet(2) == Interval(1, 3)
    assert Interval(1, 3, True, True) + FiniteSet(3) == \
        Interval(1, 3, True, False)
    X = Interval(1, 3) + FiniteSet(5)
    Y = Interval(1, 2) + FiniteSet(3)
    XandY = X.intersect(Y)
    assert 2 in X and 3 in X and 3 in XandY
    assert XandY.is_subset(X) and XandY.is_subset(Y)

    pytest.raises(TypeError, lambda: Union(1, 2, 3))

    assert X.is_iterable is False

    # issue 7843
    assert Union(S.EmptySet, FiniteSet(-sqrt(-I), sqrt(-I))) == FiniteSet(-sqrt(-I), sqrt(-I))


def test_difference():
    assert Interval(1, 3) - Interval(1, 2) == Interval(2, 3, True)
    assert Interval(1, 3) - Interval(2, 3) == Interval(1, 2, False, True)
    assert Interval(1, 3, True) - Interval(2, 3) == Interval(1, 2, True, True)
    assert Interval(1, 3, True) - Interval(2, 3, True) == \
        Interval(1, 2, True, False)
    assert Interval(0, 2) - FiniteSet(1) == \
        Union(Interval(0, 1, False, True), Interval(1, 2, True, False))

    assert FiniteSet(1, 2, 3) - FiniteSet(2) == FiniteSet(1, 3)
    assert FiniteSet('ham', 'eggs') - FiniteSet('eggs') == FiniteSet('ham')
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

    assert (0 not in S.Reals.intersect(S.Integers - FiniteSet(0)))

    assert S.EmptySet - S.Integers == S.EmptySet

    assert (S.Integers - FiniteSet(0)) - FiniteSet(1) == S.Integers - FiniteSet(0, 1)

    assert S.Reals - Union(S.Naturals, FiniteSet(pi)) == \
            Intersection(S.Reals - S.Naturals, S.Reals - FiniteSet(pi))


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

    assert (FiniteSet(5) + Interval(S.NegativeInfinity,
                                    0)).complement(S.Reals) == \
        Interval(0, 5, True, True) + Interval(5, S.Infinity, True, True)

    assert FiniteSet(1, 2, 3).complement(S.Reals) == \
        Interval(S.NegativeInfinity, 1, True, True) + \
        Interval(1, 2, True, True) + Interval(2, 3, True, True) +\
        Interval(3, S.Infinity, True, True)

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


def test_intersect():
    x = Symbol('x')
    assert Interval(0, 2).intersect(Interval(1, 2)) == Interval(1, 2)
    assert Interval(0, 2).intersect(Interval(1, 2, True)) == \
        Interval(1, 2, True)
    assert Interval(0, 2, True).intersect(Interval(1, 2)) == \
        Interval(1, 2, False, False)
    assert Interval(0, 2, True, True).intersect(Interval(1, 2)) == \
        Interval(1, 2, False, True)
    assert Interval(0, 2).intersect(Union(Interval(0, 1), Interval(2, 3))) == \
        Union(Interval(0, 1), Interval(2, 2))

    assert FiniteSet(1, 2, x).intersect(FiniteSet(x)) == FiniteSet(x)
    assert FiniteSet('ham', 'eggs').intersect(FiniteSet('ham')) == \
        FiniteSet('ham')
    assert FiniteSet(1, 2, 3, 4, 5).intersect(S.EmptySet) == S.EmptySet

    assert Interval(0, 5).intersect(FiniteSet(1, 3)) == FiniteSet(1, 3)
    assert Interval(0, 1, True, True).intersect(FiniteSet(1)) == S.EmptySet

    assert Union(Interval(0, 1), Interval(2, 3)).intersect(Interval(1, 2)) == \
        Union(Interval(1, 1), Interval(2, 2))
    assert Union(Interval(0, 1), Interval(2, 3)).intersect(Interval(0, 2)) == \
        Union(Interval(0, 1), Interval(2, 2))
    assert Union(Interval(0, 1), Interval(2, 3)).intersect(Interval(1, 2, True, True)) == \
        S.EmptySet
    assert Union(Interval(0, 1), Interval(2, 3)).intersect(S.EmptySet) == \
        S.EmptySet
    assert Union(Interval(0, 5), FiniteSet('ham')).intersect(FiniteSet(2, 3, 4, 5, 6)) == \
        FiniteSet(2, 3, 4, 5)

    # tests for the intersection alias
    assert Interval(0, 5).intersection(FiniteSet(1, 3)) == FiniteSet(1, 3)
    assert Interval(0, 1, True, True).intersection(FiniteSet(1)) == S.EmptySet

    assert Union(Interval(0, 1), Interval(2, 3)).intersection(Interval(1, 2)) == \
        Union(Interval(1, 1), Interval(2, 2))


def test_intersection():
    # iterable
    i = Intersection(FiniteSet(1, 2, 3), Interval(2, 5), evaluate=False)
    assert i.is_iterable
    assert set(i) == {Integer(2), Integer(3)}

    # challenging intervals
    x = Symbol('x', extended_real=True)
    i = Intersection(Interval(0, 3), Interval(x, 6))
    assert (5 in i) is False
    pytest.raises(TypeError, lambda: 2 in i)

    # Singleton special cases
    assert Intersection(Interval(0, 1), S.EmptySet) == S.EmptySet
    assert Intersection(S.Reals, Interval(-oo, x, True)) == Interval(-oo, x, True)

    # Products
    line = Interval(0, 5)
    i = Intersection(line**2, line**3, evaluate=False)
    assert (2, 2) not in i
    assert (2, 2, 2) not in i
    pytest.raises(ValueError, lambda: list(i))

    assert Intersection(Intersection(S.Integers, S.Naturals, evaluate=False),
                        S.Reals, evaluate=False) == \
            Intersection(S.Integers, S.Naturals, S.Reals, evaluate=False)


def test_is_disjoint():
    assert Interval(0, 2).is_disjoint(Interval(1, 2)) is False
    assert Interval(0, 2).is_disjoint(Interval(3, 4)) is True


def test_ProductSet_of_single_arg_is_arg():
    assert ProductSet(Interval(0, 1)) == Interval(0, 1)


def test_interval_subs():
    a = Symbol('a', extended_real=True)

    assert Interval(0, a).subs(a, 2) == Interval(0, 2)
    assert Interval(a, 0).subs(a, 2) == S.EmptySet


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
    assert (square.intersect(band)).measure == 20
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
    assert Interval(0, 2).contains(1) is S.true
    assert Interval(0, 2).contains(3) is S.false
    assert Interval(0, 2, True, False).contains(0) is S.false
    assert Interval(0, 2, True, False).contains(2) is S.true
    assert Interval(0, 2, False, True).contains(0) is S.true
    assert Interval(0, 2, False, True).contains(2) is S.false
    assert Interval(0, 2, True, True).contains(0) is S.false
    assert Interval(0, 2, True, True).contains(2) is S.false
    # issue sympy/sympy#10326
    assert S.Reals.contains(oo) is S.false
    assert S.Reals.contains(-oo) is S.false
    assert Interval(-oo, oo, True).contains(oo) is S.true
    assert Interval(-oo, oo).contains(-oo) is S.true

    assert FiniteSet(1, 2, 3).contains(2) is S.true
    assert FiniteSet(1, 2, Symbol('x')).contains(Symbol('x')) is S.true

    # issue 8197
    assert isinstance(FiniteSet(b).contains(-a), Contains)
    assert isinstance(FiniteSet(b).contains(a), Contains)
    assert isinstance(FiniteSet(a).contains(1), Contains)
    pytest.raises(TypeError, lambda: 1 in FiniteSet(a))

    # issue 8209
    rad1 = Integer(4)
    rad2 = Pow(2, 2, evaluate=False)
    s1 = FiniteSet(rad1)
    s2 = FiniteSet(rad2)
    assert s1 - s2 == S.EmptySet

    items = [1, 2, S.Infinity, Symbol('ham'), -1.1]
    fset = FiniteSet(*items)
    assert all(item in fset for item in items)
    assert all(fset.contains(item) is S.true for item in items)

    assert Union(Interval(0, 1), Interval(2, 5)).contains(3) is S.true
    assert Union(Interval(0, 1), Interval(2, 5)).contains(6) is S.false
    assert Union(Interval(0, 1), FiniteSet(2, 5)).contains(3) is S.false

    assert S.EmptySet.contains(1) is S.false
    assert FiniteSet(RootOf(x**3 + x - 1, 0)).contains(S.Infinity) is S.false

    assert RootOf(x**5 + x**3 + 1, 0) in S.Reals
    assert not RootOf(x**5 + x**3 + 1, 1) in S.Reals


def test_interval_symbolic():
    x = Symbol('x')
    e = Interval(0, 1)
    assert e.contains(x) == And(0 <= x, x <= 1)
    pytest.raises(TypeError, lambda: x in e)
    e = Interval(0, 1, True, True)
    assert e.contains(x) == And(0 < x, x < 1)


def test_union_contains():
    x = Symbol('x')
    i1 = Interval(0, 1)
    i2 = Interval(2, 3)
    i3 = Union(i1, i2)
    pytest.raises(TypeError, lambda: x in i3)
    e = i3.contains(x)
    assert e == Or(And(0 <= x, x <= 1), And(2 <= x, x <= 3))
    assert e.subs(x, -0.5) is false
    assert e.subs(x, 0.5) is true
    assert e.subs(x, 1.5) is false
    assert e.subs(x, 2.5) is true
    assert e.subs(x, 3.5) is false

    U = Interval(0, 2, True, True) + Interval(10, oo) + FiniteSet(-1, 2, 5, 6)
    assert all(el not in U for el in [0, 4, -oo])
    assert all(el in U for el in [2, 5, 10])


def test_is_number():
    assert Interval(0, 1).is_number is False
    assert Set().is_number is False


def test_Interval_is_left_unbounded():
    assert Interval(3, 4).is_left_unbounded is False
    assert Interval(-oo, 3).is_left_unbounded is True
    assert Interval(Float("-inf"), 3).is_left_unbounded is True


def test_Interval_is_right_unbounded():
    assert Interval(3, 4).is_right_unbounded is False
    assert Interval(3, oo).is_right_unbounded is True
    assert Interval(3, Float("+inf")).is_right_unbounded is True


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

    assert Interval(-oo, 2, right_open=False).as_relational(x) == And(Le(-oo, x), Le(x, 2))
    assert Interval(-oo, 2, right_open=True).as_relational(x) == And(Le(-oo, x), Lt(x, 2))

    assert Interval(-2, oo, left_open=False).as_relational(x) == And(Le(-2, x), Le(x, oo))
    assert Interval(-2, oo, left_open=True).as_relational(x) == And(Lt(-2, x), Le(x, oo))

    assert Interval(-oo, oo).as_relational(x) == And(Le(-oo, x), Le(x, oo))
    x = Symbol('x', extended_real=True)
    y = Symbol('y', extended_real=True)
    assert Interval(x, y).as_relational(x) == (x <= y)
    assert Interval(y, x).as_relational(x) == (y <= x)


def test_Finite_as_relational():
    x = Symbol('x')
    y = Symbol('y')

    assert FiniteSet(1, 2).as_relational(x) == Or(Eq(x, 1), Eq(x, 2))
    assert FiniteSet(y, -5).as_relational(x) == Or(Eq(x, y), Eq(x, -5))


def test_Union_as_relational():
    x = Symbol('x')
    assert (Interval(0, 1) + FiniteSet(2)).as_relational(x) == \
        Or(And(Le(0, x), Le(x, 1)), Eq(x, 2))
    assert (Interval(0, 1, True, True) + FiniteSet(1)).as_relational(x) == \
        And(Lt(0, x), Le(x, 1))


def test_Intersection_as_relational():
    x = Symbol('x')
    assert (Intersection(Interval(0, 1), FiniteSet(2),
            evaluate=False).as_relational(x)
            == And(And(Le(0, x), Le(x, 1)), Eq(x, 2)))


def test_EmptySet():
    assert S.EmptySet.as_relational(Symbol('x')) is False
    assert S.EmptySet.intersect(S.UniversalSet) == S.EmptySet
    assert S.EmptySet.boundary == S.EmptySet


def test_finite_basic():
    x = Symbol('x')
    A = FiniteSet(1, 2, 3)
    B = FiniteSet(3, 4, 5)
    AorB = Union(A, B)
    AandB = A.intersect(B)
    assert A.is_subset(AorB) and B.is_subset(AorB)
    assert AandB.is_subset(A)
    assert AandB == FiniteSet(3)

    assert A.inf == 1 and A.sup == 3
    assert AorB.inf == 1 and AorB.sup == 5
    assert FiniteSet(x, 1, 5).sup == Max(x, 5)
    assert FiniteSet(x, 1, 5).inf == Min(x, 1)

    # issue 7335
    assert FiniteSet(S.EmptySet) != S.EmptySet
    assert FiniteSet(FiniteSet(1, 2, 3)) != FiniteSet(1, 2, 3)
    assert FiniteSet((1, 2, 3)) != FiniteSet(1, 2, 3)

    # Ensure a variety of types can exist in a FiniteSet
    s = FiniteSet((1, 2), Float, A, -5, x, 'eggs', x**2, Interval)

    assert (A > B) is False
    assert (A >= B) is False
    assert (A < B) is False
    assert (A <= B) is False
    assert AorB > A and AorB > B
    assert AorB >= A and AorB >= B
    assert A >= A and A <= A
    assert A >= AandB and B >= AandB
    assert A > AandB and B > AandB


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

    assert square.complement(Interval(-oo, oo)*Interval(-oo, oo)) == Union(
        (Interval(-oo, 0, True, True) +
         Interval(1, oo, True, True))*Interval(-oo, oo),
        Interval(-oo, oo)*(Interval(-oo, 0, True, True) +
                  Interval(1, oo, True, True)))

    assert (Interval(-5, 5)**3).is_subset(Interval(-10, 10)**3)
    assert not (Interval(-10, 10)**3).is_subset(Interval(-5, 5)**3)
    assert not (Interval(-5, 5)**2).is_subset(Interval(-10, 10)**3)

    assert (Interval(.2, .5)*FiniteSet(.5)).is_subset(square)  # segment in square

    assert len(coin*coin*coin) == 8
    assert len(S.EmptySet*S.EmptySet) == 0
    assert len(S.EmptySet*coin) == 0
    pytest.raises(TypeError, lambda: len(coin*Interval(0, 2)))


def test_real():
    x = Symbol('x', extended_real=True, finite=True)

    I = Interval(0, 5)
    J = Interval(10, 20)
    A = FiniteSet(1, 2, 30, x, S.Pi)
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
    assert FiniteSet(5, 1, x, y, S.Infinity, S.NegativeInfinity).sup == \
        S.Infinity
    assert FiniteSet(5, 1, x, y, S.Infinity, S.NegativeInfinity).inf == \
        S.NegativeInfinity
    assert FiniteSet('Ham', 'Eggs').sup == Max('Ham', 'Eggs')


def test_universalset():
    U = S.UniversalSet
    x = Symbol('x')
    assert U.as_relational(x) is True
    assert U.union(Interval(2, 4)) == U

    assert U.intersect(Interval(2, 4)) == Interval(2, 4)
    assert U.measure == S.Infinity
    assert U.boundary == S.EmptySet
    assert U.contains(0) is S.true


def test_Union_of_ProductSets_shares():
    line = Interval(0, 2)
    points = FiniteSet(0, 1, 2)
    assert Union(line * line, line * points) == line * line


def test_Interval_free_symbols():
    # issue 6211
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

    # Test for Python lambda
    assert imageset(lambda x: 2*x, Interval(-2, 1)) == Interval(-4, 2)

    assert imageset(Lambda(x, a*x), Interval(0, 1)) == \
            ImageSet(Lambda(x, a*x), Interval(0, 1))

    assert imageset(Lambda(x, sin(cos(x))), Interval(0, 1)) == \
            ImageSet(Lambda(x, sin(cos(x))), Interval(0, 1))


@pytest.mark.xfail  # See: https://github.com/sympy/sympy/pull/2723#discussion_r8659826
def test_image_Intersection():
    x = Symbol('x', extended_real=True)
    y = Symbol('y', extended_real=True)
    assert imageset(x, x**2, Interval(-2, 0).intersect(Interval(x, y))) == \
           Interval(0, 4).intersect(Interval(Min(x**2, y**2), Max(x**2, y**2)))


def test_image_FiniteSet():
    x = Symbol('x', extended_real=True)
    assert imageset(x, 2*x, FiniteSet(1, 2, 3)) == FiniteSet(2, 4, 6)


def test_image_Union():
    x = Symbol('x', extended_real=True)
    assert imageset(x, x**2, Interval(-2, 0) + FiniteSet(1, 2, 3)) == \
            (Interval(0, 4) + FiniteSet(9))


def test_image_EmptySet():
    x = Symbol('x', extended_real=True)
    assert imageset(x, 2*x, S.EmptySet) == S.EmptySet


def test_issue_5724_7680():
    assert I not in S.Reals  # issue 7680
    assert Interval(-oo, oo).contains(I) is S.false


def test_boundary():
    x = Symbol('x', extended_real=True)
    y = Symbol('y', extended_real=True)
    assert FiniteSet(1).boundary == FiniteSet(1)
    assert all(Interval(0, 1, left_open, right_open).boundary == FiniteSet(0, 1)
            for left_open in (true, false) for right_open in (true, false))


def test_boundary_Union():
    assert (Interval(0, 1) + Interval(2, 3)).boundary == FiniteSet(0, 1, 2, 3)
    assert ((Interval(0, 1, False, True)
           + Interval(1, 2, True, False)).boundary == FiniteSet(0, 1, 2))

    assert (Interval(0, 1) + FiniteSet(2)).boundary == FiniteSet(0, 1, 2)
    assert Union(Interval(0, 10), Interval(5, 15), evaluate=False).boundary \
            == FiniteSet(0, 15)

    assert Union(Interval(0, 10), Interval(0, 1), evaluate=False).boundary \
            == FiniteSet(0, 10)
    assert Union(Interval(0, 10, True, True),
                 Interval(10, 15, True, True), evaluate=False).boundary \
            == FiniteSet(0, 10, 15)


@pytest.mark.xfail
def test_union_boundary_of_joining_sets():
    """ Testing the boundary of unions is a hard problem """
    assert Union(Interval(0, 10), Interval(10, 15), evaluate=False).boundary \
            == FiniteSet(0, 15)


def test_boundary_ProductSet():
    open_square = Interval(0, 1, True, True) ** 2
    assert open_square.boundary == (FiniteSet(0, 1) * Interval(0, 1)
                                  + Interval(0, 1) * FiniteSet(0, 1))

    second_square = Interval(1, 2, True, True) * Interval(0, 1, True, True)
    assert (open_square + second_square).boundary == (
        FiniteSet(0, 1) * Interval(0, 1)
        + FiniteSet(1, 2) * Interval(0, 1)
        + Interval(0, 1) * FiniteSet(0, 1)
        + Interval(1, 2) * FiniteSet(0, 1))


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


def test_issue_7841():
    pytest.raises(TypeError, lambda: x in S.Reals)


def test_Eq():
    assert Eq(Interval(0, 1), Interval(0, 1))
    assert Eq(Interval(0, 1), Interval(0, 2)) is S.false

    s1 = FiniteSet(0, 1)
    s2 = FiniteSet(1, 2)

    assert Eq(s1, s1)
    assert Eq(s1, s2) is S.false

    assert Eq(s1*s2, s1*s2)
    assert Eq(s1*s2, s2*s1) is S.false


def test_SymmetricDifference():
    assert SymmetricDifference(FiniteSet(0, 1, 2, 3, 4, 5),
            FiniteSet(2, 4, 6, 8, 10)) == FiniteSet(0, 1, 3, 5, 6, 8, 10)
    assert SymmetricDifference(FiniteSet(2, 3, 4), FiniteSet(2, 3, 4, 5)) \
            == FiniteSet(5)
    assert FiniteSet(1, 2, 3, 4, 5) ^ FiniteSet(1, 2, 5, 6) == \
            FiniteSet(3, 4, 6)
    assert Set(1, 2, 3) ^ Set(2, 3, 4) == Union(Set(1, 2, 3) - Set(2, 3, 4),
            Set(2, 3, 4) - Set(1, 2, 3))
    assert Interval(0, 4) ^ Interval(2, 5) == Union(Interval(0, 4) -
            Interval(2, 5), Interval(2, 5) - Interval(0, 4))


def test_issue_9956():
    assert Union(Interval(-oo, oo), FiniteSet(1)) == Interval(-oo, oo)
    assert Interval(-oo, oo).contains(1) is S.true
