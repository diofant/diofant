from diofant import Contains, FiniteSet, Interval, S, Symbol, false, true
from diofant.abc import x


__all__ = ()


def test_contains_basic():
    assert Contains(2, S.Integers) is true
    assert Contains(-2, S.Naturals) is false

    i = Symbol('i', integer=True)
    assert Contains(i, S.Naturals) == Contains(i, S.Naturals, evaluate=False)
    assert S.Naturals0.contains(x) == Contains(x, S.Naturals0, evaluate=False)


def test_sympyissue_6194():
    assert Contains(x, Interval(0, 1)) == (x >= 0) & (x <= 1)
    assert Contains(x, FiniteSet(0)) != false
    assert Contains(x, Interval(1, 1)) != false
    assert Contains(x, S.Integers) != false
