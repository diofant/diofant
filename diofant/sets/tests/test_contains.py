from diofant import Symbol, Contains, S, Interval, FiniteSet

from diofant.abc import x


def test_contains_basic():
    assert Contains(2, S.Integers) is S.true
    assert Contains(-2, S.Naturals) is S.false

    i = Symbol('i', integer=True)
    assert Contains(i, S.Naturals) == Contains(i, S.Naturals, evaluate=False)
    assert S.Naturals0.contains(x) == Contains(x, S.Naturals0, evaluate=False)


def test_issue_6194():
    assert Contains(x, Interval(0, 1)) == (x >= 0) & (x <= 1)
    assert Contains(x, FiniteSet(0)) != S.false
    assert Contains(x, Interval(1, 1)) != S.false
    assert Contains(x, S.Integers) != S.false
