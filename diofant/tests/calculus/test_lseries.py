from diofant import Derivative, Integer, sin
from diofant.abc import x


__all__ = ()


def test_sin():
    e = sin(x).series(x, n=None)
    assert next(e) == x
    assert next(e) == -x**3/6
    assert next(e) == x**5/120


def test_simple():
    assert list(x.series(n=None)) == [x]
    assert list(Integer(1).series(x, n=None)) == [1]
    assert [t.doit() for t in Derivative(1 + x, x).series(x, n=None)] == [1]

    # issue sympy/sympy#5183
    s = (x + 1/x).series(n=None)
    assert list(s) == [1/x, x]

    # it would be nice if all terms were grouped, but in the
    # following case that would mean that all the terms would have
    # to be known since, for example, every term has a constant in it.
    s = ((1 + x)**7).series(x, 1, n=None)
    assert [next(s) for i in range(2)] == [128, -448 + 448*x]
