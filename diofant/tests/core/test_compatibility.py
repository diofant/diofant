import pytest

from diofant import default_sort_key, ordered
from diofant.abc import x


__all__ = ()


def test_default_sort_key():
    def func(x):
        return x
    assert sorted([func, x, func], key=default_sort_key) == [func, func, x]


def test_ordered():
    # Issue sympy/sympy#7210 - this had been failing with python2/3 problems
    assert (list(ordered([{1: 3, 2: 4, 9: 10}, {1: 3}])) ==
            [{1: 3}, {1: 3, 2: 4, 9: 10}])
    # warnings should not be raised for identical items
    l = [1, 1]
    assert list(ordered(l, warn=True)) == l
    l = [[1], [2], [1]]
    assert list(ordered(l, warn=True)) == [[1], [1], [2]]
    pytest.raises(ValueError, lambda: list(ordered(['a', 'ab'],
                                                   keys=[lambda x: x[0]],
                                                   default=False,
                                                   warn=True)))
    pytest.raises(ValueError, lambda: list(ordered(['a', 'b'], default=False)))
