import pytest

from diofant import gamma, expand_func, beta, digamma, diff, conjugate
from diofant.core.function import ArgumentIndexError

from diofant.abc import x, y


def test_beta():
    assert isinstance(beta(x, y), beta)

    assert expand_func(beta(x, y)) == gamma(x)*gamma(y)/gamma(x + y)
    assert expand_func(beta(x, y) - beta(y, x)) == 0  # Symmetric
    assert expand_func(beta(x, y)) == expand_func(beta(x, y + 1) +
                                                  beta(x + 1, y)).simplify()

    assert diff(beta(x, y), x) == beta(x, y)*(digamma(x) - digamma(x + y))
    assert diff(beta(x, y), y) == beta(x, y)*(digamma(y) - digamma(x + y))
    pytest.raises(ArgumentIndexError, lambda: beta(x, y).fdiff(3))

    assert conjugate(beta(x, y)) == beta(conjugate(x), conjugate(y))
