import random

import pytest

from diofant import I, Rational
from diofant.utilities.randtest import (_randint, _randrange,
                                        random_complex_number)


__all__ = ()


def test_random_complex_number():
    random.seed(0)
    assert (random_complex_number(rational=False) ==
            2.8444218515250483 + 0.51590880588060495*I)
    random.seed(0)
    assert (random_complex_number() ==
            Rational(56888437030501, 20000000000000) +
            103181761176121*I/200000000000000)


def test_random__randrange():
    pytest.raises(ValueError, lambda: _randrange('spam'))
    pytest.raises(ValueError, lambda: _randrange([0.1, 0.5])(0))
    pytest.raises(ValueError, lambda: _randrange([1, 2, 3])(1))

    assert _randrange([1, 3, 2, 4])(2, 3) == 2


def test_random__randint():
    pytest.raises(ValueError, lambda: _randint('spam'))
    pytest.raises(ValueError, lambda: _randint([1, 3, 2, 4])(2, 1))
    pytest.raises(ValueError, lambda: _randint([3])(1, 2))

    assert _randint([1, 3, 2, 4])(2, 3) == 3
