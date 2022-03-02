from diofant import Symbol, exp
from diofant.abc import a, b, x


__all__ = ()


def test_equal():
    # pylint: disable=unneeded-not
    e1 = a + b
    e2 = 2*a*b
    e3 = a**3*b**2
    e4 = a*b + b*a
    assert not e1 == e2
    assert not e1 == e2
    assert e1 != e2
    assert e2 == e4
    assert e2 != e3
    assert not e2 == e3

    e1 = exp(x + 1/x)
    y = Symbol('x')
    e2 = exp(y + 1/y)
    assert e1 == e2
    assert not e1 != e2
    y = Symbol('y')
    e2 = exp(y + 1/y)
    assert not e1 == e2
    assert e1 != e2

    e5 = 3 + 2*x - x - x
    assert e5 == 3
    assert 3 == e5
    assert e5 != 4
    assert 4 != e5
    assert e5 != 3 + x
    assert 3 + x != e5


def test_expevalbug():
    e1 = exp(1*x)
    e3 = exp(x)
    assert e1 == e3


def test_cmp_bug1():
    class T:
        pass

    t = T()

    assert not x == t  # pylint: disable=unneeded-not
    assert x != t


def test_cmp_bug2():
    class T:
        pass

    t = T()

    assert not Symbol == t  # pylint: disable=unneeded-not
    assert Symbol != t


def test_sympyissue_4357():
    # pylint: disable=unneeded-not
    assert not Symbol == 1
    assert Symbol != 1
    assert not Symbol == 'x'
    assert Symbol != 'x'
