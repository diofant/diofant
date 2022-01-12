# this module tests that diofant works with true division turned on

from diofant import Float, Integer, Symbol


__all__ = ()


def test_truediv():
    assert 1/2 != 0
    assert Integer(1)/2 != 0


def dotest(s):
    x = Symbol('x')
    y = Symbol('y')
    l = [Integer(2), Float('1.3'), x, y, pow(x, y)*y, 5, 5.5]
    for x in l:
        for y in l:
            s(x, y)
    return True


def test_basic():
    def s(a, b):
        # pylint: disable=pointless-statement
        a
        +a
        -a
        a + b
        a - b
        a*b
        a/b
        a**b
    assert dotest(s)


def test_ibasic():
    def s(a, b):
        # pylint: disable=pointless-statement
        x = a
        x += b
        x = a
        x -= b
        x = a
        x *= b
        x = a
        x /= b
    assert dotest(s)
