from diofant import Expr, Symbol
from diofant.core.decorators import call_highest_priority


__all__ = ()


class Higher(Expr):

    _op_priority = 20.0
    result = 'high'

    is_commutative = False

    @call_highest_priority('__rmul__')
    def __mul__(self, other):
        return self.result

    @call_highest_priority('__mul__')
    def __rmul__(self, other):
        return self.result

    @call_highest_priority('__radd__')
    def __add__(self, other):
        return self.result

    @call_highest_priority('__add__')
    def __radd__(self, other):
        return self.result

    @call_highest_priority('__rsub__')
    def __sub__(self, other):
        return self.result

    @call_highest_priority('__sub__')
    def __rsub__(self, other):
        return self.result

    @call_highest_priority('__rpow__')
    def __pow__(self, other):
        return self.result

    @call_highest_priority('__pow__')
    def __rpow__(self, other):
        return self.result

    @call_highest_priority('__rtruediv__')
    def __truediv__(self, other):
        return self.result

    @call_highest_priority('__truediv__')
    def __rtruediv__(self, other):
        return self.result


class Lower(Higher):

    _op_priority = 5.0
    result = 'low'


class Lower2(Higher):
    _op_priority = 5.0
    result = 'low'

    @call_highest_priority('typo')
    def __mul__(self, other):
        return self.result


class Higher2:
    result = 'high'


def test_mul():
    x = Symbol('x')
    h = Higher()
    l = Lower()
    l2 = Lower2()
    h2 = Higher2()
    assert l*h == h*l == 'high'
    assert x*h == h*x == 'high'
    assert l*x == x*l != 'low'

    assert l2*h == 'low'
    assert l2*h2 == 'low'


def test_add():
    x = Symbol('x')
    h = Higher()
    l = Lower()
    assert l + h == h + l == 'high'
    assert x + h == h + x == 'high'
    assert l + x == x + l != 'low'


def test_sub():
    x = Symbol('x')
    h = Higher()
    l = Lower()
    assert l - h == h - l == 'high'
    assert x - h == h - x == 'high'
    assert l - x == -(x - l) != 'low'


def test_pow():
    x = Symbol('x')
    h = Higher()
    l = Lower()
    assert l**h == h**l == 'high'
    assert x**h == h**x == 'high'
    assert l**x != 'low'
    assert x**l != 'low'


def test_div():
    x = Symbol('x')
    h = Higher()
    l = Lower()
    assert l/h == h/l == 'high'
    assert x/h == h/x == 'high'
    assert l/x != 'low'
    assert x/l != 'low'
