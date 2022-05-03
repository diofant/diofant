"""This tests the basic submodule with (ideally) no reference to subclasses."""

import collections

import pytest

from diofant import (Add, Atom, Basic, Function, I, Integral, Lambda, cos,
                     default_sort_key, exp, gamma, preorder_traversal, sin)
from diofant.abc import w, x, y, z
from diofant.core.singleton import S
from diofant.core.singleton import SingletonWithManagedProperties as Singleton


__all__ = ()


b1 = Basic()
b2 = Basic(b1)
b3 = Basic(b2)
b21 = Basic(b2, b1)


def test_structure():
    assert b21.args == (b2, b1)
    assert b21.func(*b21.args) == b21
    assert bool(b1)


def test_equality():
    # pylint: disable=unneeded-not
    instances = [b1, b2, b3, b21, Basic(b1, b1, b1), Basic]
    for i, b_i in enumerate(instances):
        for j, b_j in enumerate(instances):
            assert (b_i == b_j) == (i == j)
            assert (b_i != b_j) == (i != j)

    assert Basic()
    assert Basic() != 0
    assert not Basic() == 0


def test_matches_basic():
    instances = [Basic(b1, b1, b2), Basic(b1, b2, b1), Basic(b2, b1, b1),
                 Basic(b1, b2), Basic(b2, b1), b2, b1]
    for i, b_i in enumerate(instances):
        for j, b_j in enumerate(instances):
            if i == j:
                assert b_j.match(b_i) == {}
            else:
                assert b_j.match(b_i) is None
    assert b1.match(b1) == {}


def test_has():
    assert b21.has(b1)
    assert b21.has(b3, b1)
    assert b21.has(Basic)
    assert not b1.has(b21, b3)
    assert not b21.has()


def test_subs():
    assert b21.subs({b2: b1}) == Basic(b1, b1)
    assert b21.subs({b2: b21}) == Basic(b21, b1)
    assert b3.subs({b2: b1}) == b2

    assert b21.subs([(b2, b1), (b1, b2)]) == Basic(b2, b2)

    assert b21.subs({b1: b2, b2: b1}) == Basic(b2, b2)

    pytest.raises(ValueError, lambda: b21.subs('bad arg'))
    pytest.raises(ValueError, lambda: b21.subs(b1, b2, b3))

    assert b21.subs(collections.ChainMap({b1: b2}, {b2: b1})) == Basic(b2, b2)
    assert b21.subs(collections.OrderedDict([(b2, b1), (b1, b2)])) == Basic(b2, b2)


def test_rewrite():
    assert sin(1).rewrite() == sin(1)

    f1 = sin(x) + cos(x)
    assert f1.rewrite(cos, exp) == exp(I*x)/2 + sin(x) + exp(-I*x)/2

    f2 = sin(x) + cos(y)/gamma(z)
    assert f2.rewrite(sin, exp) == -I*(exp(I*x) - exp(-I*x))/2 + cos(y)/gamma(z)


def test_atoms():
    assert b21.atoms() == set()


def test_free_symbols_empty():
    assert b21.free_symbols == set()


def test_doit():
    assert b21.doit() == b21
    assert b21.doit(deep=False) == b21


def test_S():
    assert repr(S) == 'S'


def test_xreplace():
    assert b21.xreplace({b2: b1}) == Basic(b1, b1)
    assert b21.xreplace({b2: b21}) == Basic(b21, b1)
    assert b3.xreplace({b2: b1}) == b2
    assert Basic(b1, b2).xreplace({b1: b2, b2: b1}) == Basic(b2, b1)
    assert Atom(b1).xreplace({b1: b2}) == Atom(b1)
    assert Atom(b1).xreplace({Atom(b1): b2}) == b2
    pytest.raises(TypeError, b1.xreplace)
    pytest.raises(TypeError, lambda: b1.xreplace([b1, b2]))


def test_Singleton():
    global instantiated
    instantiated = 0

    class MyNewSingleton(Basic, metaclass=Singleton):
        def __new__(cls):
            global instantiated
            instantiated += 1
            return Basic.__new__(cls)

    assert instantiated == 0
    MyNewSingleton()  # force instantiation
    assert instantiated == 1
    assert MyNewSingleton() is not Basic()
    assert MyNewSingleton() is MyNewSingleton()
    assert S.MyNewSingleton is MyNewSingleton()
    assert instantiated == 1

    class MySingletonSub(MyNewSingleton):
        pass
    assert instantiated == 1
    MySingletonSub()
    assert instantiated == 2
    assert MySingletonSub() is not MyNewSingleton()
    assert MySingletonSub() is MySingletonSub()


def test_preorder_traversal():
    expr = Basic(b21, b3)
    assert list(
        preorder_traversal(expr)) == [expr, b21, b2, b1, b1, b3, b2, b1]
    assert list(preorder_traversal(('abc', ('d', 'ef')))) == [
        ('abc', ('d', 'ef')), 'abc', ('d', 'ef'), 'd', 'ef']

    result = []
    pt = preorder_traversal(expr)
    for i in pt:
        result.append(i)
        if i == b2:
            pt.skip()
    assert result == [expr, b21, b2, b1, b3, b2]

    expr = z + w*(x + y)
    assert list(preorder_traversal([expr], keys=default_sort_key)) == \
        [[w*(x + y) + z], w*(x + y) + z, z, w*(x + y), w, x + y, x, y]
    assert list(preorder_traversal((x + y)*z, keys=True)) == \
        [z*(x + y), z, x + y, x, y]


def test_sorted_args():
    assert b21._sorted_args == b21.args
    pytest.raises(AttributeError, lambda: x._sorted_args)


def test_call():
    # See the long history of this in issues sympy/sympy#5026 and sympy/sympy#5105.

    pytest.raises(TypeError, lambda: sin(x)({x: 1, sin(x): 2}))
    pytest.raises(TypeError, lambda: sin(x)(1))

    # No effect as there are no callables
    assert sin(x).rcall(1) == sin(x)
    assert (1 + sin(x)).rcall(1) == 1 + sin(x)

    # Effect in the pressence of callables
    l = Lambda(x, 2*x)
    assert (l + x).rcall(y) == 2*y + x
    assert (x**l).rcall(2) == x**4
    # TODO UndefinedFunction does not subclass Expr
    # f = Function('f')
    # assert (2*f)(x) == 2*f(x)


def test_literal_evalf_is_number_is_zero_is_comparable():
    f = Function('f')

    # the following should not be changed without a lot of dicussion
    # `foo.is_number` should be equivalent to `not foo.free_symbols`
    # it should not attempt anything fancy; see is_zero, is_constant
    # and equals for more rigorous tests.
    assert f(1).is_number is True
    i = Integral(0, (x, x, x))
    # expressions that are symbolically 0 can be difficult to prove
    # so in case there is some easy way to know if something is 0
    # it should appear in the is_zero property for that object;
    # if is_zero is true evalf should always be able to compute that
    # zero
    assert i.evalf() == 0
    assert i.is_zero
    assert i.is_number is False
    assert i.evalf(2, strict=False) == 0

    # issue sympy/sympy#10272
    n = sin(1)**2 + cos(1)**2 - 1
    assert n.is_comparable is not True
    assert n.evalf(2, strict=False).is_comparable is not True


def test_is_evaluated():
    e = Add(x, x, evaluate=False)
    assert e.is_evaluated is False
    assert e.doit().is_evaluated is True
