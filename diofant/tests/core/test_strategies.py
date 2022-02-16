from diofant import Add, Basic, Integer
from diofant.abc import x
from diofant.core.strategies import (arguments, flatten, glom, null_safe,
                                     operator, rm_id, sort, term, unpack)


__all__ = ()


def test_rm_id():
    rmzeros = rm_id(lambda x: x == 0)
    assert rmzeros(Basic(0, 1)) == Basic(1)
    assert rmzeros(Basic(0, 0)) == Basic(0)
    assert rmzeros(Basic(2, 1)) == Basic(2, 1)


def test_glom():
    def key(x):
        return x.as_coeff_Mul()[1]

    def count(x):
        return x.as_coeff_Mul()[0]

    def newargs(cnt, arg):
        return cnt * arg

    rl = glom(key, count, newargs)

    result = rl(Add(x, -x, 3*x, 2, 3, evaluate=False))
    expected = Add(3*x, 5)
    assert set(result.args) == set(expected.args)

    result = rl(Add(*expected.args, evaluate=False))
    assert set(result.args) == set(expected.args)


def test_flatten():
    assert flatten(Basic(1, 2, Basic(3, 4))) == Basic(1, 2, 3, 4)


def test_unpack():
    assert unpack(Basic(2)) == 2
    assert unpack(Basic(2, 3)) == Basic(2, 3)


def test_sort():
    assert sort(str)(Basic(3, 1, 2)) == Basic(1, 2, 3)


def test_term():
    assert arguments(2) == ()
    assert not arguments(Integer(2))
    assert arguments(2 + x) == (2, x)
    assert operator(2 + x) == Add
    assert operator(Integer(2)) == Integer(2)
    assert term(Add, (2, x)) == 2 + x
    assert term(Integer(2), ()) == Integer(2)


def test_null_safe():
    def rl(expr):
        if expr == 1:
            return 2
    safe_rl = null_safe(rl)
    assert rl(1) == safe_rl(1)

    assert rl(3) is None
    assert safe_rl(3) == 3
