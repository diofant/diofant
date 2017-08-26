import pytest

from diofant import Integer
from diofant.core.operations import AssocOp, LatticeOp
from diofant.core.sympify import SympifyError


__all__ = ()


class MyMul(AssocOp):
    identity = Integer(1)


def test_flatten():
    assert MyMul(2, MyMul(4, 3)) == MyMul(2, 4, 3)


class join(LatticeOp):
    """Simplest possible Lattice class."""

    zero = Integer(0)
    identity = Integer(1)


def test_lattice_simple():
    assert join(join(2, 3), 4) == join(2, join(3, 4))
    assert join(2, 3) == join(3, 2)
    assert join(0, 2) == 0
    assert join(1, 2) == 2
    assert join(2, 2) == 2

    assert join(join(2, 3), 4) == join(2, 3, 4)
    assert join() == 1
    assert join(4) == 4
    assert join(1, 4, 2, 3, 1, 3, 2) == join(2, 3, 4)


def test_lattice_shortcircuit():
    pytest.raises(SympifyError, lambda: join(object))
    assert join(0, object) == 0


def test_lattice_print():
    assert str(join(5, 4, 3, 2)) == 'join(2, 3, 4, 5)'


def test_lattice_make_args():
    assert join.make_args(0) == {0}
    assert join.make_args(1) == {1}
    assert join.make_args(join(2, 3, 4)) == {Integer(2), Integer(3), Integer(4)}
