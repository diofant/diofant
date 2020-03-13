import pytest

from diofant import Integer, SympifyError
from diofant.core.operations import AssocOp, LatticeOp


__all__ = ()


class MyMul(AssocOp):
    identity = Integer(1)


def test_flatten():
    assert MyMul(2, MyMul(4, 3)) == MyMul(2, 4, 3)


class Join(LatticeOp):
    """Simplest possible Lattice class."""

    zero = Integer(0)
    identity = Integer(1)


def test_lattice_simple():
    assert Join(Join(2, 3), 4) == Join(2, Join(3, 4))
    assert Join(2, 3) == Join(3, 2)
    assert Join(0, 2) == 0
    assert Join(1, 2) == 2
    assert Join(2, 2) == 2

    assert Join(Join(2, 3), 4) == Join(2, 3, 4)
    assert Join() == 1
    assert Join(4) == 4
    assert Join(1, 4, 2, 3, 1, 3, 2) == Join(2, 3, 4)


def test_lattice_shortcircuit():
    pytest.raises(SympifyError, lambda: Join(object))
    assert Join(0, object) == 0


def test_lattice_print():
    assert str(Join(5, 4, 3, 2)) == 'Join(2, 3, 4, 5)'


def test_lattice_make_args():
    assert Join.make_args(0) == {0}
    assert Join.make_args(1) == {1}
    assert Join.make_args(Join(2, 3, 4)) == {Integer(2), Integer(3), Integer(4)}
