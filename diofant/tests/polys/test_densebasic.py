"""Tests for dense recursive polynomials' basic tools."""

from diofant import ZZ
from diofant.polys.densebasic import dmp_ground


__all__ = ()


def test_dmp_ground():
    assert dmp_ground(ZZ(0), 0) == []
    assert dmp_ground(ZZ(0), 2) == [[[]]]
    assert dmp_ground(ZZ(0), 4) == [[[[[]]]]]

    assert dmp_ground(ZZ(7), -1) == ZZ(7)
    assert dmp_ground(ZZ(1), -1) == ZZ(1)

    assert dmp_ground(ZZ(7), 0) == [ZZ(7)]
    assert dmp_ground(ZZ(7), 2) == [[[ZZ(7)]]]

    assert dmp_ground(ZZ(3), 5) == [[[[[[ZZ(3)]]]]]]
