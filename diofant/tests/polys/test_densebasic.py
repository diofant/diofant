"""Tests for dense recursive polynomials' basic tools."""

from diofant import ZZ
from diofant.polys.densebasic import dmp_convert, dmp_ground


__all__ = ()


def test_dmp_convert():
    K0, K1 = ZZ.inject('x'), ZZ

    assert dmp_convert([K0(1), K0(2)], 0, K0, K1) == [ZZ(1), ZZ(2)]
    assert dmp_convert([K1(1), K1(2)], 0, K1, K0) == [K0(1), K0(2)]

    f = [K0(1), K0(2), K0(0), K0(3)]

    assert dmp_convert(f, 0, K0, K1) == [ZZ(1), ZZ(2), ZZ(0), ZZ(3)]

    f = [K0(0), K0(1)]

    assert dmp_convert(f, 0, K0, K1) == [K1(1)]

    f = [[K0(1)], [K0(2)], [], [K0(3)]]

    assert dmp_convert(f, 1, K0, K1) == [[ZZ(1)], [ZZ(2)], [], [ZZ(3)]]

    f = [[], []]

    assert dmp_convert(f, 1, K0, K1) == [[]]


def test_dmp_ground():
    assert dmp_ground(ZZ(0), 0) == []
    assert dmp_ground(ZZ(0), 2) == [[[]]]
    assert dmp_ground(ZZ(0), 4) == [[[[[]]]]]

    assert dmp_ground(ZZ(7), -1) == ZZ(7)
    assert dmp_ground(ZZ(1), -1) == ZZ(1)

    assert dmp_ground(ZZ(7), 0) == [ZZ(7)]
    assert dmp_ground(ZZ(7), 2) == [[[ZZ(7)]]]

    assert dmp_ground(ZZ(3), 5) == [[[[[[ZZ(3)]]]]]]
