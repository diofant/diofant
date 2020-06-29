"""Tests for dense recursive polynomials' basic tools."""

from diofant import ZZ
from diofant.polys.densebasic import (dmp_convert, dmp_ground, dmp_strip,
                                      dmp_zero_p)


__all__ = ()


def test_dmp_strip():
    assert dmp_strip([], 0) == []
    assert dmp_strip([0], 0) == []
    assert dmp_strip([0, 0, 0], 0) == []

    assert dmp_strip([1], 0) == [1]
    assert dmp_strip([0, 1], 0) == [1]
    assert dmp_strip([0, 0, 0, 1], 0) == [1]

    assert dmp_strip([1, 2, 0], 0) == [1, 2, 0]
    assert dmp_strip([0, 1, 2, 0], 0) == [1, 2, 0]
    assert dmp_strip([0, 0, 0, 1, 2, 0], 0) == [1, 2, 0]

    assert dmp_strip([0, 1, 0], 0) == [1, 0]

    assert dmp_strip([0, 0, 1, 2, 3, 0], 0) == [1, 2, 3, 0]

    assert dmp_strip([0, 0, 0, 3, 0, 1], 0) == [3, 0, 1]

    assert dmp_strip([[]], 1) == [[]]
    assert dmp_strip([[], []], 1) == [[]]
    assert dmp_strip([[], [], []], 1) == [[]]

    assert dmp_strip([[[]]], 2) == [[[]]]
    assert dmp_strip([[[]], [[]]], 2) == [[[]]]
    assert dmp_strip([[[]], [[]], [[]]], 2) == [[[]]]

    assert dmp_strip([[[1]]], 2) == [[[1]]]
    assert dmp_strip([[[]], [[1]]], 2) == [[[1]]]
    assert dmp_strip([[[]], [[1]], [[]]], 2) == [[[1]], [[]]]


def test_dmp_convert():
    K0, K1 = ZZ.inject('x'), ZZ

    assert dmp_convert([K0(1), K0(2)], 0, K0, K1) == [ZZ(1), ZZ(2)]
    assert dmp_convert([K1(1), K1(2)], 0, K1, K0) == [K0(1), K0(2)]

    f = [K0(1), K0(2), K0(0), K0(3)]

    assert dmp_convert(f, 0, K0, K1) == [ZZ(1), ZZ(2), ZZ(0), ZZ(3)]

    f = [[K0(1)], [K0(2)], [], [K0(3)]]

    assert dmp_convert(f, 1, K0, K1) == [[ZZ(1)], [ZZ(2)], [], [ZZ(3)]]


def test_dmp_zero_p():
    assert dmp_zero_p([], 0) is True
    assert dmp_zero_p([[]], 1) is True

    assert dmp_zero_p([[[]]], 2) is True
    assert dmp_zero_p([[[1]]], 2) is False


def test_dmp_ground():
    assert dmp_ground(ZZ(0), 0) == []
    assert dmp_ground(ZZ(0), 2) == [[[]]]
    assert dmp_ground(ZZ(0), 4) == [[[[[]]]]]

    assert dmp_ground(ZZ(7), -1) == ZZ(7)
    assert dmp_ground(ZZ(7), 0) == [ZZ(7)]
    assert dmp_ground(ZZ(7), 2) == [[[ZZ(7)]]]
