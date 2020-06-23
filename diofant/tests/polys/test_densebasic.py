"""Tests for dense recursive polynomials' basic tools."""

import pytest

from diofant import ZZ, oo, ring
from diofant.polys.densebasic import (dmp_convert, dmp_from_dict, dmp_ground,
                                      dmp_one_p, dmp_strip, dmp_to_dict,
                                      dmp_zero_p)
from diofant.polys.specialpolys import f_polys


__all__ = ()

f_0, f_1, f_2, f_3, f_4, f_5, f_6 = [f.to_dense() for f in f_polys()]


def test_dmp_degree_in():
    R, x = ring('x', ZZ)

    assert R.dmp_degree_in(R(0), 0) == -oo
    assert R.dmp_degree_in(R(1), 0) == 0
    assert R.dmp_degree_in(x, 0) == 1
    assert R.dmp_degree_in(x**4 + 1, 0) == 4
    assert R.dmp_degree_in(x**3 + 2*x**2 + 3, 0) == 3
    assert R.dmp_degree_in(x**3 + x**2 + 2*x, 0) == 3

    R, x, y = ring('x y', ZZ)

    assert R.dmp_degree_in(R(0), 0) == -oo
    assert R.dmp_degree_in(R(1), 0) == 0
    assert R.dmp_degree_in(2*x + 1, 0) == 1
    assert R.dmp_degree_in(2*x + y**2 + 2*y + 3, 0) == 1

    pytest.raises(IndexError, lambda: R.dmp_degree_in(R(1), -5))

    R, x, y, z = ring('x y z', ZZ)

    assert R.dmp_degree_in(R(0), 0) == -oo
    assert R.dmp_degree_in(R(0), 1) == -oo
    assert R.dmp_degree_in(R(0), 2) == -oo

    assert R.dmp_degree_in(R(1), 0) == 0
    assert R.dmp_degree_in(R(1), 1) == 0
    assert R.dmp_degree_in(R(1), 2) == 0

    f = R.from_list(f_4)

    assert R.dmp_degree_in(f, 0) == 9
    assert R.dmp_degree_in(f, 1) == 12
    assert R.dmp_degree_in(f, 2) == 8

    R, x, y, z, t = ring('x y z t', ZZ)

    f = R.from_list(f_6)

    assert R.dmp_degree_in(f, 0) == 4
    assert R.dmp_degree_in(f, 1) == 4
    assert R.dmp_degree_in(f, 2) == 6
    assert R.dmp_degree_in(f, 3) == 3


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


def test_dmp_one_p():
    assert dmp_one_p([1], 0, ZZ) is True
    assert dmp_one_p([[1]], 1, ZZ) is True
    assert dmp_one_p([[[1]]], 2, ZZ) is True
    assert dmp_one_p([[[12]]], 2, ZZ) is False
    assert dmp_one_p([[[1], [1]]], 2, ZZ) is False


def test_dmp_ground():
    assert dmp_ground(ZZ(0), 0) == []
    assert dmp_ground(ZZ(0), 2) == [[[]]]
    assert dmp_ground(ZZ(0), 4) == [[[[[]]]]]

    assert dmp_ground(ZZ(7), -1) == ZZ(7)
    assert dmp_ground(ZZ(7), 0) == [ZZ(7)]
    assert dmp_ground(ZZ(7), 2) == [[[ZZ(7)]]]


def test_dmp_from_to_dict():
    assert dmp_from_dict({}, 0, ZZ) == []

    assert dmp_to_dict([], 0) == {}

    f = [3, 0, 0, 2, 0, 0, 0, 0, 8]
    h = {(8,): 3, (5,): 2, (0,): 8}

    assert dmp_from_dict(h, 0, ZZ) == f

    assert dmp_to_dict(f, 0) == h

    R,  x, y = ring('x,y', ZZ)

    f = [R(3), R(0), R(2), R(0), R(0), R(8)]
    h = {(5,): R(3), (3,): R(2), (0,): R(8)}

    assert dmp_from_dict(h, 0, R) == f

    assert dmp_to_dict(f, 0) == h

    assert dmp_to_dict([1, 0, 5, 0, 7], 0) == {(0,): 7, (2,): 5, (4,): 1}

    assert dmp_from_dict({}, 1, ZZ) == [[]]
    assert dmp_to_dict([[]], 1) == {}

    f = [[3], [], [], [2], [], [], [], [], [8]]
    g = {(8, 0): 3, (5, 0): 2, (0, 0): 8}

    assert dmp_from_dict(g, 1, ZZ) == f
    assert dmp_to_dict(f, 1) == g
