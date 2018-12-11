"""Tests for dense recursive polynomials' tools. """

import pytest

from diofant import I, Rational, sin, sqrt
from diofant.abc import x
from diofant.domains import EX, FF, QQ, ZZ
from diofant.polys.densearith import dmp_mul_ground
from diofant.polys.densebasic import (dmp_convert, dmp_normal, dmp_swap,
                                      dup_from_dict)
from diofant.polys.densetools import (dmp_clear_denoms, dmp_compose,
                                      dmp_diff_eval_in, dmp_diff_in,
                                      dmp_eval_in, dmp_eval_tail,
                                      dmp_ground_content, dmp_ground_extract,
                                      dmp_ground_monic, dmp_ground_primitive,
                                      dmp_ground_trunc, dmp_integrate_in,
                                      dmp_lift, dmp_trunc, dup_decompose,
                                      dup_mirror, dup_real_imag, dup_scale,
                                      dup_shift, dup_sign_variations,
                                      dup_transform, dup_trunc)
from diofant.polys.polyerrors import DomainError, ExactQuotientFailed
from diofant.polys.rings import ring
from diofant.polys.specialpolys import f_polys


__all__ = ()

f_0, f_1, f_2, f_3, f_4, f_5, f_6 = [ f.to_dense() for f in f_polys() ]


def test_dmp_integrate_in():
    assert dmp_integrate_in([], 1, 0, 0, QQ) == []
    assert dmp_integrate_in([], 2, 0, 0, QQ) == []

    assert dmp_integrate_in([QQ(1)], 1, 0, 0, QQ) == [QQ(1), QQ(0)]
    assert dmp_integrate_in([QQ(1)], 2, 0, 0, QQ) == [QQ(1, 2), QQ(0), QQ(0)]

    assert dmp_integrate_in([QQ(1), QQ(2), QQ(3)], 0, 0, 0, QQ) == \
        [QQ(1), QQ(2), QQ(3)]
    assert dmp_integrate_in([QQ(1), QQ(2), QQ(3)], 1, 0, 0, QQ) == \
        [QQ(1, 3), QQ(1), QQ(3), QQ(0)]
    assert dmp_integrate_in([QQ(1), QQ(2), QQ(3)], 2, 0, 0, QQ) == \
        [QQ(1, 12), QQ(1, 3), QQ(3, 2), QQ(0), QQ(0)]
    assert dmp_integrate_in([QQ(1), QQ(2), QQ(3)], 3, 0, 0, QQ) == \
        [QQ(1, 60), QQ(1, 12), QQ(1, 2), QQ(0), QQ(0), QQ(0)]

    assert dmp_integrate_in([QQ(1), QQ(2), QQ(0)], 1, 0, 0, QQ) == \
        [QQ(1, 3), QQ(1), QQ(0), QQ(0)]
    assert dmp_integrate_in([QQ(1), QQ(2), QQ(0)], 2, 0, 0, QQ) == \
        [QQ(1, 12), QQ(1, 3), QQ(0), QQ(0), QQ(0)]

    assert dmp_integrate_in(dup_from_dict({(29,): QQ(17)}, QQ), 3, 0, 0, QQ) == \
        dup_from_dict({(32,): QQ(17, 29760)}, QQ)

    assert dmp_integrate_in(dup_from_dict({(29,): QQ(17), (5,): QQ(1, 2)}, QQ), 3, 0, 0, QQ) == \
        dup_from_dict({(32,): QQ(17, 29760), (8,): QQ(1, 672)}, QQ)

    assert dmp_integrate_in([[[]]], 1, 0, 2, QQ) == [[[]]]
    assert dmp_integrate_in([[[]]], 2, 0, 2, QQ) == [[[]]]

    assert dmp_integrate_in([[[QQ(1)]]], 1, 0, 2, QQ) == [[[QQ(1)]], [[]]]
    assert dmp_integrate_in([[[QQ(1)]]], 2, 0, 2, QQ) == [[[QQ(1, 2)]], [[]], [[]]]

    assert dmp_integrate_in([[QQ(1)], [QQ(2)], [QQ(3)]], 0, 0, 1, QQ) == \
        [[QQ(1)], [QQ(2)], [QQ(3)]]
    assert dmp_integrate_in([[QQ(1)], [QQ(2)], [QQ(3)]], 1, 0, 1, QQ) == \
        [[QQ(1, 3)], [QQ(1)], [QQ(3)], []]
    assert dmp_integrate_in([[QQ(1)], [QQ(2)], [QQ(3)]], 2, 0, 1, QQ) == \
        [[QQ(1, 12)], [QQ(1, 3)], [QQ(3, 2)], [], []]
    assert dmp_integrate_in([[QQ(1)], [QQ(2)], [QQ(3)]], 3, 0, 1, QQ) == \
        [[QQ(1, 60)], [QQ(1, 12)], [QQ(1, 2)], [], [], []]

    assert dmp_integrate_in([[QQ(1)], [QQ(2), QQ(0)]], 1, 0, 1, QQ) == \
        [[QQ(1, 2)], [QQ(2), QQ(0)], []]
    assert dmp_integrate_in([[QQ(1)], [QQ(2), QQ(0)]], 2, 0, 1, QQ) == \
        [[QQ(1, 6)], [QQ(1), QQ(0)], [], []]

    f = dmp_convert(f_6, 3, ZZ, QQ)

    assert dmp_integrate_in(f, 2, 1, 3, QQ) == \
        dmp_swap(dmp_integrate_in(dmp_swap(f, 0, 1, 3, QQ), 2, 0, 3, QQ), 0, 1, 3, QQ)
    assert dmp_integrate_in(f, 3, 1, 3, QQ) == \
        dmp_swap(dmp_integrate_in(dmp_swap(f, 0, 1, 3, QQ), 3, 0, 3, QQ), 0, 1, 3, QQ)
    assert dmp_integrate_in(f, 2, 2, 3, QQ) == \
        dmp_swap(dmp_integrate_in(dmp_swap(f, 0, 2, 3, QQ), 2, 0, 3, QQ), 0, 2, 3, QQ)
    assert dmp_integrate_in(f, 3, 2, 3, QQ) == \
        dmp_swap(dmp_integrate_in(dmp_swap(f, 0, 2, 3, QQ), 3, 0, 3, QQ), 0, 2, 3, QQ)

    pytest.raises(IndexError, lambda: dmp_integrate_in(f, 2, -1, 3, QQ))
    pytest.raises(IndexError, lambda: dmp_integrate_in(f, 2, 1, 0, QQ))


def test_dmp_diff_in():
    assert dmp_diff_in([], 1, 0, 0, ZZ) == []
    assert dmp_diff_in([7], 1, 0, 0, ZZ) == []
    assert dmp_diff_in([2, 7], 1, 0, 0, ZZ) == [2]
    assert dmp_diff_in([1, 2, 1], 1, 0, 0, ZZ) == [2, 2]
    assert dmp_diff_in([1, 2, 3, 4], 1, 0, 0, ZZ) == [3, 4, 3]
    assert dmp_diff_in([1, -1, 0, 0, 2], 1, 0, 0, ZZ) == [4, -3, 0, 0]

    assert dmp_diff_in([1, 2, 3, 4], 1, 0, 0, ZZ) == [3, 4, 3]
    assert dmp_diff_in([1, 2, 3, 4], 2, 0, 0, ZZ) == [6, 4]

    f = dmp_normal([17, 34, 56, -345, 23, 76, 0, 0, 12, 3, 7], 0, ZZ)

    assert dmp_diff_in(f, 0, 0, 0, ZZ) == f
    assert dmp_diff_in(f, 2, 0, 0, ZZ) == dmp_diff_in(dmp_diff_in(f, 1, 0, 0, ZZ), 1, 0, 0, ZZ)
    assert dmp_diff_in(f, 3, 0, 0, ZZ) == dmp_diff_in(dmp_diff_in(dmp_diff_in(f, 1, 0, 0, ZZ),
                                                      1, 0, 0, ZZ), 1, 0, 0, ZZ)

    K = FF(3)
    f = dmp_normal([17, 34, 56, -345, 23, 76, 0, 0, 12, 3, 7], 0, K)

    assert dmp_diff_in(f, 1, 0, 0, K) == dmp_normal([2, 0, 1, 0, 0, 2, 0, 0, 0, 0], 0, K)
    assert dmp_diff_in(f, 2, 0, 0, K) == dmp_normal([1, 0, 0, 2, 0, 0, 0], 0, K)
    assert dmp_diff_in(f, 3, 0, 0, K) == dmp_normal([], 0, K)

    assert dmp_diff_in(f, 0, 0, 0, K) == f
    assert dmp_diff_in(f, 2, 0, 0, K) == dmp_diff_in(dmp_diff_in(f, 1, 0, 0, K), 1, 0, 0, K)
    assert dmp_diff_in(f, 3, 0, 0, K) == dmp_diff_in(dmp_diff_in(dmp_diff_in(f, 1, 0, 0, K),
                                                     1, 0, 0, K), 1, 0, 0, K)

    assert dmp_diff_in([], 1, 0, 0, ZZ) == []
    assert dmp_diff_in([[]], 1, 0, 1, ZZ) == [[]]
    assert dmp_diff_in([[[]]], 1, 0, 2, ZZ) == [[[]]]

    assert dmp_diff_in([[[1], [2]]], 1, 0, 2, ZZ) == [[[]]]

    assert dmp_diff_in([[[1]], [[]]], 1, 0, 2, ZZ) == [[[1]]]
    assert dmp_diff_in([[[3]], [[1]], [[]]], 1, 0, 2, ZZ) == [[[6]], [[1]]]

    assert dmp_diff_in([[1, 2, 3], [2, 3, 1]], 1, 0, 1, ZZ) == [[1, 2, 3]]
    assert dmp_diff_in([[1, 2, 3], [2, 3, 1]], 2, 0, 1, ZZ) == [[]]
    assert dmp_diff_in([[1, 2, 3], [2, 3, 1]], 1, 1, 1, ZZ) == [[2, 2], [4, 3]]

    assert dmp_diff_in(f_6, 0, 0, 3, ZZ) == f_6
    assert dmp_diff_in(f_6, 2, 0, 3, ZZ) == dmp_diff_in(dmp_diff_in(f_6, 1, 0, 3, ZZ), 1, 0, 3, ZZ)
    assert dmp_diff_in(f_6, 3, 0, 3, ZZ) == dmp_diff_in(
        dmp_diff_in(dmp_diff_in(f_6, 1, 0, 3, ZZ), 1, 0, 3, ZZ), 1, 0, 3, ZZ)

    f = [[1, 2, 3], [2, 3, 1]]
    assert dmp_diff_in(f, 1, 0, 1, ZZ) == [[1, 2, 3]]
    assert dmp_diff_in(f, 2, 0, 1, ZZ) == [[]]

    K = FF(23)
    F_6 = dmp_normal(f_6, 3, K)

    assert dmp_diff_in(F_6, 0, 0, 3, K) == F_6
    assert dmp_diff_in(F_6, 2, 0, 3, K) == dmp_diff_in(dmp_diff_in(F_6, 1, 0, 3, K), 1, 0, 3, K)
    assert dmp_diff_in(F_6, 3, 0, 3, K) == dmp_diff_in(
        dmp_diff_in(dmp_diff_in(F_6, 1, 0, 3, K), 1, 0, 3, K), 1, 0, 3, K)

    assert dmp_diff_in(f_6, 2, 1, 3, ZZ) == \
        dmp_swap(dmp_diff_in(dmp_swap(f_6, 0, 1, 3, ZZ), 2, 0, 3, ZZ), 0, 1, 3, ZZ)
    assert dmp_diff_in(f_6, 3, 1, 3, ZZ) == \
        dmp_swap(dmp_diff_in(dmp_swap(f_6, 0, 1, 3, ZZ), 3, 0, 3, ZZ), 0, 1, 3, ZZ)
    assert dmp_diff_in(f_6, 2, 2, 3, ZZ) == \
        dmp_swap(dmp_diff_in(dmp_swap(f_6, 0, 2, 3, ZZ), 2, 0, 3, ZZ), 0, 2, 3, ZZ)
    assert dmp_diff_in(f_6, 3, 2, 3, ZZ) == \
        dmp_swap(dmp_diff_in(dmp_swap(f_6, 0, 2, 3, ZZ), 3, 0, 3, ZZ), 0, 2, 3, ZZ)

    pytest.raises(IndexError, lambda: dmp_diff_in(f_6, 2, -1, 3, ZZ))
    pytest.raises(IndexError, lambda: dmp_diff_in(f_6, 2, 1, 0, ZZ))


def test_dmp_eval_in():
    assert dmp_eval_in([], 7, 0, 0, ZZ) == 0
    assert dmp_eval_in([1, 2], 0, 0, 0, ZZ) == 2
    assert dmp_eval_in([1, 2, 3], 7, 0, 0, ZZ) == 66
    assert dmp_eval_in([1, 2, 3], 2, 0, 0, ZZ) == 11

    assert dmp_eval_in([], 3, 0, 0, ZZ) == 0

    assert dmp_eval_in([[]], 3, 0, 1, ZZ) == []
    assert dmp_eval_in([[[]]], 3, 0, 2, ZZ) == [[]]

    assert dmp_eval_in([[1, 2]], 0, 0, 1, ZZ) == [1, 2]

    assert dmp_eval_in([[2, 3], [1, 2]], 2, 0, 1, ZZ) == [5, 8]

    assert dmp_eval_in([[[1]]], 3, 0, 2, ZZ) == [[1]]
    assert dmp_eval_in([[[1, 2]]], 3, 0, 2, ZZ) == [[1, 2]]

    assert dmp_eval_in([[3, 2], [1, 2]], 3, 0, 1, ZZ) == [10, 8]
    assert dmp_eval_in([[[3, 2]], [[1, 2]]], 3, 0, 2, ZZ) == [[10, 8]]

    assert dmp_eval_in(
        f_6, -2, 1, 3, ZZ) == dmp_eval_in(dmp_swap(f_6, 0, 1, 3, ZZ), -2, 0, 3, ZZ)
    assert dmp_eval_in(
        f_6, 7, 1, 3, ZZ) == dmp_eval_in(dmp_swap(f_6, 0, 1, 3, ZZ), 7, 0, 3, ZZ)
    assert dmp_eval_in(f_6, -2, 2, 3, ZZ) == dmp_swap(
        dmp_eval_in(dmp_swap(f_6, 0, 2, 3, ZZ), -2, 0, 3, ZZ), 0, 1, 2, ZZ)
    assert dmp_eval_in(f_6, 7, 2, 3, ZZ) == dmp_swap(
        dmp_eval_in(dmp_swap(f_6, 0, 2, 3, ZZ), 7, 0, 3, ZZ), 0, 1, 2, ZZ)

    f = [[[int(45)]], [[]], [[]], [[int(-9)], [-1], [],
                                   [int(3), int(0), int(10), int(0)]]]

    assert dmp_eval_in(f, -2, 2, 2, ZZ) == \
        [[45], [], [], [-9, -1, 0, -44]]

    pytest.raises(IndexError, lambda: dmp_eval_in(f, -2, -1, 2, ZZ))
    pytest.raises(IndexError, lambda: dmp_eval_in(f, -2, 2, 0, ZZ))


def test_dmp_eval_tail():
    assert dmp_eval_tail([[]], [1], 1, ZZ) == []
    assert dmp_eval_tail([[[]]], [1], 2, ZZ) == [[]]
    assert dmp_eval_tail([[[]]], [1, 2], 2, ZZ) == []

    assert dmp_eval_tail(f_0, [], 2, ZZ) == f_0

    assert dmp_eval_tail(f_0, [1, -17, 8], 2, ZZ) == 84496
    assert dmp_eval_tail(f_0, [-17, 8], 2, ZZ) == [-1409, 3, 85902]
    assert dmp_eval_tail(f_0, [8], 2, ZZ) == [[83, 2], [3], [302, 81, 1]]

    assert dmp_eval_tail(f_1, [-17, 8], 2, ZZ) == [-136, 15699, 9166, -27144]

    assert dmp_eval_tail(
        f_2, [-12, 3], 2, ZZ) == [-1377, 0, -702, -1224, 0, -624]
    assert dmp_eval_tail(
        f_3, [-12, 3], 2, ZZ) == [144, 82, -5181, -28872, -14868, -540]

    assert dmp_eval_tail(
        f_4, [25, -1], 2, ZZ) == [152587890625, 9765625, -59605407714843750,
                                  -3839159765625, -1562475, 9536712644531250, 610349546750, -4, 24414375000, 1562520]
    assert dmp_eval_tail(f_5, [25, -1], 2, ZZ) == [-1, -78, -2028, -17576]

    assert dmp_eval_tail(f_6, [0, 2, 4], 3, ZZ) == [5040, 0, 0, 4480]


def test_dmp_diff_eval_in():
    assert dmp_diff_eval_in(f_6, 2, 7, 1, 3, ZZ) == \
        dmp_eval_in(dmp_diff_in(dmp_swap(f_6, 0, 1, 3, ZZ), 2, 0, 3, ZZ), 7, 0, 3, ZZ)

    pytest.raises(IndexError, lambda: dmp_diff_eval_in(f_6, 2, 7, 4, 3, ZZ))


def test_dup_trunc():
    assert dup_trunc([1, 2, 3, 4, 5, 6], ZZ(3), ZZ) == [1, -1, 0, 1, -1, 0]
    assert dup_trunc([6, 5, 4, 3, 2, 1], ZZ(3), ZZ) == [-1, 1, 0, -1, 1]


def test_dmp_trunc():
    assert dmp_trunc([[]], [1, 2], 2, ZZ) == [[]]
    assert dmp_trunc([[1, 2], [1, 4, 1], [1]], [1, 2], 1, ZZ) == [[-3], [1]]


def test_dmp_ground_trunc():
    assert dmp_ground_trunc(f_0, ZZ(3), 2, ZZ) == \
        dmp_normal(
            [[[1, -1, 0], [-1]], [[]], [[1, -1, 0], [1, -1, 1], [1]]], 2, ZZ)


def test_dmp_ground_monic():
    assert dmp_ground_monic([ZZ(3), ZZ(6), ZZ(9)], 0, ZZ) == [1, 2, 3]
    assert dmp_ground_monic([QQ(3), QQ(4), QQ(2)], 0, QQ) == [1, QQ(4, 3), QQ(2, 3)]

    pytest.raises(ExactQuotientFailed, lambda: dmp_ground_monic([3, 4, 5], 0, ZZ))

    assert dmp_ground_monic([], 0, QQ) == []
    assert dmp_ground_monic([QQ(1)], 0, QQ) == [1]
    assert dmp_ground_monic([QQ(7), QQ(1), QQ(21)], 0, QQ) == [1, QQ(1, 7), 3]

    assert dmp_ground_monic([[ZZ(3)], [ZZ(6)], [ZZ(9)]], 1, ZZ) == [[1], [2], [3]]

    pytest.raises(ExactQuotientFailed, lambda: dmp_ground_monic([[3], [4], [5]], 1, ZZ))

    assert dmp_ground_monic([[]], 1, QQ) == [[]]
    assert dmp_ground_monic([[QQ(1)]], 1, QQ) == [[1]]
    assert dmp_ground_monic([[QQ(7)], [QQ(1)], [QQ(21)]], 1, QQ) == [[1], [QQ(1, 7)], [3]]


def test_dmp_ground_content():
    assert dmp_ground_content([], 0, ZZ) == 0
    assert dmp_ground_content([ZZ(+1)], 0, ZZ) == 1
    assert dmp_ground_content([ZZ(-1)], 0, ZZ) == 1
    assert dmp_ground_content([ZZ(1), ZZ(1)], 0, ZZ) == 1
    assert dmp_ground_content([ZZ(2), ZZ(2)], 0, ZZ) == 2
    assert dmp_ground_content([ZZ(1), ZZ(2), ZZ(1)], 0, ZZ) == 1
    assert dmp_ground_content([ZZ(2), ZZ(4), ZZ(2)], 0, ZZ) == 2
    assert dmp_ground_content([ZZ(6), ZZ(8), ZZ(12)], 0, ZZ) == 2
    assert dmp_ground_content([QQ(6), QQ(8), QQ(12)], 0, QQ) == 2

    assert dmp_ground_content([QQ(2, 3), QQ(4, 9)], 0, QQ) == QQ(2, 9)
    assert dmp_ground_content([QQ(2, 3), QQ(4, 5)], 0, QQ) == QQ(2, 15)

    assert dmp_ground_content([[]], 1, ZZ) == 0
    assert dmp_ground_content([[]], 1, QQ) == 0
    assert dmp_ground_content([[ZZ(+1)]], 1, ZZ) == 1
    assert dmp_ground_content([[ZZ(-1)]], 1, ZZ) == 1
    assert dmp_ground_content([[ZZ(1)], [ZZ(1)]], 1, ZZ) == 1
    assert dmp_ground_content([[ZZ(2)], [ZZ(2)]], 1, ZZ) == 2
    assert dmp_ground_content([[ZZ(1)], [ZZ(2)], [ZZ(1)]], 1, ZZ) == 1
    assert dmp_ground_content([[ZZ(2)], [ZZ(4)], [ZZ(2)]], 1, ZZ) == 2

    assert dmp_ground_content([[QQ(2, 3)], [QQ(4, 9)]], 1, QQ) == QQ(2, 9)
    assert dmp_ground_content([[QQ(2, 3)], [QQ(4, 5)]], 1, QQ) == QQ(2, 15)

    assert dmp_ground_content(f_0, 2, ZZ) == 1
    assert dmp_ground_content(dmp_mul_ground(f_0, ZZ(2), 2, ZZ), 2, ZZ) == 2

    assert dmp_ground_content(f_1, 2, ZZ) == 1
    assert dmp_ground_content(dmp_mul_ground(f_1, ZZ(3), 2, ZZ), 2, ZZ) == 3

    assert dmp_ground_content(f_2, 2, ZZ) == 1
    assert dmp_ground_content(dmp_mul_ground(f_2, ZZ(4), 2, ZZ), 2, ZZ) == 4

    assert dmp_ground_content(f_3, 2, ZZ) == 1
    assert dmp_ground_content(dmp_mul_ground(f_3, ZZ(5), 2, ZZ), 2, ZZ) == 5

    assert dmp_ground_content(f_4, 2, ZZ) == 1
    assert dmp_ground_content(dmp_mul_ground(f_4, ZZ(6), 2, ZZ), 2, ZZ) == 6

    assert dmp_ground_content(f_5, 2, ZZ) == 1
    assert dmp_ground_content(dmp_mul_ground(f_5, ZZ(7), 2, ZZ), 2, ZZ) == 7

    assert dmp_ground_content(f_6, 3, ZZ) == 1
    assert dmp_ground_content(dmp_mul_ground(f_6, ZZ(8), 3, ZZ), 3, ZZ) == 8


def test_dmp_ground_primitive():
    assert dmp_ground_primitive([], 0, ZZ) == (0, [])
    assert dmp_ground_primitive([ZZ(1)], 0, ZZ) == (1, [1])
    assert dmp_ground_primitive([ZZ(1), ZZ(1)], 0, ZZ) == (1, [1, 1])
    assert dmp_ground_primitive([ZZ(2), ZZ(2)], 0, ZZ) == (2, [1, 1])
    assert dmp_ground_primitive([ZZ(1), ZZ(2), ZZ(1)], 0, ZZ) == (1, [1, 2, 1])
    assert dmp_ground_primitive([ZZ(2), ZZ(4), ZZ(2)], 0, ZZ) == (2, [1, 2, 1])
    assert dmp_ground_primitive([ZZ(6), ZZ(8), ZZ(12)], 0, ZZ) == (2, [3, 4, 6])

    assert dmp_ground_primitive([], 0, QQ) == (0, [])
    assert dmp_ground_primitive([QQ(1)], 0, QQ) == (1, [1])
    assert dmp_ground_primitive([QQ(1), QQ(1)], 0, QQ) == (1, [1, 1])
    assert dmp_ground_primitive([QQ(2), QQ(2)], 0, QQ) == (2, [1, 1])
    assert dmp_ground_primitive([QQ(1), QQ(2), QQ(1)], 0, QQ) == (1, [1, 2, 1])
    assert dmp_ground_primitive([QQ(2), QQ(4), QQ(2)], 0, QQ) == (2, [1, 2, 1])
    assert dmp_ground_primitive([QQ(6), QQ(8), QQ(12)], 0, QQ) == (2, [3, 4, 6])

    assert dmp_ground_primitive([QQ(2, 3), QQ(4, 9)], 0, QQ) == (QQ(2, 9), [3, 2])
    assert dmp_ground_primitive([QQ(2, 3), QQ(4, 5)], 0, QQ) == (QQ(2, 15), [5, 6])

    assert dmp_ground_primitive([[]], 1, ZZ) == (0, [[]])

    assert dmp_ground_primitive(f_0, 2, ZZ) == (1, f_0)
    assert dmp_ground_primitive(dmp_mul_ground(f_0, ZZ(2), 2, ZZ), 2, ZZ) == (2, f_0)

    assert dmp_ground_primitive(f_1, 2, ZZ) == (1, f_1)
    assert dmp_ground_primitive(dmp_mul_ground(f_1, ZZ(3), 2, ZZ), 2, ZZ) == (3, f_1)

    assert dmp_ground_primitive(f_2, 2, ZZ) == (1, f_2)
    assert dmp_ground_primitive(dmp_mul_ground(f_2, ZZ(4), 2, ZZ), 2, ZZ) == (4, f_2)

    assert dmp_ground_primitive(f_3, 2, ZZ) == (1, f_3)
    assert dmp_ground_primitive(dmp_mul_ground(f_3, ZZ(5), 2, ZZ), 2, ZZ) == (5, f_3)

    assert dmp_ground_primitive(f_4, 2, ZZ) == (1, f_4)
    assert dmp_ground_primitive(dmp_mul_ground(f_4, ZZ(6), 2, ZZ), 2, ZZ) == (6, f_4)

    assert dmp_ground_primitive(f_5, 2, ZZ) == (1, f_5)
    assert dmp_ground_primitive(dmp_mul_ground(f_5, ZZ(7), 2, ZZ), 2, ZZ) == (7, f_5)

    assert dmp_ground_primitive(f_6, 3, ZZ) == (1, f_6)
    assert dmp_ground_primitive(dmp_mul_ground(f_6, ZZ(8), 3, ZZ), 3, ZZ) == (8, f_6)

    assert dmp_ground_primitive([[ZZ(2)]], 1, ZZ) == (2, [[1]])
    assert dmp_ground_primitive([[QQ(2)]], 1, QQ) == (2, [[1]])

    assert dmp_ground_primitive([[QQ(2, 3)], [QQ(4, 9)]], 1, QQ) == (QQ(2, 9), [[3], [2]])
    assert dmp_ground_primitive([[QQ(2, 3)], [QQ(4, 5)]], 1, QQ) == (QQ(2, 15), [[5], [6]])


def test_dmp_ground_extract():
    f = dmp_normal([2930944, 0, 2198208, 0, 549552, 0, 45796], 0, ZZ)
    g = dmp_normal([17585664, 0, 8792832, 0, 1099104, 0], 0, ZZ)

    F = dmp_normal([64, 0, 48, 0, 12, 0, 1], 0, ZZ)
    G = dmp_normal([384, 0, 192, 0, 24, 0], 0, ZZ)

    assert dmp_ground_extract(f, g, 0, ZZ) == (45796, F, G)

    f, g = [ZZ(6), ZZ(12), ZZ(18)], [ZZ(4), ZZ(8), ZZ(12)]

    assert dmp_ground_extract(f, g, 0, ZZ) == (2, [3, 6, 9], [2, 4, 6])

    f = dmp_normal([[2930944], [], [2198208], [], [549552], [], [45796]], 1, ZZ)
    g = dmp_normal([[17585664], [], [8792832], [], [1099104], []], 1, ZZ)

    F = dmp_normal([[64], [], [48], [], [12], [], [1]], 1, ZZ)
    G = dmp_normal([[384], [], [192], [], [24], []], 1, ZZ)

    assert dmp_ground_extract(f, g, 1, ZZ) == (45796, F, G)

    f, g = dmp_normal([1, 2], 0, ZZ), dmp_normal([3, 4], 0, ZZ)

    assert dmp_ground_extract(f, g, 0, ZZ) == (1, f, g)


def test_dup_real_imag():
    assert dup_real_imag([], ZZ) == ([[]], [[]])
    assert dup_real_imag([ZZ(1)], ZZ) == ([[1]], [[]])

    assert dup_real_imag([ZZ(1), ZZ(1)], ZZ) == ([[1], [1]], [[1, 0]])
    assert dup_real_imag([ZZ(1), ZZ(2)], ZZ) == ([[1], [2]], [[1, 0]])

    assert dup_real_imag([ZZ(1), ZZ(2), ZZ(3)], ZZ) == ([[1], [2], [-1, 0, 3]],
                                                        [[2, 0], [2, 0]])

    f = dmp_normal([1, 1, 1, 1], 0, ZZ)

    assert dup_real_imag(f, ZZ) == ([[1], [1], [-3, 0, 1], [-1, 0, 1]],
                                    [[3, 0], [2, 0], [-1, 0, 1, 0]])

    f = dmp_normal([1, 1], 0, EX)

    pytest.raises(DomainError, lambda: dup_real_imag(f, EX))

    A = QQ.algebraic_field(I)
    f = [A(1), A(I), A(0), A(-1), A(1)]

    assert dup_real_imag(f, A) == ([[1], [], [-6, -3, 0], [-1],
                                    [1, 1, 0, 0, 1]],
                                   [[4, 1], [], [-4, -3, 0, 0], [-1, 0]])

    A = QQ.algebraic_field(sqrt(2))
    f = [A(1), A(sqrt(2)), A(-1)]

    assert dup_real_imag(f, A) == ([[1], [A.unit], [-1, 0, -1]],
                                   [[2, 0], [A.unit, 0]])

    A2 = A.algebraic_field(I)
    f = [A2(1), A2(2*sqrt(2)*I), A2(I - 1)]

    assert dup_real_imag(f, A2) == ([[1], [], [-1, -2*A.unit, -1]],
                                    [[2, 2*A.unit], [1]])


def test_dup_mirror():
    assert dup_mirror([], ZZ) == []
    assert dup_mirror([1], ZZ) == [1]

    assert dup_mirror([1, 2, 3, 4, 5], ZZ) == [1, -2, 3, -4, 5]
    assert dup_mirror([1, 2, 3, 4, 5, 6], ZZ) == [-1, 2, -3, 4, -5, 6]


def test_dup_scale():
    assert dup_scale([], -1, ZZ) == []
    assert dup_scale([1], -1, ZZ) == [1]

    assert dup_scale([1, 2, 3, 4, 5], -1, ZZ) == [1, -2, 3, -4, 5]
    assert dup_scale([1, 2, 3, 4, 5], -7, ZZ) == [2401, -686, 147, -28, 5]


def test_dup_shift():
    assert dup_shift([], 1, ZZ) == []
    assert dup_shift([1], 1, ZZ) == [1]

    assert dup_shift([1, 2, 3, 4, 5], 1, ZZ) == [1, 6, 15, 20, 15]
    assert dup_shift([1, 2, 3, 4, 5], 7, ZZ) == [1, 30, 339, 1712, 3267]


def test_dup_transform():
    assert dup_transform([], [], [1, 1], ZZ) == []
    assert dup_transform([], [1], [1, 1], ZZ) == []
    assert dup_transform([], [1, 2], [1, 1], ZZ) == []

    assert dup_transform([6, -5, 4, -3, 17], [1, -3, 4], [2, -3], ZZ) == \
        [6, -82, 541, -2205, 6277, -12723, 17191, -13603, 4773]


def test_dmp_compose():
    assert dmp_compose([], [], 0, ZZ) == []
    assert dmp_compose([], [1], 0, ZZ) == []
    assert dmp_compose([], [1, 2], 0, ZZ) == []

    assert dmp_compose([1], [], 0, ZZ) == [1]

    assert dmp_compose([1, 2, 0], [], 0, ZZ) == []
    assert dmp_compose([1, 2, 1], [], 0, ZZ) == [1]

    assert dmp_compose([1, 2, 1], [1], 0, ZZ) == [4]
    assert dmp_compose([1, 2, 1], [7], 0, ZZ) == [64]

    assert dmp_compose([1, 2, 1], [1, -1], 0, ZZ) == [1, 0, 0]
    assert dmp_compose([1, 2, 1], [1, 1], 0, ZZ) == [1, 4, 4]
    assert dmp_compose([1, 2, 1], [1, 2, 1], 0, ZZ) == [1, 4, 8, 8, 4]

    assert dmp_compose([1, 2, 1], [1, 2, 1], 0, ZZ) == [1, 4, 8, 8, 4]

    assert dmp_compose([1, 1, 0], [1, -1], 0, ZZ) == [1, -1, 0]

    assert dmp_compose([[[]]], [[[]]], 2, ZZ) == [[[]]]
    assert dmp_compose([[[]]], [[[1]]], 2, ZZ) == [[[]]]
    assert dmp_compose([[[]]], [[[1]], [[2]]], 2, ZZ) == [[[]]]

    assert dmp_compose([[[1]]], [], 2, ZZ) == [[[1]]]

    assert dmp_compose([[1], [2], [ ]], [[]], 1, ZZ) == [[]]
    assert dmp_compose([[1], [2], [1]], [[]], 1, ZZ) == [[1]]

    assert dmp_compose([[1], [2], [1]], [[1]], 1, ZZ) == [[4]]
    assert dmp_compose([[1], [2], [1]], [[7]], 1, ZZ) == [[64]]

    assert dmp_compose([[1], [2], [1]], [[1], [-1]], 1, ZZ) == [[1], [ ], [ ]]
    assert dmp_compose([[1], [2], [1]], [[1], [ 1]], 1, ZZ) == [[1], [4], [4]]

    assert dmp_compose(
        [[1], [2], [1]], [[1], [2], [1]], 1, ZZ) == [[1], [4], [8], [8], [4]]


def test_dup_decompose():
    assert dup_decompose([1], ZZ) == [[1]]

    assert dup_decompose([1, 0], ZZ) == [[1, 0]]
    assert dup_decompose([1, 0, 0, 0], ZZ) == [[1, 0, 0, 0]]

    assert dup_decompose([1, 0, 0, 0, 0], ZZ) == [[1, 0, 0], [1, 0, 0]]
    assert dup_decompose(
        [1, 0, 0, 0, 0, 0, 0], ZZ) == [[1, 0, 0, 0], [1, 0, 0]]

    assert dup_decompose([7, 0, 0, 0, 1], ZZ) == [[7, 0, 1], [1, 0, 0]]
    assert dup_decompose([4, 0, 3, 0, 2], ZZ) == [[4, 3, 2], [1, 0, 0]]

    f = [1, 0, 20, 0, 150, 0, 500, 0, 625, -2, 0, -10, 9]

    assert dup_decompose(f, ZZ) == [[1, 0, 0, -2, 9], [1, 0, 5, 0]]

    f = [2, 0, 40, 0, 300, 0, 1000, 0, 1250, -4, 0, -20, 18]

    assert dup_decompose(f, ZZ) == [[2, 0, 0, -4, 18], [1, 0, 5, 0]]

    f = [1, 0, 20, -8, 150, -120, 524, -600, 865, -1034, 600, -170, 29]

    assert dup_decompose(f, ZZ) == [[1, -8, 24, -34, 29], [1, 0, 5, 0]]

    R, t = ring("t", ZZ)
    f = [6*t**2 - 42,
         48*t**2 + 96,
         144*t**2 + 648*t + 288,
         624*t**2 + 864*t + 384,
         108*t**3 + 312*t**2 + 432*t + 192]

    assert dup_decompose(f, R) == [f]


def test_dmp_lift():
    A = QQ.algebraic_field(I)
    f = [A(1), A(0), A(0), A(I), A(17*I)]

    assert dmp_lift(f, 0, A) == [1, 0, 0, 0, 0, 0, 2, 0, 578, 0, 0, 0,
                                 1, 0, -578, 0, 83521]

    pytest.raises(DomainError, lambda: dmp_lift([EX(1), EX(2)], 0, EX))


def test_dup_sign_variations():
    assert dup_sign_variations([], ZZ) == 0
    assert dup_sign_variations([1, 0], ZZ) == 0
    assert dup_sign_variations([1, 0, 2], ZZ) == 0
    assert dup_sign_variations([1, 0, 3, 0], ZZ) == 0
    assert dup_sign_variations([1, 0, 4, 0, 5], ZZ) == 0

    assert dup_sign_variations([-1, 0, 2], ZZ) == 1
    assert dup_sign_variations([-1, 0, 3, 0], ZZ) == 1
    assert dup_sign_variations([-1, 0, 4, 0, 5], ZZ) == 1

    assert dup_sign_variations([-1, -4, -5], ZZ) == 0
    assert dup_sign_variations([ 1, -4, -5], ZZ) == 1
    assert dup_sign_variations([ 1, 4, -5], ZZ) == 1
    assert dup_sign_variations([ 1, -4, 5], ZZ) == 2
    assert dup_sign_variations([-1, 4, -5], ZZ) == 2
    assert dup_sign_variations([-1, 4, 5], ZZ) == 1
    assert dup_sign_variations([-1, -4, 5], ZZ) == 1
    assert dup_sign_variations([ 1, 4, 5], ZZ) == 0

    assert dup_sign_variations([-1, 0, -4, 0, -5], ZZ) == 0
    assert dup_sign_variations([ 1, 0, -4, 0, -5], ZZ) == 1
    assert dup_sign_variations([ 1, 0, 4, 0, -5], ZZ) == 1
    assert dup_sign_variations([ 1, 0, -4, 0, 5], ZZ) == 2
    assert dup_sign_variations([-1, 0, 4, 0, -5], ZZ) == 2
    assert dup_sign_variations([-1, 0, 4, 0, 5], ZZ) == 1
    assert dup_sign_variations([-1, 0, -4, 0, 5], ZZ) == 1
    assert dup_sign_variations([ 1, 0, 4, 0, 5], ZZ) == 0


def test_dmp_clear_denoms():
    assert dmp_clear_denoms([], 0, QQ, ZZ) == (ZZ(1), [])

    assert dmp_clear_denoms([QQ(1)], 0, QQ, ZZ) == (ZZ(1), [QQ(1)])
    assert dmp_clear_denoms([QQ(7)], 0, QQ, ZZ) == (ZZ(1), [QQ(7)])

    assert dmp_clear_denoms([QQ(7, 3)], 0, QQ) == (ZZ(3), [QQ(7)])
    assert dmp_clear_denoms([QQ(7, 3)], 0, QQ, ZZ) == (ZZ(3), [QQ(7)])

    assert dmp_clear_denoms(
        [QQ(3), QQ(1), QQ(0)], 0, QQ, ZZ) == (ZZ(1), [QQ(3), QQ(1), QQ(0)])
    assert dmp_clear_denoms(
        [QQ(1), QQ(1, 2), QQ(0)], 0, QQ, ZZ) == (ZZ(2), [QQ(2), QQ(1), QQ(0)])

    assert dmp_clear_denoms([QQ(3), QQ(
        1), QQ(0)], 0, QQ, ZZ, convert=True) == (ZZ(1), [ZZ(3), ZZ(1), ZZ(0)])
    assert dmp_clear_denoms([QQ(1), QQ(
        1, 2), QQ(0)], 0, QQ, ZZ, convert=True) == (ZZ(2), [ZZ(2), ZZ(1), ZZ(0)])

    assert dmp_clear_denoms([QQ(1, 2), QQ(1, 3)], 0, QQ) == (QQ(6), [QQ(3), QQ(2)])
    assert dmp_clear_denoms([QQ(1, 2), QQ(1, 3)], 0, QQ,
                            convert=True) == (ZZ(6), [ZZ(3), ZZ(2)])

    assert dmp_clear_denoms(
        [EX(Rational(3, 2)), EX(Rational(9, 4))], 0, EX) == (EX(4), [EX(6), EX(9)])

    assert dmp_clear_denoms([EX(7)], 0, EX) == (EX(1), [EX(7)])
    assert dmp_clear_denoms([EX(sin(x)/x), EX(0)], 0, EX) == (EX(x), [EX(sin(x)), EX(0)])

    assert dmp_clear_denoms([[]], 1, QQ, ZZ) == (ZZ(1), [[]])

    assert dmp_clear_denoms([[QQ(1)]], 1, QQ, ZZ) == (ZZ(1), [[QQ(1)]])
    assert dmp_clear_denoms([[QQ(7)]], 1, QQ, ZZ) == (ZZ(1), [[QQ(7)]])

    assert dmp_clear_denoms([[QQ(7, 3)]], 1, QQ) == (ZZ(3), [[QQ(7)]])
    assert dmp_clear_denoms([[QQ(7, 3)]], 1, QQ, ZZ) == (ZZ(3), [[QQ(7)]])

    assert dmp_clear_denoms(
        [[QQ(3)], [QQ(1)], []], 1, QQ, ZZ) == (ZZ(1), [[QQ(3)], [QQ(1)], []])
    assert dmp_clear_denoms([[QQ(
        1)], [QQ(1, 2)], []], 1, QQ, ZZ) == (ZZ(2), [[QQ(2)], [QQ(1)], []])

    assert dmp_clear_denoms([QQ(3), QQ(
        1), QQ(0)], 0, QQ, ZZ, convert=True) == (ZZ(1), [ZZ(3), ZZ(1), ZZ(0)])
    assert dmp_clear_denoms([QQ(1), QQ(1, 2), QQ(
        0)], 0, QQ, ZZ, convert=True) == (ZZ(2), [ZZ(2), ZZ(1), ZZ(0)])

    assert dmp_clear_denoms([[QQ(3)], [QQ(
        1)], []], 1, QQ, ZZ, convert=True) == (ZZ(1), [[QQ(3)], [QQ(1)], []])
    assert dmp_clear_denoms([[QQ(1)], [QQ(1, 2)], []], 1, QQ, ZZ,
                            convert=True) == (ZZ(2), [[QQ(2)], [QQ(1)], []])

    assert dmp_clear_denoms(
        [[EX(Rational(3, 2))], [EX(Rational(9, 4))]], 1, EX) == (EX(4), [[EX(6)], [EX(9)]])
    assert dmp_clear_denoms([[EX(7)]], 1, EX) == (EX(1), [[EX(7)]])
    assert dmp_clear_denoms([[EX(sin(x)/x), EX(0)]], 1, EX) == (EX(x), [[EX(sin(x)), EX(0)]])
