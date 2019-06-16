"""Tests for dense recursive polynomials' basic tools. """

import random

import pytest

from diofant import oo
from diofant.domains import FF, ZZ
from diofant.polys.densebasic import (dmp_apply_pairs, dmp_convert,
                                      dmp_deflate, dmp_degree_in, dmp_eject,
                                      dmp_exclude, dmp_from_dict, dmp_ground,
                                      dmp_ground_p, dmp_include, dmp_inflate,
                                      dmp_inject, dmp_nest, dmp_normal,
                                      dmp_one, dmp_one_p, dmp_permute,
                                      dmp_raise, dmp_strip, dmp_swap,
                                      dmp_terms_gcd, dmp_to_dict, dmp_zero,
                                      dmp_zero_p, dmp_zeros, dup_inflate,
                                      dup_random, dup_reverse)
from diofant.polys.rings import ring
from diofant.polys.specialpolys import f_polys


__all__ = ()

f_0, f_1, f_2, f_3, f_4, f_5, f_6 = [f.to_dense() for f in f_polys()]


def test_dmp_LC():
    R, x = ring('x', ZZ)

    assert R.dmp_LC(0) == 0
    assert R.dmp_LC(1) == 1
    assert R.dmp_LC(2*x**3 + 3*x**2 + 4*x + 5) == 2
    assert R.dmp_LC(3*x**2 + 1) == 3
    assert R.dmp_LC(x + 2) == 1

    R, x, y = ring('x y', ZZ)
    R1 = R.drop(x)

    assert R.dmp_LC(0) == 0
    assert R.dmp_LC(2*x*y**2 + 3*x*y + 4*x + 5) == 2*R1.y**2 + 3*R1.y + 4

    R, x, y, z = ring('x y z', ZZ)
    R12 = R.drop(x)

    assert R.dmp_LC(0) == 0
    assert R.dmp_LC(2*x*y + 3*x*z + 4*x + 5) == 2*R12.y + 3*R12.z + 4


def test_dmp_TC():
    R, x = ring('x', ZZ)

    assert R.dmp_TC(0) == 0
    assert R.dmp_TC(1) == 1
    assert R.dmp_TC(x + 2) == 2
    assert R.dmp_TC(3*x**2 + 1) == 1
    assert R.dmp_TC(2*x**3 + 3*x**2 + 4*x + 5) == 5

    R, x, y = ring('x y', ZZ)

    assert R.dmp_TC(0) == 0
    assert R.dmp_TC(2*x*y**2 + 3*x*y + 4*x + 5) == 5

    R, x, y, z = ring('x y z', ZZ)

    assert R.dmp_TC(0) == 0
    assert R.dmp_TC(2*x*y + 3*x*z + 4*x + 5) == 5


def test_dmp_ground_LC():
    R, x, y = ring('x y', ZZ)

    assert R.dmp_ground_LC(0) == 0
    assert R.dmp_ground_LC(2*x*y**2 + 3*x*y + 4*x + 5) == 2

    R, x, y, z = ring('x y z', ZZ)

    assert R.dmp_ground_LC(0) == 0
    assert R.dmp_ground_LC(2*x*y + 3*x*z + 4*x + 5) == 2


def test_dmp_ground_TC():
    R, x, y = ring('x y', ZZ)

    assert R.dmp_ground_TC(0) == 0
    assert R.dmp_ground_TC(2*x*y**2 + 3*x*y + 4*x + 5) == 5

    R, x, y, z = ring('x y z', ZZ)

    assert R.dmp_ground_TC(0) == 0
    assert R.dmp_ground_TC(2*x*y + 3*x*z + 4*x + 5) == 5


def test_dmp_degree_in():
    R, x = ring('x', ZZ)

    assert R.dmp_degree_in(0, 0) == -oo
    assert R.dmp_degree_in(1, 0) == 0
    assert R.dmp_degree_in(x, 0) == 1
    assert R.dmp_degree_in(x**4 + 1, 0) == 4
    assert R.dmp_degree_in(x**3 + 2*x**2 + 3, 0) == 3
    assert R.dmp_degree_in(x**3 + x**2 + 2*x, 0) == 3

    R, x, y = ring('x y', ZZ)

    assert R.dmp_degree_in(0, 0) == -oo
    assert R.dmp_degree_in(1, 0) == 0
    assert R.dmp_degree_in(2*x + 1, 0) == 1
    assert R.dmp_degree_in(2*x + y**2 + 2*y + 3, 0) == 1

    pytest.raises(IndexError, lambda: R.dmp_degree_in(1, -5))

    R, x, y, z = ring('x y z', ZZ)

    assert R.dmp_degree_in(0, 0) == -oo
    assert R.dmp_degree_in(0, 1) == -oo
    assert R.dmp_degree_in(0, 2) == -oo

    assert R.dmp_degree_in(1, 0) == 0
    assert R.dmp_degree_in(1, 1) == 0
    assert R.dmp_degree_in(1, 2) == 0

    f = R.from_dense(f_4)

    assert R.dmp_degree_in(f, 0) == 9
    assert R.dmp_degree_in(f, 1) == 12
    assert R.dmp_degree_in(f, 2) == 8

    R, x, y, z, t = ring('x y z t', ZZ)

    f = R.from_dense(f_6)

    assert R.dmp_degree_in(f, 0) == 4
    assert R.dmp_degree_in(f, 1) == 4
    assert R.dmp_degree_in(f, 2) == 6
    assert R.dmp_degree_in(f, 3) == 3


def test_dmp_degree_list():
    R, x, y, z = ring('x y z', ZZ)

    assert R.dmp_degree_list(R.from_dense(f_0)) == (2, 2, 2)
    assert R.dmp_degree_list(R.from_dense(f_1)) == (3, 3, 3)
    assert R.dmp_degree_list(R.from_dense(f_2)) == (5, 3, 3)
    assert R.dmp_degree_list(R.from_dense(f_3)) == (5, 4, 7)
    assert R.dmp_degree_list(R.from_dense(f_4)) == (9, 12, 8)
    assert R.dmp_degree_list(R.from_dense(f_5)) == (3, 3, 3)

    R, x, y, z, t = ring('x y z t', ZZ)

    assert R.dmp_degree_list(0) == (-oo, -oo, -oo, -oo)
    assert R.dmp_degree_list(1) == (0, 0, 0, 0)

    assert R.dmp_degree_list(R.from_dense(f_6)) == (4, 4, 6, 3)


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


def test_dup_reverse():
    assert dup_reverse([1, 2, 0, 3]) == [3, 0, 2, 1]
    assert dup_reverse([1, 2, 3, 0]) == [3, 2, 1]


def test_dmp_normal():
    assert dmp_normal([0, 1.5, 2, 3], 0, ZZ) == [ZZ(1), ZZ(2), ZZ(3)]

    assert (dmp_normal([0, 0, 2, 1, 0, 11, 0], 0, ZZ) ==
            [ZZ(2), ZZ(1), ZZ(0), ZZ(11), ZZ(0)])

    assert (dmp_normal([[0], [], [0, 2, 1], [0], [11], []], 1, ZZ) ==
            [[ZZ(2), ZZ(1)], [], [ZZ(11)], []])

    F5 = FF(5)
    assert dmp_normal([5, 10, 21, -3], 0, F5) == [F5(1), F5(2)]
    F11 = FF(11)
    assert dmp_normal([11, 22, 17, 1, 0], 0, F11) == [F11(6), F11(1), F11(0)]


def test_dmp_convert():
    K0, K1 = ZZ.poly_ring('x'), ZZ

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


def test_dmp_zero():
    assert dmp_zero(0) == []
    assert dmp_zero(2) == [[[]]]


def test_dmp_one_p():
    assert dmp_one_p([1], 0, ZZ) is True
    assert dmp_one_p([[1]], 1, ZZ) is True
    assert dmp_one_p([[[1]]], 2, ZZ) is True
    assert dmp_one_p([[[12]]], 2, ZZ) is False


def test_dmp_one():
    assert dmp_one(0, ZZ) == [ZZ(1)]
    assert dmp_one(2, ZZ) == [[[ZZ(1)]]]


def test_dmp_ground_p():
    assert dmp_ground_p([], 0, 0) is True
    assert dmp_ground_p([[]], 0, 1) is True
    assert dmp_ground_p([[]], 1, 1) is False

    assert dmp_ground_p([[ZZ(1)]], 1, 1) is True
    assert dmp_ground_p([[[ZZ(2)]]], 2, 2) is True

    assert dmp_ground_p([[[ZZ(2)]]], 3, 2) is False
    assert dmp_ground_p([[[ZZ(3)], []]], 3, 2) is False

    assert dmp_ground_p([], None, 0) is True
    assert dmp_ground_p([[]], None, 1) is True

    assert dmp_ground_p([ZZ(1)], None, 0) is True
    assert dmp_ground_p([[[ZZ(1)]]], None, 2) is True

    assert dmp_ground_p([[[ZZ(3)], []]], None, 2) is False


def test_dmp_ground():
    assert dmp_ground(ZZ(0), 2) == [[[]]]

    assert dmp_ground(ZZ(7), -1) == ZZ(7)
    assert dmp_ground(ZZ(7), 0) == [ZZ(7)]
    assert dmp_ground(ZZ(7), 2) == [[[ZZ(7)]]]


def test_dmp_zeros():
    assert dmp_zeros(4, 0, ZZ) == [[], [], [], []]

    assert dmp_zeros(0, 2, ZZ) == []
    assert dmp_zeros(1, 2, ZZ) == [[[[]]]]
    assert dmp_zeros(2, 2, ZZ) == [[[[]]], [[[]]]]
    assert dmp_zeros(3, 2, ZZ) == [[[[]]], [[[]]], [[[]]]]

    assert dmp_zeros(3, -1, ZZ) == [0, 0, 0]


def test_dmp_from_to_dict():
    assert dmp_from_dict({}, 0, ZZ) == []

    assert dmp_to_dict([], 0) == {}

    f = [3, 0, 0, 2, 0, 0, 0, 0, 8]
    h = {(8,): 3, (5,): 2, (0,): 8}

    assert dmp_from_dict(h, 0, ZZ) == f

    assert dmp_to_dict(f, 0) == h

    R,  x, y = ring("x,y", ZZ)

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


def test_dmp_swap():
    f = dmp_normal([[1, 0, 0], [], [1, 0], [], [1]], 1, ZZ)
    g = dmp_normal([[1, 0, 0, 0, 0], [1, 0, 0], [1]], 1, ZZ)

    assert dmp_swap(f, 1, 1, 1, ZZ) == f

    assert dmp_swap(f, 0, 1, 1, ZZ) == g
    assert dmp_swap(g, 0, 1, 1, ZZ) == f

    pytest.raises(IndexError, lambda: dmp_swap(f, -1, -7, 1, ZZ))

    f = dmp_normal([[[2], [1, 0]], []], 2, ZZ)

    assert dmp_swap(f, 0, 1, 2, ZZ) == dmp_normal([[[2], []], [[1, 0], []]], 2, ZZ)
    assert dmp_swap(f, 1, 2, 2, ZZ) == dmp_normal([[[1], [2, 0]], [[]]], 2, ZZ)
    assert dmp_swap(f, 0, 2, 2, ZZ) == dmp_normal([[[1, 0]], [[2, 0], []]], 2, ZZ)


def test_dmp_permute():
    f = dmp_normal([[1, 0, 0], [], [1, 0], [], [1]], 1, ZZ)
    g = dmp_normal([[1, 0, 0, 0, 0], [1, 0, 0], [1]], 1, ZZ)

    assert dmp_permute(f, [0, 1], 1, ZZ) == f
    assert dmp_permute(g, [0, 1], 1, ZZ) == g

    assert dmp_permute(f, [1, 0], 1, ZZ) == g
    assert dmp_permute(g, [1, 0], 1, ZZ) == f

    f = dmp_normal([[[2], [1, 0]], []], 2, ZZ)

    assert dmp_permute(f, [1, 0, 2], 2, ZZ) == dmp_normal([[[2], []], [[1, 0], []]], 2, ZZ)
    assert dmp_permute(f, [1, 2, 0], 2, ZZ) == dmp_normal([[[1], []], [[2, 0], []]], 2, ZZ)


def test_dmp_nest():
    assert dmp_nest(ZZ(1), 2, ZZ) == [[[1]]]

    assert dmp_nest([[1]], 0, ZZ) == [[1]]
    assert dmp_nest([[1]], 1, ZZ) == [[[1]]]
    assert dmp_nest([[1]], 2, ZZ) == [[[[1]]]]


def test_dmp_raise():
    assert dmp_raise([], 2, 0, ZZ) == [[[]]]
    assert dmp_raise([[1]], 0, 1, ZZ) == [[1]]

    assert dmp_raise([[1, 2, 3], [], [2, 3]], 2, 1, ZZ) == \
        [[[[1]], [[2]], [[3]]], [[[]]], [[[2]], [[3]]]]


def test_dmp_deflate():
    assert dmp_deflate(([2],), 0, ZZ) == ((1,), ([2],))
    assert dmp_deflate(([], []), 0, ZZ) == ((1,), ([], []))

    assert dmp_deflate(([1, 2, 3],), 0, ZZ) == ((1,), ([1, 2, 3],))
    assert dmp_deflate(([1, 0, 2, 0, 3],), 0, ZZ) == ((2,), ([1, 2, 3],))

    assert dmp_deflate(([1, 0, 2, 0, 3], [2, 0, 0]), 0, ZZ) == \
        ((2,), ([1, 2, 3], [2, 0]))
    assert dmp_deflate(([1, 0, 2, 0, 3], [4, 0, 0]), 0, ZZ) == \
        ((2,), ([1, 2, 3], [4, 0]))
    assert dmp_deflate(([1, 0, 2, 0, 3], [2, 1, 0]), 0, ZZ) == \
        ((1,), ([1, 0, 2, 0, 3], [2, 1, 0]))

    assert dmp_deflate(([[]],), 1, ZZ) == \
        ((1, 1), ([[]],))
    assert dmp_deflate(([[]], [[]]), 1, ZZ) == \
        ((1, 1), ([[]], [[]]))

    assert dmp_deflate(([[1]], [[]]), 1, ZZ) == \
        ((1, 1), ([[1]], [[]]))
    assert dmp_deflate(([[1]], [[2]]), 1, ZZ) == \
        ((1, 1), ([[1]], [[2]]))
    assert dmp_deflate(([[1]], [[2, 0]]), 1, ZZ) == \
        ((1, 1), ([[1]], [[2, 0]]))

    assert dmp_deflate(([[2, 0]], [[2, 0]]), 1, ZZ) == \
        ((1, 1), ([[2, 0]], [[2, 0]]))

    assert dmp_deflate(
        ([[2]], [[2, 0, 0]]), 1, ZZ) == ((1, 2), ([[2]], [[2, 0]]))
    assert dmp_deflate(
        ([[2, 0, 0]], [[2, 0, 0]]), 1, ZZ) == ((1, 2), ([[2, 0]], [[2, 0]]))

    assert dmp_deflate(([2, 0, 0], [1, 0, 4, 0, 1]), 0, ZZ) == \
        ((2,), ([2, 0], [1, 4, 1]))

    f = [[1, 0, 0], [], [1, 0], [], [1]]
    g = [[1, 0, 1, 0], [], [1]]

    assert dmp_deflate((f,), 1, ZZ) == \
        ((2, 1), ([[1, 0, 0], [1, 0], [1]],))

    assert dmp_deflate((f, g), 1, ZZ) == \
        ((2, 1), ([[1, 0, 0], [1, 0], [1]],
                  [[1, 0, 1, 0], [1]]))


def test_dup_inflate():
    assert dup_inflate([], 17, ZZ) == []

    assert dup_inflate([1, 2, 3], 1, ZZ) == [1, 2, 3]
    assert dup_inflate([1, 2, 3], 2, ZZ) == [1, 0, 2, 0, 3]
    assert dup_inflate([1, 2, 3], 3, ZZ) == [1, 0, 0, 2, 0, 0, 3]
    assert dup_inflate([1, 2, 3], 4, ZZ) == [1, 0, 0, 0, 2, 0, 0, 0, 3]

    pytest.raises(IndexError, lambda: dup_inflate([1, 2, 3], 0, ZZ))


def test_dmp_inflate():
    assert dmp_inflate([1], (3,), 0, ZZ) == [1]

    assert dmp_inflate([[]], (3, 7), 1, ZZ) == [[]]
    assert dmp_inflate([[2]], (1, 2), 1, ZZ) == [[2]]

    assert dmp_inflate([[2, 0]], (1, 1), 1, ZZ) == [[2, 0]]
    assert dmp_inflate([[2, 0]], (1, 2), 1, ZZ) == [[2, 0, 0]]
    assert dmp_inflate([[2, 0]], (1, 3), 1, ZZ) == [[2, 0, 0, 0]]

    assert dmp_inflate([[1, 0, 0], [1], [1, 0]], (2, 1), 1, ZZ) == \
        [[1, 0, 0], [], [1], [], [1, 0]]

    pytest.raises(IndexError, lambda: dmp_inflate([[]], (-3, 7), 1, ZZ))


def test_dmp_exclude():
    assert dmp_exclude([[[]]], 2, ZZ) == ([], [[[]]], 2)
    assert dmp_exclude([[[7]]], 2, ZZ) == ([], [[[7]]], 2)

    assert dmp_exclude([1, 2, 3], 0, ZZ) == ([], [1, 2, 3], 0)
    assert dmp_exclude([[1], [2, 3]], 1, ZZ) == ([], [[1], [2, 3]], 1)

    assert dmp_exclude([[1, 2, 3]], 1, ZZ) == ([0], [1, 2, 3], 0)
    assert dmp_exclude([[1], [2], [3]], 1, ZZ) == ([1], [1, 2, 3], 0)

    assert dmp_exclude([[[1, 2, 3]]], 2, ZZ) == ([0, 1], [1, 2, 3], 0)
    assert dmp_exclude([[[1]], [[2]], [[3]]], 2, ZZ) == ([1, 2], [1, 2, 3], 0)


def test_dmp_include():
    assert dmp_include([1, 2, 3], [], 0, ZZ) == [1, 2, 3]

    assert dmp_include([1, 2, 3], [0], 0, ZZ) == [[1, 2, 3]]
    assert dmp_include([1, 2, 3], [1], 0, ZZ) == [[1], [2], [3]]

    assert dmp_include([1, 2, 3], [0, 1], 0, ZZ) == [[[1, 2, 3]]]
    assert dmp_include([1, 2, 3], [1, 2], 0, ZZ) == [[[1]], [[2]], [[3]]]


def test_dmp_inject():
    R,  x, y = ring("x,y", ZZ)

    assert dmp_inject([], 0, R) == ([[[]]], 2)
    assert dmp_inject([[]], 1, R) == ([[[[]]]], 3)

    assert dmp_inject([R(1)], 0, R) == ([[[1]]], 2)
    assert dmp_inject([[R(1)]], 1, R) == ([[[[1]]]], 3)

    assert dmp_inject([R(1), x + 2], 0, R) == ([[[1]], [[1], [2]]], 2)
    assert dmp_inject([R(1), x + 2], 0, R, front=True) == ([[[1]], [[1, 2]]], 2)

    assert dmp_inject([R(1), 2*x + 3*y + 4], 0, R) == ([[[1]], [[2], [3, 4]]], 2)

    f = [3*x**2 + 7*x*y + 5*y**2, 2*x, R(0), x*y**2 + 11]
    g = [[[3], [7, 0], [5, 0, 0]], [[2], []], [[]], [[1, 0, 0], [11]]]

    assert dmp_inject(f, 0, R) == (g, 2)


def test_dmp_eject():
    R,  x, y = ring("x,y", ZZ)

    assert dmp_eject([[[]]], 2, R) == []
    assert dmp_eject([[[[]]]], 3, R) == [[]]

    assert dmp_eject([[[1]]], 2, R) == [R(1)]
    assert dmp_eject([[[[1]]]], 3, R) == [[R(1)]]

    f = [[[1]], [[2], [3, 4]]]

    assert dmp_eject(f, 2, R) == [R(1), 2*x + 3*y + 4]
    assert dmp_eject(f, 2, R, front=True) == [R(3), x + 2*y + 4]

    f = [3*x**2 + 7*x*y + 5*y**2, 2*x, R(0), x*y**2 + 11]
    g = [[[3], [7, 0], [5, 0, 0]], [[2], []], [[]], [[1, 0, 0], [11]]]

    assert dmp_eject(g, 2, R) == f


def test_dmp_terms_gcd():
    assert dmp_terms_gcd([], 0, ZZ) == ((0,), [])
    assert dmp_terms_gcd([1, 0, 1], 0, ZZ) == ((0,), [1, 0, 1])
    assert dmp_terms_gcd([1, 0, 1, 0], 0, ZZ) == ((1,), [1, 0, 1])
    assert dmp_terms_gcd([1, 3, 1, 4, 2, 0], 0, ZZ) == ((1,), [1, 3, 1, 4, 2])
    assert dmp_terms_gcd([1, 0, 1, 0, 0], 0, ZZ) == ((2,), [1, 0, 1])

    assert dmp_terms_gcd([[]], 1, ZZ) == ((0, 0), [[]])

    assert dmp_terms_gcd([1, 0, 1, 0], 0, ZZ) == ((1,), [1, 0, 1])
    assert dmp_terms_gcd([[1], [], [1], []], 1, ZZ) == ((1, 0), [[1], [], [1]])

    assert dmp_terms_gcd([[1, 0], [], [1]], 1, ZZ) == ((0, 0), [[1, 0], [], [1]])
    assert dmp_terms_gcd([[1, 0], [1, 0, 0], [], []], 1, ZZ) == ((2, 1), [[1], [1, 0]])


def test_dmp_apply_pairs():
    def h(a, b):
        return a*b

    assert dmp_apply_pairs([1, 2, 3], [4, 5, 6], h, [], 0, ZZ) == [4, 10, 18]

    assert dmp_apply_pairs([2, 3], [4, 5, 6], h, [], 0, ZZ) == [10, 18]
    assert dmp_apply_pairs([1, 2, 3], [5, 6], h, [], 0, ZZ) == [10, 18]

    assert dmp_apply_pairs(
        [[1, 2], [3]], [[4, 5], [6]], h, [], 1, ZZ) == [[4, 10], [18]]

    assert dmp_apply_pairs(
        [[1, 2], [3]], [[4], [5, 6]], h, [], 1, ZZ) == [[8], [18]]
    assert dmp_apply_pairs(
        [[1], [2, 3]], [[4, 5], [6]], h, [], 1, ZZ) == [[5], [18]]

    def h2(x, y, z):
        return 2*x + y - z

    f = [[1], [2, 3, 4], [5]]
    g = [[3], [2, 1]]
    assert dmp_apply_pairs(f, g, h2, (1,), 1, ZZ) == [[1], [3, 5, 10], [1, 10]]
    assert dmp_apply_pairs(g, f, h2, (1,), 1, ZZ) == [[1, 2, 9], [3, 6]]

    assert dmp_apply_pairs([1, 2, 3], [3, 2, 1], h2, [1], 0, ZZ) == [4, 5, 6]


def test_dmp_slice_in():
    R, x = ring('x', ZZ)

    f = x**3 + 2*x**2 + 3*x + 4

    assert f.slice(0, 0) == 0
    assert f.slice(0, 1) == 4
    assert f.slice(0, 2) == 3*x + 4
    assert f.slice(0, 3) == 2*x**2 + 3*x + 4

    assert f.slice(0, 4) == f
    assert f.slice(0, 9) == f

    assert f.slice(1, 0) == 0
    assert f.slice(1, 1) == 0
    assert f.slice(1, 2) == 3*x
    assert f.slice(1, 3) == 2*x**2 + 3*x
    assert f.slice(1, 4) == x**3 + 2*x**2 + 3*x

    pytest.raises(IndexError, lambda: R.dmp_slice_in(f, 0, 0, -1))

    assert (x + 2).slice(0, 3) == x + 2

    R, x, y = ring('x y', ZZ)

    f = x + 2*y**2 + 3*y + 4

    assert f.slice(1, 2) == f
    assert f.slice(2, 1) == 2*y**2 + 3*y + 5


def test_dup_random():
    f = dup_random(0, -10, 10, ZZ)

    assert dmp_degree_in(f, 0, 0) == 0
    assert all(-10 <= c <= 10 for c in f)

    f = dup_random(1, -20, 20, ZZ)

    assert dmp_degree_in(f, 0, 0) == 1
    assert all(-20 <= c <= 20 for c in f)

    f = dup_random(2, -30, 30, ZZ)

    assert dmp_degree_in(f, 0, 0) == 2
    assert all(-30 <= c <= 30 for c in f)

    f = dup_random(3, -40, 40, ZZ)

    assert dmp_degree_in(f, 0, 0) == 3
    assert all(-40 <= c <= 40 for c in f)

    f = dup_random(3, -400, 400, ZZ)

    assert dmp_degree_in(f, 0, 0) == 3
    assert all(-400 <= c <= 400 for c in f)

    random.seed(11)
    assert dup_random(10, -1, 1, ZZ) == [1, 0, 0, -1, 0, 0, -1, 1, 1, 0, 1]

    for i in range(10):
        f = dup_random(3, -10, 10, ZZ, percent=50)
        assert f[0]
        assert len([c for c in f if c == 0]) == 2
