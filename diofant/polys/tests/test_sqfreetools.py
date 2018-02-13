"""Tests for square-free decomposition algorithms and related tools. """

import pytest

from diofant.core import I
from diofant.domains import FF, QQ, ZZ
from diofant.functions import sqrt
from diofant.polys.polyerrors import DomainError
from diofant.polys.rings import ring
from diofant.polys.specialpolys import f_polys


__all__ = ()

f_0, f_1, f_2, f_3, f_4, f_5, f_6 = f_polys()


def test_dup_sqf():
    R, x = ring("x", ZZ)

    assert R.dmp_sqf_part(0) == 0
    assert R(0).is_squarefree is True

    assert R.dmp_sqf_part(7) == 1
    assert R(7).is_squarefree is True

    assert R.dmp_sqf_part(2*x + 2) == x + 1
    assert (2*x + 2).is_squarefree is True

    assert R.dmp_sqf_part(x**3 + x + 1) == x**3 + x + 1
    assert (x**3 + x + 1).is_squarefree is True

    assert R.dmp_sqf_part(-x**3 + x + 1) == x**3 - x - 1
    assert (-x**3 + x + 1).is_squarefree is True

    assert R.dmp_sqf_part(2*x**3 + 3*x**2) == 2*x**2 + 3*x
    assert (2*x**3 + 3*x**2).is_squarefree is False

    assert R.dmp_sqf_part(-2*x**3 + 3*x**2) == 2*x**2 - 3*x
    assert (-2*x**3 + 3*x**2).is_squarefree is False

    assert R.dmp_sqf_part(x**3 - 3*x - 2) == x**2 - x - 2
    assert (x**3 - 3*x - 2).is_squarefree is False

    assert R.dmp_sqf_list(0) == (0, [])
    assert R.dmp_sqf_list(1) == (1, [])

    assert R.dmp_sqf_list(x) == (1, [(x, 1)])
    assert R.dmp_sqf_list(2*x**2) == (2, [(x, 2)])
    assert R.dmp_sqf_list(3*x**3) == (3, [(x, 3)])

    assert R.dmp_sqf_list(-x**5 + x**4 + x - 1) == \
        (-1, [(x**3 + x**2 + x + 1, 1), (x - 1, 2)])
    assert R.dmp_sqf_list(x**8 + 6*x**6 + 12*x**4 + 8*x**2) == \
        ( 1, [(x, 2), (x**2 + 2, 3)])

    assert R.dmp_sqf_list(2*x**2 + 4*x + 2) == (2, [(x + 1, 2)])

    R, x = ring("x", QQ)
    assert R.dmp_sqf_list(2*x**2 + 4*x + 2) == (2, [(x + 1, 2)])

    R, x = ring("x", FF(2))
    assert R.dmp_sqf_list(x**2 + 1) == (1, [(x + 1, 2)])

    R, x = ring("x", FF(3))
    assert R.dmp_sqf_list(x**10 + 2*x**7 + 2*x**4 + x) == \
        (1, [(x, 1),
             (x + 1, 3),
             (x + 2, 6)])

    R1, x = ring("x", ZZ)
    R2, y = ring("y", FF(3))

    f = x**3 + 1
    g = y**3 + 1

    assert R1.dmp_sqf_part(f) == f
    assert R2.dmp_sqf_part(g) == y + 1

    assert f.is_squarefree is True
    assert g.is_squarefree is False

    R, x, y = ring("x,y", ZZ)

    A = x**4 - 3*x**2 + 6
    D = x**6 - 5*x**4 + 5*x**2 + 4

    f, g = D, R.dmp_sub(A, R.dmp_mul(R.dmp_diff(D, 1), y))
    res = R.dmp_resultant(f, g)
    h = (4*y**2 + 1).drop(x)

    assert R.drop(x).dmp_sqf_list(res) == (45796, [(h, 3)])

    pytest.raises(DomainError, lambda: R.dmp_sqf_norm(x**2 - 1))

    Rt, t = ring("t", ZZ)
    R, x = ring("x", Rt)
    assert R.dmp_sqf_list_include(t**3*x**2) == [(t**3, 1), (x, 2)]

    K = QQ.algebraic_field(sqrt(3))
    R, x = ring("x", K)
    _, X = ring("x", QQ)
    assert R.dmp_sqf_norm(x**2 - 2) == (1, x**2 + K([QQ(-2), QQ(0)])*x + 1, X**4 - 10*X**2 + 1)


def test_dmp_sqf():
    R, x, y = ring("x,y", ZZ)
    assert R.dmp_sqf_part(0) == 0
    assert R(0).is_squarefree is True

    assert R.dmp_sqf_part(7) == 1
    assert R(7).is_squarefree is True

    assert R.dmp_sqf_list(3) == (3, [])
    assert R.dmp_sqf_list_include(3) == [(3, 1)]

    R, x, y, z = ring("x,y,z", ZZ)
    assert f_0.is_squarefree is True
    assert (f_0**2).is_squarefree is False
    assert f_1.is_squarefree is True
    assert (f_1**2).is_squarefree is False
    assert f_2.is_squarefree is True
    assert (f_2**2).is_squarefree is False
    assert f_3.is_squarefree is True
    assert (f_3**2).is_squarefree is False
    assert f_5.is_squarefree is False
    assert (f_5**2).is_squarefree is False

    assert f_4.is_squarefree is True
    assert R.dmp_sqf_part(f_4) == -f_4

    assert R.dmp_sqf_part(f_5) == x + y - z

    R, x, y, z, t = ring("x,y,z,t", ZZ)
    assert f_6.is_squarefree is True
    assert R.dmp_sqf_part(f_6) == f_6

    R, x = ring("x", ZZ)
    f = -x**5 + x**4 + x - 1

    assert R.dmp_sqf_list(f) == (-1, [(x**3 + x**2 + x + 1, 1), (x - 1, 2)])
    assert R.dmp_sqf_list_include(f) == [(-x**3 - x**2 - x - 1, 1), (x - 1, 2)]

    f = 2*x**5 + 16*x**4 + 50*x**3 + 76*x**2 + 56*x + 16

    assert R.dmp_sqf_list(f) == (2, [(x + 1, 2), (x + 2, 3)])
    assert R.dmp_sqf_list_include(f) == [(2, 1), (x + 1, 2), (x + 2, 3)]

    assert R.dmp_sqf_list(f, all=True) == (2, [(1, 1), (x + 1, 2), (x + 2, 3)])
    assert R.dmp_sqf_list_include(f, all=True) == [(2, 1), (x + 1, 2), (x + 2, 3)]

    R, x, y = ring("x,y", ZZ)
    f = -x**5 + x**4 + x - 1

    assert R.dmp_sqf_list(f) == (-1, [(x**3 + x**2 + x + 1, 1), (x - 1, 2)])
    assert R.dmp_sqf_list_include(f) == [(-x**3 - x**2 - x - 1, 1), (x - 1, 2)]

    pytest.raises(DomainError, lambda: R.dmp_sqf_norm(x**2 + y**2))

    f = -x**2 + 2*x - 1
    assert R.dmp_sqf_list_include(f) == [(-1, 1), (x - 1, 2)]

    R, x, y = ring("x,y", FF(2))
    pytest.raises(NotImplementedError, lambda: R.dmp_sqf_list(y**2 + 1))
    pytest.raises(NotImplementedError,
                  lambda: R.dmp_sqf_part(x**3 + 2*x**2*y + x*y**2))

    R, x, y = ring("x,y", QQ.algebraic_field(I))
    assert R.dmp_sqf_list(x**2 + 2*I*x - 1) == (R.one.to_dense()[0][0],
                                                [(x + I, 2)])


def test_dup_gff_list():
    R, x = ring("x", ZZ)

    f = x**5 + 2*x**4 - x**3 - 2*x**2
    assert R.dup_gff_list(f) == [(x, 1), (x + 2, 4)]

    g = x**9 - 20*x**8 + 166*x**7 - 744*x**6 + 1965*x**5 - 3132*x**4 + 2948*x**3 - 1504*x**2 + 320*x
    assert R.dup_gff_list(g) == [(x**2 - 5*x + 4, 1), (x**2 - 5*x + 4, 2), (x, 3)]

    pytest.raises(ValueError, lambda: R.dup_gff_list(0))
