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


def test_dmp_sqf():
    R, x = ring("x", ZZ)

    assert R(0).sqf_part() == 0
    assert R(0).is_squarefree is True

    assert R(7).sqf_part() == 1
    assert R(7).is_squarefree is True

    assert (2*x + 2).sqf_part() == x + 1
    assert (2*x + 2).is_squarefree is True

    assert (x**3 + x + 1).sqf_part() == x**3 + x + 1
    assert (x**3 + x + 1).is_squarefree is True

    assert (-x**3 + x + 1).sqf_part() == x**3 - x - 1
    assert (-x**3 + x + 1).is_squarefree is True

    assert (2*x**3 + 3*x**2).sqf_part() == 2*x**2 + 3*x
    assert (2*x**3 + 3*x**2).is_squarefree is False

    assert (-2*x**3 + 3*x**2).sqf_part() == 2*x**2 - 3*x
    assert (-2*x**3 + 3*x**2).is_squarefree is False

    assert (x**3 - 3*x - 2).sqf_part() == x**2 - x - 2
    assert (x**3 - 3*x - 2).is_squarefree is False

    assert R(0).sqf_list() == (0, [])
    assert R(1).sqf_list() == (1, [])

    assert x.sqf_list() == (1, [(x, 1)])
    assert (2*x**2).sqf_list() == (2, [(x, 2)])
    assert (3*x**3).sqf_list() == (3, [(x, 3)])

    assert (-x**5 + x**4 + x - 1).sqf_list() == (-1, [(x**3 + x**2 + x + 1, 1),
                                                      (x - 1, 2)])
    assert (x**8 + 6*x**6 + 12*x**4 + 8*x**2).sqf_list() == (1, [(x, 2),
                                                                 (x**2 + 2, 3)])

    assert (2*x**2 + 4*x + 2).sqf_list() == (2, [(x + 1, 2)])

    R, x = ring("x", QQ)
    assert (2*x**2 + 4*x + 2).sqf_list() == (2, [(x + 1, 2)])

    R, x = ring("x", FF(2))
    assert (x**2 + 1).sqf_list() == (1, [(x + 1, 2)])

    R, x = ring("x", FF(3))
    assert (x**10 + 2*x**7 + 2*x**4 + x).sqf_list() == (1, [(x, 1), (x + 1, 3),
                                                            (x + 2, 6)])

    R1, x = ring("x", ZZ)
    R2, y = ring("y", FF(3))

    f = x**3 + 1
    g = y**3 + 1

    assert f.sqf_part() == f
    assert g.sqf_part() == y + 1

    assert f.is_squarefree is True
    assert g.is_squarefree is False

    R, x = ring("x", FF(5))

    f = x**8 + x**7 + 3*x**6 + x**4 + 2*x**2 + 2*x + 1

    assert f.sqf_part() == x**2 + 4*x + 3

    f = 3*x**2 + 2*x + 4

    assert f.is_squarefree is True

    f = 2*x**6 + 4*x**5 + 4*x**4 + 2*x**3 + 2*x**2 + x + 4

    assert f.is_squarefree is False

    R, x = ring("x", FF(11))

    f = x**3 + 5*x**2 + 8*x + 4

    assert f.sqf_part() == x**2 + 3*x + 2

    assert R(0).is_squarefree is True
    assert R(0).sqf_list() == (0, [])

    assert R(1).is_squarefree is True
    assert R(1).sqf_list() == (1, [])

    assert (x + 1).is_squarefree is True
    assert (x + 1).sqf_list() == (1, [(x + 1, 1)])

    f = x**11 + 1

    assert f.is_squarefree is False
    assert f.sqf_list() == (1, [(x + 1, 11)])

    f = x**3 + 5*x**2 + 8*x + 4

    assert f.is_squarefree is False
    assert f.sqf_list() == (1, [(x + 1, 1), (x + 2, 2)])

    R, x = ring('x', FF(3))

    f = x**10 + 2*x**7 + 2*x**4 + x

    assert f.is_squarefree is False
    assert f.sqf_list() == (1, [(x, 1), (x + 1, 3), (x + 2, 6)])

    R, x = ring('x', FF(53))

    f = x**6 + 2*x**5 + 5*x**4 + 26*x**3 + 41*x**2 + 39*x + 38

    assert f.is_squarefree is True

    R, x = ring('x', FF(102953))

    f = x**15 + x + 1

    assert f.is_squarefree is True

    R, x, y = ring("x,y", ZZ)

    A = x**4 - 3*x**2 + 6
    D = x**6 - 5*x**4 + 5*x**2 + 4

    f, g = D, R.dmp_sub(A, R.dmp_mul(R.dmp_diff_in(D, 1, 0), y))
    res = f.resultant(g)
    h = (4*y**2 + 1).drop(x)

    assert res.sqf_list() == (45796, [(h, 3)])

    pytest.raises(DomainError, lambda: (x**2 - 1).sqf_norm())

    K = QQ.algebraic_field(sqrt(3))
    R, x = ring("x", K)
    _, X = ring("x", QQ)
    assert (x**2 - 2).sqf_norm() == (1, x**2 + K([QQ(-2), QQ(0)])*x + 1,
                                     X**4 - 10*X**2 + 1)

    R, x, y = ring("x,y", ZZ)
    assert R(0).sqf_part() == 0
    assert R(0).is_squarefree is True

    assert R(7).sqf_part() == 1
    assert R(7).is_squarefree is True

    assert R(3).sqf_list() == (3, [])

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
    assert f_4.sqf_part() == -f_4
    assert f_5.sqf_part() == x + y - z

    R, x, y, z, t = ring("x,y,z,t", ZZ)
    assert f_6.is_squarefree is True
    assert f_6.sqf_part() == f_6

    R, x = ring("x", ZZ)
    f = -x**5 + x**4 + x - 1

    assert f.sqf_list() == (-1, [(x**3 + x**2 + x + 1, 1), (x - 1, 2)])

    f = 2*x**5 + 16*x**4 + 50*x**3 + 76*x**2 + 56*x + 16

    assert f.sqf_list() == (2, [(x + 1, 2), (x + 2, 3)])

    R, x, y = ring("x,y", ZZ)
    f = -x**5 + x**4 + x - 1

    assert f.sqf_list() == (-1, [(x**3 + x**2 + x + 1, 1), (x - 1, 2)])

    pytest.raises(DomainError, lambda: (x**2 + y**2).sqf_norm())

    R, x, y = ring("x,y", FF(2))
    pytest.raises(NotImplementedError, lambda: (y**2 + 1).sqf_list())
    pytest.raises(NotImplementedError, lambda: (x**3 + 2*x**2*y + x*y**2).sqf_part())

    R, x, y = ring("x,y", QQ.algebraic_field(I))
    assert (x**2 + 2*I*x - 1).sqf_list() == (1, [(x + I, 2)])


def test_diofantissue_714():
    R, x, y, z = ring('x y z', ZZ)

    f = (x - y)*(z - 1)**2

    assert f.is_squarefree is False
    assert f.sqf_part() == (x - y)*(z - 1)
    assert f.sqf_list() == (1, [(x - y, 1), (z - 1, 2)])

    g = f
    f = (x - y)*(z - 1)
    assert (f*g).sqf_list() == (1, [(x - y, 2), (z - 1, 3)])

    assert f.is_squarefree
    assert ((x - y)*f).is_squarefree is False
