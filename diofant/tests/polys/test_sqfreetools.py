"""Tests for square-free decomposition algorithms and related tools."""

import pytest

from diofant import FF, QQ, ZZ, DomainError, I, ring, sqrt
from diofant.polys.specialpolys import f_polys


__all__ = ()

f_0, f_1, f_2, f_3, f_4, f_5, f_6 = f_polys()


def test_dmp_sqf():
    R, x = ring('x', ZZ)

    assert R(0).sqf_part() == 0
    assert R(0).is_squarefree is True

    assert R(7).sqf_part() == 1
    assert R(7).is_squarefree is True

    assert (2*x + 2).sqf_part() == x + 1
    assert (2*x + 2).is_squarefree is True

    assert (x**3 + x + 1).sqf_part() == x**3 + x + 1
    assert (x**3 + x + 1).is_squarefree is True

    assert (x**3 + 1).sqf_part() == x**3 + 1
    assert (x**3 + 1).is_squarefree is True

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

    f = x**5 - x**3 - x**2 + 1

    assert f.sqf_part() == x**4 + x**3 - x - 1
    assert f.sqf_list() == (1, [(x**3 + 2*x**2 + 2*x + 1, 1), (x - 1, 2)])

    f = 2*x**5 + 16*x**4 + 50*x**3 + 76*x**2 + 56*x + 16

    assert f.sqf_list() == (2, [(x + 1, 2), (x + 2, 3)])

    R, x = ring('x', QQ)

    assert (2*x**2 + 4*x + 2).sqf_list() == (2, [(x + 1, 2)])

    R, x = ring('x', FF(2))

    assert (x**2 + 1).sqf_list() == (1, [(x + 1, 2)])

    R, x = ring('x', FF(3))

    assert (x**3 + 1).is_squarefree is False
    assert (x**3 + 1).sqf_part() == x + 1
    assert (x**10 + 2*x**7 + 2*x**4 + x).sqf_list() == (1, [(x, 1), (x + 1, 3),
                                                            (x + 2, 6)])

    f = x**10 + 2*x**7 + 2*x**4 + x

    assert f.is_squarefree is False
    assert f.sqf_list() == (1, [(x, 1), (x + 1, 3), (x + 2, 6)])

    R, x = ring('x', FF(5))

    f = x**8 + x**7 + 3*x**6 + x**4 + 2*x**2 + 2*x + 1

    assert f.sqf_part() == x**2 + 4*x + 3

    f = 3*x**2 + 2*x + 4

    assert f.is_squarefree is True

    f = 2*x**6 + 4*x**5 + 4*x**4 + 2*x**3 + 2*x**2 + x + 4

    assert f.is_squarefree is False

    R, x = ring('x', FF(11))

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

    R, x = ring('x', FF(53))

    f = x**6 + 2*x**5 + 5*x**4 + 26*x**3 + 41*x**2 + 39*x + 38

    assert f.is_squarefree is True

    R, x = ring('x', FF(102953))

    f = x**15 + x + 1

    assert f.is_squarefree is True

    F9 = FF(9)
    R, x = ring('x', F9)

    f = x + F9(4)

    assert f.is_squarefree is True

    f = f**3

    assert f.is_squarefree is False
    assert f.sqf_list() == (1, [(x + F9(4), 3)])

    f = (x + F9(4))*(x + F9(7))

    assert f.is_squarefree is True
    assert f.sqf_list() == (1, [(f, 1)])
    assert (2*f).is_squarefree is True
    assert (2*f).sqf_list() == (2, [(f, 1)])
    assert (F9(6)*f).is_squarefree is True
    assert (F9(6)*f).sqf_list() == (F9(6), [(f, 1)])

    f *= (x + F9(7))

    assert f.is_squarefree is False
    assert f.sqf_list() == (1, [(x + F9(4), 1), (x + F9(7), 2)])

    R, x = ring('x', QQ.algebraic_field(sqrt(2)))

    assert (x**2 - 3).sqf_norm() == (1, x**2 - 2*sqrt(2)*x - 1,
                                     (x**4 - 10*x**2 + 1).set_domain(QQ))

    R, x = ring('x', QQ.algebraic_field(sqrt(3)))

    assert (x**2 - 2).sqf_norm() == (1, x**2 - 2*sqrt(3)*x + 1,
                                     (x**4 - 10*x**2 + 1).set_domain(QQ))

    R, x, y = ring('x y', FF(2))

    f = y**2 + 1

    assert f.is_squarefree is False
    assert f.sqf_part() == y + 1
    assert f.sqf_list() == (1, [(y + 1, 2)])

    f = x**3 + x*y**2

    assert f.is_squarefree is False
    assert f.sqf_part() == x*(x + y)
    assert f.sqf_list() == (1, [(x, 1), (x + y, 2)])

    f1 = 1 + x**6*y**14 + x**2*y**4 + x**7
    f2 = 1 + x**3*y**4 + x + y

    assert f1.is_squarefree
    assert f2.is_squarefree
    assert (f1**2*f2**2).sqf_list() == (1, [(f1*f2, 2)])

    R, x, y = ring('x y', FF(5))

    f = x**5*y**5 + 1

    assert f.is_squarefree is False
    assert f.sqf_list() == (1, [(x*y + 1, 5)])

    F8 = FF(8)
    R, x, y = ring('x y', F8)

    f = x**8*y**8 + 1

    assert f.is_squarefree is False
    assert f.sqf_list() == (1, [(x*y + 1, 8)])

    R, x, y = ring('x y', ZZ)

    assert R(0).sqf_part() == 0
    assert R(0).is_squarefree is True

    assert R(7).sqf_part() == 1
    assert R(7).is_squarefree is True

    assert R(3).sqf_list() == (3, [])

    f = (x + y)**2

    assert f.is_squarefree is False
    assert f.sqf_part() == x + y

    f = x**2 + y**2

    assert f.is_squarefree is True
    assert f.sqf_part() == f

    assert (x**3 + 2*x**2*y + x*y**2).sqf_part() == x**2 + x*y
    assert (x**5 + 2*x**4*y + x**3*y**2).sqf_list() == (1, [(x + y, 2), (x, 3)])

    A = x**4 - 3*x**2 + 6
    D = x**6 - 5*x**4 + 5*x**2 + 4

    f, g = D, A - D.diff(x)*y
    res = f.resultant(g)
    h = (4*y**2 + 1).drop(x)

    assert res.sqf_list() == (45796, [(h, 3)])

    pytest.raises(DomainError, lambda: (x**2 - 1).sqf_norm())

    f = -x**5 + x**4 + x - 1

    assert f.sqf_list() == (-1, [(x**3 + x**2 + x + 1, 1), (x - 1, 2)])

    pytest.raises(DomainError, lambda: (x**2 + y**2).sqf_norm())

    R, x, y = ring('x y', QQ.algebraic_field(I))

    assert (x**2 + 2*I*x - 1).sqf_list() == (1, [(x + I, 2)])

    assert (x*y + y**2).sqf_norm() == (1, x*y - I*x + y**2 - 3*I*y - 2,
                                       (x**2*y**2 + x**2 + 2*x*y**3 + 2*x*y +
                                        y**4 + 5*y**2 + 4).set_domain(QQ))

    R, x, y, z = ring('x y z', FF(2))

    f1 = x**14*z**7 + x**7*y**14 + y**7 + 1
    f2 = y**7*z**14 + z**7 + 1
    f3 = x**2*y**2*z

    f = f1**2*f2**2*f3**2

    assert f1.is_squarefree
    assert f2.is_squarefree
    assert f3.is_squarefree is False
    assert f.sqf_list() == (1, [(f1*f2*z, 2), (x*y, 4)])

    R, x, y, z = ring('x y z', FF(3))

    f = (y**2 + 1)*(x + y - 1)*(x - y + 1)**2*(x + y + z)**3

    assert f.sqf_list() == (1, [((y**2 + 1)*(x + y - 1), 1),
                                (x - y + 1, 2), (x + y + z, 3)])

    F16 = FF(16)
    R, x, y, z = ring('x y z', F16)

    f1 = F16(9)*(x**3*y**2*z + y**2*z) + F16(10)*x*y*z**3
    f2 = F16(8)*x*y*z**13 + F16(9)*(x**10*y**2*z + x**2*y**2*z)

    assert f1.is_squarefree
    assert (f1**2).sqf_list() == (F16(3), [(f1/F16(9), 2)])
    assert f2.is_squarefree
    assert (f2**2).sqf_list() == (F16(3), [(f2/F16(9), 2)])

    h = f1*f2//(y*z)**2/F16(3)

    assert (f1**2*f2**2).sqf_list() == (F16(5), [(h, 2), (y*z, 4)])

    R, x, y, z = ring('x y z', ZZ)

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

    R, x, y, z = ring('x y z', QQ)

    assert R(0).is_squarefree
    assert (x - 1).is_squarefree
    assert ((x - 1)**2).is_squarefree is False
    assert (x**2 + y**2).is_squarefree
    assert ((x + y)**2).is_squarefree is False

    R, x, y, z, t = ring('x y z t', ZZ)

    assert f_6.is_squarefree is True
    assert f_6.sqf_part() == f_6


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
