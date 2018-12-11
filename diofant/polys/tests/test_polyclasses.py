"""Tests for OO layer of several polynomial representations. """

import pytest

from diofant.domains import QQ, ZZ
from diofant.polys.polyclasses import DMP
from diofant.polys.polyerrors import ExactQuotientFailed, PolynomialError
from diofant.polys.specialpolys import f_polys


__all__ = ()

f_0, f_1, f_2, f_3, f_4, f_5, f_6 = [f.to_dense() for f in f_polys()]


def test_DMP___init__():
    f = DMP([[0], [], [0, 1, 2], [3]], ZZ)

    assert f.rep == [[1, 2], [3]]
    assert f.domain == ZZ
    assert f.lev == 1

    f = DMP([[1, 2], [3]], ZZ, 1)

    assert f.rep == [[1, 2], [3]]
    assert f.domain == ZZ
    assert f.lev == 1

    f = DMP({(1, 1): 1, (0, 0): 2}, ZZ, 1)

    assert f.rep == [[1, 0], [2]]
    assert f.domain == ZZ
    assert f.lev == 1

    assert f == DMP.from_monoms_coeffs(f.monoms(), f.coeffs(), f.lev, f.domain)


def test_DMP___eq__():
    f = DMP([[ZZ(1), ZZ(2)], [ZZ(3)]], ZZ)
    assert f == f
    assert f.eq(f)

    assert DMP([[ZZ(1), ZZ(2)], [ZZ(3)]], ZZ) == \
        DMP([[QQ(1), QQ(2)], [QQ(3)]], QQ)
    assert DMP([[QQ(1), QQ(2)], [QQ(3)]], QQ) == \
        DMP([[ZZ(1), ZZ(2)], [ZZ(3)]], ZZ)

    assert DMP([[[ZZ(1)]]], ZZ) != DMP([[ZZ(1)]], ZZ)
    assert DMP([[[ZZ(1)]]], ZZ).ne(DMP([[ZZ(1)]], ZZ))
    assert DMP([[ZZ(1)]], ZZ) != DMP([[[ZZ(1)]]], ZZ)


def test_DMP___bool__():
    assert bool(DMP([[]], ZZ)) is False
    assert bool(DMP([[1]], ZZ)) is True


def test_DMP_to_dict():
    f = DMP([[3], [], [2], [], [8]], ZZ)

    assert f.to_dict() == \
        {(4, 0): 3, (2, 0): 2, (0, 0): 8}
    assert f.to_diofant_dict() == \
        {(4, 0): ZZ.to_expr(3), (2, 0): ZZ.to_expr(2), (0, 0):
         ZZ.to_expr(8)}


def test_DMP_properties():
    assert DMP([[]], ZZ).is_zero is True
    assert DMP([[1]], ZZ).is_zero is False

    assert DMP([[1]], ZZ).is_one is True
    assert DMP([[2]], ZZ).is_one is False

    assert DMP([[1]], ZZ).is_ground is True
    assert DMP([[1], [2], [1]], ZZ).is_ground is False

    assert DMP([[1], [2, 0], [1, 0]], ZZ).is_squarefree is True
    assert DMP([[1], [2, 0], [1, 0, 0]], ZZ).is_squarefree is False

    assert DMP([[1, 2], [3]], ZZ).is_monic is True
    assert DMP([[2, 2], [3]], ZZ).is_monic is False

    assert DMP([[1, 2], [3]], ZZ).is_primitive is True
    assert DMP([[2, 4], [6]], ZZ).is_primitive is False


def test_DMP_arithmetics():
    f = DMP([[2], [2, 0]], ZZ)

    assert f*2 == DMP([[4], [4, 0]], ZZ)
    assert f//2 == DMP([[1], [1, 0]], ZZ)

    pytest.raises(ExactQuotientFailed, lambda: f.exquo_ground(3))

    f = DMP([[-5]], ZZ)
    g = DMP([[5]], ZZ)

    assert abs(f) == g

    assert -g == f

    h = DMP([[]], ZZ)

    assert f + g == h
    assert g + f == h
    assert f + 5 == h
    assert 5 + f == h

    h = DMP([[-10]], ZZ)

    assert f - g == h
    assert g - f == -h
    assert f - 5 == h
    assert 5 - f == -h

    h = DMP([[-25]], ZZ)

    assert f * g == h
    assert g * f == h
    assert f * 5 == h
    assert 5 * f == h

    h = DMP([[25]], ZZ)

    assert f**2 == h

    pytest.raises(TypeError, lambda: f**'x')

    f = DMP([[1], [], [1, 0, 0]], ZZ)
    g = DMP([[2], [-2, 0]], ZZ)

    q = DMP([[2], [2, 0]], ZZ)
    r = DMP([[8, 0, 0]], ZZ)

    assert f.pdiv(g) == (q, r)
    assert f.pquo(g) == q
    assert f.prem(g) == r

    pytest.raises(ExactQuotientFailed, lambda: f.pexquo(g))

    f = DMP([[1], [], [1, 0, 0]], ZZ)
    g = DMP([[1], [-1, 0]], ZZ)

    q = DMP([[1], [1, 0]], ZZ)
    r = DMP([[2, 0, 0]], ZZ)

    assert f.quo(g) == q

    assert divmod(f, g) == (q, r)
    assert f // g == q
    assert f // 2 == DMP([[0], [], [0, 0, 0]], QQ)
    assert f % g == r

    pytest.raises(ExactQuotientFailed, lambda: f.exquo(g))

    f = DMP([[-5]], ZZ)
    g = DMP([[5]], QQ)
    h = DMP([[]], QQ)

    assert f + g == g + f == h


def test_DMP_functionality():
    f = DMP([[1], [2, 0], [1, 0, 0]], ZZ)
    g = DMP([[1], [1, 0]], ZZ)
    h = DMP([[1]], ZZ)

    assert f.degree() == 2
    assert f.degree_list() == (2, 2)
    assert f.total_degree() == 2
    pytest.raises(TypeError, lambda: f.degree("spam"))

    assert f.LC() == ZZ(1)
    assert f.TC() == ZZ(0)
    assert f.nth(1, 1) == ZZ(2)

    pytest.raises(TypeError, lambda: f.nth(0, 'x'))

    assert f.max_norm() == 2
    assert f.l1_norm() == 4

    u = DMP([[2], [2, 0]], ZZ)

    assert f.diff(m=1, j=0) == u
    assert f.diff(m=1, j=1) == u
    pytest.raises(TypeError, lambda: f.diff(m='x', j=0))
    pytest.raises(TypeError, lambda: f.diff(j="spam"))

    u = DMP([1, 2, 1], ZZ)
    v = DMP([1, 2, 1], ZZ)

    assert f.eval(a=1, j=0) == u
    assert f.eval(a=1, j=1) == v
    pytest.raises(TypeError, lambda: f.eval(a=1, j="spam"))

    assert f.eval(1).eval(1) == ZZ(4)

    assert f.cofactors(g) == (g, g, h)
    assert f.gcd(g) == g
    assert f.lcm(g) == f

    u = DMP([[QQ(45), QQ(30), QQ(5)]], QQ)
    v = DMP([[QQ(1), QQ(2, 3), QQ(1, 9)]], QQ)

    assert u.monic() == v

    assert (4*f).content() == ZZ(4)
    assert (4*f).primitive() == (ZZ(4), f)

    f = DMP([[1], [2], [3], [4], [5], [6]], ZZ)

    assert f.trunc(3) == DMP([[1], [-1], [], [1], [-1], []], ZZ)

    f = DMP(f_4, ZZ)

    assert f.sqf_part() == -f
    assert f.sqf_list() == (ZZ(-1), [(-f, 1)])

    f = DMP([[-1], [], [], [5]], ZZ)
    g = DMP([[3, 1], [], []], ZZ)
    h = DMP([[45, 30, 5]], ZZ)

    r = DMP([675, 675, 225, 25], ZZ)

    assert f.subresultants(g) == [f, g, h]
    assert f.resultant(g) == r

    f = DMP([1, 3, 9, -13], ZZ)

    assert f.discriminant() == -11664

    f = DMP([QQ(2), QQ(0)], QQ)
    g = DMP([QQ(1), QQ(0), QQ(-16)], QQ)

    s = DMP([QQ(1, 32), QQ(0)], QQ)
    t = DMP([QQ(-1, 16)], QQ)
    h = DMP([QQ(1)], QQ)

    assert f.half_gcdex(g) == (s, h)
    assert f.gcdex(g) == (s, t, h)

    assert f.invert(g) == s

    f = DMP([[1], [2], [3]], QQ)

    pytest.raises(ValueError, lambda: f.half_gcdex(f))
    pytest.raises(ValueError, lambda: f.gcdex(f))

    pytest.raises(ValueError, lambda: f.invert(f))

    f = DMP([1, 0, 20, 0, 150, 0, 500, 0, 625, -2, 0, -10, 9], ZZ)
    g = DMP([1, 0, 0, -2, 9], ZZ)
    h = DMP([1, 0, 5, 0], ZZ)

    assert g.compose(h) == f
    assert f.decompose() == [g, h]

    f = DMP([[1], [2], [3]], QQ)

    pytest.raises(ValueError, lambda: f.decompose())
    pytest.raises(ValueError, lambda: f.sturm())

    pytest.raises(PolynomialError, lambda: f.all_coeffs())
    pytest.raises(PolynomialError, lambda: f.all_monoms())
    pytest.raises(PolynomialError, lambda: f.all_terms())

    pytest.raises(ValueError, lambda: f.shift(1))
    pytest.raises(PolynomialError, lambda: f.intervals())
    pytest.raises(PolynomialError, lambda: f.refine_root(1, 2))

    assert f.integrate() == DMP([[QQ(1, 3)], [QQ(1, 1)], [QQ(3, 1)], []], QQ)
    pytest.raises(TypeError, lambda: f.integrate(m="spam"))
    pytest.raises(TypeError, lambda: f.integrate(j="spam"))

    f = DMP([[-1], [], [], [5]], ZZ)
    g = DMP([[3, 1], [], []], QQ)

    r = DMP([675, 675, 225, 25], QQ)

    assert f.resultant(g) == r

    assert DMP([1, 2], QQ).resultant(DMP([3], ZZ)) == 3

    assert f.is_cyclotomic is False


def test_DMP_exclude():
    f = [[[[[[[[[[[[[[[[[[[[[[[[[[1]], [[]]]]]]]]]]]]]]]]]]]]]]]]]]
    J = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
         18, 19, 20, 21, 22, 24, 25]

    assert DMP(f, ZZ).exclude() == (J, DMP([1, 0], ZZ))
    assert DMP([[1], [1, 0]], ZZ).exclude() == ([], DMP([[1], [1, 0]], ZZ))


def test___hash__():
    # issue sympy/sympy#5571
    assert DMP([[1, 2], [3]], ZZ) == DMP([[int(1), int(2)], [int(3)]], ZZ)
    assert hash(DMP([[1, 2], [3]], ZZ)) == hash(DMP([[int(1), int(2)], [int(3)]], ZZ))
