"""Test sparse polynomials. """

import functools
import operator

import pytest

from diofant import oo, pi, sqrt
from diofant.abc import t, x, y, z
from diofant.core import Symbol, symbols
from diofant.domains import EX, FF, QQ, RR, ZZ
from diofant.polys.fields import field
from diofant.polys.orderings import grlex, lex
from diofant.polys.polyconfig import using
from diofant.polys.polyerrors import (CoercionFailed, ExactQuotientFailed,
                                      GeneratorsError, GeneratorsNeeded,
                                      MultivariatePolynomialError)
from diofant.polys.rings import PolyElement, PolynomialRing, ring, sring


__all__ = ()


def test_PolynomialRing___init__():
    assert len(PolynomialRing(ZZ, "x,y,z").gens) == 3
    assert len(ZZ.poly_ring(x).gens) == 1
    assert len(ZZ.poly_ring("x", "y", "z").gens) == 3
    assert len(ZZ.poly_ring(x, y, z).gens) == 3

    pytest.raises(GeneratorsNeeded, lambda: ZZ.poly_ring())
    pytest.raises(GeneratorsError, lambda: ZZ.poly_ring(0))

    assert ZZ.poly_ring(t).poly_ring("x").domain == ZZ.poly_ring(t)
    assert PolynomialRing('ZZ[t]', "x").domain == ZZ.poly_ring(t)

    pytest.raises(GeneratorsError, lambda: ZZ.poly_ring("x").poly_ring("x"))

    _lex = Symbol("lex")
    assert PolynomialRing(ZZ, "x").order == lex
    assert PolynomialRing(ZZ, "x", _lex).order == lex
    assert PolynomialRing(ZZ, "x", 'lex').order == lex

    R1 = ZZ.poly_ring("x", "y")
    R2 = ZZ.poly_ring("x", "y")
    R3 = ZZ.poly_ring("x", "y", "z")

    assert R1.x == R1.gens[0]
    assert R1.y == R1.gens[1]
    assert R1.x == R2.x
    assert R1.y == R2.y
    assert R1.x != R3.x
    assert R1.y != R3.y

    R4 = ZZ.poly_ring("gens")
    assert type(R4.gens) is tuple

    pytest.raises(GeneratorsError, lambda: PolynomialRing(ZZ, {1: 2}))
    pytest.raises(GeneratorsError, lambda: PolynomialRing(ZZ, ["x", ["y"]]))


def test_PolynomialRing___hash__():
    R, x, y, z = ring("x,y,z", QQ)
    assert hash(R)


def test_PolynomialRing___eq__():
    assert ring("x,y,z", QQ)[0] == ring("x,y,z", QQ)[0]
    assert ring("x,y,z", QQ)[0] is ring("x,y,z", QQ)[0]

    assert ring("x,y,z", QQ)[0] != ring("x,y,z", ZZ)[0]
    assert ring("x,y,z", QQ)[0] is not ring("x,y,z", ZZ)[0]

    assert ring("x,y,z", ZZ)[0] != ring("x,y,z", QQ)[0]
    assert ring("x,y,z", ZZ)[0] is not ring("x,y,z", QQ)[0]

    assert ring("x,y,z", QQ)[0] != ring("x,y", QQ)[0]
    assert ring("x,y,z", QQ)[0] is not ring("x,y", QQ)[0]

    assert ring("x,y", QQ)[0] != ring("x,y,z", QQ)[0]
    assert ring("x,y", QQ)[0] is not ring("x,y,z", QQ)[0]


def test_PolynomialRing_ring_new():
    R, x, y, z = ring("x,y,z", QQ)

    assert R.ring_new(7) == R(7)
    assert R.ring_new(7*x*y*z) == 7*x*y*z

    f = x**2 + 2*x*y + 3*x + 4*z**2 + 5*z + 6

    assert R.ring_new([[[1]], [[2], [3]], [[4, 5, 6]]]) == f
    assert R.ring_new({(2, 0, 0): 1, (1, 1, 0): 2, (1, 0, 0): 3, (0, 0, 2): 4, (0, 0, 1): 5, (0, 0, 0): 6}) == f
    assert R.ring_new([((2, 0, 0), 1), ((1, 1, 0), 2), ((1, 0, 0), 3), ((0, 0, 2), 4), ((0, 0, 1), 5), ((0, 0, 0), 6)]) == f


def test_PolynomialRing_drop():
    R,  x, y, z = ring("x,y,z", ZZ)

    assert R.drop(x) == ZZ.poly_ring("y", "z")
    assert R.drop(y) == ZZ.poly_ring("x", "z")
    assert R.drop(z) == ZZ.poly_ring("x", "y")

    assert R.drop(0) == ZZ.poly_ring("y", "z")
    assert R.drop(0).drop(0) == ZZ.poly_ring("z")
    assert R.drop(0).drop(0).drop(0) == ZZ

    assert R.drop(1) == ZZ.poly_ring("x", "z")

    assert R.drop(2) == ZZ.poly_ring("x", "y")
    assert R.drop(2).drop(1) == ZZ.poly_ring("x")
    assert R.drop(2).drop(1).drop(0) == ZZ

    pytest.raises(ValueError, lambda: R.drop(3))
    pytest.raises(ValueError, lambda: R.drop(x).drop(y))


def test_PolynomialRing_index():
    R,  x, y, z = ring("x,y,z", ZZ)

    assert R.index(0) == 0
    assert R.index(-1) == 2
    pytest.raises(ValueError, lambda: R.index(100))

    assert R.index(x) == 0
    assert R.index(y) == 1
    pytest.raises(ValueError, lambda: R.index(x + y))


def test_PolynomialRing_is_():
    R = QQ.poly_ring("x")

    assert R.is_univariate is True
    assert R.is_multivariate is False

    R = QQ.poly_ring("x", "y", "z")

    assert R.is_univariate is False
    assert R.is_multivariate is True


def test_PolynomialRing_add():
    R, x = ring("x", ZZ)
    F = [x**2 + 2*i + 3 for i in range(4)]

    assert functools.reduce(operator.add, F) == 4*x**2 + 24


def test_PolynomialRing_mul():
    R, x = ring("x", ZZ)
    F = [x**2 + 2*i + 3 for i in range(4)]

    assert functools.reduce(operator.mul, F) == (x**8 + 24*x**6 +
                                                 206*x**4 + 744*x**2 + 945)


def test_PolynomialRing_to_ground():
    R, x = ring("x", ZZ)

    pytest.raises(ValueError, lambda: R.to_ground())

    R2, x, y = ring("x,y", ZZ)
    assert R2.drop_to_ground(x) == ZZ.poly_ring("x").poly_ring("y")
    assert R2.drop_to_ground(x, y) == R2


def test_sring():
    x, y, z, t = symbols("x,y,z,t")

    R = ZZ.poly_ring("x", "y", "z")
    assert sring(x + 2*y + 3*z) == (R, R.x + 2*R.y + 3*R.z)

    R = QQ.poly_ring("x", "y", "z")
    assert sring(x + 2*y + z/3) == (R, R.x + 2*R.y + R.z/3)
    assert sring([x, 2*y, z/3]) == (R, [R.x, 2*R.y, R.z/3])

    Rt = ZZ.poly_ring("t")
    R = Rt.poly_ring("x", "y", "z")
    assert sring(x + 2*t*y + 3*t**2*z, x, y, z) == (R, R.x + 2*Rt.t*R.y + 3*Rt.t**2*R.z)

    Rt = QQ.poly_ring("t")
    R = Rt.poly_ring("x", "y", "z")
    assert sring(x + t*y/2 + t**2*z/3, x, y, z) == (R, R.x + Rt.t*R.y/2 + Rt.t**2*R.z/3)

    Rt = ZZ.frac_field("t")
    R = Rt.poly_ring("x", "y", "z")
    assert sring(x + 2*y/t + t**2*z/3, x, y, z) == (R, R.x + 2*R.y/Rt.t + Rt.t**2*R.z/3)

    R = QQ.poly_ring("x", "y")
    assert sring(x + y, domain=QQ) == (R, R.x + R.y)


def test_PolyElement___hash__():
    R,  x, y, z = ring("x,y,z", QQ)
    assert hash(x*y*z)


def test_PolyElement___eq__():
    R,  x, y = ring("x,y", ZZ)

    assert ((x*y + 5*x*y) == 6) is False
    assert ((x*y + 5*x*y) == 6*x*y) is True
    assert (6 == (x*y + 5*x*y)) is False
    assert (6*x*y == (x*y + 5*x*y)) is True

    assert ((x*y - x*y) == 0) is True
    assert (0 == (x*y - x*y)) is True

    assert ((x*y - x*y) == 1) is False
    assert (1 == (x*y - x*y)) is False

    assert ((x*y - x*y) == 1) is False
    assert (1 == (x*y - x*y)) is False

    assert ((x*y + 5*x*y) != 6) is True
    assert ((x*y + 5*x*y) != 6*x*y) is False
    assert (6 != (x*y + 5*x*y)) is True
    assert (6*x*y != (x*y + 5*x*y)) is False

    assert ((x*y - x*y) != 0) is False
    assert (0 != (x*y - x*y)) is False

    assert ((x*y - x*y) != 1) is True
    assert (1 != (x*y - x*y)) is True

    Rt, t = ring("t", ZZ)
    R,  x, y = ring("x,y", Rt)

    assert (t**3*x//x == t**3) is True
    assert (t**3*x//x == t**4) is False


def test_PolyElement__lt_le_gt_ge__():
    R, x, y = ring("x,y", ZZ)

    assert R(1) < x < x**2 < x**3
    assert R(1) <= x <= x**2 <= x**3

    assert x**3 > x**2 > x > R(1)
    assert x**3 >= x**2 >= x >= R(1)

    pytest.raises(TypeError, lambda: 1 < x)


def test_PolyElement_copy():
    R,  x, y, z = ring("x,y,z", ZZ)

    f = x*y + 3*z
    g = f.copy()

    assert f == g
    g[(1, 1, 1)] = 7
    assert f != g


def test_PolyElement_items():
    R,  x, y, z = ring("x,y,z", ZZ)

    f = x*y + 3*z
    assert list(f.items()) == [((1, 1, 0), 1), ((0, 0, 1), 3)]


def test_PolyElement_as_expr():
    R,  x, y, z = ring("x,y,z", ZZ)
    f = 3*x**2*y - x*y*z + 7*z**3 + 1

    X, Y, Z = R.symbols
    g = 3*X**2*Y - X*Y*Z + 7*Z**3 + 1

    assert f != g
    assert f.as_expr() == g

    X, Y, Z = symbols("x,y,z")
    g = 3*X**2*Y - X*Y*Z + 7*Z**3 + 1

    assert f != g
    assert f.as_expr(X, Y, Z) == g

    pytest.raises(ValueError, lambda: f.as_expr(X))


def test_PolyElement_from_expr():
    x, y, z = symbols("x,y,z")
    R, X, Y, Z = ring((x, y, z), ZZ)

    f = R.convert(1)
    assert f == 1 and isinstance(f, R.dtype)

    f = R.convert(x)
    assert f == X and isinstance(f, R.dtype)

    f = R.convert(x*y*z)
    assert f == X*Y*Z and isinstance(f, R.dtype)

    f = R.convert(x*y*z + x*y + x)
    assert f == X*Y*Z + X*Y + X and isinstance(f, R.dtype)

    f = R.convert(x**3*y*z + x**2*y**7 + 1)
    assert f == X**3*Y*Z + X**2*Y**7 + 1 and isinstance(f, R.dtype)

    pytest.raises(CoercionFailed, lambda: R.convert(1/x))
    pytest.raises(CoercionFailed, lambda: R.convert(2**x))
    pytest.raises(CoercionFailed, lambda: R.convert(7*x + sqrt(2)))

    R,  X, Y = ring((2**x, y), ZZ)
    f = R.convert(2**(2*x) + 1)
    assert f == X**2 + 1


def test_PolyElement_degree():
    R,  x, y, z = ring("x,y,z", ZZ)

    assert R(0).degree() == -oo
    assert R(1).degree() == 0
    assert (x + 1).degree() == 1
    assert (2*y**3 + z).degree() == 0
    assert (x*y**3 + z).degree() == 1
    assert (x**5*y**3 + z).degree() == 5

    assert (x**5*y**3 + z).degree(0) == 5
    assert (x**5*y**3 + z).degree(-3) == 5
    pytest.raises(ValueError, lambda: (x**5*y**3 + z).degree(100))

    assert R(0).degree(x) == -oo
    assert R(1).degree(x) == 0
    assert (x + 1).degree(x) == 1
    assert (2*y**3 + z).degree(x) == 0
    assert (x*y**3 + z).degree(x) == 1
    assert (7*x**5*y**3 + z).degree(x) == 5

    assert R(0).degree(y) == -oo
    assert R(1).degree(y) == 0
    assert (x + 1).degree(y) == 0
    assert (2*y**3 + z).degree(y) == 3
    assert (x*y**3 + z).degree(y) == 3
    assert (7*x**5*y**3 + z).degree(y) == 3

    assert R(0).degree(z) == -oo
    assert R(1).degree(z) == 0
    assert (x + 1).degree(z) == 0
    assert (2*y**3 + z).degree(z) == 1
    assert (x*y**3 + z).degree(z) == 1
    assert (7*x**5*y**3 + z).degree(z) == 1


def test_PolyElement_tail_degree():
    R,  x, y, z = ring("x,y,z", ZZ)

    assert R(0).tail_degree() == -oo
    assert R(1).tail_degree() == 0
    assert (x + 1).tail_degree() == 0
    assert (2*y**3 + x**3*z).tail_degree() == 0
    assert (x*y**3 + x**3*z).tail_degree() == 1
    assert (x**5*y**3 + x**3*z).tail_degree() == 3

    assert R(0).tail_degree(x) == -oo
    assert R(1).tail_degree(x) == 0
    assert (x + 1).tail_degree(x) == 0
    assert (2*y**3 + x**3*z).tail_degree(x) == 0
    assert (x*y**3 + x**3*z).tail_degree(x) == 1
    assert (7*x**5*y**3 + x**3*z).tail_degree(x) == 3

    assert R(0).tail_degree(y) == -oo
    assert R(1).tail_degree(y) == 0
    assert (x + 1).tail_degree(y) == 0
    assert (2*y**3 + x**3*z).tail_degree(y) == 0
    assert (x*y**3 + x**3*z).tail_degree(y) == 0
    assert (7*x**5*y**3 + x**3*z).tail_degree(y) == 0

    assert R(0).tail_degree(z) == -oo
    assert R(1).tail_degree(z) == 0
    assert (x + 1).tail_degree(z) == 0
    assert (2*y**3 + x**3*z).tail_degree(z) == 0
    assert (x*y**3 + x**3*z).tail_degree(z) == 0
    assert (7*x**5*y**3 + x**3*z).tail_degree(z) == 0


def test_PolyElement_degrees():
    R,  x, y, z = ring("x,y,z", ZZ)

    assert R(0).degrees() == (-oo, -oo, -oo)
    assert R(1).degrees() == (0, 0, 0)
    assert (x**2*y + x**3*z**2).degrees() == (3, 1, 2)


def test_PolyElement_tail_degrees():
    R,  x, y, z = ring("x,y,z", ZZ)

    assert R(0).tail_degrees() == (-oo, -oo, -oo)
    assert R(1).tail_degrees() == (0, 0, 0)
    assert (x**2*y + x**3*z**2).tail_degrees() == (2, 0, 0)


def test_PolyElement_coeff():
    R,  x, y, z = ring("x,y,z", ZZ)
    f = 3*x**2*y - x*y*z + 7*z**3 + 23

    assert f.coeff(1) == 23
    pytest.raises(ValueError, lambda: f.coeff(3))

    assert f.coeff(x) == 0
    assert f.coeff(y) == 0
    assert f.coeff(z) == 0

    assert f.coeff(x**2*y) == 3
    assert f.coeff(x*y*z) == -1
    assert f.coeff(z**3) == 7

    pytest.raises(ValueError, lambda: f.coeff(3*x**2*y))
    pytest.raises(ValueError, lambda: f.coeff(-x*y*z))
    pytest.raises(ValueError, lambda: f.coeff(7*z**3))
    pytest.raises(ValueError, lambda: f.coeff(x + y))

    f = 2*x + 3*x*y + 4*z + 5
    assert f.coeff(1) == R.domain(5)


def test_PolyElement_LC():
    R,  x, y = ring("x,y", QQ)
    assert R(0).LC == QQ(0)
    assert (x/2).LC == QQ(1, 2)
    assert (x*y/4 + x/2).LC == QQ(1, 4)


def test_PolyElement_LM():
    R,  x, y = ring("x,y", QQ)
    assert R(0).LM == (0, 0)
    assert (x/2).LM == (1, 0)
    assert (x*y/4 + x/2).LM == (1, 1)


def test_PolyElement_LT():
    R,  x, y = ring("x,y", QQ)
    assert R(0).LT == ((0, 0), QQ(0))
    assert (x/2).LT == ((1, 0), QQ(1, 2))
    assert (x*y/4 + x/2).LT == ((1, 1), QQ(1, 4))


def test_PolyElement_leading_monom():
    R,  x, y = ring("x,y", QQ)
    assert R(0).leading_monom() == 0
    assert (x/2).leading_monom() == x
    assert (x*y/4 + x/2).leading_monom() == x*y


def test_PolyElement_leading_term():
    R,  x, y = ring("x,y", QQ)
    assert R(0).leading_term() == 0
    assert (x/2).leading_term() == x/2
    assert (x*y/4 + x/2).leading_term() == x*y/4


def test_PolyElement_terms():
    R,  x, y, z = ring("x,y,z", QQ)
    terms = (x**2/3 + y**3/4 + z**4/5).terms()
    assert terms == [((2, 0, 0), QQ(1, 3)), ((0, 3, 0), QQ(1, 4)), ((0, 0, 4), QQ(1, 5))]

    R,  x, y = ring("x,y", ZZ)
    f = x*y**7 + 2*x**2*y**3

    assert f.terms() == f.terms(lex) == f.terms('lex') == [((2, 3), 2), ((1, 7), 1)]
    assert f.terms(grlex) == f.terms('grlex') == [((1, 7), 1), ((2, 3), 2)]

    R,  x, y = ring("x,y", ZZ, grlex)
    f = x*y**7 + 2*x**2*y**3

    assert f.terms() == f.terms(grlex) == f.terms('grlex') == [((1, 7), 1), ((2, 3), 2)]
    assert f.terms(lex) == f.terms('lex') == [((2, 3), 2), ((1, 7), 1)]


def test_PolyElement_monoms():
    R,  x, y, z = ring("x,y,z", QQ)
    monoms = (x**2/3 + y**3/4 + z**4/5).monoms()
    assert monoms == [(2, 0, 0), (0, 3, 0), (0, 0, 4)]

    R,  x, y = ring("x,y", ZZ)
    f = x*y**7 + 2*x**2*y**3

    assert f.monoms() == f.monoms(lex) == f.monoms('lex') == [(2, 3), (1, 7)]
    assert f.monoms(grlex) == f.monoms('grlex') == [(1, 7), (2, 3)]

    R,  x, y = ring("x,y", ZZ, grlex)
    f = x*y**7 + 2*x**2*y**3

    assert f.monoms() == f.monoms(grlex) == f.monoms('grlex') == [(1, 7), (2, 3)]
    assert f.monoms(lex) == f.monoms('lex') == [(2, 3), (1, 7)]


def test_PolyElement_coeffs():
    R,  x, y, z = ring("x,y,z", QQ)
    coeffs = (x**2/3 + y**3/4 + z**4/5).coeffs()
    assert coeffs == [QQ(1, 3), QQ(1, 4), QQ(1, 5)]

    R,  x, y = ring("x,y", ZZ)
    f = x*y**7 + 2*x**2*y**3

    assert f.coeffs() == f.coeffs(lex) == f.coeffs('lex') == [2, 1]
    assert f.coeffs(grlex) == f.coeffs('grlex') == [1, 2]

    R,  x, y = ring("x,y", ZZ, grlex)
    f = x*y**7 + 2*x**2*y**3

    assert f.coeffs() == f.coeffs(grlex) == f.coeffs('grlex') == [1, 2]
    assert f.coeffs(lex) == f.coeffs('lex') == [2, 1]


def test_PolyElement___add__():
    Rt, t = ring("t", ZZ)
    Ruv,  u, v = ring("u,v", ZZ)
    Rxyz,  x, y, z = ring("x,y,z", Ruv)

    assert dict(+x) == dict(x)

    assert dict(x + 3*y) == {(1, 0, 0): 1, (0, 1, 0): 3}

    assert dict(u + x) == dict(x + u) == {(1, 0, 0): 1, (0, 0, 0): u}
    assert dict(u + x*y) == dict(x*y + u) == {(1, 1, 0): 1, (0, 0, 0): u}
    assert dict(u + x*y + z) == dict(x*y + z + u) == {(1, 1, 0): 1, (0, 0, 1): 1, (0, 0, 0): u}

    assert dict(u*x + x) == dict(x + u*x) == {(1, 0, 0): u + 1}
    assert dict(u*x + x*y) == dict(x*y + u*x) == {(1, 1, 0): 1, (1, 0, 0): u}
    assert dict(u*x + x*y + z) == dict(x*y + z + u*x) == {(1, 1, 0): 1, (0, 0, 1): 1, (1, 0, 0): u}

    pytest.raises(TypeError, lambda: t + x)
    pytest.raises(TypeError, lambda: x + t)
    pytest.raises(TypeError, lambda: t + u)
    pytest.raises(TypeError, lambda: u + t)

    Fuv,  u, v = field("u,v", ZZ)
    Rxyz,  x, y, z = ring("x,y,z", Fuv)

    assert u + (x - u) == x
    assert dict(u + x) == dict(x + u) == {(1, 0, 0): 1, (0, 0, 0): u}

    Rxyz,  x, y, z = ring("x,y,z", EX)

    assert dict(EX(pi) + x*y*z) == dict(x*y*z + EX(pi)) == {(1, 1, 1): EX(1), (0, 0, 0): EX(pi)}

    R, x, y = ring('x, y', ZZ)

    p = x**4 + 2*y
    m = (1, 2)
    p1 = p._iadd_monom((m, 5))
    assert p == p1 and p1 == x**4 + 5*x*y**2 + 2*y

    p2 = p._iadd_monom(((0, 1), 2))
    assert p == p2 and p2 == x**4 + 5*x*y**2 + 4*y

    p3 = p._iadd_monom(((0, 1), -4))
    assert p == p3 and p3 == x**4 + 5*x*y**2


def test_PolyElement___sub__():
    Rt, t = ring("t", ZZ)
    Ruv,  u, v = ring("u,v", ZZ)
    Rxyz,  x, y, z = ring("x,y,z", Ruv)

    assert u - x == -x + u
    assert (x + u) - 2*u == x - u

    assert dict(x - 3*y) == {(1, 0, 0): 1, (0, 1, 0): -3}

    assert dict(-u + x) == dict(x - u) == {(1, 0, 0): 1, (0, 0, 0): -u}
    assert dict(-u + x*y) == dict(x*y - u) == {(1, 1, 0): 1, (0, 0, 0): -u}
    assert dict(-u + x*y + z) == dict(x*y + z - u) == {(1, 1, 0): 1, (0, 0, 1): 1, (0, 0, 0): -u}

    assert dict(-u*x + x) == dict(x - u*x) == {(1, 0, 0): -u + 1}
    assert dict(-u*x + x*y) == dict(x*y - u*x) == {(1, 1, 0): 1, (1, 0, 0): -u}
    assert dict(-u*x + x*y + z) == dict(x*y + z - u*x) == {(1, 1, 0): 1, (0, 0, 1): 1, (1, 0, 0): -u}

    pytest.raises(TypeError, lambda: t - x)
    pytest.raises(TypeError, lambda: x - t)
    pytest.raises(TypeError, lambda: t - u)
    pytest.raises(TypeError, lambda: u - t)

    Fuv,  u, v = field("u,v", ZZ)
    Rxyz,  x, y, z = ring("x,y,z", Fuv)

    assert dict(-u + x) == dict(x - u) == {(1, 0, 0): 1, (0, 0, 0): -u}

    Rxyz,  x, y, z = ring("x,y,z", EX)

    assert dict(-EX(pi) + x*y*z) == dict(x*y*z - EX(pi)) == {(1, 1, 1): EX(1), (0, 0, 0): -EX(pi)}


def test_PolyElement___mul__():
    Rt, t = ring("t", ZZ)
    Ruv,  u, v = ring("u,v", ZZ)
    Rxyz,  x, y, z = ring("x,y,z", Ruv)

    assert dict(u*x) == dict(x*u) == {(1, 0, 0): u}

    assert dict(2*u*x + z) == dict(x*2*u + z) == {(1, 0, 0): 2*u, (0, 0, 1): 1}
    assert dict(u*2*x + z) == dict(2*x*u + z) == {(1, 0, 0): 2*u, (0, 0, 1): 1}
    assert dict(2*u*x + z) == dict(x*2*u + z) == {(1, 0, 0): 2*u, (0, 0, 1): 1}
    assert dict(u*x*2 + z) == dict(x*u*2 + z) == {(1, 0, 0): 2*u, (0, 0, 1): 1}

    assert dict(2*u*x*y + z) == dict(x*y*2*u + z) == {(1, 1, 0): 2*u, (0, 0, 1): 1}
    assert dict(u*2*x*y + z) == dict(2*x*y*u + z) == {(1, 1, 0): 2*u, (0, 0, 1): 1}
    assert dict(2*u*x*y + z) == dict(x*y*2*u + z) == {(1, 1, 0): 2*u, (0, 0, 1): 1}
    assert dict(u*x*y*2 + z) == dict(x*y*u*2 + z) == {(1, 1, 0): 2*u, (0, 0, 1): 1}

    assert dict(2*u*y*x + z) == dict(y*x*2*u + z) == {(1, 1, 0): 2*u, (0, 0, 1): 1}
    assert dict(u*2*y*x + z) == dict(2*y*x*u + z) == {(1, 1, 0): 2*u, (0, 0, 1): 1}
    assert dict(2*u*y*x + z) == dict(y*x*2*u + z) == {(1, 1, 0): 2*u, (0, 0, 1): 1}
    assert dict(u*y*x*2 + z) == dict(y*x*u*2 + z) == {(1, 1, 0): 2*u, (0, 0, 1): 1}

    assert dict(3*u*(x + y) + z) == dict((x + y)*3*u + z) == {(1, 0, 0): 3*u, (0, 1, 0): 3*u, (0, 0, 1): 1}

    pytest.raises(TypeError, lambda: t*x + z)
    pytest.raises(TypeError, lambda: x*t + z)
    pytest.raises(TypeError, lambda: t*u + z)
    pytest.raises(TypeError, lambda: u*t + z)

    Fuv,  u, v = field("u,v", ZZ)
    Rxyz,  x, y, z = ring("x,y,z", Fuv)

    assert dict(u*x) == dict(x*u) == {(1, 0, 0): u}

    Rxyz,  x, y, z = ring("x,y,z", EX)

    assert dict(EX(pi)*x*y*z) == dict(x*y*z*EX(pi)) == {(1, 1, 1): EX(pi)}

    assert (x + 2*y).mul_ground(0) == Rxyz.zero

    R, x, y = ring("x,y", ZZ)
    p = x + y**2
    p1 = p._imul_num(3)
    assert p == p1 and p1 == 3*x + 3*y**2
    p2 = p._imul_num(0)
    assert p == p2 and p2 == R.zero


def test_PolyElement___floordiv__truediv__():
    R,  x, y, z = ring("x,y,z", ZZ)

    assert (2*x**2 - 4)/2 == x**2 - 2
    assert (2*x**2 - 3)/2 == x**2 - 2

    assert (x**2 - 1)//x == x
    assert (x**2 - x)//x == x - 1

    assert (x**2 - 1)//(2*x) == 0
    assert (x**2 - x)//(x - 1) == x

    R,  x, y, z = ring("x,y,z", ZZ)
    assert len((x**2/3 + y**3/4 + z**4/5).terms()) == 0

    R,  x, y, z = ring("x,y,z", QQ)
    assert len((x**2/3 + y**3/4 + z**4/5).terms()) == 3

    pytest.raises(ZeroDivisionError, lambda: x/0)

    Rt, t = ring("t", ZZ)
    Ruv,  u, v = ring("u,v", ZZ)
    Rxyz,  x, y, z = ring("x,y,z", Ruv)

    assert dict((u**2*x + u)/u) == {(1, 0, 0): u, (0, 0, 0): 1}
    pytest.raises(TypeError, lambda: u//(u**2*x + u))

    pytest.raises(TypeError, lambda: t//x)
    pytest.raises(TypeError, lambda: x//t)
    pytest.raises(TypeError, lambda: t//u)
    pytest.raises(TypeError, lambda: u//t)

    assert divmod(x, u) == (0, x)
    assert x % u == x
    assert x // u == 0
    pytest.raises(TypeError, lambda: divmod(u, x))
    pytest.raises(TypeError, lambda: u % x)
    pytest.raises(TypeError, lambda: u // x)
    pytest.raises(TypeError, lambda: divmod(u, t))
    pytest.raises(TypeError, lambda: u % t)
    pytest.raises(TypeError, lambda: u // t)
    pytest.raises(TypeError, lambda: divmod(u, sqrt(2)))
    pytest.raises(TypeError, lambda: u % sqrt(2))
    pytest.raises(TypeError, lambda: u // sqrt(2))

    R, x = ring("x", ZZ)
    f, g = x**2 + 2*x + 3, R(0)

    pytest.raises(ZeroDivisionError, lambda: divmod(f, g))
    pytest.raises(ZeroDivisionError, lambda: f % g)
    pytest.raises(ZeroDivisionError, lambda: f // g)
    pytest.raises(ZeroDivisionError, lambda: f.exquo(g))

    R,  x, y = ring("x,y", ZZ)
    f, g = x*y + 2*x + 3, R(0)

    pytest.raises(ZeroDivisionError, lambda: divmod(f, g))
    pytest.raises(ZeroDivisionError, lambda: f % g)
    pytest.raises(ZeroDivisionError, lambda: f // g)
    pytest.raises(ZeroDivisionError, lambda: f.exquo(g))

    R, x = ring("x", ZZ)

    f, g = x**2 + 1, 2*x - 4
    q, r = R(0), x**2 + 1

    assert divmod(f, g) == (q, r)
    assert f % g == r
    assert f // g == q
    pytest.raises(ExactQuotientFailed, lambda: f.exquo(g))

    assert divmod(R.zero, f) == (R.zero, R.zero)

    f, g = 3*x**3 + x**2 + x + 5, 5*x**2 - 3*x + 1
    q, r = R(0), f

    assert divmod(f, g) == (q, r)
    assert f % g == r
    assert f // g == q
    pytest.raises(ExactQuotientFailed, lambda: f.exquo(g))

    f, g = 5*x**4 + 4*x**3 + 3*x**2 + 2*x + 1, x**2 + 2*x + 3
    q, r = 5*x**2 - 6*x, 20*x + 1

    assert divmod(f, g) == (q, r)
    assert f % g == r
    assert f // g == q
    pytest.raises(ExactQuotientFailed, lambda: f.exquo(g))

    f, g = 5*x**5 + 4*x**4 + 3*x**3 + 2*x**2 + x, x**4 + 2*x**3 + 9
    q, r = 5*x - 6, 15*x**3 + 2*x**2 - 44*x + 54

    assert divmod(f, g) == (q, r)
    assert f % g == r
    assert f // g == q
    pytest.raises(ExactQuotientFailed, lambda: f.exquo(g))

    R, x = ring("x", QQ)

    f, g = x**2 + 1, 2*x - 4
    q, r = x/2 + 1, R(5)

    assert divmod(f, g) == (q, r)
    assert f % g == r
    assert f // g == q
    pytest.raises(ExactQuotientFailed, lambda: f.exquo(g))

    f, g = 3*x**3 + x**2 + x + 5, 5*x**2 - 3*x + 1
    q, r = 3*x/5 + QQ(14, 25), 52*x/25 + QQ(111, 25)

    assert divmod(f, g) == (q, r)
    assert f % g == r
    assert f // g == q
    pytest.raises(ExactQuotientFailed, lambda: f.exquo(g))

    R,  x, y = ring("x,y", ZZ)

    f, g = x**2 - y**2, x - y
    q, r = x + y, R(0)

    assert divmod(f, g) == (q, r)
    assert f % g == r
    assert f // g == q
    assert f.exquo(g) == q

    f, g = x**2 + y**2, x - y
    q, r = x + y, 2*y**2

    assert divmod(f, g) == (q, r)
    assert f % g == r
    assert f // g == q
    pytest.raises(ExactQuotientFailed, lambda: f.exquo(g))

    f, g = x**2 + y**2, -x + y
    q, r = -x - y, 2*y**2

    assert divmod(f, g) == (q, r)
    assert f % g == r
    assert f // g == q
    pytest.raises(ExactQuotientFailed, lambda: f.exquo(g))

    f, g = x**2 + y**2, 2*x - 2*y
    q, r = R(0), f

    assert divmod(f, g) == (q, r)
    assert f % g == r
    assert f // g == q
    pytest.raises(ExactQuotientFailed, lambda: f.exquo(g))

    R,  x, y = ring("x,y", QQ)

    f, g = x**2 - y**2, x - y
    q, r = x + y, R(0)

    assert divmod(f, g) == (q, r)
    assert f % g == r
    assert f // g == q
    assert f.exquo(g) == q

    f, g = x**2 + y**2, x - y
    q, r = x + y, 2*y**2

    assert divmod(f, g) == (q, r)
    assert f % g == r
    assert f // g == q
    pytest.raises(ExactQuotientFailed, lambda: f.exquo(g))

    f, g = x**2 + y**2, -x + y
    q, r = -x - y, 2*y**2

    assert divmod(f, g) == (q, r)
    assert f % g == r
    assert f // g == q
    pytest.raises(ExactQuotientFailed, lambda: f.exquo(g))

    f, g = x**2 + y**2, 2*x - 2*y
    q, r = x/2 + y/2, 2*y**2

    assert divmod(f, g) == (q, r)
    assert f % g == r
    assert f // g == q
    pytest.raises(ExactQuotientFailed, lambda: f.exquo(g))

    pytest.raises(ZeroDivisionError, lambda: f.quo_ground(0))
    pytest.raises(ZeroDivisionError, lambda: f.quo_term(((1, 1), 0)))
    pytest.raises(ZeroDivisionError, lambda: f.exquo_ground(0))

    assert R.zero.exquo_ground(2) == R.zero

    assert R.zero.quo_term(((1, 0), 1)) == R.zero
    assert g.quo_term((R.zero_monom, 2)) == x - y
    assert f.quo_term(((1, 0), 2)) == x/2


def test_PolyElement___pow__():
    R, x = ring("x", ZZ, grlex)
    f = 2*x + 3

    pytest.raises(ValueError, lambda: R.zero**0)

    assert f**0 == 1
    assert f**1 == f

    assert f**2 == f._pow_generic(2) == f._pow_multinomial(2) == 4*x**2 + 12*x + 9
    assert f**3 == f._pow_generic(3) == f._pow_multinomial(3) == 8*x**3 + 36*x**2 + 54*x + 27
    assert f**4 == f._pow_generic(4) == f._pow_multinomial(4) == 16*x**4 + 96*x**3 + 216*x**2 + 216*x + 81
    assert f**5 == f._pow_generic(5) == f._pow_multinomial(5) == 32*x**5 + 240*x**4 + 720*x**3 + 1080*x**2 + 810*x + 243

    pytest.raises(ValueError, lambda: f**-2)

    f = x**2 - 2*x + x**3 + 1
    assert f**5 == (x**15 + 5*x**14 - 25*x**12 + 5*x**11 + 71*x**10 -
                    60*x**9 - 85*x**8 + 170*x**7 - 60*x**6 - 112*x**5 +
                    170*x**4 - 115*x**3 + 45*x**2 - 10*x + 1)

    R,  x, y, z = ring("x,y,z", ZZ, grlex)
    f = x**3*y - 2*x*y**2 - 3*z + 1
    g = x**6*y**2 - 4*x**4*y**3 - 6*x**3*y*z + 2*x**3*y + 4*x**2*y**4 + 12*x*y**2*z - 4*x*y**2 + 9*z**2 - 6*z + 1

    assert f**2 == f._pow_generic(2) == f._pow_multinomial(2) == g

    R, t = ring("t", ZZ)
    f = -11200*t**4 - 2604*t**2 + 49
    g = 15735193600000000*t**16 + 14633730048000000*t**14 + 4828147466240000*t**12 \
        + 598976863027200*t**10 + 3130812416256*t**8 - 2620523775744*t**6 \
        + 92413760096*t**4 - 1225431984*t**2 + 5764801

    assert f**4 == f._pow_generic(4) == f._pow_multinomial(4) == g

    f = x**3*y - 2*x*y**2 - 3*z + x + y + 1
    assert f**4 == f._pow_generic(4)


def test_PolyElement_div():
    _, u = ring("u", ZZ)
    R, t = ring("t", ZZ)

    pytest.raises(ValueError, lambda: u.div(t))
    pytest.raises(ZeroDivisionError, lambda: R.one.div(R.zero))

    assert R.zero.div([t**2 + 1]) == ([], R.zero)

    R, x = ring("x", ZZ, grlex)

    f = x**3 - 12*x**2 - 42
    g = x - 3

    q = x**2 - 9*x - 27
    r = -123

    assert f.div([g]) == ([q], r)

    R, x = ring("x", ZZ, grlex)
    f = x**2 + 2*x + 2
    assert f.div([R(1)]) == ([f], 0)

    R, x = ring("x", QQ, grlex)
    f = x**2 + 2*x + 2
    assert f.div([R(2)]) == ([x**2/2 + x + 1], 0)

    R,  x, y = ring("x,y", ZZ, grlex)
    f = 4*x**2*y - 2*x*y + 4*x - 2*y + 8

    assert f.div([R(2)]) == ([2*x**2*y - x*y + 2*x - y + 4], 0)
    assert f.div([2*y]) == ([2*x**2 - x - 1], 4*x + 8)

    f = x - 1
    g = y - 1

    assert f.div([g]) == ([0], f)

    f = x*y**2 + 1
    G = [x*y + 1, y + 1]

    Q = [y, -1]
    r = 2

    assert f.div(G) == (Q, r)

    f = x**2*y + x*y**2 + y**2
    G = [x*y - 1, y**2 - 1]

    Q = [x + y, 1]
    r = x + y + 1

    assert f.div(G) == (Q, r)

    G = [y**2 - 1, x*y - 1]

    Q = [x + 1, x]
    r = 2*x + 1

    assert f.div(G) == (Q, r)


def test_PolyElement_monic():
    R, x = ring("x", ZZ)

    assert (2*x + 2).monic() == x + 1

    pytest.raises(ExactQuotientFailed, lambda: (2*x + 1).monic())

    R, x = ring("x", QQ)

    assert (2*x + 1).monic() == x + QQ(1, 2)


def test_PolyElement_deflate():
    R, x = ring("x", ZZ)

    assert (2*x**2).deflate(x**4 + 4*x**2 + 1) == ((2,), [2*x, x**2 + 4*x + 1])

    R,  x, y = ring("x,y", ZZ)

    assert R(0).deflate(R(0)) == ((1, 1), [0, 0])
    assert R(1).deflate(R(0)) == ((1, 1), [1, 0])
    assert R(1).deflate(R(2)) == ((1, 1), [1, 2])
    assert R(1).deflate(2*y) == ((1, 1), [1, 2*y])
    assert (2*y).deflate(2*y) == ((1, 1), [2*y, 2*y])
    assert R(2).deflate(2*y**2) == ((1, 2), [2, 2*y])
    assert (2*y**2).deflate(2*y**2) == ((1, 2), [2*y, 2*y])

    f = x**4*y**2 + x**2*y + 1
    g = x**2*y**3 + x**2*y + 1

    assert f.deflate(g) == ((2, 1), [x**2*y**2 + x*y + 1, x*y**3 + x*y + 1])


def test_PolyElement_clear_denoms():
    R,  x, y = ring("x,y", QQ)

    assert R(1).clear_denoms() == (ZZ(1), 1)
    assert R(7).clear_denoms() == (ZZ(1), 7)

    assert R(QQ(7, 3)).clear_denoms() == (3, 7)
    assert R(QQ(7, 3)).clear_denoms() == (3, 7)

    assert (3*x**2 + x).clear_denoms() == (1, 3*x**2 + x)
    assert (x**2 + x/2).clear_denoms() == (2, 2*x**2 + x)

    rQQ,  x, t = ring("x,t", QQ)
    rZZ,  X, T = ring("x,t", ZZ)

    F = [x - QQ(17824537287975195925064602467992950991718052713078834557692023531499318507213727406844943097, 413954288007559433755329699713866804710749652268151059918115348815925474842910720000)*t**7
           - QQ(4882321164854282623427463828745855894130208215961904469205260756604820743234704900167747753, 12936071500236232304854053116058337647210926633379720622441104650497671088840960000)*t**6
           - QQ(36398103304520066098365558157422127347455927422509913596393052633155821154626830576085097433, 25872143000472464609708106232116675294421853266759441244882209300995342177681920000)*t**5
           - QQ(168108082231614049052707339295479262031324376786405372698857619250210703675982492356828810819, 58212321751063045371843239022262519412449169850208742800984970927239519899784320000)*t**4
           - QQ(5694176899498574510667890423110567593477487855183144378347226247962949388653159751849449037, 1617008937529529038106756639507292205901365829172465077805138081312208886105120000)*t**3
           - QQ(154482622347268833757819824809033388503591365487934245386958884099214649755244381307907779, 60637835157357338929003373981523457721301218593967440417692678049207833228942000)*t**2
           - QQ(2452813096069528207645703151222478123259511586701148682951852876484544822947007791153163, 2425513406294293557160134959260938308852048743758697616707707121968313329157680)*t
           - QQ(34305265428126440542854669008203683099323146152358231964773310260498715579162112959703, 202126117191191129763344579938411525737670728646558134725642260164026110763140),
         t**8 + QQ(693749860237914515552, 67859264524169150569)*t**7
              + QQ(27761407182086143225024, 610733380717522355121)*t**6
              + QQ(7785127652157884044288, 67859264524169150569)*t**5
              + QQ(36567075214771261409792, 203577793572507451707)*t**4
              + QQ(36336335165196147384320, 203577793572507451707)*t**3
              + QQ(7452455676042754048000, 67859264524169150569)*t**2
              + QQ(2593331082514399232000, 67859264524169150569)*t
              + QQ(390399197427343360000, 67859264524169150569)]

    G = [3725588592068034903797967297424801242396746870413359539263038139343329273586196480000*X -
         160420835591776763325581422211936558925462474417709511019228211783493866564923546661604487873*T**7 -
         1406108495478033395547109582678806497509499966197028487131115097902188374051595011248311352864*T**6 -
         5241326875850889518164640374668786338033653548841427557880599579174438246266263602956254030352*T**5 -
         10758917262823299139373269714910672770004760114329943852726887632013485035262879510837043892416*T**4 -
         13119383576444715672578819534846747735372132018341964647712009275306635391456880068261130581248*T**3 -
         9491412317016197146080450036267011389660653495578680036574753839055748080962214787557853941760*T**2 -
         3767520915562795326943800040277726397326609797172964377014046018280260848046603967211258368000*T -
         632314652371226552085897259159210286886724229880266931574701654721512325555116066073245696000,
         610733380717522355121*T**8 +
         6243748742141230639968*T**7 +
         27761407182086143225024*T**6 +
         70066148869420956398592*T**5 +
         109701225644313784229376*T**4 +
         109009005495588442152960*T**3 +
         67072101084384786432000*T**2 +
         23339979742629593088000*T +
         3513592776846090240000]

    assert [ f.clear_denoms()[1].set_ring(rZZ) for f in F ] == G


def test_PolyElement_cofactors():
    R,  x, y = ring("x,y", ZZ)

    f, g = R(0), R(0)
    assert f.cofactors(g) == (0, 0, 0)

    f, g = R(2), R(0)
    assert f.cofactors(g) == (2, 1, 0)

    f, g = R(-2), R(0)
    assert f.cofactors(g) == (2, -1, 0)

    f, g = R(0), R(-2)
    assert f.cofactors(g) == (2, 0, -1)

    f, g = R(0), 2*x + 4
    assert f.cofactors(g) == (2*x + 4, 0, 1)

    f, g = 2*x + 4, R(0)
    assert f.cofactors(g) == (2*x + 4, 1, 0)

    f, g = R(2), R(2)
    assert f.cofactors(g) == (2, 1, 1)

    f, g = R(-2), R(2)
    assert f.cofactors(g) == (2, -1, 1)

    f, g = R(2), R(-2)
    assert f.cofactors(g) == (2, 1, -1)

    f, g = R(-2), R(-2)
    assert f.cofactors(g) == (2, -1, -1)

    f, g = x**2 + 2*x + 1, R(1)
    assert f.cofactors(g) == (1, x**2 + 2*x + 1, 1)

    f, g = x**2 + 2*x + 1, R(2)
    assert f.cofactors(g) == (1, x**2 + 2*x + 1, 2)

    f, g = 2*x**2 + 4*x + 2, R(2)
    assert f.cofactors(g) == (2, x**2 + 2*x + 1, 1)

    f, g = R(2), 2*x**2 + 4*x + 2
    assert f.cofactors(g) == (2, 1, x**2 + 2*x + 1)

    f, g = 2*x**2 + 4*x + 2, x + 1
    assert f.cofactors(g) == (x + 1, 2*x + 2, 1)

    f, g = x + 1, 2*x**2 + 4*x + 2
    assert f.cofactors(g) == (x + 1, 1, 2*x + 2)

    R, x, y, z, t = ring("x,y,z,t", ZZ)

    f, g = t**2 + 2*t + 1, 2*t + 2
    assert f.cofactors(g) == (t + 1, t + 1, 2)

    f, g = z**2*t**2 + 2*z**2*t + z**2 + z*t + z, t**2 + 2*t + 1
    h, cff, cfg = t + 1, z**2*t + z**2 + z, t + 1

    assert f.cofactors(g) == (h, cff, cfg)
    assert g.cofactors(f) == (h, cfg, cff)

    R,  x, y = ring("x,y", QQ)

    f = x**2/2 + x + QQ(1, 2)
    g = x/2 + QQ(1, 2)

    h = x + 1

    assert f.cofactors(g) == (h, g, QQ(1, 2))
    assert g.cofactors(f) == (h, QQ(1, 2), g)

    R,  x, y = ring("x,y", RR)

    f = 2.1*x*y**2 - 2.1*x*y + 2.1*x
    g = 2.1*x**3
    h = 1.0*x

    assert f.cofactors(g) == (h, f//h, g//h)
    assert g.cofactors(f) == (h, g//h, f//h)


def test_PolyElement_gcd():
    R,  x, y = ring("x,y", QQ)

    f = x**2/2 + x + QQ(1, 2)
    g = x/2 + QQ(1, 2)

    assert f.gcd(g) == x + 1

    with using(use_heu_gcd=False):
        assert f.gcd(g) == x + 1
    with using(use_heu_gcd=False, fallback_gcd_zz_method='modgcd'):
        assert f.gcd(g) == x + 1

    R, x, y = ring("x,y", QQ.algebraic_field(sqrt(2)))

    f, g = (x + sqrt(2)*y)**2, x + sqrt(2)*y

    assert f.gcd(g) == g
    with using(gcd_aa_method='modgcd'):
        assert f.gcd(g) == g


def test_PolyElement_cancel():
    R,  x, y = ring("x,y", ZZ)

    f = 2*x**3 + 4*x**2 + 2*x
    g = 3*x**2 + 3*x
    F = 2*x + 2
    G = 3

    assert f.cancel(g) == (F, G)

    assert (-f).cancel(g) == (-F, G)
    assert f.cancel(-g) == (-F, G)

    R,  x, y = ring("x,y", QQ)

    f = x**3/2 + x**2 + x/2
    g = x**2/3 + x/3
    F = 3*x + 3
    G = 2

    assert f.cancel(g) == (F, G)

    assert (-f).cancel(g) == (-F, G)
    assert f.cancel(-g) == (-F, G)

    Fx, x = field("x", ZZ)
    Rt, t = ring("t", Fx)

    f = (-x**2 - 4)/4*t
    g = t**2 + (x**2 + 2)/2

    assert f.cancel(g) == ((-x**2 - 4)*t, 4*t**2 + 2*x**2 + 4)


def test_PolyElement_max_norm():
    R,  x, y = ring("x,y", ZZ)

    assert R(0).max_norm() == 0
    assert R(1).max_norm() == 1

    assert (x**3 + 4*x**2 + 2*x + 3).max_norm() == 4


def test_PolyElement_l1_norm():
    R,  x, y = ring("x,y", ZZ)

    assert R(0).l1_norm() == 0
    assert R(1).l1_norm() == 1

    assert (x**3 + 4*x**2 + 2*x + 3).l1_norm() == 10


def test_PolyElement_diff():
    R, *X = ring("x:11", QQ)

    f = 288*X[0]**8*X[1]**6*X[4]**3*X[10]**2/5 + 8*X[0]**2*X[2]**3*X[4]**3 + 2*X[0]**2 - 2*X[1]**2

    assert f.diff(X[0]) == 2304*X[0]**7*X[1]**6*X[4]**3*X[10]**2/5 + 16*X[0]*X[2]**3*X[4]**3 + 4*X[0]
    assert f.diff(X[4]) == 864*X[0]**8*X[1]**6*X[4]**2*X[10]**2/5 + 24*X[0]**2*X[2]**3*X[4]**2
    assert f.diff(X[10]) == 576*X[0]**8*X[1]**6*X[4]**3*X[10]/5


def test_PolyElement___call__():
    R, x = ring("x", ZZ)
    f = 3*x + 1

    assert f(0) == 1
    assert f(1) == 4

    pytest.raises(ValueError, lambda: f())
    pytest.raises(ValueError, lambda: f(0, 1))

    pytest.raises(CoercionFailed, lambda: f(QQ(1, 7)))

    R,  x, y = ring("x,y", ZZ)
    f = 3*x + y**2 + 1

    assert f(0, 0) == 1
    assert f(1, 7) == 53

    Ry = R.drop(x)

    assert f(0) == Ry.y**2 + 1
    assert f(1) == Ry.y**2 + 4

    pytest.raises(ValueError, lambda: f())
    pytest.raises(ValueError, lambda: f(0, 1, 2))

    pytest.raises(CoercionFailed, lambda: f(1, QQ(1, 7)))
    pytest.raises(CoercionFailed, lambda: f(QQ(1, 7), 1))
    pytest.raises(CoercionFailed, lambda: f(QQ(1, 7), QQ(1, 7)))


def test_PolyElement_evaluate():
    R, x = ring("x", ZZ)
    f = x**3 + 4*x**2 + 2*x + 3

    r = f.evaluate(x, 0)
    assert r == 3 and not isinstance(r, PolyElement)

    pytest.raises(CoercionFailed, lambda: f.evaluate(x, QQ(1, 7)))

    R,  x, y, z = ring("x,y,z", ZZ)
    f = (x*y)**3 + 4*(x*y)**2 + 2*x*y + 3

    r = f.evaluate(x, 0)
    assert r == 3 and isinstance(r, R.drop(x).dtype)
    r = f.evaluate([(x, 0), (y, 0)])
    assert r == 3 and isinstance(r, R.drop(x, y).dtype)
    r = f.evaluate(y, 0)
    assert r == 3 and isinstance(r, R.drop(y).dtype)
    r = f.evaluate([(y, 0), (x, 0)])
    assert r == 3 and isinstance(r, R.drop(y, x).dtype)

    r = f.evaluate([(x, 0), (y, 0), (z, 0)])
    assert r == 3 and not isinstance(r, PolyElement)

    pytest.raises(CoercionFailed, lambda: f.evaluate([(x, 1), (y, QQ(1, 7))]))
    pytest.raises(CoercionFailed, lambda: f.evaluate([(x, QQ(1, 7)), (y, 1)]))
    pytest.raises(CoercionFailed, lambda: f.evaluate([(x, QQ(1, 7)), (y, QQ(1, 7))]))


def test_PolyElement_subs():
    R, x = ring("x", ZZ)
    f = x**3 + 4*x**2 + 2*x + 3

    assert f == f.subs([])
    pytest.raises(ValueError, lambda: f.subs(object()))

    r = f.subs({x: 0})
    assert r == 3 and isinstance(r, R.dtype)
    assert f.subs({(x, 0)}) == r

    pytest.raises(CoercionFailed, lambda: f.subs({x: QQ(1, 7)}))

    R,  x, y, z = ring("x,y,z", ZZ)
    f = x**3 + 4*x**2 + 2*x + 3

    r = f.subs({x: 0})
    assert r == 3 and isinstance(r, R.dtype)
    r = f.subs({x: 0, y: 0})
    assert r == 3 and isinstance(r, R.dtype)

    pytest.raises(CoercionFailed, lambda: f.subs({x: 1, y: QQ(1, 7)}))
    pytest.raises(CoercionFailed, lambda: f.subs({x: QQ(1, 7), y: 1}))
    pytest.raises(CoercionFailed, lambda: f.subs({x: QQ(1, 7), y: QQ(1, 7)}))

    f = 1 - x - y - z
    r = f.subs({y: 1})
    assert r == -x - z


def test_PolyElement_compose():
    R, x = ring("x", ZZ)
    f = x**3 + 4*x**2 + 2*x + 3

    r = f.compose(x, 0)
    assert r == 3 and isinstance(r, R.dtype)

    assert f.compose(x, x) == f
    assert f.compose(x, x**2) == x**6 + 4*x**4 + 2*x**2 + 3

    pytest.raises(CoercionFailed, lambda: f.compose(x, QQ(1, 7)))

    R,  x, y, z = ring("x,y,z", ZZ)
    f = x**3 + 4*x**2 + 2*x + 3

    r = f.compose(x, 0)
    assert r == 3 and isinstance(r, R.dtype)
    r = f.compose([(x, 0), (y, 0)])
    assert r == 3 and isinstance(r, R.dtype)
    r = f.compose({x: 0, y: 0})
    assert r == 3 and isinstance(r, R.dtype)
    pytest.raises(ValueError, lambda: f.compose("spam"))

    r = (x**3 + 4*x**2 + 2*x*y*z + 3).compose(x, y*z**2 - 1)
    q = (y*z**2 - 1)**3 + 4*(y*z**2 - 1)**2 + 2*(y*z**2 - 1)*y*z + 3
    assert r == q and isinstance(r, R.dtype)


def test_PolyElement_is_():
    R,  x, y, z = ring("x,y,z", QQ)

    assert (x - x).is_generator is False
    assert (x - x).is_ground
    assert (x - x).is_monomial
    assert (x - x).is_term

    assert (x - x + 1).is_generator is False
    assert (x - x + 1).is_ground
    assert (x - x + 1).is_monomial
    assert (x - x + 1).is_term

    assert x.is_generator
    assert x.is_ground is False
    assert x.is_monomial
    assert x.is_term

    assert (x*y).is_generator is False
    assert (x*y).is_ground is False
    assert (x*y).is_monomial
    assert (x*y).is_term

    assert (3*x).is_generator is False
    assert (3*x).is_ground is False
    assert (3*x).is_monomial is False
    assert (3*x).is_term

    assert (3*x + 1).is_generator is False
    assert (3*x + 1).is_ground is False
    assert (3*x + 1).is_monomial is False
    assert (3*x + 1).is_term is False

    assert R(0).is_zero
    assert R(1).is_zero is False

    assert R(0).is_one is False
    assert R(1).is_one

    assert (x - 1).is_monic
    assert (2*x - 1).is_monic is False

    assert (3*x + 2).is_primitive
    assert (4*x + 2).is_primitive is False

    assert (x + y + z + 1).is_linear
    assert (x*y*z + 1).is_linear is False

    assert (x*y + z + 1).is_quadratic
    assert (x*y*z + 1).is_quadratic is False

    assert R(0).is_squarefree
    assert (x - 1).is_squarefree
    assert ((x - 1)**2).is_squarefree is False
    assert (x**2 + y**2).is_squarefree
    assert ((x + y)**2).is_squarefree is False

    assert (x**2 + x + 1).is_irreducible
    assert (x**2 + 2*x + 1).is_irreducible is False

    _, t = ring("t", FF(11))

    assert (7*t + 3).is_irreducible
    assert (7*t**2 + 3*t + 1).is_irreducible is False

    _, u = ring("u", ZZ)
    f = u**16 + u**14 - u**10 - u**8 - u**6 + u**2

    assert f.is_cyclotomic is False
    assert (f + 1).is_cyclotomic

    pytest.raises(MultivariatePolynomialError, lambda: x.is_cyclotomic)

    assert f.is_positive is True
    assert f.is_negative is False
    assert f.is_nonnegative is True
    assert f.is_nonpositive is False


def test_PolyElement_drop():
    R,  x, y, z = ring("x,y,z", ZZ)

    assert R(1).drop(0).ring == ZZ.poly_ring("y", "z")
    assert R(1).drop(0).drop(0).ring == ZZ.poly_ring("z")
    assert isinstance(R(1).drop(0).drop(0).drop(0), R.dtype) is False

    pytest.raises(ValueError, lambda: z.drop(0).drop(0).drop(0))
    pytest.raises(ValueError, lambda: x.drop(0))

    f = z**2*x + 2*z*y + x*z + 1
    R2 = R.drop_to_ground(z)
    assert f.drop_to_ground(z) == z**2*R2.x + 2*z*R2.y + z*R2.x + 1

    R3 = R.drop(y, z)
    assert R3 == ZZ.poly_ring('x')
    pytest.raises(ValueError, lambda: R3.x.drop_to_ground(R3.x))


def test_PolyElement_pdiv():
    _,  x, y = ring("x,y", ZZ)

    f, g = x**2 - y**2, x - y
    q, r = x + y, 0

    assert f.pdiv(g) == (q, r)
    assert f.prem(g) == r
    assert f.pquo(g) == q
    assert f.pexquo(g) == q


def test_PolyElement_gcdex():
    _, x = ring("x", QQ)

    f, g = 2*x, x**2 - 16
    s, t, h = x/32, -QQ(1, 16), 1

    assert f.half_gcdex(g) == (s, h)
    assert f.gcdex(g) == (s, t, h)

    _, x, y = ring("x,y", QQ)

    pytest.raises(MultivariatePolynomialError, lambda: (x + y).half_gcdex(x*y))
    pytest.raises(MultivariatePolynomialError, lambda: (x + y).gcdex(x*y))


def test_PolyElement_subresultants():
    _, x = ring("x", ZZ)
    f, g, h = x**2 - 2*x + 1, x**2 - 1, 2*x - 2

    assert f.subresultants(g) == [f, g, h]


def test_PolyElement_resultant():
    _, x = ring("x", ZZ)
    f, g, h = x**2 - 2*x + 1, x**2 - 1, 0

    assert f.resultant(g) == h


def test_PolyElement_discriminant():
    _, x = ring("x", ZZ)
    f, g = x**3 + 3*x**2 + 9*x - 13, -11664

    assert f.discriminant() == g

    F, a, b, c = ring("a,b,c", ZZ)
    _, x = ring("x", F)

    f, g = a*x**2 + b*x + c, b**2 - 4*a*c

    assert f.discriminant() == g


def test_PolyElement_decompose():
    _, x = ring("x", ZZ)

    f = x**12 + 20*x**10 + 150*x**8 + 500*x**6 + 625*x**4 - 2*x**3 - 10*x + 9
    g = x**4 - 2*x + 9
    h = x**3 + 5*x

    assert g.compose(x, h) == f
    assert f.decompose() == [g, h]

    R, x, y = ring("x,y", ZZ)
    pytest.raises(MultivariatePolynomialError, lambda: (x + y).decompose())


def test_PolyElement_shift():
    _, x = ring("x", ZZ)
    assert (x**2 - 2*x + 1).shift(2) == x**2 + 2*x + 1

    R, x, y = ring("x,y", ZZ)
    pytest.raises(MultivariatePolynomialError, lambda: (x + y).shift(2))


def test_PolyElement_sturm():
    F, t = field("t", ZZ)
    _, x = ring("x", F)

    f = 1024/(15625*t**8)*x**5 - 4096/(625*t**8)*x**4 + 32/(15625*t**4)*x**3 - 128/(625*t**4)*x**2 + F(1)/62500*x - F(1)/625

    assert f.sturm() == [
        x**3 - 100*x**2 + t**4/64*x - 25*t**4/16,
        3*x**2 - 200*x + t**4/64,
        (-t**4/96 + F(20000)/9)*x + 25*t**4/18,
        (-9*t**12 - 11520000*t**8 - 3686400000000*t**4)/(576*t**8 - 245760000*t**4 + 26214400000000),
    ]

    R, x, y = ring("x,y", ZZ)
    pytest.raises(MultivariatePolynomialError, lambda: (x + y).sturm())


def test_PolyElement_sqf_norm():
    R, x = ring("x", QQ.algebraic_field(sqrt(3)))
    X = R.to_ground().x

    assert (x**2 - 2).sqf_norm() == (1, x**2 - 2*sqrt(3)*x + 1, X**4 - 10*X**2 + 1)

    R, x = ring("x", QQ.algebraic_field(sqrt(2)))
    X = R.to_ground().x

    assert (x**2 - 3).sqf_norm() == (1, x**2 - 2*sqrt(2)*x - 1, X**4 - 10*X**2 + 1)


def test_PolyElement_sqf_list():
    _, x = ring("x", ZZ)

    f = x**5 - x**3 - x**2 + 1
    g = x**3 + 2*x**2 + 2*x + 1
    h = x - 1
    p = x**4 + x**3 - x - 1

    assert f.sqf_part() == p
    assert f.sqf_list() == (1, [(g, 1), (h, 2)])


def test_PolyElement_factor_list():
    _, x = ring("x", ZZ)

    f = x**5 - x**3 - x**2 + 1

    u = x + 1
    v = x - 1
    w = x**2 + x + 1

    assert f.factor_list() == (1, [(u, 1), (v, 2), (w, 1)])


def test_PolyElement_almosteq():
    R, x, y = ring('x, y', RR)
    z = symbols('z')

    assert x.almosteq(y) is False
    assert x.almosteq(1) is False
    assert (x + 2*y).almosteq(2) is False
    assert (x + 2*y).almosteq(2*x + y) is False
    assert R.one.almosteq(2) is False
    assert R.one.almosteq(z) is False
