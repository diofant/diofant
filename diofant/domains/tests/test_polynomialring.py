"""Tests for the PolynomialRing classes. """

import pytest

from diofant import sqrt
from diofant.abc import x, y
from diofant.domains import QQ, ZZ
from diofant.polys.orderings import build_product_order
from diofant.polys.polyerrors import CoercionFailed, GeneratorsNeeded


__all__ = ()


ALG = QQ.algebraic_field(sqrt(2), sqrt(3))


def test_build_order():
    R = QQ.poly_ring(x, y, order=build_product_order((("lex", x),
                                                      ("ilex", y)), (x, y)))
    assert R.order((1, 5)) == ((1,), (-5,))


def test_globalring():
    Qxy = QQ.frac_field(x, y)
    R = QQ.poly_ring(x, y)
    X = R.convert(x)
    Y = R.convert(y)

    assert x in R
    assert 1/x not in R
    assert 1/(1 + x) not in R
    assert Y in R
    assert X.ring == R.ring
    assert X * (Y**2 + 1) == R.convert(x * (y**2 + 1))
    assert X * Y == R.convert(x * y)
    assert X + Y == R.convert(x + y)
    assert X - Y == R.convert(x - y)
    assert X + 1 == R.convert(x + 1)
    assert X**2 / X == X

    assert R.from_PolynomialRing(ZZ.poly_ring(x, y).convert(x), ZZ.poly_ring(x, y)) == X
    assert R.from_FractionField(Qxy.convert(x), Qxy) == X


def test_localring():
    Qxy = QQ.frac_field(x, y)
    R = QQ.poly_ring(x, y, order="ilex")
    X = R.convert(x)
    Y = R.convert(y)

    assert x in R
    assert 1/x not in R
    assert Y in R
    assert X.ring == R.ring
    assert X + Y == R.convert(x + y)
    assert X - Y == R.convert(x - y)
    assert X + 1 == R.convert(x + 1)
    assert X**2 / X == X

    assert R.from_PolynomialRing(ZZ.poly_ring(x, y).convert(x), ZZ.poly_ring(x, y)) == X
    assert R.from_FractionField(Qxy.convert(x), Qxy) == X


def test_conversion():
    L = QQ.poly_ring(x, y, order="ilex")
    G = QQ.poly_ring(x, y)

    assert L.convert(x) == L.convert(G.convert(x), G)
    assert G.convert(x) == G.convert(L.convert(x), L)
    pytest.raises(CoercionFailed, lambda: G.convert(L.convert(1/(1 + x)), L))

    R = ALG.poly_ring(x, y)
    assert R.convert(ALG.new(1), ALG) == R.new(1)
    pytest.raises(CoercionFailed,
                  lambda: R.convert(ALG.new(1), QQ.algebraic_field(sqrt(2))))


def test_units():
    R = QQ.poly_ring(x)
    assert R.is_unit(R.convert(1))
    assert not R.is_unit(R.convert(x))
    assert not R.is_unit(R.convert(1 + x))

    R = QQ.poly_ring(x, order='ilex')
    assert R.is_unit(R.convert(1))
    assert not R.is_unit(R.convert(x))

    R = ZZ.poly_ring(x)
    assert R.is_unit(R.convert(1))
    assert not R.is_unit(R.convert(2))
    assert not R.is_unit(R.convert(x))
    assert not R.is_unit(R.convert(1 + x))


def test_poly_frac():
    pytest.raises(GeneratorsNeeded, lambda: QQ.poly_ring())
    pytest.raises(GeneratorsNeeded, lambda: QQ.frac_field())


def test_methods():
    R = QQ.poly_ring(x)
    X = R.convert(x)

    assert R.is_nonpositive(-X) is True
    assert R.is_positive(X) is True

    assert R.gcdex(X**3 - X, X**2) == (-1, X, X)

    assert R.factorial(3) == 6

    F = QQ.frac_field(y)
    Y = F.convert(y)
    assert F.is_nonpositive(-Y) is True
    assert F.is_nonnegative(Y) is True

    assert F.factorial(3) == 6
