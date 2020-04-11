"""Tests for the PolynomialRing classes."""

import pytest

from diofant import FF, QQ, ZZ, CoercionFailed, GeneratorsNeeded, sqrt
from diofant.abc import x, y
from diofant.polys.orderings import build_product_order


__all__ = ()


ALG = QQ.algebraic_field(sqrt(2), sqrt(3))


def test_build_order():
    R = QQ.poly_ring(x, y, order=build_product_order((('lex', x),
                                                      ('ilex', y)), (x, y)))
    assert R.order((1, 5)) == ((1,), (-5,))


def test_globalring():
    Qxy = QQ.inject(x, y).field
    R = QQ.inject(x, y)
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
    assert X**2 // X == X

    assert R.convert(ZZ.inject(x, y).convert(x), ZZ.inject(x, y)) == X
    assert R.convert(Qxy.convert(x), Qxy) == X


def test_localring():
    Qxy = QQ.inject(x, y).field
    R = QQ.poly_ring(x, y, order='ilex')
    X = R.convert(x)
    Y = R.convert(y)

    assert x in R
    assert 1/x not in R
    assert Y in R
    assert X.ring == R.ring
    assert X + Y == R.convert(x + y)
    assert X - Y == R.convert(x - y)
    assert X + 1 == R.convert(x + 1)
    assert X**2 // X == X

    assert R.convert(ZZ.inject(x, y).convert(x), ZZ.inject(x, y)) == X
    assert R.convert(Qxy.convert(x), Qxy) == X


def test_conversion():
    L = QQ.poly_ring(x, y, order='ilex')
    G = QQ.inject(x, y)

    assert L.convert(x) == L.convert(G.convert(x), G)
    assert G.convert(x) == G.convert(L.convert(x), L)
    pytest.raises(CoercionFailed, lambda: G.convert(L.convert(1/(1 + x)), L))

    R = ALG.inject(x, y)
    assert R.convert(ALG(1), ALG) == R(1)
    pytest.raises(CoercionFailed,
                  lambda: R.convert(ALG(1), QQ.algebraic_field(sqrt(2))))

    R = R.drop(y)
    pytest.raises(CoercionFailed, lambda: R.convert(G(y), R))

    pytest.raises(CoercionFailed, lambda: R.convert(FF(8)(2)))


def test_units():
    R = QQ.inject(x)
    assert R.convert(1) == R.one
    assert R.convert(x) != R.one
    assert R.convert(1 + x) != R.one

    R = QQ.poly_ring(x, order='ilex')
    assert R.convert(1) == R.one
    assert R.convert(x) != R.one

    R = ZZ.inject(x)
    assert R.convert(1) == R.one
    assert R.convert(2) != R.one
    assert R.convert(x) != R.one
    assert R.convert(1 + x) != R.one


def test_poly_frac():
    pytest.raises(GeneratorsNeeded, lambda: QQ.inject())
    pytest.raises(GeneratorsNeeded, lambda: QQ.inject().field)


def test_methods():
    R = QQ.inject(x)
    X = R.convert(x)

    assert R.is_normal(-X) is False
    assert R.is_normal(+X) is True

    assert R.gcdex(X**3 - X, X**2) == (-1, X, X)

    F = QQ.inject(y).field
    Y = F.convert(y)
    assert F.is_normal(-Y) is False
    assert F.is_normal(+Y) is True
