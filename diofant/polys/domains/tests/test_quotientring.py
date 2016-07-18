"""Tests for quotient rings."""

import pytest

from diofant import QQ, ZZ
from diofant.polys.polyerrors import NotReversible

from diofant.abc import x, y


def test_QuotientRingElement():
    R = QQ.old_poly_ring(x)/[x**10]
    X = R.convert(x)

    assert X*(X + 1) == R.convert(x**2 + x)
    assert X*x == R.convert(x**2)
    assert x*X == R.convert(x**2)
    assert X + x == R.convert(2*x)
    assert x + X == 2*X
    assert X**2 == R.convert(x**2)
    assert 1/(1 - X) == R.convert(sum(x**i for i in range(10)))
    assert X**10 == R.zero
    assert X != x

    pytest.raises(NotReversible, lambda: 1/X)


def test_QuotientRing():
    I = QQ.old_poly_ring(x).ideal(x**2 + 1)
    R = QQ.old_poly_ring(x)/I

    assert R == QQ.old_poly_ring(x)/[x**2 + 1]
    assert R == QQ.old_poly_ring(x)/QQ.old_poly_ring(x).ideal(x**2 + 1)
    assert R != QQ.old_poly_ring(x)

    assert R.convert(1)/x == -x + I
    assert -1 + I == x**2 + I
    assert R.convert(ZZ(1), ZZ) == 1 + I
    assert R.convert(R.convert(x), R) == R.convert(x)

    X = R.convert(x)
    Y = QQ.old_poly_ring(x).convert(x)
    assert -1 + I == X**2 + I
    assert -1 + I == Y**2 + I
    assert R.to_diofant(X) == x

    pytest.raises(ValueError, lambda: QQ.old_poly_ring(x)/QQ.old_poly_ring(x, y).ideal(x))

    R = QQ.old_poly_ring(x, order="ilex")
    I = R.ideal(x)
    assert R.convert(1) + I == (R/I).convert(1)
