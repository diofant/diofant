from diofant import MatrixSymbol
from diofant.abc import n
from diofant.matrices.expressions.factorizations import (LofCholesky, eig, lu,
                                                         qr, svd)


__all__ = ()

X = MatrixSymbol('X', n, n)


def test_LU():
    L, U = lu(X)
    assert L.shape == U.shape == X.shape


def test_Cholesky():
    L = LofCholesky(X)
    assert L.shape == X.shape


def test_QR():
    Q_, R = qr(X)
    assert Q_.shape == R.shape == X.shape


def test_svd():
    U, S, V = svd(X)
    assert U.shape == S.shape == V.shape == X.shape


def test_eig():
    E, V = eig(X)
    assert E.shape == V.shape == X.shape
