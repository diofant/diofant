from diofant.matrices.expressions.factorizations import (lu, LofCholesky,
                                                         qr, svd, eig)
from diofant import Symbol, MatrixSymbol

__all__ = ()

n = Symbol('n')
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
