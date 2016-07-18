from diofant.matrices.expressions.factorizations import lu, LofCholesky, qr, svd
from diofant import Symbol, MatrixSymbol


n = Symbol('n')
X = MatrixSymbol('X', n, n)


def test_LU():
    L, U = lu(X)
    assert L.shape == U.shape == X.shape


def test_Cholesky():
    L = LofCholesky(X)


def test_QR():
    Q_, R = qr(X)
    assert Q_.shape == R.shape == X.shape


def test_svd():
    U, S, V = svd(X)
    assert U.shape == S.shape == V.shape == X.shape
