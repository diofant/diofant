from diofant.abc import n
from diofant.matrices.expressions import MatrixSymbol
from diofant.matrices.expressions.diagonal import DiagonalMatrix, DiagonalOf


__all__ = ()

x = MatrixSymbol('x', n, 1)
X = MatrixSymbol('X', n, n)
D = DiagonalMatrix(x)
d = DiagonalOf(X)


def test_DiagonalMatrix():
    assert D.shape == (n, n)
    assert D[1, 2] == 0
    assert D[1, 1] == x[1, 0]


def test_DiagonalOf():
    assert d.shape == (n, 1)
    assert d[2, 0] == X[2, 2]
