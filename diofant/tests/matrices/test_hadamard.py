import pytest

from diofant import (HadamardProduct, MatrixSymbol, ShapeError,
                     hadamard_product, symbols)
from diofant.abc import k, m, n


__all__ = ()

Z = MatrixSymbol('Z', n, n)
A = MatrixSymbol('A', n, m)
B = MatrixSymbol('B', n, m)
C = MatrixSymbol('C', m, k)


def test_HadamardProduct():
    assert HadamardProduct(A, B, A).shape == A.shape

    pytest.raises(ShapeError, lambda: HadamardProduct(A, B.T))
    pytest.raises(TypeError, lambda: HadamardProduct(A, n))
    pytest.raises(TypeError, lambda: HadamardProduct(A, 1))

    assert HadamardProduct(A, 2*B, -A)[1, 1] == \
        -2 * A[1, 1] * B[1, 1] * A[1, 1]

    mix = HadamardProduct(Z*A, B)*C
    assert mix.shape == (n, k)

    assert set(HadamardProduct(A, B, A).T.args) == {A.T, A.T, B.T}


def test_HadamardProduct_isnt_commutative():
    assert HadamardProduct(A, B) != HadamardProduct(B, A)


def test_mixed_indexing():
    X = MatrixSymbol('X', 2, 2)
    Y = MatrixSymbol('Y', 2, 2)
    Z = MatrixSymbol('Z', 2, 2)

    assert (X*HadamardProduct(Y, Z))[0, 0] == \
        X[0, 0]*Y[0, 0]*Z[0, 0] + X[0, 1]*Y[1, 0]*Z[1, 0]


def test_canonicalize():
    X = MatrixSymbol('X', 2, 2)
    expr = HadamardProduct(X, check=False)
    assert isinstance(expr, HadamardProduct)
    expr2 = expr.doit()  # unpack is called
    assert isinstance(expr2, MatrixSymbol)


def test_hadamard():
    m, n, p = symbols('m, n, p', integer=True)
    A = MatrixSymbol('A', m, n)
    B = MatrixSymbol('B', m, n)
    C = MatrixSymbol('C', m, p)
    with pytest.raises(TypeError):
        hadamard_product()
    assert hadamard_product(A) == A
    assert isinstance(hadamard_product(A, B), HadamardProduct)
    assert hadamard_product(A, B).doit() == hadamard_product(A, B)
    with pytest.raises(ShapeError):
        hadamard_product(A, C)
