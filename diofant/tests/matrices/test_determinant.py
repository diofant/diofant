import pytest

from diofant import (Determinant, Identity, Matrix, MatrixExpr, MatrixSymbol,
                     ShapeError, Transpose, ZeroMatrix, det, eye, symbols)


__all__ = ()

n = symbols('n', integer=True)
A = MatrixSymbol('A', n, n)
B = MatrixSymbol('B', n, n)
C = MatrixSymbol('C', 3, 4)


def test_det():
    assert isinstance(Determinant(A), Determinant)
    assert not isinstance(Determinant(A), MatrixExpr)
    pytest.raises(ShapeError, lambda: Determinant(C))
    assert det(eye(3)) == 1
    assert det(Matrix(3, 3, [1, 3, 2, 4, 1, 3, 2, 5, 2])) == 17
    A / det(A)  # Make sure this is possible

    pytest.raises(TypeError, lambda: Determinant(1))

    assert Determinant(A).arg is A


def test_eval_determinant():
    assert det(Matrix()) == 1
    assert det(Identity(n)) == 1
    assert det(ZeroMatrix(n, n)) == 0
    assert det(Transpose(A)) == det(A)
