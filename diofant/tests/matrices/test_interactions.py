"""
We have a few different kind of Matrices
Matrix, ImmutableMatrix, MatrixExpr

Here we test the extent to which they cooperate
"""

import pytest

from diofant import (Identity, ImmutableMatrix, MatAdd, Matrix, MatrixExpr,
                     MatrixSymbol, eye)
from diofant.abc import a
from diofant.matrices.matrices import classof


__all__ = ()

SM = MatrixSymbol('X', 3, 3)
MM = Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
IM = ImmutableMatrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
meye = eye(3)
imeye = ImmutableMatrix(eye(3))
ideye = Identity(3)


def test_IM_MM():
    assert isinstance(MM + IM, Matrix)
    assert isinstance(IM + MM, ImmutableMatrix)
    assert isinstance(2*IM + MM, ImmutableMatrix)
    assert MM.equals(IM)
    assert MM.equals(2) is False


def test_ME_MM():
    assert isinstance(Identity(3) + MM, MatrixExpr)
    assert isinstance(SM + MM, MatAdd)
    assert isinstance(MM + SM, MatAdd)
    assert (Identity(3) + MM)[1, 1] == 6


def test_equality():
    a, b, c = Identity(3), eye(3), ImmutableMatrix(eye(3))
    for x in [a, b, c]:
        for y in [a, b, c]:
            assert x.equals(y)


def test_matrix_symbol_MM():
    X = MatrixSymbol('X', 3, 3)
    Y = eye(3) + X
    assert Y[1, 1] == 1 + X[1, 1]


def test_matrix_symbol_vector_matrix_multiplication():
    # issue sympy/sympy#11090
    SV = MatrixSymbol('V', 3, 1)
    A = MM * SV
    B = IM * SV
    assert A == B
    C = (SV.T * MM.T).T
    assert B == C
    D = (SV.T * IM.T).T
    assert C == D


def test_indexing_interactions():
    assert (a * IM)[1, 1] == 5*a
    assert (SM + IM)[1, 1] == SM[1, 1] + IM[1, 1]
    assert (SM * IM)[1, 1] == SM[1, 0]*IM[0, 1] + SM[1, 1]*IM[1, 1] + \
        SM[1, 2]*IM[2, 1]


def test_classof():
    A = Matrix(3, 3, range(9))
    B = ImmutableMatrix(3, 3, range(9))
    C = MatrixSymbol('C', 3, 3)
    assert classof(A, A) == Matrix
    assert classof(B, B) == ImmutableMatrix
    assert classof(A, B) == ImmutableMatrix
    assert classof(B, A) == ImmutableMatrix
    pytest.raises(TypeError, lambda: classof(A, C))
