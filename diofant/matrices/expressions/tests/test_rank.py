import pytest

from diofant.core import symbols
from diofant.matrices import Matrix, eye
from diofant.matrices.expressions import (Identity, MatrixExpr, MatrixSymbol,
                                          Rank, ZeroMatrix, rank)


__all__ = ()

n = symbols('n', integer=True)
A = MatrixSymbol('A', n, n)


def test_rank():
    assert isinstance(Rank(A), Rank)
    assert not isinstance(Rank(A), MatrixExpr)
    assert rank(eye(3)) == 3

    pytest.raises(TypeError, lambda: Rank(1))

    assert Rank(A).arg is A


def test_eval_rank():
    assert rank(Identity(n)) == n
    assert rank(ZeroMatrix(n, n)) == 0
    assert rank(Matrix([[1, 2], [3, 4]])) == 2
