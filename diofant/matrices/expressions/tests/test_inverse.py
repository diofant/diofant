import pytest

from diofant.core import symbols
from diofant.matrices import Identity, ShapeError, det, eye
from diofant.matrices.expressions import Inverse, MatrixSymbol


__all__ = ()

n, m, l = symbols('n m l', integer=True)
A = MatrixSymbol('A', n, m)
B = MatrixSymbol('B', m, l)
C = MatrixSymbol('C', n, n)
D = MatrixSymbol('D', n, n)
E = MatrixSymbol('E', m, n)


def test_inverse():
    pytest.raises(ShapeError, lambda: Inverse(A))
    pytest.raises(ShapeError, lambda: Inverse(A*B))
    pytest.raises(TypeError, lambda: Inverse(1))

    assert Inverse(C).shape == (n, n)
    assert Inverse(A*E).shape == (n, n)
    assert Inverse(E*A).shape == (m, m)
    assert Inverse(C).inverse() == C
    assert isinstance(Inverse(Inverse(C)), Inverse)

    assert C.inverse().inverse() == C

    assert C.inverse()*C == Identity(C.rows)

    assert Identity(n).inverse() == Identity(n)
    assert (3*Identity(n)).inverse() == Identity(n)/3

    # Simplifies Muls if possible (i.e. submatrices are square)
    assert (C*D).inverse() == D.inverse()*C.inverse()
    # But still works when not possible
    assert isinstance((A*E).inverse(), Inverse)

    assert Inverse(eye(3)).doit() == eye(3)
    assert Inverse(eye(3)).doit(deep=False) == eye(3)

    assert det(Inverse(C)) == 1/det(C)
