from diofant import (Adjoint, Integer, Matrix, MatrixSymbol, Transpose,
                     adjoint, conjugate, eye, symbols, trace, transpose)


__all__ = ()

n, m, l, k, p = symbols('n m l k p', integer=True)
A = MatrixSymbol('A', n, m)
B = MatrixSymbol('B', m, l)
C = MatrixSymbol('C', n, n)
D = MatrixSymbol('D', n, n)


def test_adjoint():
    Sq = MatrixSymbol('Sq', n, n)

    assert Adjoint(A).shape == (m, n)
    assert Adjoint(A*B).shape == (l, n)
    assert adjoint(Adjoint(A)) == A
    assert isinstance(Adjoint(Adjoint(A)), Adjoint)

    assert conjugate(Adjoint(A)) == Transpose(A) == Adjoint(A).conjugate()
    assert transpose(Adjoint(A)) == Adjoint(Transpose(A)) == Transpose(A).adjoint()

    assert Adjoint(eye(3)).doit() == Adjoint(eye(3)).doit(deep=False) == eye(3)

    assert Adjoint(Integer(5)).doit() == Integer(5)

    assert Adjoint(Matrix([[1, 2], [3, 4]])).doit() == Matrix([[1, 3], [2, 4]])

    assert adjoint(trace(Sq)) == conjugate(trace(Sq))
    assert trace(adjoint(Sq)) == conjugate(trace(Sq))

    assert Adjoint(Sq)[0, 1] == conjugate(Sq[1, 0])

    assert Adjoint(A*B).doit() == Adjoint(B) * Adjoint(A)
    assert Adjoint(C + D).doit() == Adjoint(C) + Adjoint(D)
