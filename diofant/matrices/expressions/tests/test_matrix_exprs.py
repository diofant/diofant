import pytest

from diofant.abc import t, x
from diofant.core import Add, Eq, Integer, Mul, Rational, symbols
from diofant.external import import_module
from diofant.functions import cos, sin, sqrt, transpose
from diofant.matrices import (Adjoint, Identity, ImmutableMatrix, Inverse,
                              MatAdd, MatMul, MatPow, Matrix, MatrixExpr,
                              MatrixSymbol, ShapeError, Transpose, ZeroMatrix)
from diofant.matrices.expressions.matexpr import MatrixElement
from diofant.matrices.expressions.slice import MatrixSlice
from diofant.simplify import simplify


numpy = import_module('numpy')

__all__ = ()

n, m, l, k, p = symbols('n m l k p', integer=True)
A = MatrixSymbol('A', n, m)
B = MatrixSymbol('B', m, l)
C = MatrixSymbol('C', n, n)
D = MatrixSymbol('D', n, n)
E = MatrixSymbol('E', m, n)


def test_shape():
    assert A.shape == (n, m)
    assert (A*B).shape == (n, l)
    pytest.raises(ShapeError, lambda: B*A)


def test_matexpr():
    assert (x*A).shape == A.shape
    assert (x*A).__class__ == MatMul
    assert 2*A - A - A == ZeroMatrix(*A.shape)
    assert (A*B).shape == (n, l)


def test_subs():
    A = MatrixSymbol('A', n, m)
    B = MatrixSymbol('B', m, l)
    C = MatrixSymbol('C', m, l)

    assert A.subs(n, m).shape == (m, m)

    assert (A*B).subs(B, C) == A*C

    assert (A*B).subs(l, n).is_square


def test_ZeroMatrix():
    A = MatrixSymbol('A', n, m)
    Z = ZeroMatrix(n, m)

    assert A + Z == A
    assert A*Z.T == ZeroMatrix(n, n)
    assert Z*A.T == ZeroMatrix(n, n)
    assert A - A == ZeroMatrix(*A.shape)

    assert not Z

    assert transpose(Z) == ZeroMatrix(m, n)
    assert Z.conjugate() == Z

    assert ZeroMatrix(n, n)**0 == Identity(n)
    with pytest.raises(ShapeError):
        Z**0
    with pytest.raises(ShapeError):
        Z**2


def test_ZeroMatrix_doit():
    Znn = ZeroMatrix(Add(n, n, evaluate=False), n)
    assert isinstance(Znn.rows, Add)
    assert Znn.doit() == ZeroMatrix(2*n, n)
    assert isinstance(Znn.doit().rows, Mul)


def test_Identity():
    A = MatrixSymbol('A', n, m)
    In = Identity(n)
    Im = Identity(m)

    assert A*Im == A
    assert In*A == A

    assert transpose(In) == In
    assert In.inverse() == In
    assert In.conjugate() == In


def test_Identity_doit():
    Inn = Identity(Add(n, n, evaluate=False))
    assert isinstance(Inn.rows, Add)
    assert Inn.doit() == Identity(2*n)
    assert isinstance(Inn.doit().rows, Mul)


def test_addition():
    A = MatrixSymbol('A', n, m)
    B = MatrixSymbol('B', n, m)

    assert isinstance(A + B, MatAdd)
    assert (A + B).shape == A.shape
    assert isinstance(A - A + 2*B, MatMul)

    pytest.raises(ShapeError, lambda: A + B.T)
    pytest.raises(TypeError, lambda: A + 1)
    pytest.raises(TypeError, lambda: 5 + A)
    pytest.raises(TypeError, lambda: 5 - A)

    assert A + ZeroMatrix(n, m) - A == ZeroMatrix(n, m)
    with pytest.raises(TypeError):
        ZeroMatrix(n, m) + Integer(0)


def test_multiplication():
    A = MatrixSymbol('A', n, m)
    B = MatrixSymbol('B', m, l)
    C = MatrixSymbol('C', n, n)

    assert (2*A*B).shape == (n, l)

    assert (A*0*B) == ZeroMatrix(n, l)

    pytest.raises(ShapeError, lambda: B*A)
    assert (2*A).shape == A.shape

    assert A * ZeroMatrix(m, m) * B == ZeroMatrix(n, l)

    assert C * Identity(n) * C.inverse() == Identity(n)

    assert B/2 == Rational(1, 2)*B
    pytest.raises(NotImplementedError, lambda: 2/B)

    A = MatrixSymbol('A', n, n)
    B = MatrixSymbol('B', n, n)
    assert Identity(n) * (A + B) == A + B


def test_MatPow():
    A = MatrixSymbol('A', n, n)

    AA = MatPow(A, 2)
    assert AA.exp == 2
    assert AA.base == A
    assert (A**n).exp == n

    assert A**0 == Identity(n)
    assert A**1 == A
    assert A**2 == AA
    assert A**-1 == Inverse(A)
    assert A**Rational(1, 2) == sqrt(A)
    pytest.raises(ShapeError, lambda: MatrixSymbol('B', 3, 2)**2)


def test_MatrixSymbol():
    X = MatrixSymbol('X', n, m)
    assert X.shape == (n, m)
    pytest.raises(TypeError, lambda: MatrixSymbol('X', n, m)(t))  # issue sympy/sympy#5855
    assert X.doit() == X.doit(deep=False) == X
    assert X.canonicalize() == X


def test_dense_conversion():
    X = MatrixSymbol('X', 2, 2)
    assert ImmutableMatrix(X) == ImmutableMatrix(2, 2, lambda i, j: X[i, j])
    assert Matrix(X) == Matrix(2, 2, lambda i, j: X[i, j])


def test_free_symbols():
    assert (C*D).free_symbols == {C, D}


def test_zero_matmul():
    assert isinstance(0 * MatrixSymbol('X', 2, 2), MatrixExpr)


def test_matadd_simplify():
    A = MatrixSymbol('A', 1, 1)
    assert simplify(MatAdd(A, ImmutableMatrix([[sin(x)**2 + cos(x)**2]]))) == \
        MatAdd(A, ImmutableMatrix([[1]]))


def test_matmul_simplify():
    A = MatrixSymbol('A', 1, 1)
    assert simplify(MatMul(A, ImmutableMatrix([[sin(x)**2 + cos(x)**2]]))) == \
        MatMul(A, ImmutableMatrix([[1]]))


def test_invariants():
    A = MatrixSymbol('A', n, m)
    B = MatrixSymbol('B', m, l)
    X = MatrixSymbol('X', n, n)
    objs = [Identity(n), ZeroMatrix(m, n), MatMul(A, B), MatAdd(A, A),
            Transpose(A), Adjoint(A), Inverse(X), MatPow(X, 2), MatPow(X, -1),
            MatPow(X, 0)]
    for obj in objs:
        assert obj == obj.__class__(*obj.args)


def test_indexing():
    A = MatrixSymbol('A', n, m)
    A[1, 2]
    A[l, k]
    A[l+1, k+1]

    assert A[:] == MatrixSlice(A, (0, n, 1), (0, m, 1))


def test_single_indexing():
    A = MatrixSymbol('A', 2, 3)
    assert A[1] == A[0, 1]
    assert A[3] == A[1, 0]
    assert list(A[:2, :2]) == [A[0, 0], A[0, 1], A[1, 0], A[1, 1]]
    pytest.raises(IndexError, lambda: A[6])
    pytest.raises(IndexError, lambda: A[n])
    B = MatrixSymbol('B', n, m)
    pytest.raises(IndexError, lambda: B[1])


def test_MatrixElement_diff():
    assert (A[3, 0]*A[0, 0]).diff(A[0, 0]) == A[3, 0]


def test_MatrixElement_doit():
    u = MatrixSymbol('u', 2, 1)
    v = ImmutableMatrix([3, 5])
    assert u[0, 0].subs(u, v).doit() == v[0, 0]
    assert u[0, 0].subs(u, v).doit(deep=False) == v[0, 0]


def test_identity_powers():
    M = Identity(n)
    assert MatPow(M, 3).doit() == M**3
    assert M**n == M
    assert MatPow(M, 0).doit() == M**2
    assert M**-2 == M
    assert MatPow(M, -2).doit() == M**0
    N = Identity(3)
    assert MatPow(N, 2).doit() == N**n
    assert MatPow(N, 3).doit() == N
    assert MatPow(N, -2).doit() == N**4
    assert MatPow(N, 2).doit() == N**0


def test_Zero_power():
    z1 = ZeroMatrix(n, n)
    assert z1**4 == z1
    pytest.raises(ValueError, lambda: z1**-2)
    assert z1**0 == Identity(n)
    assert MatPow(z1, 2).doit() == z1**2
    pytest.raises(ValueError, lambda: MatPow(z1, -2).doit())
    z2 = ZeroMatrix(3, 3)
    assert MatPow(z2, 4).doit() == z2**4
    pytest.raises(ValueError, lambda: z2**-3)
    assert z2**3 == MatPow(z2, 3).doit()
    assert z2**0 == Identity(3)
    pytest.raises(ValueError, lambda: MatPow(z2, -1).doit())


def test_diofantissue_469():
    A = MatrixSymbol("A", n, n)
    B = MatrixSymbol("B", n, n)
    expr = Eq(A, B)
    assert simplify(expr) == expr


@pytest.mark.skipif(numpy is None, reason="Couldn't import numpy.")
def test_array_coeersion():
    A = MatrixSymbol('A', 2, 2)
    assert numpy.array(A)[1][0] == MatrixElement(A, 1, 0)
