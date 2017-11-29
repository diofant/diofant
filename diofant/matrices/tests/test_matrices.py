import collections

import pytest

from diofant import (Abs, Basic, E, Float, Function, I, Integer, Max, Min, N,
                     Poly, Pow, PurePoly, Rational, StrPrinter, Symbol, cos,
                     exp, oo, pi, simplify, sin, sqrt, sstr, symbols, sympify,
                     trigsimp)
from diofant.abc import a, b, c, d, k, n, x, y, z
from diofant.core.compatibility import iterable
from diofant.external import import_module
from diofant.matrices import (GramSchmidt, ImmutableMatrix,
                              ImmutableSparseMatrix, Matrix, SparseMatrix,
                              casoratian, diag, eye, hessian,
                              matrix_multiply_elementwise, ones, randMatrix,
                              rot_axis1, rot_axis2, rot_axis3, vandermonde,
                              wronskian, zeros)
from diofant.matrices.matrices import (DeferredVector, MatrixError,
                                       NonSquareMatrixError, ShapeError,
                                       mgamma)
from diofant.utilities.iterables import capture, flatten


__all__ = ()

# don't re-order this list
classes = (Matrix, SparseMatrix, ImmutableMatrix, ImmutableSparseMatrix)

numpy = import_module('numpy')


def test_args():
    for c, cls in enumerate(classes):
        m = cls.zeros(3, 2)
        # all should give back the same type of arguments, e.g. ints for shape
        assert m.shape == (3, 2) and all(type(i) is int for i in m.shape)
        assert m.rows == 3 and type(m.rows) is int
        assert m.cols == 2 and type(m.cols) is int
        if not c % 2:
            assert type(m._mat) is list
        else:
            assert type(m._smat) is dict


def test_division():
    v = Matrix(1, 2, [x, y])
    assert v.__truediv__(z) == Matrix(1, 2, [x/z, y/z])
    assert v.__truediv__(z) == Matrix(1, 2, [x/z, y/z])
    assert v/z == Matrix(1, 2, [x/z, y/z])


def test_sum():
    m = Matrix([[1, 2, 3], [x, y, x], [2*y, -50, z*x]])
    assert m + m == Matrix([[2, 4, 6], [2*x, 2*y, 2*x], [4*y, -100, 2*z*x]])
    n = Matrix(1, 2, [1, 2])
    pytest.raises(ShapeError, lambda: m + n)


def test_addition():
    a = Matrix((
        (1, 2),
        (3, 1),
    ))

    b = Matrix((
        (1, 2),
        (3, 0),
    ))

    assert a + b == a.add(b) == Matrix([[2, 4], [6, 1]])


def test_fancy_index_matrix():
    for M in (Matrix, SparseMatrix):
        a = M(3, 3, range(9))
        assert a == a[:, :]
        assert a[1, :] == Matrix(1, 3, [3, 4, 5])
        assert a[:, 1] == Matrix([1, 4, 7])
        assert a[[0, 1], :] == Matrix([[0, 1, 2], [3, 4, 5]])
        assert a[[0, 1], 2] == a[[0, 1], [2]]
        assert a[2, [0, 1]] == a[[2], [0, 1]]
        assert a[:, [0, 1]] == Matrix([[0, 1], [3, 4], [6, 7]])
        assert a[0, 0] == 0
        assert a[0:2, :] == Matrix([[0, 1, 2], [3, 4, 5]])
        assert a[:, 0:2] == Matrix([[0, 1], [3, 4], [6, 7]])
        assert a[::2, 1] == a[[0, 2], 1]
        assert a[1, ::2] == a[1, [0, 2]]
        a = M(3, 3, range(9))
        assert a[[0, 2, 1, 2, 1], :] == Matrix([
            [0, 1, 2],
            [6, 7, 8],
            [3, 4, 5],
            [6, 7, 8],
            [3, 4, 5]])
        assert a[:, [0, 2, 1, 2, 1]] == Matrix([
            [0, 2, 1, 2, 1],
            [3, 5, 4, 5, 4],
            [6, 8, 7, 8, 7]])

    a = SparseMatrix.zeros(3)
    a[1, 2] = 2
    a[0, 1] = 3
    a[2, 0] = 4
    assert a.extract([1, 1], [2]) == Matrix([
        [2],
        [2]])
    assert a.extract([1, 0], [2, 2, 2]) == Matrix([
        [2, 2, 2],
        [0, 0, 0]])
    assert a.extract([1, 0, 1, 2], [2, 0, 1, 0]) == Matrix([
        [2, 0, 0, 0],
        [0, 0, 3, 0],
        [2, 0, 0, 0],
        [0, 4, 0, 4]])


def test_multiplication():
    a = Matrix((
        (1, 2),
        (3, 1),
        (0, 6),
    ))

    b = Matrix((
        (1, 2),
        (3, 0),
    ))

    c = a*b
    assert c == a.multiply(b)
    assert c[0, 0] == 7
    assert c[0, 1] == 2
    assert c[1, 0] == 6
    assert c[1, 1] == 6
    assert c[2, 0] == 18
    assert c[2, 1] == 0

    h = matrix_multiply_elementwise(a, c)
    assert h == a.multiply_elementwise(c)
    assert h[0, 0] == 7
    assert h[0, 1] == 4
    assert h[1, 0] == 18
    assert h[1, 1] == 6
    assert h[2, 0] == 0
    assert h[2, 1] == 0
    pytest.raises(ShapeError, lambda: matrix_multiply_elementwise(a, b))

    c = b * x
    assert isinstance(c, Matrix)
    assert c[0, 0] == x
    assert c[0, 1] == 2*x
    assert c[1, 0] == 3*x
    assert c[1, 1] == 0

    c2 = x * b
    assert c == c2

    c = 5 * b
    assert isinstance(c, Matrix)
    assert c[0, 0] == 5
    assert c[0, 1] == 2*5
    assert c[1, 0] == 3*5
    assert c[1, 1] == 0


def test_power():
    pytest.raises(NonSquareMatrixError, lambda: Matrix((1, 2))**2)

    R = Rational
    A = Matrix([[2, 3], [4, 5]])
    assert (A**-3)[:] == [R(-269)/8, R(153)/8, R(51)/2, R(-29)/2]
    assert (A**5)[:] == [6140, 8097, 10796, 14237]
    A = Matrix([[2, 1, 3], [4, 2, 4], [6, 12, 1]])
    assert (A**3)[:] == [290, 262, 251, 448, 440, 368, 702, 954, 433]
    assert A**0 == eye(3)
    assert A**1 == A
    assert (Matrix([[2]]) ** 100)[0, 0] == 2**100
    assert eye(2)**10000000 == eye(2)
    assert Matrix([[1, 2], [3, 4]])**Integer(2) == Matrix([[7, 10], [15, 22]])

    A = Matrix([[33, 24], [48, 57]])
    assert (A**Rational(1, 2))[:] == [5, 2, 4, 7]
    A = Matrix([[0, 4], [-1, 5]])
    assert (A**Rational(1, 2))**2 == A

    assert Matrix([[1, 0], [1, 1]])**Rational(1, 2) == Matrix([[1, 0], [Rational(1, 2), 1]])
    assert Matrix([[1, a], [0, 1]])**n == Matrix([[1, a*n], [0, 1]])
    assert Matrix([[b, a], [0, b]])**n == Matrix([[b**n, a*b**(n - 1)*n], [0, b**n]])
    assert Matrix([[a, 1, 0], [0, a, 1], [0, 0, a]])**n == Matrix([
        [a**n, a**(n - 1)*n, a**(n - 2)*(n - 1)*n/2],
        [0, a**n, a**(n - 1)*n],
        [0, 0, a**n]])
    assert Matrix([[a, 1, 0], [0, a, 0], [0, 0, b]])**n == Matrix([
        [a**n, a**(n - 1)*n, 0],
        [0, a**n, 0],
        [0, 0, b**n]])


def test_creation():
    pytest.raises(ValueError, lambda: Matrix(5, 5, range(20)))
    pytest.raises(IndexError, lambda: Matrix((1, 2))[2])
    with pytest.raises(IndexError):
        Matrix((1, 2))[1:2] = 5
    with pytest.raises(IndexError):
        Matrix((1, 2))[3] = 5

    assert Matrix() == Matrix([]) == Matrix([[]]) == Matrix(0, 0, [])

    a = Matrix([[x, 0], [0, 0]])
    m = a
    assert m.cols == m.rows
    assert m.cols == 2
    assert m[:] == [x, 0, 0, 0]

    b = Matrix(2, 2, [x, 0, 0, 0])
    m = b
    assert m.cols == m.rows
    assert m.cols == 2
    assert m[:] == [x, 0, 0, 0]

    assert a == b

    assert Matrix(b) == b

    c = Matrix((
        Matrix((
            (1, 2, 3),
            (4, 5, 6)
        )),
        (7, 8, 9)
    ))
    assert c.cols == 3
    assert c.rows == 3
    assert c[:] == [1, 2, 3, 4, 5, 6, 7, 8, 9]

    assert Matrix(eye(2)) == eye(2)
    assert ImmutableMatrix(ImmutableMatrix(eye(2))) == ImmutableMatrix(eye(2))
    assert ImmutableMatrix(c) == c.as_immutable()
    assert Matrix(ImmutableMatrix(c)) == ImmutableMatrix(c).as_mutable()

    assert c is not Matrix(c)


def test_tolist():
    lst = [[1, Rational(1, 2), x*y, 0], [x, y, z, x**2], [y, -1, z*x, 3]]
    m = Matrix(lst)
    assert m.tolist() == lst


def test_as_mutable():
    assert zeros(0, 3).as_mutable() == zeros(0, 3)
    assert zeros(0, 3).as_immutable() == ImmutableMatrix(zeros(0, 3))
    assert zeros(3, 0).as_immutable() == ImmutableMatrix(zeros(3, 0))


def test_determinant():

    for M in [Matrix(), Matrix([[1]])]:
        assert (
            M.det() ==
            M.det_bareis() ==
            M.berkowitz_det() ==
            M.det_LU_decomposition() ==
            1)

    M = Matrix(( (-3,  2),
                 ( 8, -5) ))

    assert M.det(method="bareis") == -1
    assert M.det(method="berkowitz") == -1

    M = Matrix(( (x,   1),
                 (y, 2*y) ))

    assert M.det(method="bareis") == 2*x*y - y
    assert M.det(method="berkowitz") == 2*x*y - y

    M = Matrix(( (1, 1, 1),
                 (1, 2, 3),
                 (1, 3, 6) ))

    assert M.det(method="bareis") == 1
    assert M.det(method="berkowitz") == 1

    M = Matrix(( ( 3, -2,  0, 5),
                 (-2,  1, -2, 2),
                 ( 0, -2,  5, 0),
                 ( 5,  0,  3, 4) ))

    assert M.det(method="bareis") == -289
    assert M.det(method="berkowitz") == -289

    M = Matrix(( ( 1,  2,  3,  4),
                 ( 5,  6,  7,  8),
                 ( 9, 10, 11, 12),
                 (13, 14, 15, 16) ))

    assert M.det(method="bareis") == 0
    assert M.det(method="berkowitz") == 0

    M = Matrix(( (3, 2, 0, 0, 0),
                 (0, 3, 2, 0, 0),
                 (0, 0, 3, 2, 0),
                 (0, 0, 0, 3, 2),
                 (2, 0, 0, 0, 3) ))

    assert M.det(method="bareis") == 275
    assert M.det(method="berkowitz") == 275

    M = Matrix(( (1, 0,  1,  2, 12),
                 (2, 0,  1,  1,  4),
                 (2, 1,  1, -1,  3),
                 (3, 2, -1,  1,  8),
                 (1, 1,  1,  0,  6) ))

    assert M.det(method="bareis") == -55
    assert M.det(method="berkowitz") == -55

    M = Matrix(( (-5,  2,  3,  4,  5),
                 ( 1, -4,  3,  4,  5),
                 ( 1,  2, -3,  4,  5),
                 ( 1,  2,  3, -2,  5),
                 ( 1,  2,  3,  4, -1) ))

    assert M.det(method="bareis") == 11664
    assert M.det(method="berkowitz") == 11664

    M = Matrix(( ( 2,  7, -1, 3, 2),
                 ( 0,  0,  1, 0, 1),
                 (-2,  0,  7, 0, 2),
                 (-3, -2,  4, 5, 3),
                 ( 1,  0,  0, 0, 1) ))

    assert M.det(method="bareis") == 123
    assert M.det(method="berkowitz") == 123

    M = Matrix(( (x, y, z),
                 (1, 0, 0),
                 (y, z, x) ))

    assert M.det(method="bareis") == z**2 - x*y
    assert M.det(method="berkowitz") == z**2 - x*y


def test_det_LU_decomposition():

    for M in [Matrix(), Matrix([[1]])]:
        assert M.det(method="det_LU") == 1

    M = Matrix(( (-3,  2),
                 ( 8, -5) ))

    assert M.det(method="det_LU") == -1

    M = Matrix(( (x,   1),
                 (y, 2*y) ))

    assert M.det(method="det_LU") == 2*x*y - y

    M = Matrix(( (1, 1, 1),
                 (1, 2, 3),
                 (1, 3, 6) ))

    assert M.det(method="det_LU") == 1

    M = Matrix(( ( 3, -2,  0, 5),
                 (-2,  1, -2, 2),
                 ( 0, -2,  5, 0),
                 ( 5,  0,  3, 4) ))

    assert M.det(method="det_LU") == -289

    M = Matrix(( (3, 2, 0, 0, 0),
                 (0, 3, 2, 0, 0),
                 (0, 0, 3, 2, 0),
                 (0, 0, 0, 3, 2),
                 (2, 0, 0, 0, 3) ))

    assert M.det(method="det_LU") == 275

    M = Matrix(( (1, 0,  1,  2, 12),
                 (2, 0,  1,  1,  4),
                 (2, 1,  1, -1,  3),
                 (3, 2, -1,  1,  8),
                 (1, 1,  1,  0,  6) ))

    assert M.det(method="det_LU") == -55

    M = Matrix(( (-5,  2,  3,  4,  5),
                 ( 1, -4,  3,  4,  5),
                 ( 1,  2, -3,  4,  5),
                 ( 1,  2,  3, -2,  5),
                 ( 1,  2,  3,  4, -1) ))

    assert M.det(method="det_LU") == 11664

    M = Matrix(( ( 2,  7, -1, 3, 2),
                 ( 0,  0,  1, 0, 1),
                 (-2,  0,  7, 0, 2),
                 (-3, -2,  4, 5, 3),
                 ( 1,  0,  0, 0, 1) ))

    assert M.det(method="det_LU") == 123

    M = Matrix(( (x, y, z),
                 (1, 0, 0),
                 (y, z, x) ))

    assert M.det(method="det_LU") == z**2 - x*y


def test_berkowitz_minors():
    B = Matrix(2, 2, [1, 2, 2, 1])

    assert B.berkowitz_minors() == (1, -3)


def test_slicing():
    m0 = eye(4)
    assert m0[:3, :3] == eye(3)
    assert m0[2:4, 0:2] == zeros(2)

    m1 = Matrix(3, 3, lambda i, j: i + j)
    assert m1[0, :] == Matrix(1, 3, (0, 1, 2))
    assert m1[1:3, 1] == Matrix(2, 1, (2, 3))

    m2 = Matrix([[0, 1, 2, 3], [4, 5, 6, 7], [8, 9, 10, 11], [12, 13, 14, 15]])
    assert m2[:, -1] == Matrix(4, 1, [3, 7, 11, 15])
    assert m2[-2:, :] == Matrix([[8, 9, 10, 11], [12, 13, 14, 15]])


def test_submatrix_assignment():
    m = zeros(4)
    m[2:4, 2:4] = eye(2)
    assert m == Matrix(((0, 0, 0, 0),
                        (0, 0, 0, 0),
                        (0, 0, 1, 0),
                        (0, 0, 0, 1)))
    m[:2, :2] = eye(2)
    assert m == eye(4)
    m[:, 0] = Matrix(4, 1, (1, 2, 3, 4))
    assert m == Matrix(((1, 0, 0, 0),
                        (2, 1, 0, 0),
                        (3, 0, 1, 0),
                        (4, 0, 0, 1)))
    m[:, :] = zeros(4)
    assert m == zeros(4)
    m[:, :] = [(1, 2, 3, 4), (5, 6, 7, 8), (9, 10, 11, 12), (13, 14, 15, 16)]
    assert m == Matrix(((1, 2, 3, 4),
                        (5, 6, 7, 8),
                        (9, 10, 11, 12),
                        (13, 14, 15, 16)))
    m[:2, 0] = [0, 0]
    assert m == Matrix(((0, 2, 3, 4),
                        (0, 6, 7, 8),
                        (9, 10, 11, 12),
                        (13, 14, 15, 16)))

    with pytest.raises(ValueError):
        m[:, 1] = object()


def test_extract():
    m = Matrix(4, 3, lambda i, j: i*3 + j)
    assert m.extract([0, 1, 3], [0, 1]) == Matrix(3, 2, [0, 1, 3, 4, 9, 10])
    assert m.extract([0, 3], [0, 0, 2]) == Matrix(2, 3, [0, 0, 2, 9, 9, 11])
    assert m.extract(range(4), range(3)) == m
    pytest.raises(IndexError, lambda: m.extract([4], [0]))
    pytest.raises(IndexError, lambda: m.extract([0], [3]))


def test_reshape():
    m0 = eye(3)
    assert m0.reshape(1, 9) == Matrix(1, 9, (1, 0, 0, 0, 1, 0, 0, 0, 1))
    m1 = Matrix(3, 4, lambda i, j: i + j)
    assert m1.reshape(
        4, 3) == Matrix(((0, 1, 2), (3, 1, 2), (3, 4, 2), (3, 4, 5)))
    assert m1.reshape(2, 6) == Matrix(((0, 1, 2, 3, 1, 2), (3, 4, 2, 3, 4, 5)))


def test_applyfunc():
    m0 = eye(3)
    assert m0.applyfunc(lambda x: 2*x) == eye(3)*2
    assert m0.applyfunc(lambda x: 0) == zeros(3)


def test_expand():
    m0 = Matrix([[x*(x + y), 2], [((x + y)*y)*x, x*(y + x*(x + y))]])
    # Test if expand() returns a matrix
    m1 = m0.expand()
    assert m1 == Matrix(
        [[x*y + x**2, 2], [x*y**2 + y*x**2, x*y + y*x**2 + x**3]])

    a = Symbol('a', extended_real=True)

    assert Matrix([exp(I*a)]).expand(complex=True) == \
        Matrix([cos(a) + I*sin(a)])

    assert Matrix([[0, 1, 2], [0, 0, -1], [0, 0, 0]]).exp() == Matrix([
        [1, 1, Rational(3, 2)],
        [0, 1, -1],
        [0, 0, 1]]
    )


def test_random():
    M = randMatrix(3, 3)
    M = randMatrix(3, 3, seed=3)
    M = randMatrix(3, 4, 0, 150)
    M = randMatrix(3, symmetric=True, percent=50)
    M = randMatrix(3, symmetric=True)
    S = M.copy()
    S.simplify()
    assert S == M  # doesn't fail when elements are Numbers, not int


def test_LUdecomp():
    testmat = Matrix([[0, 2, 5, 3],
                      [3, 3, 7, 4],
                      [8, 4, 0, 2],
                      [-2, 6, 3, 4]])
    L, U, p = testmat.LUdecomposition()
    assert L.is_lower
    assert U.is_upper
    assert (L*U).permuteBkwd(p) - testmat == zeros(4)

    testmat = Matrix([[6, -2, 7, 4],
                      [0, 3, 6, 7],
                      [1, -2, 7, 4],
                      [-9, 2, 6, 3]])
    L, U, p = testmat.LUdecomposition()
    assert L.is_lower
    assert U.is_upper
    assert (L*U).permuteBkwd(p) - testmat == zeros(4)

    M = Matrix(((1, x, 1), (2, y, 0), (y, 0, z)))
    L, U, p = M.LUdecomposition()
    assert L.is_lower
    assert U.is_upper
    assert (L*U).permuteBkwd(p) - M == zeros(3)

    mL = Matrix((
        (1, 0, 0),
        (2, 3, 0),
    ))
    assert mL.is_lower is True
    assert mL.is_upper is False
    mU = Matrix((
        (1, 2, 3),
        (0, 4, 5),
    ))
    assert mU.is_lower is False
    assert mU.is_upper is True

    # test FF LUdecomp
    M = Matrix([[1, 3, 3],
                [3, 2, 6],
                [3, 2, 2]])
    P, L, Dee, U = M.LUdecompositionFF()
    assert P*M == L*Dee.inv()*U

    M = Matrix([[1,  2, 3,  4],
                [3, -1, 2,  3],
                [3,  1, 3, -2],
                [6, -1, 0,  2]])
    P, L, Dee, U = M.LUdecompositionFF()
    assert P*M == L*Dee.inv()*U

    M = Matrix([[0, 0, 1],
                [2, 3, 0],
                [3, 1, 4]])
    P, L, Dee, U = M.LUdecompositionFF()
    assert P*M == L*Dee.inv()*U

    pytest.raises(ValueError,
                  lambda: Matrix([[0, 2], [0, 0]]).LUdecompositionFF())


def test_LUsolve():
    A = Matrix([[2, 3, 5],
                [3, 6, 2],
                [8, 3, 6]])
    x = Matrix(3, 1, [3, 7, 5])
    b = A*x
    soln = A.LUsolve(b)
    assert soln == x
    A = Matrix([[0, -1, 2],
                [5, 10, 7],
                [8,  3, 4]])
    x = Matrix(3, 1, [-1, 2, 5])
    b = A*x
    soln = A.LUsolve(b)
    assert soln == x


def test_QRsolve():
    A = Matrix([[2, 3, 5],
                [3, 6, 2],
                [8, 3, 6]])
    x = Matrix(3, 1, [3, 7, 5])
    b = A*x
    soln = A.QRsolve(b)
    assert soln == x
    x = Matrix([[1, 2], [3, 4], [5, 6]])
    b = A*x
    soln = A.QRsolve(b)
    assert soln == x

    A = Matrix([[0, -1, 2],
                [5, 10, 7],
                [8,  3, 4]])
    x = Matrix(3, 1, [-1, 2, 5])
    b = A*x
    soln = A.QRsolve(b)
    assert soln == x
    x = Matrix([[7, 8], [9, 10], [11, 12]])
    b = A*x
    soln = A.QRsolve(b)
    assert soln == x


def test_inverse():
    A = eye(4)
    assert A.inv() == eye(4)
    assert A.inv(method="LU") == eye(4)
    assert A.inv(method="ADJ") == eye(4)
    A = Matrix([[2, 3, 5],
                [3, 6, 2],
                [8, 3, 6]])
    Ainv = A.inv()
    assert A*Ainv == eye(3)
    assert A.inv(method="LU") == Ainv
    assert A.inv(method="ADJ") == Ainv

    pytest.raises(ValueError, lambda: A.inv(method="SPAM"))

    # test that immutability is not a problem
    cls = ImmutableMatrix
    m = cls([[48, 49, 31],
             [ 9, 71, 94],
             [59, 28, 65]])
    assert all(type(m.inv(s)) is cls for s in 'GE ADJ LU'.split())
    cls = ImmutableSparseMatrix
    m = cls([[48, 49, 31],
             [ 9, 71, 94],
             [59, 28, 65]])
    assert all(type(m.inv(s)) is cls for s in 'CH LDL'.split())

    A = ImmutableSparseMatrix([[20682, 12468,  9899, 12111],
                               [12468, 18618, 14513,  9401],
                               [ 9899, 14513, 12770,  6823],
                               [12111,  9401,  6823,  7894]])
    Ainv = A.inv('CH')
    assert A.inv('LDL') == Ainv
    assert A*Ainv == Ainv*A
    assert Matrix(Ainv*A) == eye(4)


def test_matrix_inverse_mod():
    A = Matrix(2, 1, [1, 0])
    pytest.raises(NonSquareMatrixError, lambda: A.inv_mod(2))
    A = Matrix(2, 2, [1, 0, 0, 0])
    pytest.raises(ValueError, lambda: A.inv_mod(2))
    A = Matrix(2, 2, [1, 2, 3, 4])
    Ai = Matrix(2, 2, [1, 1, 0, 1])
    assert A.inv_mod(3) == Ai
    A = Matrix(2, 2, [1, 0, 0, 1])
    assert A.inv_mod(2) == A


def test_util():
    R = Rational

    v1 = Matrix(1, 3, [1, 2, 3])
    v2 = Matrix(1, 3, [3, 4, 5])
    assert v1.norm() == sqrt(14)
    assert v1.project(v2) == Matrix(1, 3, [R(39)/25, R(52)/25, R(13)/5])
    assert Matrix.zeros(1, 2) == Matrix(1, 2, [0, 0])
    assert ones(1, 2) == Matrix(1, 2, [1, 1])
    assert v1.copy() == v1
    # cofactor
    assert eye(3) == eye(3).cofactorMatrix()
    test = Matrix([[1, 3, 2], [2, 6, 3], [2, 3, 6]])
    assert test.cofactorMatrix() == \
        Matrix([[27, -6, -6], [-12, 2, 3], [-3, 1, 0]])
    test = Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    assert test.cofactorMatrix() == \
        Matrix([[-3, 6, -3], [6, -12, 6], [-3, 6, -3]])

    printer = StrPrinter()
    assert Matrix().table(printer) == '[]'


def test_jacobian_hessian():
    L = Matrix(1, 2, [x**2*y, 2*y**2 + x*y])
    syms = [x, y]
    assert L.jacobian(syms) == Matrix([[2*x*y, x**2], [y, 4*y + x]])

    L = Matrix(1, 2, [x, x**2*y**3])
    assert L.jacobian(syms) == Matrix([[1, 0], [2*x*y**3, x**2*3*y**2]])

    f = x**2*y
    syms = [x, y]
    assert hessian(f, syms) == Matrix([[2*y, 2*x], [2*x, 0]])

    f = x**2*y**3
    assert hessian(f, syms) == \
        Matrix([[2*y**3, 6*x*y**2], [6*x*y**2, 6*x**2*y]])

    f = z + x*y**2
    g = x**2 + 2*y**3
    ans = Matrix([[0,   2*y],
                  [2*y, 2*x]])
    assert ans == hessian(f, Matrix([x, y]))
    assert ans == hessian(f, Matrix([x, y]).T)
    assert hessian(f, (y, x), [g]) == Matrix([
        [     0, 6*y**2, 2*x],
        [6*y**2,    2*x, 2*y],
        [   2*x,    2*y,   0]])

    pytest.raises(ValueError, lambda: hessian(Basic(), Matrix([x, y])))
    pytest.raises(ValueError, lambda: hessian(f, Matrix([x, y]),
                                              Matrix([Basic(), y])))


def test_QR():
    A = Matrix([[1, 2], [2, 3]])
    Q, S = A.QRdecomposition()
    R = Rational
    assert Q == Matrix([
        [  5**R(-1, 2),  (R(2)/5)*(R(1)/5)**R(-1, 2)],
        [2*5**R(-1, 2), (-R(1)/5)*(R(1)/5)**R(-1, 2)]])
    assert S == Matrix([[5**R(1, 2), 8*5**R(-1, 2)], [0, (R(1)/5)**R(1, 2)]])
    assert Q*S == A
    assert Q.T * Q == eye(2)

    A = Matrix([[1, 1, 1], [1, 1, 3], [2, 3, 4]])
    Q, R = A.QRdecomposition()
    assert Q.T * Q == eye(Q.cols)
    assert R.is_upper
    assert A == Q*R


def test_QR_non_square():
    A = Matrix([[9, 0, 26], [12, 0, -7], [0, 4, 4], [0, -3, -3]])
    Q, R = A.QRdecomposition()
    assert Q.T * Q == eye(Q.cols)
    assert R.is_upper
    assert A == Q*R

    A = Matrix([[1, -1, 4], [1, 4, -2], [1, 4, 2], [1, -1, 0]])
    Q, R = A.QRdecomposition()
    assert Q.T * Q == eye(Q.cols)
    assert R.is_upper
    assert A == Q*R


def test_nullspace():
    # first test reduced row-ech form
    R = Rational

    M = Matrix([[5, 7, 2,  1],
                [1, 6, 2, -1]])
    out, tmp = M.rref()
    assert out == Matrix([[1, 0, -R(2)/23, R(13)/23],
                          [0, 1,  R(8)/23, R(-6)/23]])

    M = Matrix([[-5, -1,  4, -3, -1],
                [ 1, -1, -1,  1,  0],
                [-1,  0,  0,  0,  0],
                [ 4,  1, -4,  3,  1],
                [-2,  0,  2, -2, -1]])
    assert M*M.nullspace()[0] == Matrix(5, 1, [0]*5)

    M = Matrix([[ 1,  3, 0,  2,  6, 3, 1],
                [-2, -6, 0, -2, -8, 3, 1],
                [ 3,  9, 0,  0,  6, 6, 2],
                [-1, -3, 0,  1,  0, 9, 3]])
    out, tmp = M.rref()
    assert out == Matrix([[1, 3, 0, 0, 2, 0, 0],
                          [0, 0, 0, 1, 2, 0, 0],
                          [0, 0, 0, 0, 0, 1, R(1)/3],
                          [0, 0, 0, 0, 0, 0, 0]])

    # now check the vectors
    basis = M.nullspace()
    assert basis[0] == Matrix([-3, 1, 0, 0, 0, 0, 0])
    assert basis[1] == Matrix([0, 0, 1, 0, 0, 0, 0])
    assert basis[2] == Matrix([-2, 0, 0, -2, 1, 0, 0])
    assert basis[3] == Matrix([0, 0, 0, 0, 0, R(-1)/3, 1])
    assert basis == M.nullspace(simplify=True)

    # issue sympy/sympy#4797; just see that we can do it when rows > cols
    M = Matrix([[1, 2], [2, 4], [3, 6]])
    assert M.nullspace()


def test_wronskian():
    assert wronskian([cos(x), sin(x)], x) == cos(x)**2 + sin(x)**2
    assert wronskian([exp(x), exp(2*x)], x) == exp(3*x)
    assert wronskian([exp(x), x], x) == exp(x) - x*exp(x)
    assert wronskian([1, x, x**2], x) == 2
    w1 = -6*exp(x)*sin(x)*x + 6*cos(x)*exp(x)*x**2 - 6*exp(x)*cos(x)*x - \
        exp(x)*cos(x)*x**3 + exp(x)*sin(x)*x**3
    assert wronskian([exp(x), cos(x), x**3], x).expand() == w1
    assert wronskian([exp(x), cos(x), x**3], x, method='berkowitz').expand() \
        == w1
    w2 = -x**3*cos(x)**2 - x**3*sin(x)**2 - 6*x*cos(x)**2 - 6*x*sin(x)**2
    assert wronskian([sin(x), cos(x), x**3], x).expand() == w2
    assert wronskian([sin(x), cos(x), x**3], x, method='berkowitz').expand() \
        == w2
    assert wronskian([], x) == 1


def test_eigen():
    R = Rational

    assert eye(3).charpoly(x) == Poly((x - 1)**3, x)
    assert eye(3).charpoly(y) == Poly((y - 1)**3, y)

    M = Matrix([[1, 0, 0],
                [0, 1, 0],
                [0, 0, 1]])

    assert M.eigenvals(multiple=False) == {1: 3}

    assert M.eigenvects() == (
        [(1, 3, [Matrix([1, 0, 0]),
                 Matrix([0, 1, 0]),
                 Matrix([0, 0, 1])])])

    M = Matrix([[0, 1, 1],
                [1, 0, 0],
                [1, 1, 1]])

    assert M.eigenvals() == {2: 1, -1: 1, 0: 1}

    assert M.eigenvects() == (
        [
            (-1, 1, [Matrix([-1, 1, 0])]),
            ( 0, 1, [Matrix([0, -1, 1])]),
            ( 2, 1, [Matrix([R(2, 3), R(1, 3), 1])])
        ])

    M = Matrix([[x, 0],
                [0, 1]])

    assert M.eigenvals() == {x: 1, 1: 1}

    M = Matrix([[1, -1],
                [1,  3]])
    assert M.eigenvects() == ([(2, 2, [Matrix(2, 1, [-1, 1])])])

    M = Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    a = R(15, 2)
    b = 3*33**R(1, 2)
    c = R(13, 2)
    d = (R(33, 8) + 3*b/8)
    e = (R(33, 8) - 3*b/8)

    def NS(e, n):
        return str(N(e, n))
    r = [
        (a - b/2, 1, [Matrix([(12 + 24/(c - b/2))/((c - b/2)*e) + 3/(c - b/2),
                              (6 + 12/(c - b/2))/e, 1])]),
        (      0, 1, [Matrix([1, -2, 1])]),
        (a + b/2, 1, [Matrix([(12 + 24/(c + b/2))/((c + b/2)*d) + 3/(c + b/2),
                              (6 + 12/(c + b/2))/d, 1])]),
    ]
    r1 = [(NS(r[i][0], 2), NS(r[i][1], 2),
           [NS(j, 2) for j in r[i][2][0]]) for i in range(len(r))]
    r = M.eigenvects()
    r2 = [(NS(r[i][0], 2), NS(r[i][1], 2),
           [NS(j, 2) for j in r[i][2][0]]) for i in range(len(r))]
    assert sorted(r1) == sorted(r2)

    eps = Symbol('eps', extended_real=True)

    M = Matrix([[abs(eps), I*eps    ],
                [-I*eps,   abs(eps) ]])

    assert M.eigenvects() == (
        [
            ( 0, 1, [Matrix([[-I*eps/abs(eps)], [1]])]),
            ( 2*abs(eps), 1, [ Matrix([[I*eps/abs(eps)], [1]]) ] ),
        ])

    M = Matrix(3, 3, [1, 2, 0, 0, 3, 0, 2, -4, 2])
    M._eigenvects = M.eigenvects(simplify=False)
    assert max(i.q for i in M._eigenvects[0][2][0]) > 1
    M._eigenvects = M.eigenvects(simplify=True)
    assert max(i.q for i in M._eigenvects[0][2][0]) == 1
    M = Matrix([[Rational(1, 4), 1], [1, 1]])
    assert M.eigenvects(simplify=True) == [
        (Rational(5, 8) + sqrt(73)/8, 1, [Matrix([[8/(3 + sqrt(73))], [1]])]),
        (-sqrt(73)/8 + Rational(5, 8), 1, [Matrix([[8/(-sqrt(73) + 3)], [1]])])]
    assert M.eigenvects(simplify=False) == [
        (Rational(5, 8) + sqrt(73)/8, 1,
         [Matrix([[-1/(-sqrt(73)/8 + Rational(-3, 8))], [1]])]),
        (-sqrt(73)/8 + Rational(5, 8), 1,
         [Matrix([[-1/(Rational(-3, 8) + sqrt(73)/8)], [1]])]),
    ]

    m = Matrix([[1, .6, .6], [.6, .9, .9], [.9, .6, .6]])
    evals = {-sqrt(385)/20 + Rational(5, 4): 1, sqrt(385)/20 + Rational(5, 4): 1, 0: 1}
    assert m.eigenvals() == evals
    nevals = list(sorted(m.eigenvals(rational=False)))
    sevals = list(sorted(evals))
    assert all(abs(nevals[i] - sevals[i]) < 1e-9 for i in range(len(nevals)))

    # issue sympy/sympy#10719
    assert Matrix([]).eigenvals() == {}
    assert Matrix([]).eigenvects() == []


def test_subs():
    assert Matrix([[1, x], [x, 4]]).subs(x, 5) == Matrix([[1, 5], [5, 4]])
    assert Matrix([[x, 2], [x + y, 4]]).subs([[x, -1], [y, -2]]) == \
        Matrix([[-1, 2], [-3, 4]])
    assert Matrix([[x, 2], [x + y, 4]]).subs([(x, -1), (y, -2)]) == \
        Matrix([[-1, 2], [-3, 4]])
    assert Matrix([[x, 2], [x + y, 4]]).subs({x: -1, y: -2}) == \
        Matrix([[-1, 2], [-3, 4]])
    assert Matrix([x*y]).subs({x: y - 1, y: x - 1}, simultaneous=True) == \
        Matrix([(x - 1)*(y - 1)])

    for cls in classes:
        assert Matrix([[2, 0], [0, 2]]) == cls.eye(2).subs(1, 2)


def test_xreplace():
    assert Matrix([[1, x], [x, 4]]).xreplace({x: 5}) == \
        Matrix([[1, 5], [5, 4]])
    assert Matrix([[x, 2], [x + y, 4]]).xreplace({x: -1, y: -2}) == \
        Matrix([[-1, 2], [-3, 4]])
    for cls in classes:
        assert Matrix([[2, 0], [0, 2]]) == cls.eye(2).xreplace({1: 2})


def test_simplify():
    f = Function('f')

    m = Matrix([[1, x], [x + 1/x, x - 1]])
    m = m.row_join(eye(m.cols))
    raw = m.rref(simplify=lambda x: x)[0]
    assert raw != m.rref(simplify=True)[0]

    M = Matrix([[            1/x + 1/y,                 (x + x*y) / x  ],
                [ (f(x) + y*f(x))/f(x), 2 * (1/n - cos(n * pi)/n) / pi ]])
    M.simplify()
    assert M == Matrix([[ (x + y)/(x * y),                        1 + y ],
                        [           1 + y, 2*((1 - 1*cos(pi*n))/(pi*n)) ]])
    eq = (1 + x)**2
    M = Matrix([[eq]])
    M.simplify()
    assert M == Matrix([[eq]])
    M.simplify(ratio=oo) == M
    assert M == Matrix([[eq.simplify(ratio=oo)]])


def test_transpose():
    M = Matrix([[1, 2, 3, 4, 5, 6, 7, 8, 9, 0],
                [1, 2, 3, 4, 5, 6, 7, 8, 9, 0]])
    assert M.T == Matrix( [ [1, 1],
                            [2, 2],
                            [3, 3],
                            [4, 4],
                            [5, 5],
                            [6, 6],
                            [7, 7],
                            [8, 8],
                            [9, 9],
                            [0, 0] ])
    assert M.T.T == M
    assert M.T == M.transpose()


def test_conjugate():
    M = Matrix([[0, I, 5],
                [1, 2, 0]])

    assert M.T == Matrix([[0, 1],
                          [I, 2],
                          [5, 0]])

    assert M.C == Matrix([[0, -I, 5],
                          [1,  2, 0]])
    assert M.C == M.conjugate()

    assert M.H == M.T.C
    assert M.H == Matrix([[ 0, 1],
                          [-I, 2],
                          [ 5, 0]])


def test_conj_dirac():
    pytest.raises(AttributeError, lambda: eye(3).D)

    M = Matrix([[1, I, I, I],
                [0, 1, I, I],
                [0, 0, 1, I],
                [0, 0, 0, 1]])

    assert M.D == Matrix([[ 1,  0,  0,  0],
                          [-I,  1,  0,  0],
                          [-I, -I, -1,  0],
                          [-I, -I,  I, -1]])


def test_trace():
    M = Matrix([[1, 0, 0],
                [0, 5, 0],
                [0, 0, 8]])
    assert M.trace() == 14


def test_shape():
    M = Matrix([[x, 0, 0],
                [0, y, 0]])
    assert M.shape == (2, 3)


def test_col_row_op():
    M = Matrix([[x, 0, 0],
                [0, y, 0]])
    M.row_op(1, lambda r, j: r + j + 1)
    assert M == Matrix([[x,     0, 0],
                        [1, y + 2, 3]])

    M.col_op(0, lambda c, j: c + y**j)
    assert M == Matrix([[x + 1,     0, 0],
                        [1 + y, y + 2, 3]])

    # neither row nor slice give copies that allow the original matrix to
    # be changed
    assert M.row(0) == Matrix([[x + 1, 0, 0]])
    r1 = M.row(0)
    r1[0] = 42
    assert M[0, 0] == x + 1
    r1 = M[0, :-1]  # also testing negative slice
    r1[0] = 42
    assert M[0, 0] == x + 1
    c1 = M.col(0)
    assert c1 == Matrix([x + 1, 1 + y])
    c1[0] = 0
    assert M[0, 0] == x + 1
    c1 = M[:, 0]
    c1[0] = 42
    assert M[0, 0] == x + 1


def test_zip_row_op():
    for cls in classes[:2]:  # XXX: immutable matrices don't support row ops
        M = cls.eye(3)
        M.zip_row_op(1, 0, lambda v, u: v + 2*u)
        assert M == cls([[1, 0, 0],
                         [2, 1, 0],
                         [0, 0, 1]])

        M = cls.eye(3)*2
        M[0, 1] = -1
        M.zip_row_op(1, 0, lambda v, u: v + 2*u)
        assert M == cls([[2, -1, 0],
                         [4,  0, 0],
                         [0,  0, 2]])


def test_sympyissue_3950():
    m = Matrix([1, 2, 3])
    a = Matrix([1, 2, 3])
    b = Matrix([2, 2, 3])
    assert not (m in [])
    assert not (m in [1])
    assert m != 1
    assert m == a
    assert m != b


def test_sympyissue_3981():
    class Index1:
        def __index__(self):
            return 1

    class Index2:
        def __index__(self):
            return 2
    index1 = Index1()
    index2 = Index2()

    m = Matrix([1, 2, 3])

    assert m[index2] == 3

    m[index2] = 5
    assert m[2] == 5

    m = Matrix([[1, 2, 3], [4, 5, 6]])
    assert m[index1, index2] == 6
    assert m[1, index2] == 6
    assert m[index1, 2] == 6

    m[index1, index2] = 4
    assert m[1, 2] == 4
    m[1, index2] = 6
    assert m[1, 2] == 6
    m[index1, 2] = 8
    assert m[1, 2] == 8


def test_evalf():
    a = Matrix([sqrt(5), 6])
    assert all(a.evalf()[i] == a[i].evalf() for i in range(2))
    assert all(a.evalf(2)[i] == a[i].evalf(2) for i in range(2))
    assert all(a.n(2)[i] == a[i].n(2) for i in range(2))


def test_is_symbolic():
    a = Matrix([[x, x], [x, x]])
    assert a.is_symbolic() is True
    a = Matrix([[1, 2, 3, 4], [5, 6, 7, 8]])
    assert a.is_symbolic() is False
    a = Matrix([[1, 2, 3, 4], [5, 6, x, 8]])
    assert a.is_symbolic() is True
    a = Matrix([[1, x, 3]])
    assert a.is_symbolic() is True
    a = Matrix([[1, 2, 3]])
    assert a.is_symbolic() is False
    a = Matrix([[1], [x], [3]])
    assert a.is_symbolic() is True
    a = Matrix([[1], [2], [3]])
    assert a.is_symbolic() is False


def test_is_upper():
    a = Matrix([[1, 2, 3]])
    assert a.is_upper is True
    a = Matrix([[1], [2], [3]])
    assert a.is_upper is False


def test_is_lower():
    a = Matrix([[1, 2, 3]])
    assert a.is_lower is False
    a = Matrix([[1], [2], [3]])
    assert a.is_lower is True


def test_is_nilpotent():
    a = Matrix(4, 4, [0, 2, 1, 6, 0, 0, 1, 2, 0, 0, 0, 3, 0, 0, 0, 0])
    assert a.is_nilpotent()
    a = Matrix([[1, 0], [0, 1]])
    assert not a.is_nilpotent()
    a = Matrix([])
    assert a.is_nilpotent()


def test_zeros_ones_fill():
    n, m = 3, 5

    a = zeros(n, m)
    a.fill( 5 )

    b = 5 * ones(n, m)

    assert a == b
    assert a.rows == b.rows == 3
    assert a.cols == b.cols == 5
    assert a.shape == b.shape == (3, 5)
    assert zeros(2) == zeros(2, 2)
    assert ones(2) == ones(2, 2)
    assert zeros(2, 3) == Matrix(2, 3, [0]*6)
    assert ones(2, 3) == Matrix(2, 3, [1]*6)


def test_empty_zeros():
    a = zeros(0)
    assert a == Matrix()
    a = zeros(0, 2)
    assert a.rows == 0
    assert a.cols == 2
    a = zeros(2, 0)
    assert a.rows == 2
    assert a.cols == 0


def test_sympyissue_3749():
    a = Matrix([[x**2, x*y], [x*sin(y), x*cos(y)]])
    assert a.diff(x) == Matrix([[2*x, y], [sin(y), cos(y)]])
    assert Matrix([
        [x, -x, x**2],
        [exp(x), 1/x - exp(-x), x + 1/x]]).limit(x, oo) == \
        Matrix([[oo, -oo, oo], [oo, 0, oo]])
    assert Matrix([
        [(exp(x) - 1)/x, 2*x + y*x, x**x ],
        [1/x, abs(x), abs(sin(x + 1))]]).limit(x, 0) == \
        Matrix([[1, 0, 1], [oo, 0, sin(1)]])
    assert a.integrate(x) == Matrix([
        [Rational(1, 3)*x**3, y*x**2/2],
        [x**2*sin(y)/2, x**2*cos(y)/2]])


def test_inv_iszerofunc():
    A = eye(4)
    A.col_swap(0, 1)
    for method in "GE", "LU":
        assert A.inv(method=method, iszerofunc=lambda x: x == 0) == \
            A.inv(method="ADJ")


def test_jacobian_metrics():
    rho, phi = symbols("rho,phi")
    X = Matrix([rho*cos(phi), rho*sin(phi)])
    Y = Matrix([rho, phi])
    J = X.jacobian(Y)
    assert J == X.jacobian(Y.T)
    assert J == (X.T).jacobian(Y)
    assert J == (X.T).jacobian(Y.T)
    g = J.T*eye(J.shape[0])*J
    g = g.applyfunc(trigsimp)
    assert g == Matrix([[1, 0], [0, rho**2]])


def test_jacobian2():
    rho, phi = symbols("rho,phi")
    X = Matrix([rho*cos(phi), rho*sin(phi), rho**2])
    Y = Matrix([rho, phi])
    J = Matrix([
        [cos(phi), -rho*sin(phi)],
        [sin(phi),  rho*cos(phi)],
        [   2*rho,             0],
    ])
    assert X.jacobian(Y) == J


def test_sympyissue_4564():
    X = Matrix([exp(x + y + z), exp(x + y + z), exp(x + y + z)])
    Y = Matrix([x, y, z])
    for i in range(1, 3):
        for j in range(1, 3):
            X_slice = X[:i, :]
            Y_slice = Y[:j, :]
            J = X_slice.jacobian(Y_slice)
            assert J.rows == i
            assert J.cols == j
            for k in range(j):
                assert J[:, k] == X_slice


def test_nonvectorJacobian():
    X = Matrix([[exp(x + y + z), exp(x + y + z)],
                [exp(x + y + z), exp(x + y + z)]])
    pytest.raises(TypeError, lambda: X.jacobian(Matrix([x, y, z])))
    X = X[0, :]
    Y = Matrix([[x, y], [x, z]])
    pytest.raises(TypeError, lambda: X.jacobian(Y))
    pytest.raises(TypeError, lambda: X.jacobian(Matrix([ [x, y], [x, z] ])))


def test_vec():
    m = Matrix([[1, 3], [2, 4]])
    m_vec = m.vec()
    assert m_vec.cols == 1
    for i in range(4):
        assert m_vec[i] == i + 1


def test_vech():
    m = Matrix([[1, 2], [2, 3]])
    m_vech = m.vech()
    assert m_vech.cols == 1
    for i in range(3):
        assert m_vech[i] == i + 1
    m_vech = m.vech(diagonal=False)
    assert m_vech[0] == 2

    m = Matrix([[1, x*(x + y)], [y*x + x**2, 1]])
    m_vech = m.vech(diagonal=False)
    assert m_vech[0] == x*(x + y)

    m = Matrix([[1, x*(x + y)], [y*x, 1]])
    m_vech = m.vech(diagonal=False, check_symmetry=False)
    assert m_vech[0] == y*x


def test_vech_errors():
    m = Matrix([[1, 3]])
    pytest.raises(ShapeError, lambda: m.vech())
    m = Matrix([[1, 3], [2, 4]])
    pytest.raises(ValueError, lambda: m.vech())
    pytest.raises(ShapeError, lambda: Matrix([ [1, 3] ]).vech())
    pytest.raises(ValueError, lambda: Matrix([ [1, 3], [2, 4] ]).vech())


def test_diag():
    a = Matrix([[1, 2], [2, 3]])
    b = Matrix([[3, x], [y, 3]])
    c = Matrix([[3, x, 3], [y, 3, z], [x, y, z]])
    assert diag(a, b, b) == Matrix([
        [1, 2, 0, 0, 0, 0],
        [2, 3, 0, 0, 0, 0],
        [0, 0, 3, x, 0, 0],
        [0, 0, y, 3, 0, 0],
        [0, 0, 0, 0, 3, x],
        [0, 0, 0, 0, y, 3],
    ])
    assert diag(a, b, c) == Matrix([
        [1, 2, 0, 0, 0, 0, 0],
        [2, 3, 0, 0, 0, 0, 0],
        [0, 0, 3, x, 0, 0, 0],
        [0, 0, y, 3, 0, 0, 0],
        [0, 0, 0, 0, 3, x, 3],
        [0, 0, 0, 0, y, 3, z],
        [0, 0, 0, 0, x, y, z],
    ])
    assert diag(a, c, b) == Matrix([
        [1, 2, 0, 0, 0, 0, 0],
        [2, 3, 0, 0, 0, 0, 0],
        [0, 0, 3, x, 3, 0, 0],
        [0, 0, y, 3, z, 0, 0],
        [0, 0, x, y, z, 0, 0],
        [0, 0, 0, 0, 0, 3, x],
        [0, 0, 0, 0, 0, y, 3],
    ])
    a = Matrix([x, y, z])
    b = Matrix([[1, 2], [3, 4]])
    c = Matrix([[5, 6]])
    assert diag(a, 7, b, c) == Matrix([
        [x, 0, 0, 0, 0, 0],
        [y, 0, 0, 0, 0, 0],
        [z, 0, 0, 0, 0, 0],
        [0, 7, 0, 0, 0, 0],
        [0, 0, 1, 2, 0, 0],
        [0, 0, 3, 4, 0, 0],
        [0, 0, 0, 0, 5, 6],
    ])
    assert diag(1, [2, 3], [[4, 5]]) == Matrix([
        [1, 0, 0, 0],
        [0, 2, 0, 0],
        [0, 3, 0, 0],
        [0, 0, 4, 5]])

    pytest.raises(ValueError, lambda: diag(1, 2, 3, spam=123))


def test_get_diag_blocks1():
    a = Matrix([[1, 2], [2, 3]])
    b = Matrix([[3, x], [y, 3]])
    c = Matrix([[3, x, 3], [y, 3, z], [x, y, z]])
    assert a.get_diag_blocks() == [a]
    assert b.get_diag_blocks() == [b]
    assert c.get_diag_blocks() == [c]


def test_get_diag_blocks2():
    a = Matrix([[1, 2], [2, 3]])
    b = Matrix([[3, x], [y, 3]])
    c = Matrix([[3, x, 3], [y, 3, z], [x, y, z]])
    assert diag(a, b, b).get_diag_blocks() == [a, b, b]
    assert diag(a, b, c).get_diag_blocks() == [a, b, c]
    assert diag(a, c, b).get_diag_blocks() == [a, c, b]
    assert diag(c, c, b).get_diag_blocks() == [c, c, b]


def test_inv_block():
    a = Matrix([[1, 2], [2, 3]])
    b = Matrix([[3, x], [y, 3]])
    c = Matrix([[3, x, 3], [y, 3, z], [x, y, z]])
    A = diag(a, b, b)
    assert A.inv(try_block_diag=True) == diag(a.inv(), b.inv(), b.inv())
    A = diag(a, b, c)
    assert A.inv(try_block_diag=True) == diag(a.inv(), b.inv(), c.inv())
    A = diag(a, c, b)
    assert A.inv(try_block_diag=True) == diag(a.inv(), c.inv(), b.inv())
    A = diag(a, a, b, a, c, a)
    assert A.inv(try_block_diag=True) == diag(
        a.inv(), a.inv(), b.inv(), a.inv(), c.inv(), a.inv())
    assert A.inv(try_block_diag=True, method="ADJ") == diag(
        a.inv(method="ADJ"), a.inv(method="ADJ"), b.inv(method="ADJ"),
        a.inv(method="ADJ"), c.inv(method="ADJ"), a.inv(method="ADJ"))


def test_creation_args():
    """
    Check that matrix dimensions can be specified using any reasonable type
    (see issue sympy/sympy#4614).
    """
    pytest.raises(ValueError, lambda: zeros(3, -1))
    pytest.raises(TypeError, lambda: zeros(1, 2, 3, 4))
    assert zeros(int(3)) == zeros(3)
    assert zeros(Integer(3)) == zeros(3)
    assert zeros(3.) == zeros(3)
    assert eye(int(3)) == eye(3)
    assert eye(Integer(3)) == eye(3)
    assert eye(3.) == eye(3)
    assert ones(int(3), Integer(4)) == ones(3, 4)
    pytest.raises(TypeError, lambda: Matrix(5))
    pytest.raises(TypeError, lambda: Matrix(1, 2))


def test_diagonal_symmetrical():
    m = Matrix(2, 2, [0, 1, 1, 0])
    assert not m.is_diagonal()
    assert m.is_symmetric()
    assert m.is_symmetric(simplify=False)

    m = Matrix(2, 2, [1, 0, 0, 1])
    assert m.is_diagonal()

    m = diag(1, 2, 3)
    assert m.is_diagonal()
    assert m.is_symmetric()

    m = Matrix(3, 3, [1, 0, 0, 0, 2, 0, 0, 0, 3])
    assert m == diag(1, 2, 3)

    m = Matrix(2, 3, zeros(2, 3))
    assert not m.is_symmetric()
    assert m.is_diagonal()

    m = Matrix(((5, 0), (0, 6), (0, 0)))
    assert m.is_diagonal()

    m = Matrix(((5, 0, 0), (0, 6, 0)))
    assert m.is_diagonal()

    m = Matrix(3, 3, [1, x**2 + 2*x + 1, y, (x + 1)**2, 2, 0, y, 0, 3])
    assert m.is_symmetric()
    assert not m.is_symmetric(simplify=False)
    assert m.expand().is_symmetric(simplify=False)


def test_diagonalization():
    m = Matrix(3, 2, [-3, 1, -3, 20, 3, 10])
    assert not m.is_diagonalizable()
    assert not m.is_symmetric()
    pytest.raises(NonSquareMatrixError, lambda: m.diagonalize())

    # diagonalizable
    m = diag(1, 2, 3)
    (P, D) = m.diagonalize()
    assert P == eye(3)
    assert D == m

    m = Matrix(2, 2, [0, 1, 1, 0])
    assert m.is_symmetric()
    assert m.is_diagonalizable()
    (P, D) = m.diagonalize()
    assert P.inv() * m * P == D

    m = Matrix(2, 2, [1, 0, 0, 3])
    assert m.is_symmetric()
    assert m.is_diagonalizable()
    (P, D) = m.diagonalize()
    assert P.inv() * m * P == D
    assert P == eye(2)
    assert D == m

    m = Matrix(2, 2, [1, 1, 0, 0])
    assert m.is_diagonalizable()
    (P, D) = m.diagonalize()
    assert P.inv() * m * P == D

    m = Matrix(3, 3, [1, 2, 0, 0, 3, 0, 2, -4, 2])
    assert m.is_diagonalizable()
    (P, D) = m.diagonalize()
    assert P.inv() * m * P == D
    for i in P:
        assert i.as_numer_denom()[1] == 1

    m = Matrix(2, 2, [1, 0, 0, 0])
    assert m.is_diagonal()
    assert m.is_diagonalizable()
    (P, D) = m.diagonalize()
    assert P.inv() * m * P == D
    assert P == Matrix([[0, 1], [1, 0]])

    # diagonalizable, complex only
    m = Matrix(2, 2, [0, 1, -1, 0])
    assert not m.is_diagonalizable(True)
    pytest.raises(MatrixError, lambda: m.diagonalize(True))
    assert m.is_diagonalizable()
    (P, D) = m.diagonalize()
    assert P.inv() * m * P == D

    # not diagonalizable
    m = Matrix(2, 2, [0, 1, 0, 0])
    assert not m.is_diagonalizable()
    pytest.raises(MatrixError, lambda: m.diagonalize())

    m = Matrix(3, 3, [-3, 1, -3, 20, 3, 10, 2, -2, 4])
    assert not m.is_diagonalizable()
    pytest.raises(MatrixError, lambda: m.diagonalize())

    # symbolic
    m = Matrix(2, 2, [a, c, c, b])
    assert m.is_symmetric()
    assert m.is_diagonalizable()


@pytest.mark.xfail
def test_eigen_vects():
    m = Matrix(2, 2, [1, 0, 0, I])
    pytest.raises(NotImplementedError, lambda: m.is_diagonalizable(True))
    # !!! bug because of eigenvects() or roots(x**2 + (-1 - I)*x + I, x)
    # see issue sympy/sympy#5292
    assert not m.is_diagonalizable(True)
    pytest.raises(MatrixError, lambda: m.diagonalize(True))
    (P, D) = m.diagonalize(True)


def test_jordan_form():

    m = Matrix(3, 2, [-3, 1, -3, 20, 3, 10])
    pytest.raises(NonSquareMatrixError, lambda: m.jordan_form())

    # diagonalizable
    m = Matrix(3, 3, [7, -12, 6, 10, -19, 10, 12, -24, 13])
    Jmust = Matrix(3, 3, [-1, 0, 0, 0, 1, 0, 0, 0, 1])
    P, J = m.jordan_form()
    assert Jmust == J
    assert Jmust == m.diagonalize()[1]

    # m = Matrix(3, 3, [0, 6, 3, 1, 3, 1, -2, 2, 1])
    # m.jordan_form()  # very long
    # m.jordan_form()  #

    # diagonalizable, complex only

    # Jordan cells
    # complexity: one of eigenvalues is zero
    m = Matrix(3, 3, [0, 1, 0, -4, 4, 0, -2, 1, 2])
    # The blocks are ordered according to the value of their eigenvalues,
    # in order to make the matrix compatible with .diagonalize()
    Jmust = Matrix(3, 3, [2, 1, 0, 0, 2, 0, 0, 0, 2])
    P, J = m.jordan_form()
    assert Jmust == J
    P, Jcells = m.jordan_cells()
    # same here see 1456ff
    assert Jcells[1] == Matrix(1, 1, [2])
    assert Jcells[0] == Matrix(2, 2, [2, 1, 0, 2])

    # complexity: all of eigenvalues are equal
    m = Matrix(3, 3, [2, 6, -15, 1, 1, -5, 1, 2, -6])
    # Jmust = Matrix(3, 3, [-1, 0, 0, 0, -1, 1, 0, 0, -1])
    # same here see 1456ff
    Jmust = Matrix(3, 3, [-1, 1, 0, 0, -1, 0, 0, 0, -1])
    P, J = m.jordan_form()
    assert Jmust == J

    # complexity: two of eigenvalues are zero
    m = Matrix(3, 3, [4, -5, 2, 5, -7, 3, 6, -9, 4])
    Jmust = Matrix(3, 3, [0, 1, 0, 0, 0, 0, 0, 0, 1])
    P, J = m.jordan_form()
    assert Jmust == J

    m = Matrix(4, 4, [6, 5, -2, -3, -3, -1, 3, 3, 2, 1, -2, -3, -1, 1, 5, 5])
    Jmust = Matrix(4, 4, [2, 1, 0, 0,
                          0, 2, 0, 0,
                          0, 0, 2, 1,
                          0, 0, 0, 2])
    P, J = m.jordan_form()
    assert Jmust == J

    m = Matrix(4, 4, [6, 2, -8, -6, -3, 2, 9, 6, 2, -2, -8, -6, -1, 0, 3, 4])
    # Jmust = Matrix(4, 4, [2, 0, 0, 0, 0, 2, 1, 0, 0, 0, 2, 0, 0, 0, 0, -2])
    # same here see 1456ff
    Jmust = Matrix(4, 4, [-2, 0, 0, 0,
                          0, 2, 1, 0,
                          0, 0, 2, 0,
                          0, 0, 0, 2])
    P, J = m.jordan_form()
    assert Jmust == J

    m = Matrix(4, 4, [5, 4, 2, 1, 0, 1, -1, -1, -1, -1, 3, 0, 1, 1, -1, 2])
    assert not m.is_diagonalizable()
    Jmust = Matrix(4, 4, [1, 0, 0, 0, 0, 2, 0, 0, 0, 0, 4, 1, 0, 0, 0, 4])
    P, J = m.jordan_form()
    assert Jmust == J

    # the following tests are new and include (some) test the cases where the old
    # algorithm failed due to the fact that the block structure can
    # *NOT* be determined  from algebraic and geometric multiplicity alone
    # This can be seen most easily when one lets compute the J.c.f. of a matrix that
    # is in J.c.f already.
    m = Matrix(4, 4, [2, 1, 0, 0,
                      0, 2, 1, 0,
                      0, 0, 2, 0,
                      0, 0, 0, 2])
    P, J = m.jordan_form()
    assert m == J

    m = Matrix(4, 4, [2, 1, 0, 0,
                      0, 2, 0, 0,
                      0, 0, 2, 1,
                      0, 0, 0, 2])
    P, J = m.jordan_form()
    assert m == J


def test_jordan_form_complex_sympyissue_9274():
    A = Matrix([[ 2,  4,  1,  0],
                [-4,  2,  0,  1],
                [ 0,  0,  2,  4],
                [ 0,  0, -4,  2]])
    p = 2 - 4*I
    q = 2 + 4*I
    Jmust1 = Matrix([[p, 1, 0, 0],
                     [0, p, 0, 0],
                     [0, 0, q, 1],
                     [0, 0, 0, q]])
    Jmust2 = Matrix([[q, 1, 0, 0],
                     [0, q, 0, 0],
                     [0, 0, p, 1],
                     [0, 0, 0, p]])
    P, J = A.jordan_form()
    assert J == Jmust1 or J == Jmust2
    assert simplify(P*J*P.inv()) == A


def test_sympyissue_10220():
    # two non-orthogonal Jordan blocks with eigenvalue 1
    M = Matrix([[1, 0, 0, 1],
                [0, 1, 1, 0],
                [0, 0, 1, 1],
                [0, 0, 0, 1]])
    P, C = M.jordan_cells()
    assert P == Matrix([[0, 1, 0, 1],
                        [1, 0, 0, 0],
                        [0, 1, 0, 0],
                        [0, 0, 1, 0]])
    assert len(C) == 2


def test_Matrix_berkowitz_charpoly():
    UA, K_i, K_w = symbols('UA K_i K_w')

    A = Matrix([[-K_i - UA + K_i**2/(K_i + K_w),       K_i*K_w/(K_i + K_w)],
                [           K_i*K_w/(K_i + K_w), -K_w + K_w**2/(K_i + K_w)]])

    charpoly = A.berkowitz_charpoly(x)

    assert charpoly == \
        Poly(x**2 + (K_i*UA + K_w*UA + 2*K_i*K_w)/(K_i + K_w)*x +
             K_i*K_w*UA/(K_i + K_w), x, domain='ZZ(K_i,K_w,UA)')

    assert type(charpoly) is PurePoly

    A = Matrix([[1, 3], [2, 0]])

    assert A.charpoly() == A.charpoly(x) == PurePoly(x**2 - x - 6)


def test_exp():
    m = Matrix([[3, 4], [0, -2]])
    m_exp = Matrix([[exp(3), -4*exp(-2)/5 + 4*exp(3)/5], [0, exp(-2)]])
    assert m.exp() == m_exp
    assert exp(m.as_immutable()) == m_exp

    m = Matrix([[1, 0], [0, 1]])
    assert m.exp() == Matrix([[E, 0], [0, E]])
    assert exp(m.as_immutable()) == Matrix([[E, 0], [0, E]])


def test_has():
    A = Matrix(((x, y), (2, 3)))
    assert A.has(x)
    assert not A.has(z)
    assert A.has(Symbol)

    A = A.subs(x, 2)
    assert not A.has(x)


def test_errors():
    pytest.raises(ValueError, lambda: Matrix([[1, 2], [1]]))
    pytest.raises(TypeError, lambda: Matrix(2, 2, 1))
    pytest.raises(IndexError, lambda: Matrix([[1, 2]])[1.2, 5])
    pytest.raises(IndexError, lambda: Matrix([[1, 2]])[1, 5.2])
    pytest.raises(ValueError, lambda: randMatrix(3, c=4, symmetric=True))
    pytest.raises(ValueError, lambda: Matrix([1, 2]).reshape(4, 6))
    pytest.raises(ShapeError,
                  lambda: Matrix([[1, 2], [3, 4]]).copyin_matrix([1, 0], Matrix([1, 2])))
    pytest.raises(TypeError, lambda: Matrix([[1, 2], [3, 4]]).copyin_list([0,
                                                                           1], set()))
    pytest.raises(NonSquareMatrixError, lambda: Matrix([[1, 2, 3], [2, 3, 0]]).inv())
    pytest.raises(ShapeError,
                  lambda: Matrix(1, 2, [1, 2]).row_join(Matrix([[1, 2], [3, 4]])))
    pytest.raises(
        ShapeError, lambda: Matrix([1, 2]).col_join(Matrix([[1, 2], [3, 4]])))
    pytest.raises(ShapeError, lambda: Matrix([1]).row_insert(1, Matrix([[1,
                                                                         2], [3, 4]])))
    pytest.raises(ShapeError, lambda: Matrix([1]).col_insert(1, Matrix([[1,
                                                                         2], [3, 4]])))
    pytest.raises(NonSquareMatrixError, lambda: Matrix([1, 2]).trace())
    pytest.raises(TypeError, lambda: Matrix([1]).applyfunc(1))
    pytest.raises(ShapeError, lambda: Matrix([1]).LUsolve(Matrix([[1, 2], [3, 4]])))
    pytest.raises(MatrixError, lambda: Matrix([[1, 2, 3],
                                               [4, 5, 6],
                                               [7, 8, 9]]).QRdecomposition())
    pytest.raises(MatrixError, lambda: Matrix(1, 2, [1, 2]).QRdecomposition())
    pytest.raises(
        NonSquareMatrixError, lambda: Matrix([1, 2]).LUdecomposition_Simple())
    pytest.raises(ValueError, lambda: Matrix([[1, 2], [3, 4]]).minorEntry(4, 5))
    pytest.raises(ValueError, lambda: Matrix([[1, 2], [3, 4]]).minorMatrix(4, 5))
    pytest.raises(TypeError, lambda: Matrix([1, 2, 3]).cross(1))
    pytest.raises(TypeError, lambda: Matrix([1, 2, 3]).dot(1))
    pytest.raises(ShapeError, lambda: Matrix([1, 2, 3]).dot(Matrix([1, 2])))
    pytest.raises(ShapeError, lambda: Matrix([1, 2]).dot([]))
    pytest.raises(TypeError, lambda: Matrix([1, 2]).dot('a'))
    pytest.raises(NonSquareMatrixError, lambda: Matrix([1, 2, 3]).exp())
    pytest.raises(ShapeError, lambda: Matrix([[1, 2], [3, 4]]).normalized())
    pytest.raises(ValueError, lambda: Matrix([1, 2]).inv(method='not a method'))
    pytest.raises(NonSquareMatrixError, lambda: Matrix([1, 2]).inverse_GE())
    pytest.raises(ValueError, lambda: Matrix([[1, 2], [1, 2]]).inverse_GE())
    pytest.raises(NonSquareMatrixError, lambda: Matrix([1, 2]).inverse_ADJ())
    pytest.raises(ValueError, lambda: Matrix([[1, 2], [1, 2]]).inverse_ADJ())
    pytest.raises(NonSquareMatrixError, lambda: Matrix([1, 2]).inverse_LU())
    pytest.raises(NonSquareMatrixError, lambda: Matrix([1, 2]).is_nilpotent())
    pytest.raises(NonSquareMatrixError, lambda: Matrix([1, 2]).det())
    pytest.raises(ValueError,
                  lambda: Matrix([[1, 2], [3, 4]]).det(method='Not a real method'))
    pytest.raises(NonSquareMatrixError, lambda: Matrix([1, 2]).det_bareis())
    pytest.raises(NonSquareMatrixError, lambda: Matrix([1, 2]).berkowitz())
    pytest.raises(NonSquareMatrixError, lambda: Matrix([1, 2]).berkowitz_det())
    pytest.raises(ValueError,
                  lambda: hessian(Matrix([[1, 2], [3, 4]]), Matrix([[1, 2], [2, 1]])))
    pytest.raises(ValueError, lambda: hessian(Matrix([[1, 2], [3, 4]]), []))
    pytest.raises(ValueError, lambda: hessian(x**2, 'a'))
    pytest.raises(ValueError,
                  lambda: Matrix([[5, 10, 7],
                                  [0, -1, 2],
                                  [8, 3, 4]]).LUdecomposition_Simple(iszerofunc=lambda x: abs(x) <= 4))
    pytest.raises(TypeError,
                  lambda: Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])**(0.5))
    pytest.raises(IndexError, lambda: eye(3)[5, 2])
    pytest.raises(IndexError, lambda: eye(3)[2, 5])
    M = Matrix(((1, 2, 3, 4), (5, 6, 7, 8), (9, 10, 11, 12), (13, 14, 15, 16)))
    pytest.raises(ValueError, lambda: M.det('method=LU_decomposition()'))
    pytest.raises(NonSquareMatrixError, lambda: ones(2, 3).det_LU_decomposition())
    with pytest.raises(TypeError):
        M[(1,)] = 1
    pytest.raises(AttributeError, lambda: zeros(0).jordan_cells())


def test_len():
    assert len(Matrix()) == 0
    assert len(Matrix([[1, 2]])) == len(Matrix([[1], [2]])) == 2
    assert len(Matrix(0, 2, lambda i, j: 0)) == \
        len(Matrix(2, 0, lambda i, j: 0)) == 0
    assert len(Matrix([[0, 1, 2], [3, 4, 5]])) == 6
    assert Matrix([1]) == Matrix([[1]])
    assert not Matrix()
    assert Matrix() == Matrix([])


def test_integrate():
    A = Matrix(((1, 4, x), (y, 2, 4), (10, 5, x**2)))
    assert A.integrate(x) == \
        Matrix(((x, 4*x, x**2/2), (x*y, 2*x, 4*x), (10*x, 5*x, x**3/3)))
    assert A.integrate(y) == \
        Matrix(((y, 4*y, x*y), (y**2/2, 2*y, 4*y), (10*y, 5*y, y*x**2)))


def test_limit():
    A = Matrix(((1, 4, sin(x)/x), (y, 2, 4), (10, 5, x**2 + 1)))
    assert A.limit(x, 0) == Matrix(((1, 4, 1), (y, 2, 4), (10, 5, 1)))


def test_diff():
    A = Matrix(((1, 4, x), (y, 2, 4), (10, 5, x**2 + 1)))
    assert A.diff(x) == Matrix(((0, 0, 1), (0, 0, 0), (0, 0, 2*x)))
    assert A.diff(y) == Matrix(((0, 0, 0), (1, 0, 0), (0, 0, 0)))


def test_getattr():
    A = Matrix(((1, 4, x), (y, 2, 4), (10, 5, x**2 + 1)))
    pytest.raises(AttributeError, lambda: A.nonexistantattribute)
    assert getattr(A, 'diff')(x) == Matrix(((0, 0, 1), (0, 0, 0), (0, 0, 2*x)))


def test_hessenberg():
    A = Matrix([[3, 4, 1], [2, 4, 5], [0, 1, 2]])
    assert A.is_upper_hessenberg
    A = A.T
    assert A.is_lower_hessenberg
    A[0, -1] = 1
    assert A.is_lower_hessenberg is False

    A = Matrix([[3, 4, 1], [2, 4, 5], [3, 1, 2]])
    assert not A.is_upper_hessenberg


def test_cholesky():
    pytest.raises(NonSquareMatrixError, lambda: Matrix((1, 2)).cholesky())
    pytest.raises(ValueError, lambda: Matrix(((1, 2), (3, 4))).cholesky())
    A = Matrix(((25, 15, -5), (15, 18, 0), (-5, 0, 11)))
    assert A.cholesky() * A.cholesky().T == A
    assert A.cholesky().is_lower
    assert A.cholesky() == Matrix([[5, 0, 0], [3, 3, 0], [-1, 1, 3]])


def test_LDLdecomposition():
    pytest.raises(NonSquareMatrixError, lambda: Matrix((1, 2)).LDLdecomposition())
    pytest.raises(ValueError, lambda: Matrix(((1, 2), (3, 4))).LDLdecomposition())
    A = Matrix(((25, 15, -5), (15, 18, 0), (-5, 0, 11)))
    L, D = A.LDLdecomposition()
    assert L * D * L.T == A
    assert L.is_lower
    assert L == Matrix([[1, 0, 0], [ Rational(3, 5), 1, 0], [Rational(-1, 5), Rational(1, 3), 1]])
    assert D.is_diagonal()
    assert D == Matrix([[25, 0, 0], [0, 9, 0], [0, 0, 9]])


def test_cholesky_solve():
    A = Matrix([[2, 3, 5],
                [3, 6, 2],
                [8, 3, 6]])
    x = Matrix(3, 1, [3, 7, 5])
    b = A*x
    soln = A.cholesky_solve(b)
    assert soln == x
    A = Matrix([[0, -1, 2],
                [5, 10, 7],
                [8,  3, 4]])
    x = Matrix(3, 1, [-1, 2, 5])
    b = A*x
    soln = A.cholesky_solve(b)
    assert soln == x

    A = Matrix([[-12, -3, -8],
                [ -3, -4, -6],
                [ -8, -6, -6]])
    x = Matrix(3, 1, [Rational(-6, 83), Rational(15, 83), Rational(-7, 83)])
    assert A.cholesky_solve(A*x) == x


def test_LDLsolve():
    A = Matrix([[2, 3, 5],
                [3, 6, 2],
                [8, 3, 6]])
    x = Matrix(3, 1, [3, 7, 5])
    b = A*x
    soln = A.LDLsolve(b)
    assert soln == x
    A = Matrix([[0, -1, 2],
                [5, 10, 7],
                [8,  3, 4]])
    x = Matrix(3, 1, [-1, 2, 5])
    b = A*x
    soln = A.LDLsolve(b)
    assert soln == x


def test_lower_triangular_solve():

    pytest.raises(NonSquareMatrixError,
                  lambda: Matrix([1, 0]).lower_triangular_solve(Matrix([0, 1])))
    pytest.raises(ShapeError,
                  lambda: Matrix([[1, 0], [0, 1]]).lower_triangular_solve(Matrix([1])))
    pytest.raises(ValueError,
                  lambda: Matrix([[2, 1], [1, 2]]).lower_triangular_solve(
                      Matrix([[1, 0], [0, 1]])))

    A = Matrix([[1, 0], [0, 1]])
    B = Matrix([[x, y], [y, x]])
    C = Matrix([[4, 8], [2, 9]])
    S = Matrix([[0, 0], [0, 1]])

    assert A.lower_triangular_solve(B) == B
    assert A.lower_triangular_solve(C) == C
    pytest.raises(ValueError, lambda: S.lower_triangular_solve(B))


def test_upper_triangular_solve():

    pytest.raises(NonSquareMatrixError,
                  lambda: Matrix([1, 0]).upper_triangular_solve(Matrix([0, 1])))
    pytest.raises(TypeError,
                  lambda: Matrix([[1, 0], [0, 1]]).upper_triangular_solve(Matrix([1])))
    pytest.raises(TypeError,
                  lambda: Matrix([[2, 1], [1, 2]]).upper_triangular_solve(
                      Matrix([[1, 0], [0, 1]])))

    A = Matrix([[1, 0], [0, 1]])
    B = Matrix([[x, y], [y, x]])
    C = Matrix([[2, 4], [3, 8]])
    S = Matrix([[0, 0], [0, 1]])

    assert A.upper_triangular_solve(B) == B
    assert A.upper_triangular_solve(C) == C
    pytest.raises(ValueError, lambda: S.upper_triangular_solve(B))


def test_diagonal_solve():
    pytest.raises(TypeError, lambda: ones(2).diagonal_solve(Matrix([1])))
    pytest.raises(TypeError, lambda: eye(2).diagonal_solve(Matrix([1])))
    A = Matrix([[1, 0], [0, 1]])*2
    B = Matrix([[x, y], [y, x]])
    assert A.diagonal_solve(B) == B/2


def test_matrix_norm():
    # Vector Tests
    # Test columns and symbols
    x = Symbol('x', real=True)
    v = Matrix([cos(x), sin(x)])
    assert trigsimp(v.norm(2)) == 1
    assert v.norm(10) == Pow(cos(x)**10 + sin(x)**10, Rational(1, 10))

    # Test Rows
    A = Matrix([[5, Rational(3, 2)]])
    assert A.norm() == Pow(25 + Rational(9, 4), Rational(1, 2))
    assert A.norm(oo) == max(A._mat)
    assert A.norm(-oo) == min(A._mat)

    # Matrix Tests
    # Intuitive test
    A = Matrix([[1, 1], [1, 1]])
    assert A.norm(2) == 2
    assert A.norm(-2) == 0
    assert A.norm('frobenius') == 2
    assert eye(10).norm(2) == eye(10).norm(-2) == 1

    # Test with Symbols and more complex entries
    A = Matrix([[3, y, y], [x, Rational(1, 2), -pi]])
    assert (A.norm('fro')
            == sqrt(Rational(37, 4) + 2*abs(y)**2 + pi**2 + x**2))

    # Check non-square
    A = Matrix([[1, 2, -3], [4, 5, Rational(13, 2)]])
    assert A.norm(2) == sqrt(Rational(389, 8) + sqrt(78665)/8)
    assert A.norm(-2) == Integer(0)
    assert A.norm('frobenius') == sqrt(389)/2

    # Test properties of matrix norms
    # https//en.wikipedia.org/wiki/Matrix_norm#Definition
    # Two matrices
    A = Matrix([[1, 2], [3, 4]])
    B = Matrix([[5, 5], [-2, 2]])
    C = Matrix([[0, -I], [I, 0]])
    D = Matrix([[1, 0], [0, -1]])
    L = [A, B, C, D]
    alpha = Symbol('alpha', extended_real=True)

    for order in ['fro', 2, -2]:
        # Zero Check
        assert zeros(3).norm(order) == Integer(0)
        # Check Triangle Inequality for all Pairs of Matrices
        for X in L:
            for Y in L:
                assert simplify(X.norm(order) + Y.norm(order) >=
                                (X + Y).norm(order))
        # Scalar multiplication linearity
        for M in [A, B, C, D]:
            if order in [2, -2]:
                # Abs is causing tests to fail when Abs(alpha) is inside a Max
                # or Min. The tests produce mathematically true statements that
                # are too complex to be simplified well.
                continue
            try:
                assert ((alpha*M).norm(order) ==
                        abs(alpha) * M.norm(order))
            except NotImplementedError:
                pass  # Some Norms fail on symbolic matrices due to Max issue

    # Test Properties of Vector Norms
    # https//en.wikipedia.org/wiki/Vector_norm
    # Two column vectors
    a = Matrix([1, 1 - 1*I, -3])
    b = Matrix([Rational(1, 2), 1*I, 1])
    c = Matrix([-1, -1, -1])
    d = Matrix([3, 2, I])
    e = Matrix([Integer(1e2), Rational(1, 1e2), 1])
    L = [a, b, c, d, e]
    alpha = Symbol('alpha', extended_real=True)

    for order in [1, 2, -1, -2, oo, -oo, pi]:
        # Zero Check
        if order > 0:
            assert Matrix([0, 0, 0]).norm(order) == Integer(0)
        # Triangle inequality on all pairs
        if order >= 1:  # Triangle InEq holds only for these norms
            for v in L:
                for w in L:
                    assert simplify(v.norm(order) + w.norm(order) >=
                                    (v + w).norm(order))
        # Linear to scalar multiplication
        if order in [1, 2, -1, -2, oo, -oo]:
            for vec in L:
                try:
                    assert simplify((alpha*v).norm(order) -
                                    (abs(alpha) * v.norm(order))) == 0
                except NotImplementedError:
                    pass  # Some Norms fail on symbolics due to Max issue


def test_singular_values():
    x = Symbol('x', extended_real=True)

    A = Matrix([[0, 1*I], [2, 0]])
    assert A.singular_values() == [2, 1]

    A = eye(3)
    A[1, 1] = x
    A[2, 2] = 5
    vals = A.singular_values()
    assert 1 in vals and 5 in vals and abs(x) in vals

    A = Matrix([[sin(x), cos(x)], [-cos(x), sin(x)]])
    vals = [sv.trigsimp() for sv in A.singular_values()]
    assert vals == [Integer(1), Integer(1)]


def test_condition_number():
    x = Symbol('x', extended_real=True)
    A = eye(3)
    A[0, 0] = 10
    A[2, 2] = Rational(1, 10)
    assert A.condition_number() == 100

    A[1, 1] = x
    assert A.condition_number() == Max(10, Abs(x)) / Min(Rational(1, 10), Abs(x))

    M = Matrix([[cos(x), sin(x)], [-sin(x), cos(x)]])
    Mc = M.condition_number()
    assert all(Float(1.).epsilon_eq(Mc.subs(x, val).evalf()) for val in
               [Rational(1, 5), Rational(1, 2), Rational(1, 10), pi/2, pi, 7*pi/4 ])

    assert Matrix([]).condition_number() == 0  # issue sympy/sympy#10782


def test_equality():
    A = Matrix(((1, 2, 3), (4, 5, 6), (7, 8, 9)))
    B = Matrix(((9, 8, 7), (6, 5, 4), (3, 2, 1)))
    assert A == A[:, :]
    assert not A != A[:, :]
    assert not A == B
    assert A != B
    assert A != 10
    assert not A == 10

    # A SparseMatrix can be equal to a Matrix
    C = SparseMatrix(((1, 0, 0), (0, 1, 0), (0, 0, 1)))
    D = Matrix(((1, 0, 0), (0, 1, 0), (0, 0, 1)))
    assert C == D
    assert not C != D


def test_col_join():
    assert eye(3).col_join(Matrix([[7, 7, 7]])) == \
        Matrix([[1, 0, 0],
                [0, 1, 0],
                [0, 0, 1],
                [7, 7, 7]])


def test_row_insert():
    r4 = Matrix([[4, 4, 4]])
    for i in range(-4, 5):
        l = [1, 0, 0]
        l.insert(i, 4)
        assert flatten(eye(3).row_insert(i, r4).col(0).tolist()) == l


def test_col_insert():
    c4 = Matrix([4, 4, 4])
    for i in range(-4, 5):
        l = [0, 0, 0]
        l.insert(i, 4)
        assert flatten(zeros(3).col_insert(i, c4).row(0).tolist()) == l


def test_normalized():
    assert Matrix([3, 4]).normalized() == \
        Matrix([Rational(3, 5), Rational(4, 5)])


def test_print_nonzero():
    assert capture(lambda: eye(3).print_nonzero()) == \
        '[X  ]\n[ X ]\n[  X]\n'
    assert capture(lambda: eye(3).print_nonzero('.')) == \
        '[.  ]\n[ . ]\n[  .]\n'


def test_zeros_eye():
    assert Matrix.eye(3) == eye(3)
    assert Matrix.zeros(3) == zeros(3)
    assert ones(3, 4) == Matrix(3, 4, [1]*12)

    i = Matrix([[1, 0], [0, 1]])
    z = Matrix([[0, 0], [0, 0]])
    for cls in classes:
        m = cls.eye(2)
        assert i == m  # but m == i will fail if m is immutable
        assert i == eye(2, cls=cls)
        assert type(m) == cls
        m = cls.zeros(2)
        assert z == m
        assert z == zeros(2, cls=cls)
        assert type(m) == cls


def test_vandermonde():
    assert vandermonde(1) == Matrix([1])
    assert vandermonde(2, [x, y]) == Matrix([[1, x], [1, y]])
    assert vandermonde(3) == Matrix([[1, Symbol('C0'), Symbol('C0')**2],
                                     [1, Symbol('C1'), Symbol('C1')**2],
                                     [1, Symbol('C2'), Symbol('C2')**2]])


def test_is_zero():
    assert Matrix().is_zero
    assert Matrix([[0, 0], [0, 0]]).is_zero
    assert zeros(3, 4).is_zero
    assert not eye(3).is_zero
    assert Matrix([[x, 0], [0, 0]]).is_zero is None
    assert SparseMatrix([[x, 0], [0, 0]]).is_zero is None
    assert ImmutableMatrix([[x, 0], [0, 0]]).is_zero is None
    assert ImmutableSparseMatrix([[x, 0], [0, 0]]).is_zero is None
    assert Matrix([[x, 1], [0, 0]]).is_zero is False
    a = Symbol('a', nonzero=True)
    assert Matrix([[a, 0], [0, 0]]).is_zero is False


def test_rotation_matrices():
    # This tests the rotation matrices by rotating about an axis and back.
    theta = pi/3
    r3_plus = rot_axis3(theta)
    r3_minus = rot_axis3(-theta)
    r2_plus = rot_axis2(theta)
    r2_minus = rot_axis2(-theta)
    r1_plus = rot_axis1(theta)
    r1_minus = rot_axis1(-theta)
    assert r3_minus*r3_plus*eye(3) == eye(3)
    assert r2_minus*r2_plus*eye(3) == eye(3)
    assert r1_minus*r1_plus*eye(3) == eye(3)

    # Check the correctness of the trace of the rotation matrix
    assert r1_plus.trace() == 1 + 2*cos(theta)
    assert r2_plus.trace() == 1 + 2*cos(theta)
    assert r3_plus.trace() == 1 + 2*cos(theta)

    # Check that a rotation with zero angle doesn't change anything.
    assert rot_axis1(0) == eye(3)
    assert rot_axis2(0) == eye(3)
    assert rot_axis3(0) == eye(3)


def test_DeferredVector():
    assert str(DeferredVector("vector")[4]) == "vector[4]"
    assert repr(DeferredVector("vector")) == "DeferredVector('vector')"
    assert sympify(DeferredVector("d")) == DeferredVector("d")
    pytest.raises(IndexError, lambda: DeferredVector('X')[-1])


def test_DeferredVector_not_iterable():
    assert not iterable(DeferredVector('X'))


def test_DeferredVector_Matrix():
    pytest.raises(TypeError, lambda: Matrix(DeferredVector("V")))


def test_GramSchmidt():
    R = Rational
    m1 = Matrix(1, 2, [1, 2])
    m2 = Matrix(1, 2, [2, 3])
    assert GramSchmidt([m1, m2]) == \
        [Matrix(1, 2, [1, 2]), Matrix(1, 2, [R(2)/5, R(-1)/5])]
    assert GramSchmidt([m1.T, m2.T]) == \
        [Matrix(2, 1, [1, 2]), Matrix(2, 1, [R(2)/5, R(-1)/5])]
    # from wikipedia
    assert GramSchmidt([Matrix([3, 1]), Matrix([2, 2])], True) == [
        Matrix([3*sqrt(10)/10, sqrt(10)/10]),
        Matrix([-sqrt(10)/10, 3*sqrt(10)/10])]

    pytest.raises(ValueError, lambda: GramSchmidt([Matrix([1, 2]),
                                                   Matrix([2, 4])]))


def test_casoratian():
    assert casoratian([1, 2, 3, 4], 1) == 0
    assert casoratian([1, 2, 3, 4], 1, zero=False) == 0


def test_zero_dimension_multiply():
    assert (Matrix()*zeros(0, 3)).shape == (0, 3)
    assert zeros(3, 0)*zeros(0, 3) == zeros(3, 3)
    assert zeros(0, 3)*zeros(3, 0) == Matrix()


def test_slice_sympyissue_5983():
    m = Matrix(2, 2, range(4))
    assert m[1, :] == Matrix([[2, 3]])
    assert m[-1, :] == Matrix([[2, 3]])
    assert m[:, 1] == Matrix([[1, 3]]).T
    assert m[:, -1] == Matrix([[1, 3]]).T
    pytest.raises(IndexError, lambda: m[2, :])
    pytest.raises(IndexError, lambda: m[2, 2])


def test_slice_sympyissue_6500():
    assert zeros(0, 3)[:, -1].shape == (0, 1)
    assert zeros(3, 0)[0, :] == Matrix(1, 0, [])


def test_copyin():
    s = zeros(3, 3)
    s[3] = 1
    assert s[:, 0] == Matrix([0, 1, 0])
    assert s[3] == 1
    assert s[3: 4] == [1]
    s[1, 1] = 42
    assert s[1, 1] == 42
    assert s[1, 1:] == Matrix([[42, 0]])
    s[1, 1:] = Matrix([[5, 6]])
    assert s[1, :] == Matrix([[1, 5, 6]])
    s[1, 1:] = [[42, 43]]
    assert s[1, :] == Matrix([[1, 42, 43]])
    s[0, 0] = 17
    assert s[:, :1] == Matrix([17, 1, 0])
    s[0, 0] = [1, 1, 1]
    assert s[:, 0] == Matrix([1, 1, 1])
    s[0, 0] = Matrix([1, 1, 1])
    assert s[:, 0] == Matrix([1, 1, 1])
    s[0, 0] = SparseMatrix([1, 1, 1])
    assert s[:, 0] == Matrix([1, 1, 1])


def test_invertible_check():
    assert Matrix([[1, 2], [1, 2]]).rref() == (Matrix([[1, 2], [0, 0]]), [0])
    pytest.raises(ValueError, lambda: Matrix([[1, 2], [1, 2]]).inv())
    m = Matrix([
        [-1, -1,  0],
        [ x,  1,  1],
        [ 1,  x, -1],
    ])
    assert len(m.rref()[1]) == 2
    assert m.rref()[0] != eye(3)
    pytest.raises(ValueError, lambda: m.inv(method="ADJ"))
    pytest.raises(ValueError, lambda: m.inv(method="GE"))
    pytest.raises(ValueError, lambda: m.inv(method="LU"))


@pytest.mark.xfail
def test_sympyissue_3959():
    e = x*y
    assert e.subs(x, Matrix([3, 5, 3])) == Matrix([3, 5, 3])*y


def test_sympyissue_5964():
    assert str(Matrix([[1, 2], [3, 4]])) == 'Matrix([\n[1, 2],\n[3, 4]])'


def test_sympyissue_7604():
    assert sstr(Matrix([[x, 2*y], [y**2, x + 3]])) == \
        'Matrix([\n[   x,   2*y],\n[y**2, x + 3]])'


def test_is_Identity():
    assert eye(3).is_Identity
    assert eye(3).as_immutable().is_Identity
    assert not zeros(3).is_Identity
    assert not ones(3).is_Identity
    # issue sympy/sympy#6242
    assert not Matrix([[1, 0, 0]]).is_Identity
    # issue sympy/sympy#8854
    assert SparseMatrix(3, 3, {(0, 0): 1, (1, 1): 1, (2, 2): 1}).is_Identity
    assert not SparseMatrix(2, 3, range(6)).is_Identity
    assert not SparseMatrix(3, 3, {(0, 0): 1, (1, 1): 1}).is_Identity
    assert not SparseMatrix(3, 3, {(0, 0): 1, (1, 1): 1, (2, 2): 1, (0, 1): 2, (0, 2): 3}).is_Identity


def test_dot():
    assert ones(1, 3).dot(ones(3, 1)) == 3
    assert ones(1, 3).dot([1, 1, 1]) == 3
    assert ones(2, 3).dot([1, 1]) == [2, 2, 2]
    pytest.raises(ShapeError, lambda: ones(1, 3).dot([1, 1]))
    pytest.raises(ShapeError, lambda: ones(1, 3).dot(ones(2, 1)))


def test_dual():
    B_x, B_y, B_z, E_x, E_y, E_z = symbols(
        'B_x B_y B_z E_x E_y E_z', extended_real=True)
    F = Matrix((
        (   0,  E_x,  E_y,  E_z),
        (-E_x,    0,  B_z, -B_y),
        (-E_y, -B_z,    0,  B_x),
        (-E_z,  B_y, -B_x,    0)
    ))
    Fd = Matrix((
        (  0, -B_x, -B_y, -B_z),
        (B_x,    0,  E_z, -E_y),
        (B_y, -E_z,    0,  E_x),
        (B_z,  E_y, -E_x,    0)
    ))
    assert F.dual().equals(Fd)
    assert eye(3).dual().equals(zeros(3))
    assert F.dual().dual().equals(-F)


def test_anti_symmetric():
    assert Matrix([1, 2]).is_anti_symmetric() is False
    m = Matrix(3, 3, [0, x**2 + 2*x + 1, y, -(x + 1)**2, 0, x*y, -y, -x*y, 0])
    assert m.is_anti_symmetric() is True
    assert m.is_anti_symmetric(simplify=False) is False
    assert m.is_anti_symmetric(simplify=lambda x: x) is False

    # tweak to fail
    m[2, 1] = -m[2, 1]
    assert m.is_anti_symmetric() is False
    # untweak
    m[2, 1] = -m[2, 1]

    m = m.expand()
    assert m.is_anti_symmetric(simplify=False) is True
    m[0, 0] = 1
    assert m.is_anti_symmetric() is False


def test_normalize_sort_diogonalization():
    A = Matrix(((1, 2), (2, 1)))
    P, Q = A.diagonalize(normalize=True)
    assert P*P.T == P.T*P == eye(P.cols)
    P, Q = A.diagonalize(normalize=True, sort=True)
    assert P*P.T == P.T*P == eye(P.cols)
    assert P*Q*P.inv() == A


def test_sympyissue_5321():
    pytest.raises(ValueError, lambda: Matrix([[1, 2, 3], Matrix(0, 1, [])]))


def test_sympyissue_5320():
    assert Matrix.hstack(eye(2), 2*eye(2)) == Matrix([
        [1, 0, 2, 0],
        [0, 1, 0, 2]
    ])
    assert Matrix.vstack(eye(2), 2*eye(2)) == Matrix([
        [1, 0],
        [0, 1],
        [2, 0],
        [0, 2]
    ])
    cls = SparseMatrix
    assert cls.hstack(cls(eye(2)), cls(2*eye(2))) == Matrix([
        [1, 0, 2, 0],
        [0, 1, 0, 2]
    ])


def test_cross():
    a = [1, 2, 3]
    b = [3, 4, 5]
    col = Matrix([-2, 4, -2])
    row = col.T

    def test(M, ans):
        assert ans == M
        assert type(M) == cls
    for cls in classes:
        A = cls(a)
        B = cls(b)
        test(A.cross(B), col)
        test(A.cross(B.T), col)
        test(A.T.cross(B.T), row)
        test(A.T.cross(B), row)
    pytest.raises(ShapeError, lambda:
                  Matrix(1, 2, [1, 1]).cross(Matrix(1, 2, [1, 1])))


def test_hash():
    for cls in classes[-2:]:
        s = {cls.eye(1), cls.eye(1)}
        assert len(s) == 1 and s.pop() == cls.eye(1)
    # issue sympy/sympy#3979
    for cls in classes[:2]:
        assert not isinstance(cls.eye(1), collections.Hashable)


@pytest.mark.xfail
def test_sympyissue_3979():
    # when this passes, delete this and change the [1:2]
    # to [:2] in the test_hash above for issue sympy/sympy#3979
    cls = classes[0]
    pytest.raises(AttributeError, lambda: hash(cls.eye(1)))


def test_adjoint():
    dat = [[0, I], [1, 0]]
    ans = Matrix([[0, 1], [-I, 0]])
    for cls in classes:
        assert ans == cls(dat).adjoint()


def test_simplify_immutable():
    assert simplify(ImmutableMatrix([[sin(x)**2 + cos(x)**2]])) == \
        ImmutableMatrix([[1]])


def test_rank():
    m = Matrix([[1, 2], [x, 1 - 1/x]])
    assert m.rank() == 2
    n = Matrix(3, 3, range(1, 10))
    assert n.rank() == 2
    p = zeros(3)
    assert p.rank() == 0


def test_replace():
    F, G = symbols('F, G', cls=Function)
    K = Matrix(2, 2, lambda i, j: G(i+j))
    M = Matrix(2, 2, lambda i, j: F(i+j))
    N = M.replace(F, G)
    assert N == K

    K = Matrix(2, 2, [G(0), G(1), G(1), G(2)])
    M = Matrix(2, 2, lambda i, j: F(i+j))
    N = M.replace(F, G)
    assert N == K


def test_atoms():
    m = Matrix([[1, 2], [x, 1 - 1/x]])
    assert m.atoms() == {Integer(1), Integer(2), Integer(-1), x}
    assert m.atoms(Symbol) == {x}


def test_pinv():
    # Pseudoinverse of an invertible matrix is the inverse.
    A1 = Matrix([[a, b], [c, d]])
    assert simplify(A1.pinv()) == simplify(A1.inv())
    # Test the four properties of the pseudoinverse for various matrices.
    As = [Matrix([[13, 104], [2212, 3], [-3, 5]]),
          Matrix([[1, 7, 9], [11, 17, 19]]),
          Matrix([a, b])]
    for A in As:
        A_pinv = A.pinv()
        AAp = A * A_pinv
        ApA = A_pinv * A
        assert simplify(AAp * A) == A
        assert simplify(ApA * A_pinv) == A_pinv
        assert AAp.H == AAp
        assert ApA.H == ApA

    assert zeros(3).pinv() == zeros(3)


def test_pinv_solve():
    # Fully determined system (unique result, identical to other solvers).
    A = Matrix([[1, 5], [7, 9]])
    B = Matrix([12, 13])
    assert A.pinv_solve(B) == A.cholesky_solve(B)
    assert A.pinv_solve(B) == A.LDLsolve(B)
    assert A.pinv_solve(B) == Matrix([Rational(-43, 26), Rational(71, 26)])
    assert A * A.pinv() * B == B
    # Fully determined, with two-dimensional B matrix.
    B = Matrix([[12, 13, 14], [15, 16, 17]])
    assert A.pinv_solve(B) == A.cholesky_solve(B)
    assert A.pinv_solve(B) == A.LDLsolve(B)
    assert A.pinv_solve(B) == Matrix([[-33, -37, -41], [69, 75, 81]]) / 26
    assert A * A.pinv() * B == B
    # Underdetermined system (infinite results).
    A = Matrix([[1, 0, 1], [0, 1, 1]])
    B = Matrix([5, 7])
    solution = A.pinv_solve(B, Matrix([x, y, z]))
    assert solution == Matrix([[x/3 + y/3 - z/3 + 1],
                               [x/3 + y/3 - z/3 + 3],
                               [-x/3 - y/3 + z/3 + 4]])
    assert A * A.pinv() * B == B
    # Overdetermined system (least squares results).
    A = Matrix([[1, 0], [0, 0], [0, 1]])
    B = Matrix([3, 2, 1])
    assert A.pinv_solve(B) == Matrix([3, 1])
    # Proof the solution is not exact.
    assert A * A.pinv() * B != B


@pytest.mark.xfail
def test_pinv_rank_deficient():
    # Test the four properties of the pseudoinverse for various matrices.
    As = [Matrix([[1, 1, 1], [2, 2, 2]]),
          Matrix([[1, 0], [0, 0]])]
    for A in As:
        A_pinv = A.pinv()
        AAp = A * A_pinv
        ApA = A_pinv * A
        assert simplify(AAp * A) == A
        assert simplify(ApA * A_pinv) == A_pinv
        assert AAp.H == AAp
        assert ApA.H == ApA
    # Test solving with rank-deficient matrices.
    A = Matrix([[1, 0], [0, 0]])
    # Exact, non-unique solution.
    B = Matrix([3, 0])
    solution = A.pinv_solve(B)
    w1 = solution.atoms(Symbol).pop()
    assert w1.name == 'w1_0'
    assert solution == Matrix([3, w1])
    assert A * A.pinv() * B == B
    # Least squares, non-unique solution.
    B = Matrix([3, 1])
    solution = A.pinv_solve(B)
    w1 = solution.atoms(Symbol).pop()
    assert w1.name == 'w1_0'
    assert solution == Matrix([3, w1])
    assert A * A.pinv() * B != B


def test_sympyissue_7201():
    assert ones(0, 1) + ones(0, 1) == Matrix(0, 1, [])
    assert ones(1, 0) + ones(1, 0) == Matrix(1, 0, [])


def test_free_symbols():
    for M in ImmutableMatrix, ImmutableSparseMatrix, Matrix, SparseMatrix:
        assert M([[x], [0]]).free_symbols == {x}


@pytest.mark.skipif(numpy is None, reason="no numpy")
def test_from_ndarray():
    """See issue sympy/sympy#7465."""
    assert Matrix(numpy.array([1, 2, 3])) == Matrix([1, 2, 3])
    assert Matrix(numpy.array([[1, 2, 3]])) == Matrix([[1, 2, 3]])
    assert Matrix(numpy.array([[1, 2, 3], [4, 5, 6]])) == \
        Matrix([[1, 2, 3], [4, 5, 6]])
    assert Matrix(numpy.array([x, y, z])) == Matrix([x, y, z])
    pytest.raises(NotImplementedError, lambda: Matrix(numpy.array([[
        [1, 2], [3, 4]], [[5, 6], [7, 8]]])))


def test_hermitian():
    a = Matrix([[1, I], [-I, 1]])
    assert a.is_hermitian
    a[0, 0] = 2*I
    assert a.is_hermitian is False
    a[0, 0] = x
    assert a.is_hermitian is None
    a[0, 1] = a[1, 0]*I
    assert a.is_hermitian is False


def test_sympyissue_9457_9467_9876():
    # for row_del(index)
    M = Matrix([[1, 2, 3], [2, 3, 4], [3, 4, 5]])
    M.row_del(1)
    assert M == Matrix([[1, 2, 3], [3, 4, 5]])
    N = Matrix([[1, 2, 3], [2, 3, 4], [3, 4, 5]])
    N.row_del(-2)
    assert N == Matrix([[1, 2, 3], [3, 4, 5]])
    O = Matrix([[1, 2, 3], [5, 6, 7], [9, 10, 11]])
    O.row_del(-1)
    assert O == Matrix([[1, 2, 3], [5, 6, 7]])
    P = Matrix([[1, 2, 3], [2, 3, 4], [3, 4, 5]])
    pytest.raises(IndexError, lambda: P.row_del(10))
    Q = Matrix([[1, 2, 3], [2, 3, 4], [3, 4, 5]])
    pytest.raises(IndexError, lambda: Q.row_del(-10))

    # for col_del(index)
    M = Matrix([[1, 2, 3], [2, 3, 4], [3, 4, 5]])
    M.col_del(1)
    assert M == Matrix([[1, 3], [2, 4], [3, 5]])
    N = Matrix([[1, 2, 3], [2, 3, 4], [3, 4, 5]])
    N.col_del(-2)
    assert N == Matrix([[1, 3], [2, 4], [3, 5]])
    P = Matrix([[1, 2, 3], [2, 3, 4], [3, 4, 5]])
    pytest.raises(IndexError, lambda: P.col_del(10))
    Q = Matrix([[1, 2, 3], [2, 3, 4], [3, 4, 5]])
    pytest.raises(IndexError, lambda: Q.col_del(-10))


def test_sympyissue_9422():
    x, y = symbols('x y', commutative=False)
    M = eye(2)
    M1 = Matrix(2, 2, [x, y, y, z])
    assert y*x*M != x*y*M
    assert b*a*M == a*b*M
    assert x*M1 != M1*x
    assert a*M1 == M1*a
    assert y*x*M == Matrix([[y*x, 0], [0, y*x]])


def test_sympyissue_9480():
    m = Matrix([[-5 + 5*sqrt(2), -5],
                [-5*sqrt(2)/2 + 5, -5*sqrt(2)/2]])
    assert m.rank() == 1


def test_diofantissue_288():
    m = Matrix([[-exp(I*k)*I/(4*k) + Rational(1, 2) + exp(-I*k)*I/(4*k),
                 exp(I*k)*I/(4*k) + Rational(1, 2) - exp(-I*k)*I/(4*k),
                 exp(I*k)/4 + Rational(1, 2) + exp(-I*k)/4,
                 -exp(I*k)/4 - Rational(1, 2) - exp(-I*k)/4],
                [exp(I*k)*I/(4*k) + Rational(1, 2) - exp(-I*k)*I/(4*k),
                 -exp(I*k)*I/(4*k) + Rational(1, 2) + exp(-I*k)*I/(4*k),
                 -exp(I*k)/4 - Rational(1, 2) - exp(-I*k)/4,
                 exp(I*k)/4 + Rational(1, 2) + exp(-I*k)/4],
                [exp(I*k)/4 + Rational(1, 2) + exp(-I*k)/4,
                 -exp(I*k)/4 - Rational(1, 2) - exp(-I*k)/4,
                 exp(I*k)*I*k/4 - exp(-I*k)*I*k/4,
                 -exp(I*k)*I*k/4 + exp(-I*k)*I*k/4],
                [-exp(I*k)/4 - 1/2 - exp(-I*k)/4,
                 exp(I*k)/4 + 1/2 + exp(-I*k)/4,
                 -exp(I*k)*I*k/4 + exp(-I*k)*I*k/4,
                 exp(I*k)*I*k/4 - exp(-I*k)*I*k/4]])
    assert m.det() == 0
    assert m.rank() != 4


def test_sympyissue_11434():
    ax, ay, bx, by, cx, cy, dx, dy, ex, ey, t0, t1 = \
        symbols('a_x a_y b_x b_y c_x c_y d_x d_y e_x e_y t_0 t_1')
    M = Matrix([[ax, ay, ax*t0, ay*t0, 0],
                [bx, by, bx*t0, by*t0, 0],
                [cx, cy, cx*t0, cy*t0, 1],
                [dx, dy, dx*t0, dy*t0, 1],
                [ex, ey, 2*ex*t1 - ex*t0, 2*ey*t1 - ey*t0, 0]])
    assert M.det() == 0
    assert M.rank() == 4


def test_solve_least_squares():
    M = Matrix([[1, 2], [2, 3], [3, 4]])
    r = M.solve_least_squares(Matrix([8, 14, 18]))
    ans = Matrix([[Rational(5, 3)], [Rational(10, 3)]])
    assert r == ans
    r = M.solve_least_squares(Matrix([8, 14, 18]), method=None)
    assert r == ans


def test_mgamma():
    assert mgamma(0) == Matrix(((1, 0, 0, 0),
                                (0, 1, 0, 0),
                                (0, 0, -1, 0),
                                (0, 0, 0, -1)))
    assert mgamma(1) == Matrix(((0, 0, 0, 1),
                                (0, 0, 1, 0),
                                (0, -1, 0, 0),
                                (-1, 0, 0, 0)))
    assert mgamma(2) == Matrix(((0, 0, 0, -I),
                                (0, 0, I, 0),
                                (0, I, 0, 0),
                                (-I, 0, 0, 0)))
    assert mgamma(3) == Matrix(((0, 0, 1, 0),
                                (0, 0, 0, -1),
                                (-1, 0, 0, 0),
                                (0, 1, 0, 0)))
    assert mgamma(5) == Matrix(((0, 0, 1, 0),
                                (0, 0, 0, 1),
                                (1, 0, 0, 0),
                                (0, 1, 0, 0)))
    assert mgamma(0, lower=True) == mgamma(0)
    assert mgamma(1, lower=True) == -mgamma(1)
    pytest.raises(IndexError, lambda: mgamma(4))


def test_sympyissue_10770():
    M = Matrix([])
    a = ['col_insert', 'row_join'], Matrix([9, 6, 3])
    b = ['row_insert', 'col_join'], a[1].T
    c = ['row_insert', 'col_insert'], Matrix([[1, 2], [3, 4]])
    for ops, m in (a, b, c):
        for op in ops:
            f = getattr(M, op)
            new = f(m) if 'join' in op else f(42, m)
            assert new == m and id(new) != id(m)
