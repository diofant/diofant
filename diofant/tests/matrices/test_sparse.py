import pytest

from diofant import (I, Matrix, MutableDenseMatrix, MutableSparseMatrix,
                     PurePoly, Rational, ShapeError, SparseMatrix, eye, ones,
                     zeros)
from diofant.abc import x, y, z


__all__ = ()


def test_sparse_matrix():
    def sparse_eye(n):
        return SparseMatrix.eye(n)

    def sparse_zeros(n):
        return SparseMatrix.zeros(n)

    # creation args
    pytest.raises(TypeError, lambda: SparseMatrix(1, 2))
    pytest.raises(ValueError, lambda: SparseMatrix(2, 2, (1, 3, 4, 5, 6)))

    a = SparseMatrix((
        (1, 0),
        (0, 1)))
    assert SparseMatrix(a) == a

    a = MutableSparseMatrix([])
    b = MutableDenseMatrix([1, 2])
    assert a.row_join(b) == b
    assert a.col_join(b) == b
    assert type(a.row_join(b)) == type(a)
    assert type(a.col_join(b)) == type(a)

    # test element assignment
    a = SparseMatrix((
        (1, 0),
        (0, 1)))

    a[3] = 4
    assert a[1, 1] == 4
    a[3] = 1

    a[0, 0] = 2
    assert a == SparseMatrix((
        (2, 0),
        (0, 1)))
    a[1, 0] = 5
    assert a == SparseMatrix((
        (2, 0),
        (5, 1)))
    a[1, 1] = 0
    assert a == SparseMatrix((
        (2, 0),
        (5, 0)))
    assert a._smat == {(0, 0): 2, (1, 0): 5}

    # test_multiplication
    a = SparseMatrix((
        (1, 2),
        (3, 1),
        (0, 6)))

    b = SparseMatrix((
        (1, 2),
        (3, 0)))

    c = a*b
    assert c[0, 0] == 7
    assert c[0, 1] == 2
    assert c[1, 0] == 6
    assert c[1, 1] == 6
    assert c[2, 0] == 18
    assert c[2, 1] == 0

    c = b * x
    assert isinstance(c, SparseMatrix)
    assert c[0, 0] == x
    assert c[0, 1] == 2*x
    assert c[1, 0] == 3*x
    assert c[1, 1] == 0

    c = 5 * b
    assert isinstance(c, SparseMatrix)
    assert c[0, 0] == 5
    assert c[0, 1] == 2*5
    assert c[1, 0] == 3*5
    assert c[1, 1] == 0

    # test_power
    A = SparseMatrix([[2, 3], [4, 5]])
    assert (A**5)[:] == [6140, 8097, 10796, 14237]
    A = SparseMatrix([[2, 1, 3], [4, 2, 4], [6, 12, 1]])
    assert (A**3)[:] == [290, 262, 251, 448, 440, 368, 702, 954, 433]

    # test_creation
    a = SparseMatrix([[x, 0], [0, 0]])
    m = a
    assert m.cols == m.rows
    assert m.cols == 2
    assert m[:] == [x, 0, 0, 0]
    b = SparseMatrix(2, 2, [x, 0, 0, 0])
    m = b
    assert m.cols == m.rows
    assert m.cols == 2
    assert m[:] == [x, 0, 0, 0]

    assert a == b
    S = sparse_eye(3)
    del S[1, :]
    assert S == SparseMatrix([
                             [1, 0, 0],
                             [0, 0, 1]])
    S = sparse_eye(3)
    del S[:, 1]
    assert S == SparseMatrix([
                             [1, 0],
                             [0, 0],
                             [0, 1]])
    S = SparseMatrix.eye(3)
    S[2, 1] = 2
    S.col_swap(1, 0)
    assert S == SparseMatrix([[0, 1, 0],
                              [1, 0, 0],
                              [2, 0, 1]])
    S.row_swap(0, 1)
    assert S == SparseMatrix([[1, 0, 0],
                              [0, 1, 0],
                              [2, 0, 1]])
    S.col_swap(0, 1)
    assert S == SparseMatrix([[0, 1, 0],
                              [1, 0, 0],
                              [0, 2, 1]])
    S.row_swap(0, 2)
    assert S == SparseMatrix([[0, 2, 1],
                              [1, 0, 0],
                              [0, 1, 0]])

    S.col_swap(0, 2)
    assert S == SparseMatrix([[1, 2, 0],
                              [0, 0, 1],
                              [0, 1, 0]])

    a = SparseMatrix(1, 2, [1, 2])
    b = a.copy()
    c = a.copy()
    assert a[0] == 1
    del a[0, :]
    assert a == SparseMatrix(0, 2, [])
    del b[:, 1]
    assert b == SparseMatrix(1, 1, [1])

    # test_determinant
    assert SparseMatrix(1, 1, [0]).det() == 0
    assert SparseMatrix([[1]]).det() == 1
    assert SparseMatrix(((-3, 2), (8, -5))).det() == -1
    assert SparseMatrix(((x, 1), (y, 2*y))).det() == 2*x*y - y
    assert SparseMatrix(((1, 1, 1),
                         (1, 2, 3),
                         (1, 3, 6))).det() == 1
    assert SparseMatrix((( 3, -2,  0, 5),
                         (-2,  1, -2, 2),
                         ( 0, -2,  5, 0),
                         ( 5,  0,  3, 4))).det() == -289
    assert SparseMatrix((( 1,  2,  3,  4),
                         ( 5,  6,  7,  8),
                         ( 9, 10, 11, 12),
                         (13, 14, 15, 16))).det() == 0
    assert SparseMatrix(((3, 2, 0, 0, 0),
                         (0, 3, 2, 0, 0),
                         (0, 0, 3, 2, 0),
                         (0, 0, 0, 3, 2),
                         (2, 0, 0, 0, 3))).det() == 275
    assert SparseMatrix(((1, 0,  1,  2, 12),
                         (2, 0,  1,  1,  4),
                         (2, 1,  1, -1,  3),
                         (3, 2, -1,  1,  8),
                         (1, 1,  1,  0,  6))).det() == -55
    assert SparseMatrix(((-5,  2,  3,  4,  5),
                         ( 1, -4,  3,  4,  5),
                         ( 1,  2, -3,  4,  5),
                         ( 1,  2,  3, -2,  5),
                         ( 1,  2,  3,  4, -1))).det() == 11664
    assert SparseMatrix((( 2,  7, -1, 3, 2),
                         ( 0,  0,  1, 0, 1),
                         (-2,  0,  7, 0, 2),
                         (-3, -2,  4, 5, 3),
                         ( 1,  0,  0, 0, 1))).det() == 123

    # test_slicing
    m0 = sparse_eye(4)
    assert m0[:3, :3] == sparse_eye(3)
    assert m0[2:4, 0:2] == sparse_zeros(2)

    m1 = SparseMatrix(3, 3, lambda i, j: i + j)
    assert m1[0, :] == SparseMatrix(1, 3, (0, 1, 2))
    assert m1[1:3, 1] == SparseMatrix(2, 1, (2, 3))

    m2 = SparseMatrix(
        [[0, 1, 2, 3], [4, 5, 6, 7], [8, 9, 10, 11], [12, 13, 14, 15]])
    assert m2[:, -1] == SparseMatrix(4, 1, [3, 7, 11, 15])
    assert m2[-2:, :] == SparseMatrix([[8, 9, 10, 11], [12, 13, 14, 15]])

    assert SparseMatrix([[1, 2], [3, 4]])[[1], [1]] == Matrix([[4]])

    # test_submatrix_assignment
    m = sparse_zeros(4)
    m[2:4, 2:4] = sparse_eye(2)
    assert m == SparseMatrix([(0, 0, 0, 0),
                              (0, 0, 0, 0),
                              (0, 0, 1, 0),
                              (0, 0, 0, 1)])
    assert len(m._smat) == 2
    m[:2, :2] = sparse_eye(2)
    assert m == sparse_eye(4)
    m[:, 0] = SparseMatrix(4, 1, (1, 2, 3, 4))
    assert m == SparseMatrix([(1, 0, 0, 0),
                              (2, 1, 0, 0),
                              (3, 0, 1, 0),
                              (4, 0, 0, 1)])
    m[:, :] = sparse_zeros(4)
    assert m == sparse_zeros(4)
    m[:, :] = ((1, 2, 3, 4), (5, 6, 7, 8), (9, 10, 11, 12), (13, 14, 15, 16))
    assert m == SparseMatrix((( 1,  2,  3,  4),
                              ( 5,  6,  7,  8),
                              ( 9, 10, 11, 12),
                              (13, 14, 15, 16)))
    m[:2, 0] = [0, 0]
    assert m == SparseMatrix((( 0,  2,  3,  4),
                              ( 0,  6,  7,  8),
                              ( 9, 10, 11, 12),
                              (13, 14, 15, 16)))

    # test_reshape
    m0 = sparse_eye(3)
    assert m0.reshape(1, 9) == SparseMatrix(1, 9, (1, 0, 0, 0, 1, 0, 0, 0, 1))
    m1 = SparseMatrix(3, 4, lambda i, j: i + j)
    assert m1.reshape(4, 3) == \
        SparseMatrix([(0, 1, 2), (3, 1, 2), (3, 4, 2), (3, 4, 5)])
    assert m1.reshape(2, 6) == \
        SparseMatrix([(0, 1, 2, 3, 1, 2), (3, 4, 2, 3, 4, 5)])

    # test_applyfunc
    m0 = sparse_eye(3)
    assert m0.applyfunc(lambda x: 2*x) == sparse_eye(3)*2
    assert m0.applyfunc(lambda x: 0) == sparse_zeros(3)

    # test_LUdecomp
    testmat = SparseMatrix([[ 0, 2, 5, 3],
                            [ 3, 3, 7, 4],
                            [ 8, 4, 0, 2],
                            [-2, 6, 3, 4]])
    L, U, p = testmat.LUdecomposition()
    assert L.is_lower
    assert U.is_upper
    assert (L*U).permuteBkwd(p) - testmat == sparse_zeros(4)

    testmat = SparseMatrix([[ 6, -2, 7, 4],
                            [ 0,  3, 6, 7],
                            [ 1, -2, 7, 4],
                            [-9,  2, 6, 3]])
    L, U, p = testmat.LUdecomposition()
    assert L.is_lower
    assert U.is_upper
    assert (L*U).permuteBkwd(p) - testmat == sparse_zeros(4)

    M = Matrix(((1, x, 1), (2, y, 0), (y, 0, z)))
    L, U, p = M.LUdecomposition()
    assert L.is_lower
    assert U.is_upper
    assert (L*U).permuteBkwd(p) - M == sparse_zeros(3)

    # test_LUsolve
    A = SparseMatrix([[2, 3, 5],
                      [3, 6, 2],
                      [8, 3, 6]])
    B = SparseMatrix(3, 1, [3, 7, 5])
    b = A*B
    soln = A.LUsolve(b)
    assert soln == B
    A = SparseMatrix([[0, -1, 2],
                      [5, 10, 7],
                      [8,  3, 4]])
    B = SparseMatrix(3, 1, [-1, 2, 5])
    b = A*B
    soln = A.LUsolve(b)
    assert soln == B

    # test_inverse
    A = sparse_eye(4)
    assert A.inv() == sparse_eye(4)
    assert A.inv(method='CH') == sparse_eye(4)
    assert A.inv(method='LDL') == sparse_eye(4)

    A = SparseMatrix([[2, 3, 5],
                      [3, 6, 2],
                      [7, 2, 6]])
    Ainv = SparseMatrix(Matrix(A).inv())
    assert A*Ainv == sparse_eye(3)
    assert A.inv(method='CH') == Ainv
    assert A.inv(method='LDL') == Ainv

    A = SparseMatrix([[2, 3, 5],
                      [3, 6, 2],
                      [5, 2, 6]])
    Ainv = SparseMatrix(Matrix(A).inv())
    assert A*Ainv == sparse_eye(3)
    assert A.inv(method='CH') == Ainv
    assert A.inv(method='LDL') == Ainv

    # test_cross
    v1 = Matrix(1, 3, [1, 2, 3])
    v2 = Matrix(1, 3, [3, 4, 5])
    assert v1.cross(v2) == Matrix(1, 3, [-2, 4, -2])
    assert v1.norm(2)**2 == 14

    # conjugate
    a = SparseMatrix(((1, 2 + I), (3, 4)))
    assert a.C == SparseMatrix([
        [1, 2 - I],
        [3,     4]])

    # mul
    assert a*Matrix(2, 2, [1, 0, 0, 1]) == a
    assert a + Matrix(2, 2, [1, 1, 1, 1]) == SparseMatrix([
        [2, 3 + I],
        [4,     5]])
    assert a*0 == Matrix([[0, 0], [0, 0]])

    # col join
    assert a.col_join(sparse_eye(2)) == SparseMatrix([
        [1, 2 + I],
        [3,     4],
        [1,     0],
        [0,     1]])
    A = SparseMatrix(ones(3))
    B = eye(3)
    assert A.col_join(B) == Matrix([[1, 1, 1], [1, 1, 1], [1, 1, 1],
                                    [1, 0, 0], [0, 1, 0], [0, 0, 1]])

    # row join
    A = SparseMatrix(((1, 0, 1), (0, 1, 0), (1, 1, 0)))
    B = Matrix(((1, 0, 0), (0, 1, 0), (0, 0, 1)))
    assert A.row_join(B) == Matrix([[1, 0, 1, 1, 0, 0],
                                    [0, 1, 0, 0, 1, 0],
                                    [1, 1, 0, 0, 0, 1]])

    # symmetric
    assert not a.is_symmetric(simplify=False)
    assert sparse_eye(3).is_symmetric(simplify=False)

    # test_cofactor
    assert sparse_eye(3) == sparse_eye(3).cofactorMatrix()
    test = SparseMatrix([[1, 3, 2], [2, 6, 3], [2, 3, 6]])
    assert test.cofactorMatrix() == \
        SparseMatrix([[27, -6, -6], [-12, 2, 3], [-3, 1, 0]])
    test = SparseMatrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    assert test.cofactorMatrix() == \
        SparseMatrix([[-3, 6, -3], [6, -12, 6], [-3, 6, -3]])

    # test_jacobian
    L = SparseMatrix(1, 2, [x**2*y, 2*y**2 + x*y])
    syms = [x, y]
    assert L.jacobian(syms) == Matrix([[2*x*y, x**2], [y, 4*y + x]])

    L = SparseMatrix(1, 2, [x, x**2*y**3])
    assert L.jacobian(syms) == SparseMatrix([[1, 0], [2*x*y**3, x**2*3*y**2]])

    # test_QR
    A = Matrix([[1, 2], [2, 3]])
    Q, S = A.QRdecomposition()
    R = Rational
    assert Q == Matrix([
        [  5**R(-1, 2),  (R(2)/5)*(R(1)/5)**R(-1, 2)],
        [2*5**R(-1, 2), (-R(1)/5)*(R(1)/5)**R(-1, 2)]])
    assert S == Matrix([
        [5**R(1, 2),     8*5**R(-1, 2)],
        [         0, (R(1)/5)**R(1, 2)]])
    assert Q*S == A
    assert Q.T * Q == sparse_eye(2)

    R = Rational
    # test nullspace
    # first test reduced row-ech form

    M = SparseMatrix([[5, 7, 2, 1],
                      [1, 6, 2, -1]])
    out, _ = M.rref()
    assert out == Matrix([[1, 0, -R(2)/23, R(13)/23],
                          [0, 1,  R(8)/23, R(-6)/23]])

    M = SparseMatrix([[ 1,  3, 0,  2,  6, 3, 1],
                      [-2, -6, 0, -2, -8, 3, 1],
                      [ 3,  9, 0,  0,  6, 6, 2],
                      [-1, -3, 0,  1,  0, 9, 3]])

    out, _ = M.rref()
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

    # test eigen
    sparse_eye3 = sparse_eye(3)
    assert sparse_eye3.charpoly(x) == PurePoly((x - 1)**3)
    assert sparse_eye3.charpoly(y) == PurePoly((y - 1)**3)

    # test values
    M = Matrix([( 0, 1, -1),
                ( 1, 1,  0),
                (-1, 0,  1)])
    vals = M.eigenvals()
    assert sorted(vals) == [-1, 1, 2]

    R = Rational
    M = Matrix([[1, 0, 0],
                [0, 1, 0],
                [0, 0, 1]])
    assert M.eigenvects() == [(1, 3, [
        Matrix([1, 0, 0]),
        Matrix([0, 1, 0]),
        Matrix([0, 0, 1])])]
    M = Matrix([[5, 0, 2],
                [3, 2, 0],
                [0, 0, 1]])
    assert M.eigenvects() == [(1, 1, [Matrix([R(-1)/2, R(3)/2, 1])]),
                              (2, 1, [Matrix([0, 1, 0])]),
                              (5, 1, [Matrix([1, 1, 0])])]

    assert M.zeros(3, 5) == SparseMatrix(3, 5, {})
    A = SparseMatrix(10, 10, {(0, 0): 18, (0, 9): 12, (1, 4): 18, (2, 7): 16, (3, 9): 12, (4, 2): 19, (5, 7): 16, (6, 2): 12, (9, 7): 18})
    assert A.row_list() == [(0, 0, 18), (0, 9, 12), (1, 4, 18), (2, 7, 16), (3, 9, 12), (4, 2, 19), (5, 7, 16), (6, 2, 12), (9, 7, 18)]
    assert A.col_list() == [(0, 0, 18), (4, 2, 19), (6, 2, 12), (1, 4, 18), (2, 7, 16), (5, 7, 16), (9, 7, 18), (0, 9, 12), (3, 9, 12)]
    assert SparseMatrix.eye(2).nnz() == 2

    M = SparseMatrix.eye(3)*2
    M[1, 0] = -1
    M.col_op(1, lambda v, i: v + 2*M[i, 0])
    assert M == Matrix([[ 2, 4, 0], [-1, 0, 0], [ 0, 0, 2]])

    M = SparseMatrix.zeros(3)
    M.fill(1)
    assert M == ones(3)

    assert SparseMatrix(ones(0, 3)).tolist() == []


def test_eq():
    A = SparseMatrix(((1, 2), (3, 4)))
    assert A != 1
    assert A != zeros(2, 1)


def test_transpose():
    assert SparseMatrix(((1, 2), (3, 4))).transpose() == \
        SparseMatrix(((1, 3), (2, 4)))


def test_trace():
    assert SparseMatrix(((1, 2), (3, 4))).trace() == 5
    assert SparseMatrix(((0, 0), (0, 4))).trace() == 4


def test_CL_RL():
    assert SparseMatrix(((1, 2), (3, 4))).row_list() == \
        [(0, 0, 1), (0, 1, 2), (1, 0, 3), (1, 1, 4)]
    assert SparseMatrix(((1, 2), (3, 4))).col_list() == \
        [(0, 0, 1), (1, 0, 3), (0, 1, 2), (1, 1, 4)]


def test_add():
    assert SparseMatrix(((1, 0), (0, 1))) + SparseMatrix(((0, 1), (1, 0))) == \
        SparseMatrix(((1, 1), (1, 1)))
    a = SparseMatrix(100, 100, lambda i, j: int(j != 0 and i % j == 0))
    b = SparseMatrix(100, 100, lambda i, j: int(i != 0 and j % i == 0))
    assert len(a._smat) + len(b._smat) - len((a + b)._smat) > 0


def test_errors():
    pytest.raises(ValueError, lambda: SparseMatrix(1.4, 2, lambda i, j: 0))
    pytest.raises(ValueError, lambda: SparseMatrix(2, 2, 1))
    pytest.raises(TypeError, lambda: SparseMatrix([1, 2, 3], [1, 2]))
    pytest.raises(ValueError, lambda: SparseMatrix([[1, 2], [3, 4]])[(1, 2, 3)])
    pytest.raises(IndexError, lambda: SparseMatrix([[1, 2], [3, 4]])[5])
    pytest.raises(ValueError, lambda: SparseMatrix([[1, 2], [3, 4]])[1, 2, 3])
    pytest.raises(TypeError,
                  lambda: SparseMatrix([[1, 2],
                                        [3, 4]]).copyin_list([0, 1], set()))
    pytest.raises(IndexError, lambda: SparseMatrix([[1, 2], [3, 4]])[1, 2])
    pytest.raises(TypeError, lambda: SparseMatrix([1, 2, 3]).cross(1))
    pytest.raises(IndexError, lambda: SparseMatrix(1, 2, [1, 2])[3])
    pytest.raises(ShapeError,
                  lambda: SparseMatrix(1, 2,
                                       [1, 2]) + SparseMatrix(2, 1, [2, 1]))
    pytest.raises(IndexError, lambda: SparseMatrix([1, 2, 3])[3, 0])
    pytest.raises(TypeError, lambda: SparseMatrix([1, 2, 3]).applyfunc(1))
    pytest.raises(ValueError, lambda: SparseMatrix([1, 2, 3]).reshape(2, 2))
    pytest.raises(ValueError,
                  lambda: SparseMatrix([[2, 3], [4, 1]]).cholesky())
    pytest.raises(ValueError,
                  lambda: SparseMatrix([[2, 3], [4, 1]]).LDLdecomposition())
    pytest.raises(ValueError, lambda: SparseMatrix([[2, 3], [4, 1]]).add(1))
    pytest.raises(ShapeError,
                  lambda: SparseMatrix([[1, 2],
                                        [3, 4]]).row_join(Matrix([[1, 2]])))
    pytest.raises(ShapeError,
                  lambda: SparseMatrix([[1, 2],
                                        [3, 4]]).col_join(Matrix([1, 2])))
    pytest.raises(ShapeError,
                  lambda: SparseMatrix([[1, 2],
                                        [3, 4]]).copyin_matrix([1, 0],
                                                               Matrix([1, 2])))


def test_len():
    assert not SparseMatrix()
    assert SparseMatrix() == SparseMatrix([])
    assert SparseMatrix() == SparseMatrix([[]])


def test_sparse_zeros_sparse_eye():
    assert SparseMatrix.eye(3) == eye(3, cls=SparseMatrix)
    assert len(SparseMatrix.eye(3)._smat) == 3
    assert SparseMatrix.zeros(3) == zeros(3, cls=SparseMatrix)
    assert len(SparseMatrix.zeros(3)._smat) == 0


def test_copyin():
    s = SparseMatrix(3, 3, {})
    s[1, 0] = 1
    assert s[:, 0] == SparseMatrix(Matrix([0, 1, 0]))
    assert s[3] == 1
    assert s[3: 4] == [1]
    s[1, 1] = 42
    assert s[1, 1] == 42
    assert s[1, 1:] == SparseMatrix([[42, 0]])
    s[1, 1:] = Matrix([[5, 6]])
    assert s[1, :] == SparseMatrix([[1, 5, 6]])
    s[1, 1:] = [[42, 43]]
    assert s[1, :] == SparseMatrix([[1, 42, 43]])
    s[0, 0] = 17
    assert s[:, :1] == SparseMatrix([17, 1, 0])
    s[0, 0] = [1, 1, 1]
    assert s[:, 0] == SparseMatrix([1, 1, 1])
    s[0, 0] = Matrix([1, 1, 1])
    assert s[:, 0] == SparseMatrix([1, 1, 1])
    s[0, 0] = SparseMatrix([1, 1, 1])
    assert s[:, 0] == SparseMatrix([1, 1, 1])


def test_sparse_solve():
    A = SparseMatrix(((25, 15, -5), (15, 18, 0), (-5, 0, 11)))
    assert A.cholesky() == Matrix([
        [ 5, 0, 0],
        [ 3, 3, 0],
        [-1, 1, 3]])
    assert A.cholesky() * A.cholesky().T == Matrix([
        [25, 15, -5],
        [15, 18, 0],
        [-5, 0, 11]])

    A = SparseMatrix(((25, 15, -5), (15, 18, 0), (-5, 0, 11)))
    L, D = A.LDLdecomposition()
    assert 15*L == Matrix([
        [15, 0, 0],
        [ 9, 15, 0],
        [-3, 5, 15]])
    assert D == Matrix([
        [25, 0, 0],
        [ 0, 9, 0],
        [ 0, 0, 9]])
    assert L * D * L.T == A

    A = SparseMatrix(((3, 0, 2), (0, 0, 1), (1, 2, 0)))
    assert A.inv() * A == SparseMatrix(eye(3))

    A = SparseMatrix([
        [ 2, -1, 0],
        [-1, 2, -1],
        [ 0, 0, 2]])
    ans = SparseMatrix([
        [Rational(2, 3), Rational(1, 3), Rational(1, 6)],
        [Rational(1, 3), Rational(2, 3), Rational(1, 3)],
        [     0,      0, Rational(1, 2)]])
    assert A.inv(method='CH') == ans
    assert A.inv(method='LDL') == ans
    assert A * ans == SparseMatrix(eye(3))

    s = A.solve(A[:, 0], 'LDL')
    assert A*s == A[:, 0]
    s = A.solve(A[:, 0], 'CH')
    assert A*s == A[:, 0]
    A = A.col_join(A)
    s = A.solve_least_squares(A[:, 0], 'CH')
    assert A*s == A[:, 0]
    s = A.solve_least_squares(A[:, 0], 'LDL')
    assert A*s == A[:, 0]

    pytest.raises(ValueError, lambda: SparseMatrix([[1, 0, 1],
                                                    [0, 0, 1]]).solve([1, 1]))
    pytest.raises(ValueError, lambda: SparseMatrix([[1, 0], [0, 0],
                                                    [2, 1]]).solve([1, 1, 1]))


def test_hermitian():
    a = SparseMatrix([[0, I], [-I, 0]])
    assert a.is_hermitian
    a = SparseMatrix([[1, I], [-I, 1]])
    assert a.is_hermitian
    a[0, 0] = 2*I
    assert a.is_hermitian is False
    a[0, 0] = x
    assert a.is_hermitian is None
    a[0, 1] = a[1, 0]*I
    assert a.is_hermitian is False


def test_fill():
    a = SparseMatrix([[0, I], [-I, 0]])
    a.fill(0)
    assert a == Matrix([[0, 0], [0, 0]])
