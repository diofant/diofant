import collections
import functools
import math
from types import FunctionType

from ..core import (Add, Atom, Basic, Dummy, Expr, Float, I, Integer, Pow,
                    Symbol, count_ops, oo, symbols)
from ..core.compatibility import as_int, is_sequence
from ..core.logic import fuzzy_and
from ..core.sympify import sympify
from ..functions import Max, Min, exp, factorial, sqrt
from ..polys import PurePoly, cancel, gcd, roots
from ..printing.defaults import DefaultPrinting
from ..simplify import nsimplify, signsimp
from ..simplify import simplify as _simplify
from ..utilities import default_sort_key, flatten


def _iszero(x):
    """Returns True if x is zero."""
    r = x.equals(0)
    if r is None:
        raise NotImplementedError(f'Zero-decision problem for {x}')
    return r


class MatrixError(Exception):
    """Generic matrix error."""


class ShapeError(ValueError, MatrixError):
    """Wrong matrix shape."""


class NonSquareMatrixError(ShapeError):
    """Raised when a square matrix is expected."""


class MatrixBase(DefaultPrinting):
    """Base class for matrices."""

    # Added just for numpy compatibility
    __array_priority__ = 11

    is_Matrix = True
    is_Identity: bool | None = None
    _class_priority = 3
    _sympify = staticmethod(sympify)

    __hash__ = None  # type: ignore[assignment]

    @classmethod
    def _handle_creation_inputs(cls, *args, **kwargs):
        """Return the number of rows, cols and flat matrix elements.

        Examples
        ========

        Matrix can be constructed as follows:

        * from a nested list of iterables

        >>> Matrix(((1, 2 + I), (3, 4)))
        Matrix([
        [1, 2 + I],
        [3,     4]])

        * from un-nested iterable (interpreted as a column)

        >>> Matrix([1, 2])
        Matrix([
        [1],
        [2]])

        * from un-nested iterable with dimensions

        >>> Matrix(1, 2, [1, 2])
        Matrix([[1, 2]])

        * from no arguments (a 0 x 0 matrix)

        >>> Matrix()
        Matrix(0, 0, [])

        * from a rule

        >>> Matrix(2, 2, lambda i, j: i/(j + 1))
        Matrix([
        [0,   0],
        [1, 1/2]])

        """
        from .sparse import SparseMatrixBase

        flat_list = None

        if len(args) == 1:
            # Matrix(SparseMatrix(...))
            if isinstance(args[0], SparseMatrixBase):
                return args[0].rows, args[0].cols, flatten(args[0].tolist())

            # Matrix(Matrix(...))
            elif isinstance(args[0], MatrixBase):
                return args[0].rows, args[0].cols, args[0]._mat

            # Matrix(MatrixSymbol('X', 2, 2))
            elif isinstance(args[0], Basic) and args[0].is_Matrix:
                return args[0].rows, args[0].cols, args[0].as_explicit()._mat

            # Matrix(numpy.ones((2, 2)))
            elif hasattr(args[0], '__array__'):
                # NumPy array or matrix or some other object that implements
                # __array__. So let's first use this method to get a
                # numpy.array() and then make a python list out of it.
                arr = args[0].__array__()
                if len(arr.shape) == 2:
                    rows, cols = arr.shape[0], arr.shape[1]
                    flat_list = [cls._sympify(i) for i in arr.ravel()]
                    return rows, cols, flat_list
                elif len(arr.shape) == 1:
                    rows, cols = arr.shape[0], 1
                    flat_list = [Integer(0)]*rows
                    for i, a in enumerate(arr):
                        flat_list[i] = cls._sympify(a)
                    return rows, cols, flat_list
                else:
                    raise NotImplementedError(
                        'Diofant supports just 1D and 2D matrices')

            # Matrix([1, 2, 3]) or Matrix([[1, 2], [3, 4]])
            elif is_sequence(args[0]):
                in_mat = []
                ncol = set()
                for row in args[0]:
                    if isinstance(row, MatrixBase):
                        in_mat.extend(row.tolist())
                        assert row.cols or row.rows
                        ncol.add(row.cols)
                    else:
                        in_mat.append(row)
                        try:
                            ncol.add(len(row))
                        except TypeError:
                            ncol.add(1)
                if len(ncol) > 1:
                    raise ValueError(f'Got rows of variable lengths: {sorted(ncol)}')
                cols = ncol.pop() if ncol else 0
                rows = len(in_mat) if cols else 0
                if rows:
                    if not is_sequence(in_mat[0]):
                        cols = 1
                        flat_list = [cls._sympify(i) for i in in_mat]
                        return rows, cols, flat_list
                flat_list = []
                for j in range(rows):
                    for i in range(cols):
                        flat_list.append(cls._sympify(in_mat[j][i]))

        elif len(args) == 3:
            rows = as_int(args[0])
            cols = as_int(args[1])

            # Matrix(2, 2, lambda i, j: i+j)
            if len(args) == 3 and isinstance(args[2], collections.abc.Callable):
                op = args[2]
                flat_list = []
                for i in range(rows):
                    flat_list.extend(
                        [cls._sympify(op(cls._sympify(i), cls._sympify(j)))
                         for j in range(cols)])

            # Matrix(2, 2, [1, 2, 3, 4])
            elif len(args) == 3 and is_sequence(args[2]):
                flat_list = args[2]
                if len(flat_list) != rows*cols:
                    raise ValueError('List length should be equal to rows*columns')
                flat_list = [cls._sympify(i) for i in flat_list]

        # Matrix()
        elif len(args) == 0:
            # Empty Matrix
            rows = cols = 0
            flat_list = []

        if flat_list is None:
            raise TypeError('Data type not understood')

        return rows, cols, flat_list

    def _setitem(self, key, value):
        """Helper to set value at location given by key.

        Examples
        ========

        >>> m = Matrix(((1, 2+I), (3, 4)))
        >>> m
        Matrix([
        [1, 2 + I],
        [3,     4]])
        >>> m[1, 0] = 9
        >>> m
        Matrix([
        [1, 2 + I],
        [9,     4]])
        >>> m[1, 0] = [[0, 1]]

        To replace row r you assign to position r*m where m
        is the number of columns:

        >>> M = zeros(4)
        >>> m = M.cols
        >>> M[3*m] = ones(1, m)*2
        >>> M
        Matrix([
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        [2, 2, 2, 2]])

        And to replace column c you can assign to position c:

        >>> M[2] = ones(m, 1)*4
        >>> M
        Matrix([
        [0, 0, 4, 0],
        [0, 0, 4, 0],
        [0, 0, 4, 0],
        [2, 2, 4, 2]])

        """
        from . import Matrix

        is_slice = isinstance(key, slice)
        i, j = key = self.key2ij(key)
        is_mat = isinstance(value, MatrixBase)
        if type(i) is slice or type(j) is slice:
            if is_mat:
                self.copyin_matrix(key, value)
                return
            if not isinstance(value, Expr) and is_sequence(value):
                self.copyin_list(key, value)
                return
            raise ValueError(f'unexpected value: {value}')
        if not is_mat and not isinstance(value, Basic) and is_sequence(value):
            value = Matrix(value)
            is_mat = True
        if is_mat:
            assert not is_slice
            key = (slice(i, i + value.rows), slice(j, j + value.cols))
            self.copyin_matrix(key, value)
        else:
            return i, j, self._sympify(value)

    def copy(self):
        """Returns the copy of a matrix."""
        return self._new(self.rows, self.cols, self._mat)

    def trace(self):
        """Returns the trace of a matrix."""
        if not self.is_square:
            raise NonSquareMatrixError()
        return self._eval_trace()

    def inv(self, method=None, **kwargs):
        r"""Returns the inverse of the matrix.

        Parameters
        ==========

        method : {'GE', 'LU', 'ADJ', 'CH', 'LDL'} or None
            Selects algorithm for inversion.  For dense matrices
            available {'GE', 'LU', 'ADJ'}, default is 'GE'.  For
            sparse: {'CH', 'LDL'}, default is 'LDL'.

        Raises
        ======

        ValueError
            If the determinant of the matrix is zero.

        See Also
        ========

        inverse_LU
        inverse_GE
        inverse_ADJ

        """
        if not self.is_square:
            raise NonSquareMatrixError()
        if method is not None:
            kwargs['method'] = method
        return self._eval_inverse(**kwargs)

    def inv_mod(self, m):
        r"""
        Returns the inverse of the matrix `K` (mod `m`), if it exists.

        Method to find the matrix inverse of `K` (mod `m`) implemented in this function:

        * Compute `\mathrm{adj}(K) = \mathrm{cof}(K)^t`, the adjoint matrix of `K`.

        * Compute `r = 1/\mathrm{det}(K) \pmod m`.

        * `K^{-1} = r\cdot \mathrm{adj}(K) \pmod m`.

        Examples
        ========

        >>> A = Matrix(2, 2, [1, 2, 3, 4])
        >>> A.inv_mod(5)
        Matrix([
        [3, 1],
        [4, 2]])
        >>> A.inv_mod(3)
        Matrix([
        [1, 1],
        [0, 1]])

        """
        from ..ntheory import totient
        if not self.is_square:
            raise NonSquareMatrixError()
        N = self.cols
        phi = totient(m)
        det_K = self.det()
        if gcd(det_K, m) != 1:
            raise ValueError(f'Matrix is not invertible (mod {m:d})')
        det_inv = pow(int(det_K), int(phi - 1), int(m))
        K_adj = self.cofactorMatrix().transpose()
        K_inv = self.__class__(N, N, [det_inv*K_adj[i, j] % m for i in range(N) for j in range(N)])
        return K_inv

    def transpose(self):
        """Matrix transposition."""
        return self._eval_transpose()

    T = property(transpose, None, None, 'Matrix transposition.')

    def conjugate(self):
        """By-element conjugation."""
        return self._eval_conjugate()

    C = property(conjugate, None, None, 'By-element conjugation.')

    def adjoint(self):
        """Conjugate transpose or Hermitian conjugation."""
        return self.T.C

    @property
    def H(self):
        """Return Hermite conjugate.

        Examples
        ========

        >>> m = Matrix((0, 1 + I, 2, 3))
        >>> m
        Matrix([
        [    0],
        [1 + I],
        [    2],
        [    3]])
        >>> m.H
        Matrix([[0, 1 - I, 2, 3]])

        See Also
        ========

        conjugate: By-element conjugation
        D: Dirac conjugation

        """
        return self.T.C

    @property
    def D(self):
        """Return Dirac conjugate (if self.rows == 4).

        Examples
        ========

        >>> m = Matrix((0, 1 + I, 2, 3))
        >>> m.D
        Matrix([[0, 1 - I, -2, -3]])
        >>> m = (eye(4) + I*eye(4))
        >>> m[0, 3] = 2
        >>> m.D
        Matrix([
        [1 - I,     0,      0,      0],
        [    0, 1 - I,      0,      0],
        [    0,     0, -1 + I,      0],
        [    2,     0,      0, -1 + I]])

        If the matrix does not have 4 rows an AttributeError will be raised
        because this property is only defined for matrices with 4 rows.

        >>> Matrix(eye(2)).D
        Traceback (most recent call last):
        ...
        AttributeError: Matrix has no attribute D.

        See Also
        ========

        conjugate: By-element conjugation
        H: Hermite conjugation

        """
        if self.rows != 4:
            raise AttributeError
        return self.H*mgamma(0)

    def __array__(self):
        from .dense import matrix2numpy
        return matrix2numpy(self)

    def __len__(self):
        """Return the number of elements of self.

        Implemented mainly so bool(Matrix()) == False.

        """
        return int(self.rows*self.cols)

    @property
    def shape(self):
        """The shape (dimensions) of the matrix as the 2-tuple (rows, cols).

        Examples
        ========

        >>> M = zeros(2, 3)
        >>> M.shape
        (2, 3)
        >>> M.rows
        2
        >>> M.cols
        3

        """
        return self.rows, self.cols

    def __sub__(self, a):
        return self + (-a)

    def __mul__(self, other):
        """Return self*other where other is either a scalar or a matrix
        of compatible dimensions.

        Examples
        ========

        >>> A = Matrix([[1, 2, 3], [4, 5, 6]])
        >>> 2*A == A*2 == Matrix([[2, 4, 6], [8, 10, 12]])
        True
        >>> B = Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        >>> A*B
        Matrix([
        [30, 36, 42],
        [66, 81, 96]])
        >>> B*A
        Traceback (most recent call last):
        ...
        ShapeError: Matrices size mismatch.
        >>>

        See Also
        ========

        matrix_multiply_elementwise

        """
        if getattr(other, 'is_Matrix', False):
            A = self
            B = other
            if A.cols != B.rows:
                raise ShapeError('Matrices size mismatch.')
            if A.cols == 0:
                return classof(A, B)._new(A.rows, B.cols, lambda i, j: 0)
            try:
                blst = B.T.tolist()
            except AttributeError:
                return NotImplemented
            alst = A.tolist()
            return classof(A, B)._new(A.rows, B.cols, lambda i, j:
                                      functools.reduce(lambda k, l: k + l,
                                                       [a_ik * b_kj for a_ik, b_kj in zip(alst[i], blst[j])]))
        else:
            return self._new(self.rows, self.cols,
                             [i*other for i in self._mat])

    def __rmul__(self, a):
        if getattr(a, 'is_Matrix', False):
            return self._new(a)*self
        return self._new(self.rows, self.cols, [a*i for i in self._mat])

    def __pow__(self, num):
        from ..functions import binomial
        from . import MutableMatrix, diag, eye

        if not self.is_square:
            raise NonSquareMatrixError()
        if isinstance(num, (Integer, int)):
            n = int(num)
            if n < 0:
                return self.inv()**-n   # A**-2 = (A**-1)**2
            a = eye(self.cols)
            s = self
            while n:
                if n % 2:
                    a *= s
                    n -= 1
                if not n:
                    break
                s *= s
                n //= 2
            return self._new(a)
        elif isinstance(num, Expr):

            def jordan_cell_power(jc, n):
                N = jc.shape[0]
                l = jc[0, 0]
                for i in range(N):
                    for j in range(N - i):
                        bn = binomial(n, i)
                        if isinstance(bn, binomial):
                            bn = bn._eval_expand_func()
                        jc[j, i + j] = l**(n - i)*bn

            jordan_cells, P = self.jordan_cells()
            # Make sure jordan_cells matrices are mutable:
            jordan_cells = [MutableMatrix(j) for j in jordan_cells]
            for j in jordan_cells:
                jordan_cell_power(j, num)
            return self._new(P*diag(*jordan_cells)*P.inv())
        else:
            raise TypeError(
                'Only Diofant expressions or int objects are supported as exponent for matrices')

    def __add__(self, other):
        """Return self + other, raising ShapeError if shapes don't match."""
        if getattr(other, 'is_Matrix', False):
            A = self
            B = other
            if A.shape != B.shape:
                raise ShapeError('Matrix size mismatch.')
            alst = A.tolist()
            blst = B.tolist()
            ret = [Integer(0)]*A.rows
            for i in range(A.shape[0]):
                ret[i] = [j + k for j, k in zip(alst[i], blst[i])]
            rv = classof(A, B)._new(ret)
            if 0 in A.shape:
                rv = rv.reshape(*A.shape)
            return rv
        raise TypeError(f'cannot add matrix and {type(other)}')

    def __radd__(self, other):
        return self + other

    def __truediv__(self, other):
        return self*(Integer(1)/other)

    def __neg__(self):
        return -1*self

    def multiply(self, b):
        """Returns self*b

        See Also
        ========

        dot
        cross
        multiply_elementwise

        """
        return self*b

    def add(self, b):
        """Return self + b."""
        return self + b

    def table(self, printer, rowstart='[', rowend=']', rowsep='\n',
              colsep=', ', align='right'):
        r"""
        String form of Matrix as a table.

        ``printer`` is the printer to use for on the elements (generally
        something like StrPrinter())

        ``rowstart`` is the string used to start each row (by default '[').

        ``rowend`` is the string used to end each row (by default ']').

        ``rowsep`` is the string used to separate rows (by default a newline).

        ``colsep`` is the string used to separate columns (by default ', ').

        ``align`` defines how the elements are aligned. Must be one of 'left',
        'right', or 'center'.  You can also use '<', '>', and '^' to mean the
        same thing, respectively.

        This is used by the string printer for Matrix.

        Examples
        ========

        >>> from diofant.printing.str import StrPrinter
        >>> M = Matrix([[1, 2], [-33, 4]])
        >>> printer = StrPrinter()
        >>> M.table(printer)
        '[  1, 2]\n[-33, 4]'
        >>> print(M.table(printer))
        [  1, 2]
        [-33, 4]
        >>> print(M.table(printer, rowsep=',\n'))
        [  1, 2],
        [-33, 4]
        >>> print(M.table(printer, colsep=' '))
        [  1 2]
        [-33 4]
        >>> print(M.table(printer, align='center'))
        [ 1 , 2]
        [-33, 4]
        >>> print(M.table(printer, rowstart='{', rowend='}'))
        {  1, 2}
        {-33, 4}

        """
        # Handle zero dimensions:
        if self.rows == 0 or self.cols == 0:
            return '[]'
        # Build table of string representations of the elements
        res = []
        # Track per-column max lengths for pretty alignment
        maxlen = [0] * self.cols
        for i in range(self.rows):
            res.append([])
            for j in range(self.cols):
                s = printer._print(self[i, j])
                res[-1].append(s)
                maxlen[j] = max(len(s), maxlen[j])
        # Patch strings together
        align = {
            'left': 'ljust',
            'right': 'rjust',
            'center': 'center',
            '<': 'ljust',
            '>': 'rjust',
            '^': 'center',
        }[align]
        for i, row in enumerate(res):
            for j, elem in enumerate(row):
                row[j] = getattr(elem, align)(maxlen[j])
            res[i] = rowstart + colsep.join(row) + rowend
        return rowsep.join(res)

    def _format_str(self, printer):
        if self.rows == 0 or self.cols == 0:
            return f'Matrix({self.rows}, {self.cols}, [])'
        if self.rows == 1:
            return 'Matrix([%s])' % self.table(printer, rowsep=',\n')  # noqa: SFS101
        return 'Matrix([\n%s])' % self.table(printer, rowsep=',\n')  # noqa: SFS101

    def _repr_pretty_(self, p, cycle):
        from ..printing import pretty
        p.text(pretty(self))

    def _repr_latex_(self):
        from ..printing import latex
        return latex(self, mode='equation')

    def cholesky(self):
        """Returns the Cholesky decomposition L of a matrix A
        such that L * L.T = A

        A must be a square, symmetric, positive-definite
        and non-singular matrix.

        Examples
        ========

        >>> A = Matrix(((25, 15, -5), (15, 18, 0), (-5, 0, 11)))
        >>> A.cholesky()
        Matrix([
        [ 5, 0, 0],
        [ 3, 3, 0],
        [-1, 1, 3]])
        >>> A.cholesky() * A.cholesky().T
        Matrix([
        [25, 15, -5],
        [15, 18,  0],
        [-5,  0, 11]])

        See Also
        ========

        LDLdecomposition
        LUdecomposition
        QRdecomposition

        """
        if not self.is_square:
            raise NonSquareMatrixError('Matrix must be square.')
        if not self.is_symmetric():
            raise ValueError('Matrix must be symmetric.')
        return self._cholesky()

    def LDLdecomposition(self):
        """Returns the LDL Decomposition (L, D) of matrix A,
        such that L * D * L.T == A
        This method eliminates the use of square root.
        Further this ensures that all the diagonal entries of L are 1.
        A must be a square, symmetric, positive-definite
        and non-singular matrix.

        Examples
        ========

        >>> A = Matrix(((25, 15, -5), (15, 18, 0), (-5, 0, 11)))
        >>> L, D = A.LDLdecomposition()
        >>> L
        Matrix([
        [   1,   0, 0],
        [ 3/5,   1, 0],
        [-1/5, 1/3, 1]])
        >>> D
        Matrix([
        [25, 0, 0],
        [ 0, 9, 0],
        [ 0, 0, 9]])
        >>> L * D * L.T * A.inv() == eye(A.rows)
        True

        See Also
        ========

        cholesky
        LUdecomposition
        QRdecomposition

        """
        if not self.is_square:
            raise NonSquareMatrixError('Matrix must be square.')
        if not self.is_symmetric():
            raise ValueError('Matrix must be symmetric.')
        return self._LDLdecomposition()

    def lower_triangular_solve(self, rhs):
        """Solves Ax = B, where A is a lower triangular matrix.

        See Also
        ========

        upper_triangular_solve
        cholesky_solve
        diagonal_solve
        LDLsolve
        LUsolve
        QRsolve
        pinv_solve

        """
        if not self.is_square:
            raise NonSquareMatrixError('Matrix must be square.')
        if rhs.rows != self.rows:
            raise ShapeError('Matrices size mismatch.')
        if not self.is_lower:
            raise ValueError('Matrix must be lower triangular.')
        return self._lower_triangular_solve(rhs)

    def upper_triangular_solve(self, rhs):
        """Solves Ax = B, where A is an upper triangular matrix.

        See Also
        ========

        lower_triangular_solve
        cholesky_solve
        diagonal_solve
        LDLsolve
        LUsolve
        QRsolve
        pinv_solve

        """
        if not self.is_square:
            raise NonSquareMatrixError('Matrix must be square.')
        if rhs.rows != self.rows:
            raise TypeError('Matrix size mismatch.')
        if not self.is_upper:
            raise TypeError('Matrix is not upper triangular.')
        return self._upper_triangular_solve(rhs)

    def cholesky_solve(self, rhs):
        """Solves Ax = B using Cholesky decomposition,
        for a general square non-singular matrix.
        For a non-square matrix with rows > cols,
        the least squares solution is returned.

        See Also
        ========

        lower_triangular_solve
        upper_triangular_solve
        diagonal_solve
        LDLsolve
        LUsolve
        QRsolve
        pinv_solve

        """
        if self.is_symmetric():
            L = self._cholesky()
        elif self.rows >= self.cols:
            L = (self.T*self)._cholesky()
            rhs = self.T*rhs
        else:
            raise NotImplementedError('Under-determined System.')
        Y = L._lower_triangular_solve(rhs)
        return (L.T)._upper_triangular_solve(Y)

    def diagonal_solve(self, rhs):
        """Solves Ax = B efficiently, where A is a diagonal Matrix,
        with non-zero diagonal entries.

        Examples
        ========

        >>> A = eye(2)*2
        >>> B = Matrix([[1, 2], [3, 4]])
        >>> A.diagonal_solve(B) == B/2
        True

        See Also
        ========

        lower_triangular_solve
        upper_triangular_solve
        cholesky_solve
        LDLsolve
        LUsolve
        QRsolve
        pinv_solve

        """
        if not self.is_diagonal():
            raise TypeError('Matrix should be diagonal')
        if rhs.rows != self.rows:
            raise TypeError('Size mis-match')
        return self._diagonal_solve(rhs)

    def LDLsolve(self, rhs):
        """Solves Ax = B using LDL decomposition,
        for a general square and non-singular matrix.

        For a non-square matrix with rows > cols,
        the least squares solution is returned.

        Examples
        ========

        >>> A = eye(2)*2
        >>> B = Matrix([[1, 2], [3, 4]])
        >>> A.LDLsolve(B) == B/2
        True

        See Also
        ========

        LDLdecomposition
        lower_triangular_solve
        upper_triangular_solve
        cholesky_solve
        diagonal_solve
        LUsolve
        QRsolve
        pinv_solve

        """
        if self.is_symmetric():
            L, D = self.LDLdecomposition()
        elif self.rows >= self.cols:
            L, D = (self.T*self).LDLdecomposition()
            rhs = self.T*rhs
        else:
            raise NotImplementedError('Under-determined System.')
        Y = L._lower_triangular_solve(rhs)
        Z = D._diagonal_solve(Y)
        return (L.T)._upper_triangular_solve(Z)

    def solve_least_squares(self, rhs, method='CH'):
        """Return the least-square fit to the data.

        By default the cholesky_solve routine is used (method='CH'); other
        methods of matrix inversion can be used.

        Examples
        ========

        >>> A = Matrix([1, 2, 3])
        >>> B = Matrix([2, 3, 4])
        >>> S = Matrix(A.row_join(B))
        >>> S
        Matrix([
        [1, 2],
        [2, 3],
        [3, 4]])

        If each line of S represent coefficients of Ax + By
        and x and y are [2, 3] then S*xy is:

        >>> r = S*Matrix([2, 3])
        >>> r
        Matrix([
        [ 8],
        [13],
        [18]])

        But let's add 1 to the middle value and then solve for the
        least-squares value of xy:

        >>> xy = S.solve_least_squares(Matrix([8, 14, 18]))
        >>> xy
        Matrix([
        [ 5/3],
        [10/3]])

        The error is given by S*xy - r:

        >>> S*xy - r
        Matrix([
        [1/3],
        [1/3],
        [1/3]])
        >>> _.norm().evalf(2)
        0.58

        If a different xy is used, the norm will be higher:

        >>> xy += ones(2, 1)/10
        >>> (S*xy - r).norm().evalf(2)
        1.5

        See Also
        ========

        inv

        """
        if method == 'CH':
            return self.cholesky_solve(rhs)
        t = self.T
        return (t*self).inv(method=method)*t*rhs

    def extract(self, rowsList, colsList):
        """Return a submatrix by specifying a list of rows and columns.
        Negative indices can be given. All indices must be in the range
        -n <= i < n where n is the number of rows or columns.

        Examples
        ========

        >>> m = Matrix(4, 3, range(12))
        >>> m
        Matrix([
        [0,  1,  2],
        [3,  4,  5],
        [6,  7,  8],
        [9, 10, 11]])
        >>> m.extract([0, 1, 3], [0, 1])
        Matrix([
        [0,  1],
        [3,  4],
        [9, 10]])

        Rows or columns can be repeated:

        >>> m.extract([0, 0, 1], [-1])
        Matrix([
        [2],
        [2],
        [5]])

        Every other row can be taken by using range to provide the indices:

        >>> m.extract(range(0, m.rows, 2), [-1])
        Matrix([
        [2],
        [8]])

        """
        cols = self.cols
        flat_list = self._mat
        rowsList = [a2idx(k, self.rows) for k in rowsList]
        colsList = [a2idx(k, self.cols) for k in colsList]
        return self._new(len(rowsList), len(colsList),
                         lambda i, j: flat_list[rowsList[i]*cols + colsList[j]])

    def key2bounds(self, keys):
        """Converts a key with potentially mixed types of keys (integer and slice)
        into a tuple of ranges and raises an error if any index is out of self's
        range.

        See Also
        ========

        key2ij

        """
        islice, jslice = (isinstance(k, slice) for k in keys)
        if islice:
            assert self.rows
            rlo, rhi = keys[0].indices(self.rows)[:2]
        else:
            rlo = a2idx(keys[0], self.rows)
            rhi = rlo + 1
        if jslice:
            assert self.cols
            clo, chi = keys[1].indices(self.cols)[:2]
        else:
            clo = a2idx(keys[1], self.cols)
            chi = clo + 1
        return rlo, rhi, clo, chi

    def key2ij(self, key):
        """Converts key into canonical form, converting integers or indexable
        items into valid integers for self's range or returning slices
        unchanged.

        See Also
        ========

        key2bounds

        """
        if is_sequence(key):
            if not len(key) == 2:
                raise TypeError('key must be a sequence of length 2')
            return [a2idx(i, n) if not isinstance(i, slice) else i
                    for i, n in zip(key, self.shape)]
        elif isinstance(key, slice):
            return key.indices(len(self))[:2]
        else:
            return divmod(a2idx(key, len(self)), self.cols)

    def evalf(self, dps=15, **options):
        """Apply evalf() to each element of self."""
        return self.applyfunc(lambda i: i.evalf(dps, **options))

    def atoms(self, *types):
        """Returns the atoms that form the current object.

        Examples
        ========

        >>> Matrix([[x]])
        Matrix([[x]])
        >>> _.atoms()
        {x}

        """
        if types:
            types = tuple(t if isinstance(t, type) else type(t) for t in types)
        else:
            types = Atom,
        result = set()
        for i in self:
            result.update(i.atoms(*types))
        return result

    @property
    def free_symbols(self):
        """Returns the free symbols within the matrix.

        Examples
        ========

        >>> Matrix([[x], [1]]).free_symbols
        {x}

        """
        return set().union(*[i.free_symbols for i in self])

    def subs(self, *args, **kwargs):  # should mirror core.basic.subs
        """Return a new matrix with subs applied to each entry.

        Examples
        ========

        >>> SparseMatrix(1, 1, [x])
        Matrix([[x]])
        >>> _.subs({x: y})
        Matrix([[y]])
        >>> Matrix(_).subs({y: x})
        Matrix([[x]])

        """
        return self.applyfunc(lambda x: x.subs(*args, **kwargs))

    def xreplace(self, rule):  # should mirror core.basic.xreplace
        """Return a new matrix with xreplace applied to each entry.

        Examples
        ========

        >>> SparseMatrix(1, 1, [x])
        Matrix([[x]])
        >>> _.xreplace({x: y})
        Matrix([[y]])
        >>> Matrix(_).xreplace({y: x})
        Matrix([[x]])

        """
        return self.applyfunc(lambda x: x.xreplace(rule))

    def expand(self, deep=True, modulus=None, power_base=True, power_exp=True,
               mul=True, log=True, multinomial=True, basic=True, **hints):
        """Apply core.function.expand to each entry of the matrix.

        Examples
        ========

        >>> Matrix(1, 1, [x*(x+1)])
        Matrix([[x*(x + 1)]])
        >>> _.expand()
        Matrix([[x**2 + x]])

        """
        return self.applyfunc(lambda x: x.expand(
                              deep, modulus, power_base, power_exp, mul, log, multinomial, basic,
                              **hints))

    def simplify(self, ratio=1.7, measure=count_ops):
        """Apply simplify to each element of the matrix.

        Examples
        ========

        >>> SparseMatrix(1, 1, [x*sin(y)**2 + x*cos(y)**2])
        Matrix([[x*sin(y)**2 + x*cos(y)**2]])
        >>> _.simplify()
        Matrix([[x]])

        """
        return self.applyfunc(lambda x: x.simplify(ratio, measure))
    _eval_simplify = simplify

    def print_nonzero(self, symb='X'):
        """Shows location of non-zero entries for fast shape lookup.

        Examples
        ========

        >>> m = Matrix(2, 3, lambda i, j: i*3+j)
        >>> m
        Matrix([
        [0, 1, 2],
        [3, 4, 5]])
        >>> m.print_nonzero()
        [ XX]
        [XXX]
        >>> m = eye(4)
        >>> m.print_nonzero('x')
        [x   ]
        [ x  ]
        [  x ]
        [   x]

        """
        s = []
        for i in range(self.rows):
            line = []
            for j in range(self.cols):
                if self[i, j] == 0:
                    line.append(' ')
                else:
                    line.append(str(symb))
            s.append(f"[{''.join(line)}]")
        print('\n'.join(s))

    def LUsolve(self, rhs, iszerofunc=_iszero):
        """Solve the linear system Ax = rhs for x where A = self.

        This is for symbolic matrices, for real or complex ones use
        mpmath.lu_solve or mpmath.qr_solve.

        See Also
        ========

        lower_triangular_solve
        upper_triangular_solve
        cholesky_solve
        diagonal_solve
        LDLsolve
        QRsolve
        pinv_solve
        LUdecomposition

        """
        if rhs.rows != self.rows:
            raise ShapeError('`self` and `rhs` must have the same number of rows.')

        A, perm = self.LUdecomposition_Simple(iszerofunc=_iszero)
        n = self.rows
        b = rhs.permuteFwd(perm).as_mutable()
        # forward substitution, all diag entries are scaled to 1
        for i in range(n):
            for j in range(i):
                scale = A[i, j]
                b.zip_row_op(i, j, lambda x, y: x - y*scale)
        # backward substitution
        for i in range(n - 1, -1, -1):
            for j in range(i + 1, n):
                scale = A[i, j]
                b.zip_row_op(i, j, lambda x, y: x - y*scale)
            scale = A[i, i]
            b.row_op(i, lambda x, _: x/scale)
        return rhs.__class__(b)

    def LUdecomposition(self, iszerofunc=_iszero):
        """Returns the decomposition LU and the row swaps p.

        Examples
        ========

        >>> a = Matrix([[4, 3], [6, 3]])
        >>> L, U, _ = a.LUdecomposition()
        >>> L
        Matrix([
        [  1, 0],
        [3/2, 1]])
        >>> U
        Matrix([
        [4,    3],
        [0, -3/2]])

        See Also
        ========

        cholesky
        LDLdecomposition
        QRdecomposition
        LUdecomposition_Simple
        LUdecompositionFF
        LUsolve

        """
        combined, p = self.LUdecomposition_Simple(iszerofunc=_iszero)
        L = self.zeros(self.rows)
        U = self.zeros(self.rows)
        for i in range(self.rows):
            for j in range(self.rows):
                if i > j:
                    L[i, j] = combined[i, j]
                else:
                    if i == j:
                        L[i, i] = 1
                    U[i, j] = combined[i, j]
        return L, U, p

    def LUdecomposition_Simple(self, iszerofunc=_iszero):
        """Returns A comprised of L, U (L's diag entries are 1) and
        p which is the list of the row swaps (in order).

        See Also
        ========

        LUdecomposition
        LUdecompositionFF
        LUsolve

        """
        if not self.is_square:
            raise NonSquareMatrixError('A Matrix must be square to apply LUdecomposition_Simple().')
        n = self.rows
        A = self.as_mutable()
        p = []
        # factorization
        for j in range(n):
            for i in range(j):
                for k in range(i):
                    A[i, j] = A[i, j] - A[i, k]*A[k, j]
            pivot = -1
            for i in range(j, n):
                for k in range(j):
                    A[i, j] = A[i, j] - A[i, k]*A[k, j]
                # find the first non-zero pivot, includes any expression
                if pivot == -1 and not iszerofunc(A[i, j]):
                    pivot = i
            if pivot < 0:
                # this result is based on iszerofunc's analysis of the possible pivots, so even though
                # the element may not be strictly zero, the supplied iszerofunc's evaluation gave True
                raise ValueError('No nonzero pivot found; inversion failed.')
            if pivot != j:  # row must be swapped
                A.row_swap(pivot, j)
                p.append([pivot, j])
            scale = 1 / A[j, j]
            for i in range(j + 1, n):
                A[i, j] = A[i, j]*scale
        return A, p

    def LUdecompositionFF(self):
        """Compute a fraction-free LU decomposition.

        Returns 4 matrices P, L, D, U such that PA = L D**-1 U.
        If the elements of the matrix belong to some integral domain I, then all
        elements of L, D and U are guaranteed to belong to I.

        **Reference**
            - W. Zhou & D.J. Jeffrey, "Fraction-free matrix factors: new forms
              for LU and QR factors". Frontiers in Computer Science in China,
              Vol 2, no. 1, pp. 67-80, 2008.

        See Also
        ========

        LUdecomposition
        LUdecomposition_Simple
        LUsolve

        """
        from . import SparseMatrix
        zeros = SparseMatrix.zeros
        eye = SparseMatrix.eye

        n, m = self.rows, self.cols
        U, L, P = self.as_mutable(), eye(n), eye(n)
        DD = zeros(n, n)
        oldpivot = 1

        for k in range(n - 1):
            if U[k, k] == 0:
                for kpivot in range(k + 1, n):
                    if U[kpivot, k]:
                        break
                else:
                    raise ValueError('Matrix is not full rank')
                U[k, k:], U[kpivot, k:] = U[kpivot, k:], U[k, k:]
                L[k, :k], L[kpivot, :k] = L[kpivot, :k], L[k, :k]
                P[k, :], P[kpivot, :] = P[kpivot, :], P[k, :]
            L[k, k] = Ukk = U[k, k]
            DD[k, k] = oldpivot*Ukk
            for i in range(k + 1, n):
                L[i, k] = Uik = U[i, k]
                for j in range(k + 1, m):
                    U[i, j] = (Ukk*U[i, j] - U[k, j]*Uik) / oldpivot
                U[i, k] = 0
            oldpivot = Ukk
        DD[n - 1, n - 1] = oldpivot
        return P, L, DD, U

    def cofactorMatrix(self, method='berkowitz'):
        """Return a matrix containing the cofactor of each element.

        See Also
        ========

        cofactor
        minorEntry
        minorMatrix
        adjugate

        """
        out = self._new(self.rows, self.cols, lambda i, j:
                        self.cofactor(i, j, method))
        return out

    def minorEntry(self, i, j, method='berkowitz'):
        """Calculate the minor of an element.

        See Also
        ========

        minorMatrix
        cofactor
        cofactorMatrix

        """
        if not 0 <= i < self.rows or not 0 <= j < self.cols:
            raise ValueError('`i` and `j` must satisfy 0 <= i < `self.rows` ' +
                             f'({self.rows:d})' + f'and 0 <= j < `self.cols` ({self.cols:d}).')
        return self.minorMatrix(i, j).det(method)

    def minorMatrix(self, i, j):
        """Creates the minor matrix of a given element.

        See Also
        ========

        minorEntry
        cofactor
        cofactorMatrix

        """
        if not 0 <= i < self.rows or not 0 <= j < self.cols:
            raise ValueError('`i` and `j` must satisfy 0 <= i < `self.rows` ' +
                             f'({self.rows:d})' + f'and 0 <= j < `self.cols` ({self.cols:d}).')
        M = self.as_mutable()
        del M[i, :]
        del M[:, j]
        return self._new(M)

    def cofactor(self, i, j, method='berkowitz'):
        """Calculate the cofactor of an element.

        See Also
        ========

        cofactorMatrix
        minorEntry
        minorMatrix

        """
        if (i + j) % 2 == 0:
            return self.minorEntry(i, j, method)
        else:
            return -1*self.minorEntry(i, j, method)

    def jacobian(self, X):
        """Calculates the Jacobian matrix (derivative of a vectorial function).

        Parameters
        ==========

        self : vector of expressions representing functions f_i(x_1, ..., x_n).
        X : set of x_i's in order, it can be a list or a Matrix

        Both self and X can be a row or a column matrix in any order
        (i.e., jacobian() should always work).

        Examples
        ========

        >>> from diofant.abc import phi, rho
        >>> X = Matrix([rho*cos(phi), rho*sin(phi), rho**2])
        >>> Y = Matrix([rho, phi])
        >>> X.jacobian(Y)
        Matrix([
        [cos(phi), -rho*sin(phi)],
        [sin(phi),  rho*cos(phi)],
        [   2*rho,             0]])
        >>> X = Matrix([rho*cos(phi), rho*sin(phi)])
        >>> X.jacobian(Y)
        Matrix([
        [cos(phi), -rho*sin(phi)],
        [sin(phi),  rho*cos(phi)]])

        See Also
        ========

        diofant.matrices.dense.hessian
        diofant.matrices.dense.wronskian

        """
        if not isinstance(X, MatrixBase):
            X = self._new(X)
        # Both X and self can be a row or a column matrix, so we need to make
        # sure all valid combinations work, but everything else fails:
        if self.shape[0] == 1:
            m = self.shape[1]
        elif self.shape[1] == 1:
            m = self.shape[0]
        else:
            raise TypeError('self must be a row or a column matrix')
        if X.shape[0] == 1:
            n = X.shape[1]
        elif X.shape[1] == 1:
            n = X.shape[0]
        else:
            raise TypeError('X must be a row or a column matrix')

        # m is the number of functions and n is the number of variables
        # computing the Jacobian is now easy:
        return self._new(m, n, lambda j, i: self[j].diff(X[i]))

    def QRdecomposition(self):
        """Return Q, R where A = Q*R, Q is orthogonal and R is upper triangular.

        Examples
        ========

        This is the example from wikipedia:

        >>> A = Matrix([[12, -51, 4], [6, 167, -68], [-4, 24, -41]])
        >>> Q, R = A.QRdecomposition()
        >>> Q
        Matrix([
        [ 6/7, -69/175, -58/175],
        [ 3/7, 158/175,   6/175],
        [-2/7,    6/35,  -33/35]])
        >>> R
        Matrix([
        [14,  21, -14],
        [ 0, 175, -70],
        [ 0,   0,  35]])
        >>> A == Q*R
        True

        QR factorization of an identity matrix:

        >>> A = Matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        >>> Q, R = A.QRdecomposition()
        >>> Q
        Matrix([
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1]])
        >>> R
        Matrix([
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1]])

        See Also
        ========

        cholesky
        LDLdecomposition
        LUdecomposition
        QRsolve

        """
        cls = self.__class__
        mat = self.as_mutable()

        if not mat.rows >= mat.cols:
            raise MatrixError(
                'The number of rows must be greater than columns')
        n = mat.rows
        m = mat.cols
        rank = n
        row_reduced = mat.rref()[0]
        for i in range(row_reduced.rows):
            if row_reduced[i, :].norm() == 0:
                rank -= 1
        if not rank == mat.cols:
            raise MatrixError('The rank of the matrix must match the columns')
        Q, R = mat.zeros(n, m), mat.zeros(m)
        for j in range(m):      # for each column vector
            tmp = mat[:, j]     # take original v
            for i in range(j):
                # subtract the project of mat on new vector
                tmp -= Q[:, i]*mat[:, j].dot(Q[:, i])
                tmp.expand()
            # normalize it
            R[j, j] = tmp.norm()
            Q[:, j] = tmp / R[j, j]
            if Q[:, j].norm() != 1:
                raise NotImplementedError("Couldn't normalize the "
                                          f'vector {j:d}.')
            for i in range(j):
                R[i, j] = Q[:, i].dot(mat[:, j])
        return cls(Q), cls(R)

    def QRsolve(self, b):
        """Solve the linear system 'Ax = b'.

        'self' is the matrix 'A', the method argument is the vector
        'b'.  The method returns the solution vector 'x'.  If 'b' is a
        matrix, the system is solved for each column of 'b' and the
        return value is a matrix of the same shape as 'b'.

        This method is slower (approximately by a factor of 2) but
        more stable for floating-point arithmetic than the LUsolve method.
        However, LUsolve usually uses an exact arithmetic, so you don't need
        to use QRsolve.

        This is mainly for educational purposes and symbolic matrices, for real
        (or complex) matrices use mpmath.qr_solve.

        See Also
        ========

        lower_triangular_solve
        upper_triangular_solve
        cholesky_solve
        diagonal_solve
        LDLsolve
        LUsolve
        pinv_solve
        QRdecomposition

        """
        Q, R = self.as_mutable().QRdecomposition()
        y = Q.T*b

        # back substitution to solve R*x = y:
        # We build up the result "backwards" in the vector 'x' and reverse it
        # only in the end.
        x = []
        n = R.rows
        for j in range(n - 1, -1, -1):
            tmp = y[j, :]
            for k in range(j + 1, n):
                tmp -= R[j, k]*x[n - 1 - k]
            x.append(tmp / R[j, j])
        return self._new([row._mat for row in reversed(x)])

    def cross(self, b):
        """Return the cross product of `self` and `b` relaxing the condition
        of compatible dimensions: if each has 3 elements, a matrix of the
        same type and shape as `self` will be returned. If `b` has the same
        shape as `self` then common identities for the cross product (like
        `a x b = - b x a`) will hold.

        See Also
        ========

        dot
        multiply
        multiply_elementwise

        """
        if not is_sequence(b):
            raise TypeError(f'`b` must be an ordered iterable or Matrix, not {type(b)}.')
        if not self.rows * self.cols == b.rows * b.cols == 3:
            raise ShapeError('Dimensions incorrect for cross product.')
        return self._new(self.rows, self.cols, ((self[1]*b[2] - self[2]*b[1]),
                                                (self[2]*b[0] - self[0]*b[2]),
                                                (self[0]*b[1] - self[1]*b[0])))

    def dot(self, b):
        """Return the dot product of Matrix self and b relaxing the condition
        of compatible dimensions: if either the number of rows or columns are
        the same as the length of b then the dot product is returned. If self
        is a row or column vector, a scalar is returned. Otherwise, a list
        of results is returned (and in that case the number of columns in self
        must match the length of b).

        Examples
        ========

        >>> M = Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        >>> v = [1, 1, 1]
        >>> M[0, :].dot(v)
        6
        >>> M[:, 0].dot(v)
        12
        >>> M.dot(v)
        [6, 15, 24]

        See Also
        ========

        cross
        multiply
        multiply_elementwise

        """
        from . import Matrix

        if not isinstance(b, MatrixBase):
            if is_sequence(b):
                if len(b) != self.cols and len(b) != self.rows:
                    raise ShapeError('Dimensions incorrect for dot product.')
                return self.dot(Matrix(b))
            else:
                raise TypeError(f'`b` must be an ordered iterable or Matrix, not {type(b)}')

        mat = self
        if mat.cols == b.rows:
            if b.cols != 1:
                mat = mat.T
                b = b.T
            prod = flatten((mat*b).tolist())
            if len(prod) == 1:
                return prod[0]
            return prod
        if mat.cols == b.cols:
            return mat.dot(b.T)
        elif mat.rows == b.rows:
            return mat.T.dot(b)
        else:
            raise ShapeError('Dimensions incorrect for dot product.')

    def multiply_elementwise(self, b):
        """Return the Hadamard product (elementwise product) of A and B

        Examples
        ========

        >>> A = Matrix([[0, 1, 2], [3, 4, 5]])
        >>> B = Matrix([[1, 10, 100], [100, 10, 1]])
        >>> A.multiply_elementwise(B)
        Matrix([
        [  0, 10, 200],
        [300, 40,   5]])

        See Also
        ========

        cross
        dot
        multiply

        """
        from . import matrix_multiply_elementwise

        return matrix_multiply_elementwise(self, b)

    def values(self):
        """Return non-zero values of self."""
        return [i for i in flatten(self.tolist()) if not i.is_zero]

    def norm(self, ord=None):
        """Return the Norm of a Matrix or Vector.
        In the simplest case this is the geometric size of the vector
        Other norms can be specified by the ord parameter


        =====  ============================  ==========================
        ord    norm for matrices             norm for vectors
        =====  ============================  ==========================
        None   Frobenius norm                2-norm
        'fro'  Frobenius norm                - does not exist
        inf    --                            max(abs(x))
        -inf   --                            min(abs(x))
        1      --                            as below
        -1     --                            as below
        2      2-norm (largest sing. value)  as below
        -2     smallest singular value       as below
        other  - does not exist              sum(abs(x)**ord)**(1./ord)
        =====  ============================  ==========================

        Examples
        ========

        >>> x = Symbol('x', real=True)
        >>> v = Matrix([cos(x), sin(x)])
        >>> trigsimp(v.norm())
        1
        >>> v.norm(10)
        (sin(x)**10 + cos(x)**10)**(1/10)
        >>> A = Matrix([[1, 1], [1, 1]])
        >>> A.norm(2)  # Spectral norm (max of |Ax|/|x| under 2-vector-norm)
        2
        >>> A.norm(-2)  # Inverse spectral norm (smallest singular value)
        0
        >>> A.norm()  # Frobenius Norm
        2
        >>> Matrix([1, -2]).norm(oo)
        2
        >>> Matrix([-1, 2]).norm(-oo)
        1

        See Also
        ========

        normalized

        """
        # Row or Column Vector Norms
        vals = list(self.values()) or [0]
        if self.rows == 1 or self.cols == 1:
            if ord == 2 or ord is None:  # Common case sqrt(<x, x>)
                return sqrt(Add(*(abs(i)**2 for i in vals)))

            elif ord == 1:  # sum(abs(x))
                return Add(*(abs(i) for i in vals))

            elif ord == oo:  # max(abs(x))
                return Max(*[abs(i) for i in vals])

            elif ord == -oo:  # min(abs(x))
                return Min(*[abs(i) for i in vals])

            # Otherwise generalize the 2-norm, Sum(x_i**ord)**(1/ord)
            # Note that while useful this is not mathematically a norm
            return Pow(Add(*(abs(i)**ord for i in vals)), Integer(1) / ord)

        # Matrix Norms
        else:
            if ord == 2:  # Spectral Norm
                # Maximum singular value
                return Max(*self.singular_values())

            elif ord == -2:
                # Minimum singular value
                return Min(*self.singular_values())

            elif (ord is None or isinstance(ord, str) and ord.lower() in
                    ['f', 'fro', 'frobenius', 'vector']):
                # Reshape as vector and send back to norm function
                return self.vec().norm(ord=2)

            else:
                raise NotImplementedError('Matrix Norms under development')

    def normalized(self):
        """Return the normalized version of ``self``.

        See Also
        ========

        norm

        """
        if self.rows != 1 and self.cols != 1:
            raise ShapeError('A Matrix must be a vector to normalize.')
        norm = self.norm()
        out = self.applyfunc(lambda i: i / norm)
        return out

    def project(self, v):
        """Return the projection of ``self`` onto the line containing ``v``.

        Examples
        ========

        >>> V = Matrix([sqrt(3)/2, Rational(1, 2)])
        >>> x = Matrix([[1, 0]])
        >>> V.project(x)
        Matrix([[sqrt(3)/2, 0]])
        >>> V.project(-x)
        Matrix([[sqrt(3)/2, 0]])

        """
        return v*(self.dot(v) / v.dot(v))

    def permuteBkwd(self, perm):
        """Permute the rows of the matrix with the given permutation in reverse.

        Examples
        ========

        >>> M = eye(3)
        >>> M.permuteBkwd([[0, 1], [0, 2]])
        Matrix([
        [0, 1, 0],
        [0, 0, 1],
        [1, 0, 0]])

        See Also
        ========

        permuteFwd

        """
        copy = self.copy()
        for i in range(len(perm) - 1, -1, -1):
            copy.row_swap(perm[i][0], perm[i][1])
        return copy

    def permuteFwd(self, perm):
        """Permute the rows of the matrix with the given permutation.

        Examples
        ========

        >>> M = eye(3)
        >>> M.permuteFwd([[0, 1], [0, 2]])
        Matrix([
        [0, 0, 1],
        [1, 0, 0],
        [0, 1, 0]])

        See Also
        ========

        permuteBkwd

        """
        copy = self.copy()
        for p in perm:
            copy.row_swap(p[0], p[1])
        return copy

    def exp(self):
        """Return the exponentiation of a square matrix."""
        if not self.is_square:
            raise NonSquareMatrixError('Exponentiation is valid only '
                                       'for square matrices')
        cells, P = self.jordan_cells()

        def _jblock_exponential(b):
            # This function computes the matrix exponential for one single Jordan block
            nr = b.rows
            l = b[0, 0]
            if nr == 1:
                res = exp(l)
            else:
                from . import eye

                # extract the diagonal part
                d = b[0, 0]*eye(nr)
                # and the nilpotent part
                n = b-d
                # compute its exponential
                nex = eye(nr)
                for i in range(1, nr):
                    nex = nex+n**i/factorial(i)
                # combine the two parts
                res = exp(b[0, 0])*nex
            return res

        blocks = list(map(_jblock_exponential, cells))
        from . import diag
        eJ = diag(* blocks)
        # n = self.rows
        ret = P*eJ*P.inv()
        return type(self)(ret)

    @property
    def is_square(self):
        """Checks if a matrix is square.

        A matrix is square if the number of rows equals the number of columns.
        The empty matrix is square by definition, since the number of rows and
        the number of columns are both zero.

        Examples
        ========

        >>> a = Matrix([[1, 2, 3], [4, 5, 6]])
        >>> b = Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        >>> c = Matrix([])
        >>> a.is_square
        False
        >>> b.is_square
        True
        >>> c.is_square
        True

        """
        return self.rows == self.cols

    @property
    def is_zero(self):
        """Checks if a matrix is a zero matrix.

        A matrix is zero if every element is zero.  A matrix need not be square
        to be considered zero.  The empty matrix is zero by the principle of
        vacuous truth.  For a matrix that may or may not be zero (e.g.
        contains a symbol), this will be None

        Examples
        ========

        >>> a = Matrix([[0, 0], [0, 0]])
        >>> b = zeros(3, 4)
        >>> c = Matrix([[0, 1], [0, 0]])
        >>> d = Matrix([])
        >>> e = Matrix([[x, 0], [0, 0]])
        >>> a.is_zero
        True
        >>> b.is_zero
        True
        >>> c.is_zero
        False
        >>> d.is_zero
        True
        >>> e.is_zero

        """
        if any(i.is_nonzero for i in self):
            return False
        if not any(i.is_zero is None for i in self):
            return True

    def is_nilpotent(self):
        """Checks if a matrix is nilpotent.

        A matrix B is nilpotent if for some integer k, B**k is
        a zero matrix.

        Examples
        ========

        >>> a = Matrix([[0, 0, 0], [1, 0, 0], [1, 1, 0]])
        >>> a.is_nilpotent()
        True

        >>> a = Matrix([[1, 0, 1], [1, 0, 0], [1, 1, 0]])
        >>> a.is_nilpotent()
        False

        """
        if not self:
            return True
        if not self.is_square:
            raise NonSquareMatrixError(
                'Nilpotency is valid only for square matrices')
        x = Dummy('x')
        if self.charpoly(x).args[0] == x**self.rows:
            return True
        return False

    @property
    def is_upper(self):
        """Check if matrix is an upper triangular matrix. True can be returned
        even if the matrix is not square.

        Examples
        ========

        >>> m = Matrix(2, 2, [1, 0, 0, 1])
        >>> m
        Matrix([
        [1, 0],
        [0, 1]])
        >>> m.is_upper
        True

        >>> m = Matrix(4, 3, [5, 1, 9, 0, 4, 6, 0, 0, 5, 0, 0, 0])
        >>> m
        Matrix([
        [5, 1, 9],
        [0, 4, 6],
        [0, 0, 5],
        [0, 0, 0]])
        >>> m.is_upper
        True

        >>> m = Matrix(2, 3, [4, 2, 5, 6, 1, 1])
        >>> m
        Matrix([
        [4, 2, 5],
        [6, 1, 1]])
        >>> m.is_upper
        False

        See Also
        ========

        is_lower
        is_diagonal
        is_upper_hessenberg

        """
        return all(self[i, j].is_zero
                   for i in range(1, self.rows)
                   for j in range(min(i, self.cols)))

    @property
    def is_lower(self):
        """Check if matrix is a lower triangular matrix. True can be returned
        even if the matrix is not square.

        Examples
        ========

        >>> m = Matrix(2, 2, [1, 0, 0, 1])
        >>> m
        Matrix([
        [1, 0],
        [0, 1]])
        >>> m.is_lower
        True

        >>> m = Matrix(4, 3, [0, 0, 0, 2, 0, 0, 1, 4, 0, 6, 6, 5])
        >>> m
        Matrix([
        [0, 0, 0],
        [2, 0, 0],
        [1, 4, 0],
        [6, 6, 5]])
        >>> m.is_lower
        True

        >>> m = Matrix(2, 2, [x**2 + y, y**2 + x, 0, x + y])
        >>> m
        Matrix([
        [x**2 + y, x + y**2],
        [       0,    x + y]])
        >>> m.is_lower
        False

        See Also
        ========

        is_upper
        is_diagonal
        is_lower_hessenberg

        """
        return all(self[i, j].is_zero
                   for i in range(self.rows)
                   for j in range(i + 1, self.cols))

    @property
    def is_hermitian(self):
        """Checks if the matrix is Hermitian.

        In a Hermitian matrix element i,j is the complex conjugate of
        element j,i.

        Examples
        ========

        >>> a = Matrix([[1, I], [-I, 1]])
        >>> a
        Matrix([
        [ 1, I],
        [-I, 1]])
        >>> a.is_hermitian
        True
        >>> a[0, 0] = 2*I
        >>> a.is_hermitian
        False
        >>> a[0, 0] = x
        >>> a.is_hermitian
        >>> a[0, 1] = a[1, 0]*I
        >>> a.is_hermitian
        False

        """
        def cond():
            yield self.is_square
            yield fuzzy_and(
                self[i, i].is_extended_real for i in range(self.rows))
            yield fuzzy_and(
                (self[i, j] - self[j, i].conjugate()).is_zero
                for i in range(self.rows)
                for j in range(i + 1, self.cols))
        return fuzzy_and(i for i in cond())

    @property
    def is_upper_hessenberg(self):
        """Checks if the matrix is the upper-Hessenberg form.

        The upper hessenberg matrix has zero entries
        below the first subdiagonal.

        Examples
        ========

        >>> a = Matrix([[1, 4, 2, 3], [3, 4, 1, 7], [0, 2, 3, 4], [0, 0, 1, 3]])
        >>> a
        Matrix([
        [1, 4, 2, 3],
        [3, 4, 1, 7],
        [0, 2, 3, 4],
        [0, 0, 1, 3]])
        >>> a.is_upper_hessenberg
        True

        See Also
        ========

        is_lower_hessenberg
        is_upper

        """
        return all(self[i, j].is_zero
                   for i in range(2, self.rows)
                   for j in range(min(self.cols, i - 1)))

    @property
    def is_lower_hessenberg(self):
        r"""Checks if the matrix is in the lower-Hessenberg form.

        The lower hessenberg matrix has zero entries
        above the first superdiagonal.

        Examples
        ========

        >>> a = Matrix([[1, 2, 0, 0], [5, 2, 3, 0], [3, 4, 3, 7], [5, 6, 1, 1]])
        >>> a
        Matrix([
        [1, 2, 0, 0],
        [5, 2, 3, 0],
        [3, 4, 3, 7],
        [5, 6, 1, 1]])
        >>> a.is_lower_hessenberg
        True

        See Also
        ========

        is_upper_hessenberg
        is_lower

        """
        return all(self[i, j].is_zero
                   for i in range(self.rows)
                   for j in range(i + 2, self.cols))

    def is_symbolic(self):
        """Checks if any elements contain Symbols.

        Examples
        ========

        >>> M = Matrix([[x, y], [1, 0]])
        >>> M.is_symbolic()
        True

        """
        return any(element.has(Symbol) for element in self.values())

    def is_symmetric(self, simplify=True):
        """Check if matrix is symmetric matrix,
        that is square matrix and is equal to its transpose.

        By default, simplifications occur before testing symmetry.
        They can be skipped using 'simplify=False'; while speeding things a bit,
        this may however induce false negatives.

        Examples
        ========

        >>> m = Matrix(2, 2, [0, 1, 1, 2])
        >>> m
        Matrix([
        [0, 1],
        [1, 2]])
        >>> m.is_symmetric()
        True

        >>> m = Matrix(2, 2, [0, 1, 2, 0])
        >>> m
        Matrix([
        [0, 1],
        [2, 0]])
        >>> m.is_symmetric()
        False

        >>> m = Matrix(2, 3, [0, 0, 0, 0, 0, 0])
        >>> m
        Matrix([
        [0, 0, 0],
        [0, 0, 0]])
        >>> m.is_symmetric()
        False

        >>> m = Matrix(3, 3, [1, x**2 + 2*x + 1, y, (x + 1)**2, 2, 0, y, 0, 3])
        >>> m
        Matrix([
        [         1, x**2 + 2*x + 1, y],
        [(x + 1)**2,              2, 0],
        [         y,              0, 3]])
        >>> m.is_symmetric()
        True

        If the matrix is already simplified, you may speed-up is_symmetric()
        test by using 'simplify=False'.

        >>> m.is_symmetric(simplify=False)
        False
        >>> m1 = m.expand()
        >>> m1.is_symmetric(simplify=False)
        True

        """
        if not self.is_square:
            return False
        if simplify:
            delta = self - self.transpose()
            delta.simplify()
            return delta.equals(self.zeros(self.rows, self.cols))
        else:
            return self == self.transpose()

    def is_anti_symmetric(self, simplify=True):
        """Check if matrix M is an antisymmetric matrix,
        that is, M is a square matrix with all M[i, j] == -M[j, i].

        When ``simplify=True`` (default), the sum M[i, j] + M[j, i] is
        simplified before testing to see if it is zero. By default,
        the Diofant simplify function is used. To use a custom function
        set simplify to a function that accepts a single argument which
        returns a simplified expression. To skip simplification, set
        simplify to False but note that although this will be faster,
        it may induce false negatives.

        Examples
        ========

        >>> m = Matrix(2, 2, [0, 1, -1, 0])
        >>> m
        Matrix([
        [ 0, 1],
        [-1, 0]])
        >>> m.is_anti_symmetric()
        True
        >>> m = Matrix(2, 3, [0, 0, x, -y, 0, 0])
        >>> m
        Matrix([
        [ 0, 0, x],
        [-y, 0, 0]])
        >>> m.is_anti_symmetric()
        False

        >>> m = Matrix(3, 3, [0, x**2 + 2*x + 1, y,
        ...                   -(x + 1)**2, 0, x*y,
        ...                   -y, -x*y, 0])

        Simplification of matrix elements is done by default so even
        though two elements which should be equal and opposite wouldn't
        pass an equality test, the matrix is still reported as
        anti-symmetric:

        >>> m[0, 1] == -m[1, 0]
        False
        >>> m.is_anti_symmetric()
        True

        If 'simplify=False' is used for the case when a Matrix is already
        simplified, this will speed things up. Here, we see that without
        simplification the matrix does not appear anti-symmetric:

        >>> m.is_anti_symmetric(simplify=False)
        False

        But if the matrix were already expanded, then it would appear
        anti-symmetric and simplification in the is_anti_symmetric routine
        is not needed:

        >>> m = m.expand()
        >>> m.is_anti_symmetric(simplify=False)
        True

        """
        # accept custom simplification
        simpfunc = simplify if isinstance(simplify, FunctionType) else \
            _simplify if simplify else False

        if not self.is_square:
            return False
        n = self.rows
        if simplify:
            for i in range(n):
                # diagonal
                if not simpfunc(self[i, i]).is_zero:
                    return False
                # others
                for j in range(i + 1, n):
                    diff = self[i, j] + self[j, i]
                    if not simpfunc(diff).is_zero:
                        return False
            return True
        else:
            for i in range(n):
                for j in range(i, n):
                    if self[i, j] != -self[j, i]:
                        return False
            return True

    def is_diagonal(self):
        """Check if matrix is diagonal,
        that is matrix in which the entries outside the main diagonal are all zero.

        Examples
        ========

        >>> m = Matrix(2, 2, [1, 0, 0, 2])
        >>> m
        Matrix([
        [1, 0],
        [0, 2]])
        >>> m.is_diagonal()
        True

        >>> m = Matrix(2, 2, [1, 1, 0, 2])
        >>> m
        Matrix([
        [1, 1],
        [0, 2]])
        >>> m.is_diagonal()
        False

        >>> m = diag(1, 2, 3)
        >>> m
        Matrix([
        [1, 0, 0],
        [0, 2, 0],
        [0, 0, 3]])
        >>> m.is_diagonal()
        True

        See Also
        ========

        is_lower
        is_upper
        is_diagonalizable
        diagonalize

        """
        for i in range(self.rows):
            for j in range(self.cols):
                if i != j and self[i, j]:
                    return False
        return True

    def det(self, method='bareiss'):
        """Computes the matrix determinant using the method "method".

        Possible values for "method":
          bareiss ... det_bareis
          berkowitz ... berkowitz_det
          det_LU ... det_LU_decomposition

        See Also
        ========

        det_bareiss
        berkowitz_det
        det_LU_decomposition

        """
        # if methods were made internal and all determinant calculations
        # passed through here, then these lines could be factored out of
        # the method routines
        if not self.is_square:
            raise NonSquareMatrixError()
        if not self:
            return Integer(1)
        if method == 'bareiss':
            return self.det_bareiss()
        elif method == 'berkowitz':
            return self.berkowitz_det()
        elif method == 'det_LU':
            return self.det_LU_decomposition()
        else:
            raise ValueError(f"Determinant method '{method}' unrecognized")

    def det_bareiss(self):
        """Compute matrix determinant using Bareiss' fraction-free
        algorithm which is an extension of the well known Gaussian
        elimination method. This approach is best suited for dense
        symbolic matrices and will result in a determinant with
        minimal number of fractions. It means that less term
        rewriting is needed on resulting formulae.

        TODO: Implement algorithm for sparse matrices (SFF),
        Hong R. Lee, B.David Saunders, Fraction Free Gaussian Elimination
        for Sparse Matrices, In Journal of Symbolic Computation, Volume 19,
        Issue 5, 1995, Pages 393-402, ISSN 0747-7171,
        https://www.sciencedirect.com/science/article/pii/S074771718571022X.

        See Also
        ========

        det
        berkowitz_det

        """
        if not self.is_square:
            raise NonSquareMatrixError()
        if not self:
            return Integer(1)

        M, n = self.copy().as_mutable(), self.rows

        if n == 1:
            det = M[0, 0]
        elif n == 2:
            det = M[0, 0]*M[1, 1] - M[0, 1]*M[1, 0]
        elif n == 3:
            det = (M[0, 0]*M[1, 1]*M[2, 2] + M[0, 1]*M[1, 2]*M[2, 0] + M[0, 2]*M[1, 0]*M[2, 1]) - \
                  (M[0, 2]*M[1, 1]*M[2, 0] + M[0, 0]*M[1, 2]*M[2, 1] + M[0, 1]*M[1, 0]*M[2, 2])
        else:
            sign = 1  # track current sign in case of column swap

            for k in range(n - 1):
                # look for a pivot in the current column
                # and assume det == 0 if none is found
                if M[k, k] == 0:
                    for i in range(k + 1, n):
                        if M[i, k]:
                            M.row_swap(i, k)
                            sign *= -1
                            break
                    else:
                        return Integer(0)

                # proceed with Bareiss' fraction-free (FF)
                # form of Gaussian elimination algorithm
                for i in range(k + 1, n):
                    for j in range(k + 1, n):
                        D = M[k, k]*M[i, j] - M[i, k]*M[k, j]

                        if k > 0:
                            D /= M[k - 1, k - 1]

                        if D.is_Atom:
                            M[i, j] = D
                        else:
                            M[i, j] = cancel(D)

            det = sign*M[n - 1, n - 1]

        return det.expand()

    def det_LU_decomposition(self):
        """Compute matrix determinant using LU decomposition

        Note that this method fails if the LU decomposition itself
        fails. In particular, if the matrix has no inverse this method
        will fail.

        TODO: Implement algorithm for sparse matrices (SFF),
        Hong R. Lee, B.David Saunders, Fraction Free Gaussian Elimination
        for Sparse Matrices, In Journal of Symbolic Computation, Volume 19,
        Issue 5, 1995, Pages 393-402, ISSN 0747-7171,
        https://www.sciencedirect.com/science/article/pii/S074771718571022X.

        See Also
        ========

        det
        det_bareiss
        berkowitz_det

        """
        if not self.is_square:
            raise NonSquareMatrixError()
        if not self:
            return Integer(1)

        M, n = self.copy(), self.rows
        p, prod = [], 1
        l, u, p = M.LUdecomposition()
        if len(p) % 2:
            prod = -1

        for k in range(n):
            prod = prod*u[k, k]*l[k, k]

        return prod.expand()

    def adjugate(self, method='berkowitz'):
        """Returns the adjugate matrix.

        Adjugate matrix is the transpose of the cofactor matrix.

        https://en.wikipedia.org/wiki/Adjugate

        See Also
        ========

        cofactorMatrix
        transpose
        berkowitz

        """
        return self.cofactorMatrix(method).T

    def inverse_LU(self, iszerofunc=_iszero):
        """Calculates the inverse using LU decomposition.

        See Also
        ========

        diofant.matrices.matrices.MatrixBase.inv
        inverse_GE
        inverse_ADJ

        """
        if not self.is_square:
            raise NonSquareMatrixError()

        ok = self.rref(simplify=True)[0]
        if any(iszerofunc(ok[j, j]) for j in range(ok.rows)):
            raise ValueError('Matrix det == 0; not invertible.')

        return self.LUsolve(self.eye(self.rows), iszerofunc=_iszero)

    def inverse_GE(self, iszerofunc=_iszero):
        """Calculates the inverse using Gaussian elimination.

        See Also
        ========

        diofant.matrices.matrices.MatrixBase.inv
        inverse_LU
        inverse_ADJ

        """
        from . import Matrix

        if not self.is_square:
            raise NonSquareMatrixError('A Matrix must be square to invert.')

        big = Matrix.hstack(self.as_mutable(), Matrix.eye(self.rows))
        red = big.rref(iszerofunc=iszerofunc, simplify=True)[0]
        if any(iszerofunc(red[j, j]) for j in range(red.rows)):
            raise ValueError('Matrix det == 0; not invertible.')

        return self._new(red[:, big.rows:])

    def inverse_ADJ(self, iszerofunc=_iszero):
        """Calculates the inverse using the adjugate matrix and a determinant.

        See Also
        ========

        diofant.matrices.matrices.MatrixBase.inv
        inverse_LU
        inverse_GE

        """
        if not self.is_square:
            raise NonSquareMatrixError('A Matrix must be square to invert.')

        d = self.berkowitz_det()
        zero = d.equals(0)
        if zero:
            raise ValueError('Matrix det == 0; not invertible.')

        return self.adjugate() / d

    def rref(self, iszerofunc=_iszero, simplify=False):
        """Return reduced row-echelon form of matrix and indices of pivot vars.

        To simplify elements before finding nonzero pivots set simplify=True
        (to use the default Diofant simplify function) or pass a custom
        simplify function.

        Examples
        ========

        >>> m = Matrix([[1, 2], [x, 1 - 1/x]])
        >>> m.rref()
        (Matrix([
        [1, 0],
        [0, 1]]), [0, 1])

        """
        simpfunc = simplify if isinstance(
            simplify, FunctionType) else _simplify
        # pivot: index of next row to contain a pivot
        pivot, r = 0, self.as_mutable()
        # pivotlist: indices of pivot variables (non-free)
        pivotlist = []
        for i in range(r.cols):
            if pivot == r.rows:
                break
            if simplify:
                r[pivot, i] = simpfunc(r[pivot, i])
            if iszerofunc(r[pivot, i]):
                for k in range(pivot, r.rows):
                    if simplify and k > pivot:
                        r[k, i] = simpfunc(r[k, i])
                    if not iszerofunc(r[k, i]):
                        r.row_swap(pivot, k)
                        break
                else:
                    continue
            scale = r[pivot, i]
            r.row_op(pivot, lambda x, _: x / scale)
            for j in range(r.rows):
                if j == pivot:
                    continue
                scale = r[j, i]
                r.zip_row_op(j, pivot, lambda x, y: x - y*scale)
            pivotlist.append(i)
            pivot += 1
        return self._new(r), pivotlist

    def rank(self, iszerofunc=_iszero, simplify=False):
        """
        Returns the rank of a matrix

        >>> m = Matrix([[1, 2], [x, 1 - 1/x]])
        >>> m.rank()
        2
        >>> n = Matrix(3, 3, range(1, 10))
        >>> n.rank()
        2

        """
        row_reduced = self.rref(iszerofunc=iszerofunc, simplify=simplify)
        rank = len(row_reduced[-1])
        return rank

    def nullspace(self, simplify=False, iszerofunc=_iszero):
        """Returns list of vectors (Matrix objects) that span nullspace of self."""
        from . import zeros

        simpfunc = simplify if isinstance(
            simplify, FunctionType) else _simplify
        reduced, pivots = self.rref(simplify=simpfunc, iszerofunc=iszerofunc)

        basis = []
        # create a set of vectors for the basis
        for i in range(self.cols - len(pivots)):
            basis.append(zeros(self.cols, 1))
        # contains the variable index to which the vector corresponds
        basiskey, cur = [-1]*len(basis), 0
        for i in range(self.cols):
            if i not in pivots:
                basiskey[cur] = i
                cur += 1
        for i in range(self.cols):
            if i not in pivots:  # free var, just set vector's ith place to 1
                basis[basiskey.index(i)][i, 0] = 1
            else:               # add negative of nonpivot entry to corr vector
                for j in range(i + 1, self.cols):
                    line = pivots.index(i)
                    v = reduced[line, j]
                    if simplify:
                        v = simpfunc(v)
                    if v:
                        if j in pivots:
                            # XXX: Is this the correct error?
                            raise NotImplementedError("Couldn't compute the "
                                                      'nullspace of `self`.')
                        basis[basiskey.index(j)][i, 0] = -v
        return [self._new(b) for b in basis]

    def berkowitz(self):
        """The Berkowitz algorithm.

           Given N x N matrix with symbolic content, compute efficiently
           coefficients of characteristic polynomials of 'self' and all
           its square sub-matrices composed by removing both i-th row
           and column, without division in the ground domain.

           This method is particularly useful for computing determinant,
           principal minors and characteristic polynomial, when 'self'
           has complicated coefficients e.g. polynomials. Semi-direct
           usage of this algorithm is also important in computing
           efficiently sub-resultant PRS.

           Assuming that M is a square matrix of dimension N x N and
           I is N x N identity matrix,  then the following following
           definition of characteristic polynomial is begin used:

                          charpoly(M) = det(t*I - M)

           As a consequence, all polynomials generated by Berkowitz
           algorithm are monic.

           >>> M = Matrix([[x, y, z], [1, 0, 0], [y, z, x]])

           >>> p, q, r = M.berkowitz()

           >>> p  # 1 x 1 M's sub-matrix
           (1, -x)

           >>> q  # 2 x 2 M's sub-matrix
           (1, -x, -y)

           >>> r  # 3 x 3 M's sub-matrix
           (1, -2*x, x**2 - y*z - y, x*y - z**2)

           For more information on the implemented algorithm refer to:

           [1] S.J. Berkowitz, On computing the determinant in small
               parallel time using a small number of processors, ACM,
               Information Processing Letters 18, 1984, pp. 147-150

           [2] M. Keber, Division-Free computation of sub-resultants
               using Bezout matrices, Tech. Report MPI-I-2006-1-006,
               Saarbrücken, 2006

        See Also
        ========

        berkowitz_det
        berkowitz_minors
        berkowitz_charpoly
        berkowitz_eigenvals

        """
        from . import zeros

        if not self.is_square:
            raise NonSquareMatrixError()

        A, N = self, self.rows
        transforms = [0]*(N - 1)

        for n in range(N, 1, -1):
            T, k = zeros(n + 1, n), n - 1

            R, C = -A[k, :k], A[:k, k]
            A, a = A[:k, :k], -A[k, k]

            items = [C]

            for i in range(n - 2):
                items.append(A*items[i])

            for i, B in enumerate(items):
                items[i] = (R*B)[0, 0]

            items = [Integer(1), a] + items

            for i in range(n):
                T[i:, i] = items[:n - i + 1]

            transforms[k - 1] = T

        polys = [self._new([Integer(1), -A[0, 0]])]

        for i, T in enumerate(transforms):
            polys.append(T*polys[i])

        return tuple(map(tuple, polys))

    def berkowitz_det(self):
        """Computes determinant using Berkowitz method.

        See Also
        ========

        det
        berkowitz

        """
        if not self.is_square:
            raise NonSquareMatrixError()
        if not self:
            return Integer(1)
        poly = self.berkowitz()[-1]
        sign = (-1)**(len(poly) - 1)
        return sign*poly[-1]

    def berkowitz_minors(self):
        """Computes principal minors using Berkowitz method.

        See Also
        ========

        berkowitz

        """
        sign, minors = Integer(-1), []

        for poly in self.berkowitz():
            minors.append(sign*poly[-1])
            sign = -sign

        return tuple(minors)

    def berkowitz_charpoly(self, x=Dummy('lambda'), simplify=_simplify):
        """Computes characteristic polynomial minors using Berkowitz method.

        A PurePoly is returned so using different variables for ``x`` does
        not affect the comparison or the polynomials:

        Examples
        ========

        >>> A = Matrix([[1, 3], [2, 0]])
        >>> A.berkowitz_charpoly(x) == A.berkowitz_charpoly(y)
        True

        Specifying ``x`` is optional; a Dummy with name ``lambda`` is used by
        default (which looks good when pretty-printed in unicode):

        >>> A.berkowitz_charpoly().as_expr()
        _lambda**2 - _lambda - 6

        Be sure your provided ``x`` doesn't clash with existing symbols:

        >>> A = Matrix([[1, 2], [x, 0]])
        >>> A.charpoly().as_expr()
        -2*x + _lambda**2 - _lambda
        >>> A.charpoly(x).as_expr()
        Traceback (most recent call last):
        ...
        GeneratorsError: polynomial ring and it's ground domain share generators
        >>> A.charpoly(y).as_expr()
        -2*x + y**2 - y

        See Also
        ========

        berkowitz

        """
        return PurePoly(list(map(simplify, reversed(self.berkowitz()[-1]))), x)

    charpoly = berkowitz_charpoly

    def berkowitz_eigenvals(self, **flags):
        """Computes eigenvalues of a Matrix using Berkowitz method.

        See Also
        ========

        berkowitz

        """
        return roots(self.berkowitz_charpoly(Dummy('x')), **flags)

    def eigenvals(self, **flags):
        """Return eigen values using the berkowitz_eigenvals routine.

        Since the roots routine doesn't always work well with Floats,
        they will be replaced with Rationals before calling that
        routine. If this is not desired, set flag ``rational`` to False.

        """
        # roots doesn't like Floats, so replace them with Rationals
        # unless the nsimplify flag indicates that this has already
        # been done, e.g. in eigenvects
        mat = self
        if not mat:
            return {}
        if flags.pop('rational', True):
            if any(v.has(Float) for v in mat):
                mat = mat._new(mat.rows, mat.cols,
                               [nsimplify(v, rational=True) for v in mat])

        flags.pop('simplify', None)  # pop unsupported flag
        return mat.berkowitz_eigenvals(**flags)

    def eigenvects(self, **flags):
        """Return list of triples (eigenval, multiplicity, basis).

        The flag ``simplify`` has two effects:
            1) if bool(simplify) is True, as_content_primitive()
            will be used to tidy up normalization artifacts;
            2) if nullspace needs simplification to compute the
            basis, the simplify flag will be passed on to the
            nullspace routine which will interpret it there.

        If the matrix contains any Floats, they will be changed to Rationals
        for computation purposes, but the answers will be returned after being
        evaluated with evalf. If it is desired to removed small imaginary
        portions during the evalf step, pass a value for the ``chop`` flag.

        """
        from . import eye

        primitive = bool(flags.get('simplify', False))
        chop = flags.pop('chop', False)

        flags.pop('multiple', None)  # remove this if it's there

        # roots doesn't like Floats, so replace them with Rationals
        float = False
        mat = self
        if any(v.has(Float) for v in self):
            float = True
            mat = mat._new(mat.rows, mat.cols, [nsimplify(
                v, rational=True) for v in mat])
            flags['rational'] = False  # to tell eigenvals not to do this

        out, vlist = [], mat.eigenvals(**flags)
        vlist = list(vlist.items())
        vlist.sort(key=default_sort_key)
        flags.pop('rational', None)

        for r, k in vlist:
            tmp = mat.as_mutable() - eye(mat.rows)*r
            basis = tmp.nullspace()
            # whether tmp.is_symbolic() is True or False, it is possible that
            # the basis will come back as [] in which case simplification is
            # necessary.

            if primitive:
                # the relationship A*e = lambda*e will still hold if we change the
                # eigenvector; so if simplify is True we tidy up any normalization
                # artifacts with as_content_primtive (default) and remove any pure Integer
                # denominators.
                l = 1
                for i, b in enumerate(basis[0]):
                    c, p = signsimp(b).as_content_primitive()
                    if c != 1:
                        b = c*p
                        l = math.lcm(l, c.denominator)
                    basis[0][i] = b
                if l != 1:
                    basis[0] *= l
            if float:
                out.append((r.evalf(chop=chop), k, [
                           mat._new(b).evalf(chop=chop) for b in basis]))
            else:
                out.append((r, k, [mat._new(b) for b in basis]))
        return out

    def singular_values(self):
        """Compute the singular values of a Matrix

        Examples
        ========

        >>> x = Symbol('x', real=True)
        >>> A = Matrix([[0, 1, 0], [0, x, 0], [-1, 0, 0]])
        >>> A.singular_values()
        [sqrt(x**2 + 1), 1, 0]

        See Also
        ========

        condition_number

        """
        mat = self.as_mutable()
        # Compute eigenvalues of A.H A
        valmultpairs = (mat.H*mat).eigenvals()

        # Expands result from eigenvals into a simple list
        vals = []
        for k, v in valmultpairs.items():
            vals += [sqrt(k)]*v  # dangerous! same k in several spots!
        # sort them in descending order
        vals.sort(reverse=True, key=default_sort_key)

        return vals

    def condition_number(self):
        """Returns the condition number of a matrix.

        This is the maximum singular value divided by the minimum singular value

        Examples
        ========

        >>> A = Matrix([[1, 0, 0], [0, 10, 0], [0, 0, Rational(1, 10)]])
        >>> A.condition_number()
        100

        See Also
        ========

        singular_values

        """
        if not self:
            return Integer(0)
        singularvalues = self.singular_values()
        return Max(*singularvalues) / Min(*singularvalues)

    def integrate(self, *args):
        """Integrate each element of the matrix.

        Examples
        ========

        >>> M = Matrix([[x, y], [1, 0]])
        >>> M.integrate(x)
        Matrix([
        [x**2/2, x*y],
        [     x,   0]])
        >>> M.integrate((x, 0, 2))
        Matrix([
        [2, 2*y],
        [2,   0]])

        See Also
        ========

        limit
        diff

        """
        return self._new(self.rows, self.cols,
                         lambda i, j: self[i, j].integrate(*args))

    def limit(self, *args):
        """Calculate the limit of each element in the matrix.

        Examples
        ========

        >>> M = Matrix([[x, y], [1, 0]])
        >>> M.limit(x, 2)
        Matrix([
        [2, y],
        [1, 0]])

        See Also
        ========

        integrate
        diff

        """
        return self._new(self.rows, self.cols,
                         lambda i, j: self[i, j].limit(*args))

    def diff(self, *args):
        """Calculate the derivative of each element in the matrix.

        Examples
        ========

        >>> M = Matrix([[x, y], [1, 0]])
        >>> M.diff(x)
        Matrix([
        [1, 0],
        [0, 0]])

        See Also
        ========

        integrate
        limit

        """
        return self._new(self.rows, self.cols,
                         lambda i, j: self[i, j].diff(*args))

    def vec(self):
        """Return the Matrix converted into a one column matrix by stacking columns

        Examples
        ========

        >>> m = Matrix([[1, 3], [2, 4]])
        >>> m
        Matrix([
        [1, 3],
        [2, 4]])
        >>> m.vec()
        Matrix([
        [1],
        [2],
        [3],
        [4]])

        See Also
        ========

        vech

        """
        return self.T.reshape(len(self), 1)

    def vech(self, diagonal=True, check_symmetry=True):
        """Return the unique elements of a symmetric Matrix as a one column matrix
        by stacking the elements in the lower triangle.

        Arguments:
        diagonal -- include the diagonal cells of self or not
        check_symmetry -- checks symmetry of self but not completely reliably

        Examples
        ========

        >>> m = Matrix([[1, 2], [2, 3]])
        >>> m
        Matrix([
        [1, 2],
        [2, 3]])
        >>> m.vech()
        Matrix([
        [1],
        [2],
        [3]])
        >>> m.vech(diagonal=False)
        Matrix([[2]])

        See Also
        ========

        vec

        """
        from . import zeros

        c = self.cols
        if c != self.rows:
            raise ShapeError('Matrix must be square')
        if check_symmetry:
            self.simplify()
            if self != self.transpose():
                raise ValueError('Matrix appears to be asymmetric; consider check_symmetry=False')
        count = 0
        if diagonal:
            v = zeros(c*(c + 1) // 2, 1)
            for j in range(c):
                for i in range(j, c):
                    v[count] = self[i, j]
                    count += 1
        else:
            v = zeros(c*(c - 1) // 2, 1)
            for j in range(c):
                for i in range(j + 1, c):
                    v[count] = self[i, j]
                    count += 1
        return v

    def get_diag_blocks(self):
        """Obtains the square sub-matrices on the main diagonal of a square matrix.

        Useful for inverting symbolic matrices or solving systems of
        linear equations which may be decoupled by having a block diagonal
        structure.

        Examples
        ========

        >>> A = Matrix([[1, 3, 0, 0], [y, z*z, 0, 0], [0, 0, x, 0], [0, 0, 0, 0]])
        >>> a1, a2, a3 = A.get_diag_blocks()
        >>> a1
        Matrix([
        [1,    3],
        [y, z**2]])
        >>> a2
        Matrix([[x]])
        >>> a3
        Matrix([[0]])

        """
        sub_blocks = []

        def recurse_sub_blocks(M):
            i = 1
            while i <= M.rows:
                if i == 1:
                    to_the_right = M[0, i:]
                    to_the_bottom = M[i:, 0]
                else:
                    to_the_right = M[:i, i:]
                    to_the_bottom = M[i:, :i]
                if not any(to_the_right) and not any(to_the_bottom):
                    sub_blocks.append(M[:i, :i])
                    if M.shape != M[:i, :i].shape:
                        return recurse_sub_blocks(M[i:, i:])
                i += 1
        recurse_sub_blocks(self)
        return sub_blocks

    def diagonalize(self, reals_only=False, sort=False, normalize=False):
        """
        Return (P, D), where D is diagonal and

            D = P^-1 * M * P

        where M is current matrix.

        Examples
        ========

        >>> m = Matrix(3, 3, [1, 2, 0, 0, 3, 0, 2, -4, 2])
        >>> m
        Matrix([
        [1,  2, 0],
        [0,  3, 0],
        [2, -4, 2]])
        >>> (P, D) = m.diagonalize()
        >>> D
        Matrix([
        [1, 0, 0],
        [0, 2, 0],
        [0, 0, 3]])
        >>> P
        Matrix([
        [-1, 0, -1],
        [ 0, 0, -1],
        [ 2, 1,  2]])
        >>> P.inv() * m * P
        Matrix([
        [1, 0, 0],
        [0, 2, 0],
        [0, 0, 3]])

        See Also
        ========

        is_diagonal
        is_diagonalizable

        """
        from . import diag

        if not self.is_square:
            raise NonSquareMatrixError()
        if not self.is_diagonalizable(reals_only, False):
            self._diagonalize_clear_subproducts()
            raise MatrixError('Matrix is not diagonalizable')
        assert self._eigenvects is not None
        if sort:
            self._eigenvects.sort(key=default_sort_key)
            self._eigenvects.reverse()
        diagvals = []
        P = self._new(self.rows, 0, [])
        for eigenval, multiplicity, vects in self._eigenvects:
            for k in range(multiplicity):
                diagvals.append(eigenval)
                vec = vects[k]
                if normalize:
                    vec = vec / vec.norm()
                P = P.col_insert(P.cols, vec)
        D = diag(*diagvals)
        self._diagonalize_clear_subproducts()
        return P, D

    def is_diagonalizable(self, reals_only=False, clear_subproducts=True):
        """Check if matrix is diagonalizable.

        If reals_only==True then check that diagonalized matrix consists of the only not complex values.

        Some subproducts could be used further in other methods to avoid double calculations,
        By default (if clear_subproducts==True) they will be deleted.

        Examples
        ========

        >>> m = Matrix(3, 3, [1, 2, 0, 0, 3, 0, 2, -4, 2])
        >>> m
        Matrix([
        [1,  2, 0],
        [0,  3, 0],
        [2, -4, 2]])
        >>> m.is_diagonalizable()
        True
        >>> m = Matrix(2, 2, [0, 1, 0, 0])
        >>> m
        Matrix([
        [0, 1],
        [0, 0]])
        >>> m.is_diagonalizable()
        False
        >>> m = Matrix(2, 2, [0, 1, -1, 0])
        >>> m
        Matrix([
        [ 0, 1],
        [-1, 0]])
        >>> m.is_diagonalizable()
        True
        >>> m.is_diagonalizable(True)
        False

        See Also
        ========

        is_diagonal
        diagonalize

        """
        # pylint: disable=attribute-defined-outside-init
        if not self.is_square:
            return False
        res = False
        self._is_symbolic = self.is_symbolic()
        self._is_symmetric = self.is_symmetric()
        self._eigenvects = None
        self._eigenvects = self.eigenvects(simplify=True)
        all_iscorrect = True
        for eigenval, multiplicity, vects in self._eigenvects:
            if len(vects) != multiplicity:
                all_iscorrect = False
                break
            if reals_only and not eigenval.is_extended_real:
                all_iscorrect = False
                break
        res = all_iscorrect
        if clear_subproducts:
            self._diagonalize_clear_subproducts()
        return res

    def _diagonalize_clear_subproducts(self):
        del self._is_symbolic
        del self._is_symmetric
        del self._eigenvects

    def jordan_cell(self, eigenval, n):
        n = int(n)
        from . import MutableMatrix
        out = MutableMatrix.zeros(n)
        for i in range(n-1):
            out[i, i] = eigenval
            out[i, i+1] = 1
        out[n-1, n-1] = eigenval
        return type(self)(out)

    def _jordan_block_structure(self):
        # To every eigenvalue may belong `i` blocks with size s(i)
        # and a chain of generalized eigenvectors
        # which will be determined by the following computations:
        # for every eigenvalue we will add a dictionary
        # containing, for all blocks, the blocksizes and the attached chain vectors
        # that will eventually be used to form the transformation P
        jordan_block_structures = {}
        _eigenvects = self.eigenvects()
        ev = self.eigenvals()
        if len(ev) == 0:
            raise AttributeError('could not compute the eigenvalues')
        for eigenval, multiplicity, vects in _eigenvects:
            l_jordan_chains = {}
            geometrical = len(vects)
            assert geometrical != 0
            if geometrical == multiplicity:
                # The Jordan chains have all length 1 and consist of only one vector
                # which is the eigenvector of course
                chains = []
                for v in vects:
                    chain = [v]
                    chains.append(chain)
                l_jordan_chains[1] = chains
                jordan_block_structures[eigenval] = l_jordan_chains
            else:
                # Up to now we know nothing about the sizes of the blocks of our Jordan matrix.
                # Note that knowledge of algebraic and geometrical multiplicity
                # will *NOT* be sufficient to determine this structure.
                # The blocksize `s` could be defined as the minimal `k` where
                # `kernel(self-lI)^k = kernel(self-lI)^(k+1)`
                # The extreme case would be that k = (multiplicity-geometrical+1)
                # but the blocks could be smaller.

                # Consider for instance the following matrix

                # [2 1 0 0]
                # [0 2 1 0]
                # [0 0 2 0]
                # [0 0 0 2]

                # which coincides with it own Jordan canonical form.
                # It has only one eigenvalue l=2 of (algebraic) multiplicity=4.
                # It has two eigenvectors, one belonging to the last row (blocksize 1)
                # and one being the last part of a jordan chain of length 3 (blocksize of the first block).

                # Note again that it is not not possible to obtain this from the algebraic and geometrical
                # multiplicity alone. This only gives us an upper limit for the dimension of one of
                # the subspaces (blocksize of according jordan block) given by
                # max=(multiplicity-geometrical+1) which is reached for our matrix
                # but not for

                # [2 1 0 0]
                # [0 2 0 0]
                # [0 0 2 1]
                # [0 0 0 2]

                # although multiplicity=4 and geometrical=2 are the same for this matrix.

                from . import MutableMatrix
                I = MutableMatrix.eye(self.rows)
                l = eigenval
                M = (self-l*I)

                # We will store the matrices `(self-l*I)^k` for further computations
                # for convenience only we store `Ms[0]=(sefl-lI)^0=I`
                # so the index is the same as the power for all further Ms entries
                # We also store the vectors that span these kernels (Ns[0] = [])
                # and also their dimensions `a_s`
                # this is mainly done for debugging since the number of blocks of a given size
                # can be computed from the a_s, in order to check our result which is obtained simpler
                # by counting the number of Jordan chains for `a` given `s`
                # `a_0` is `dim(Kernel(Ms[0]) = dim (Kernel(I)) = 0` since `I` is regular

                l_jordan_chains = {}
                Ms = [I]
                Ns = [[]]
                a = [0]
                smax = 0
                M_new = Ms[-1]*M
                Ns_new = M_new.nullspace()
                a_new = len(Ns_new)
                Ms.append(M_new)
                Ns.append(Ns_new)
                while a_new > a[-1]:  # as long as the nullspaces increase compute further powers
                    a.append(a_new)
                    M_new = Ms[-1]*M
                    Ns_new = M_new.nullspace()
                    a_new = len(Ns_new)
                    Ms.append(M_new)
                    Ns.append(Ns_new)
                    smax += 1

                # We now have `Ms[-1]=((self-l*I)**s)=Z=0`.
                # We also know the size of the biggest Jordan block
                # associated with `l` to be `s`.
                # Now let us proceed with the computation of the associate part of the transformation matrix `P`.
                # We already know the kernel (=nullspace)  `K_l` of (self-lI) which consists of the
                # eigenvectors belonging to eigenvalue `l`.
                # The dimension of this space is the geometric multiplicity of eigenvalue `l`.
                # For every eigenvector ev out of `K_l`, there exists a subspace that is
                # spanned by the Jordan chain of ev. The dimension of this subspace is
                # represented by the length `s` of the Jordan block.
                # The chain itself is given by `{e_0,..,e_s-1}` where:
                # `e_k+1 =(self-lI)e_k (*)`
                # and
                # `e_s-1=ev`
                # So it would be possible to start with the already known `ev` and work backwards until one
                # reaches `e_0`. Unfortunately this can not be done by simply solving system (*) since its matrix
                # is singular (by definition of the eigenspaces).
                # This approach would force us a choose in every step the degree of freedom undetermined
                # by (*). This is difficult to implement with computer algebra systems and also quite inefficient.
                # We therefore reformulate the problem in terms of nullspaces.
                # To do so we start from the other end and choose `e0`'s out of
                # `E=Kernel(self-lI)^s / Kernel(self-lI)^(s-1)`
                # Note that `Kernel(self-lI)^s = Kernel(Z) = V` (the whole vector space).
                # So in the first step `s=smax` this restriction turns out to actually restrict nothing at all
                # and the only remaining condition is to choose vectors in `Kernel(self-lI)^(s-1)`.
                # Subsequently we compute `e_1=(self-lI)e_0`, `e_2=(self-lI)*e_1` and so on.
                # The subspace `E` can have a dimension larger than one.
                # That means that we have more than one Jordan block of size `s` for the eigenvalue `l`
                # and as many Jordan chains (this is the case in the second example).
                # In this case we start as many Jordan chains and have as many blocks of size `s` in the jcf.
                # We now have all the Jordan blocks of size `s` but there might be others attached to the same
                # eigenvalue that are smaller.
                # So we will do the same procedure also for `s-1` and so on until 1 (the lowest possible order
                # where the Jordan chain is of length 1 and just represented by the eigenvector).

                for s in reversed(range(1, smax+1)):
                    S = Ms[s]
                    # We want the vectors in `Kernel((self-lI)^s)`,
                    # but without those in `Kernel(self-lI)^s-1`
                    # so we will add their adjoints as additional equations
                    # to the system formed by `S` to get the orthogonal
                    # complement.
                    # (`S` will no longer be quadratic.)

                    exclude_vectors = Ns[s-1]
                    for k in range(a[s-1]):
                        S = S.col_join((exclude_vectors[k]).adjoint())

                    # We also want to exclude the vectors
                    # in the chains for the bigger blocks
                    # that we have already computed (if there are any).
                    # (That is why we start with the biggest s).

                    # Since Jordan blocks are not orthogonal in general
                    # (in the original space), only those chain vectors
                    # that are on level s (index `s-1` in a chain)
                    # are added.

                    for chain_list in l_jordan_chains.values():
                        for chain in chain_list:
                            S = S.col_join(chain[s-1].adjoint())

                    e0s = S.nullspace()
                    # Determine the number of chain leaders
                    # for blocks of size `s`.
                    n_e0 = len(e0s)
                    s_chains = []
                    # s_cells=[]
                    for i in range(n_e0):
                        chain = [e0s[i]]
                        for k in range(1, s):
                            v = M*chain[k-1]
                            chain.append(v)

                        # We want the chain leader appear as the last of the block.
                        chain.reverse()
                        s_chains.append(chain)
                    l_jordan_chains[s] = s_chains
            jordan_block_structures[eigenval] = l_jordan_chains
        return jordan_block_structures

    def jordan_form(self, calc_transformation=True):
        r"""Return Jordan form J of current matrix.

        Also (if calc_transformation=True) the transformation P such that

        .. math::
            J = P^{-1} \cdot M \cdot P

        and the jordan blocks forming J will be calculated.

        Examples
        ========

        >>> m = Matrix([[+6,  5, -2, -3],
        ...             [-3, -1,  3,  3],
        ...             [+2,  1, -2, -3],
        ...             [-1,  1,  5,  5]])
        >>> J, P = m.jordan_form()
        >>> J
        Matrix([
        [2, 1, 0, 0],
        [0, 2, 0, 0],
        [0, 0, 2, 1],
        [0, 0, 0, 2]])

        See Also
        ========

        jordan_cells

        """
        from . import diag
        res = self.jordan_cells(calc_transformation=calc_transformation)
        if calc_transformation:
            Jcells, P = res
        else:
            Jcells = res
        J = diag(*Jcells)
        J = type(self)(J)
        if calc_transformation:
            return J, P
        else:
            return J

    def jordan_cells(self, calc_transformation=True):
        r"""Return a list of Jordan cells of current matrix.
        This list shape Jordan matrix J.

        If calc_transformation=False, then transformation P such that

        .. math::
            J = P^{-1} \cdot M \cdot P

        will not be calculated.

        Examples
        ========

        >>> m = Matrix([[+6,  5, -2, -3],
        ...             [-3, -1,  3,  3],
        ...             [+2,  1, -2, -3],
        ...             [-1,  1,  5,  5]])

        >>> Jcells, P = m.jordan_cells()
        >>> Jcells[0]
        Matrix([
        [2, 1],
        [0, 2]])
        >>> Jcells[1]
        Matrix([
        [2, 1],
        [0, 2]])

        See Also
        ========

        jordan_form

        """
        from . import MutableMatrix

        n = self.rows
        Jcells = []
        Pcols_new = []
        jordan_block_structures = self._jordan_block_structure()

        # Order according to default_sort_key, this makes sure the order is the same as in .diagonalize():
        for eigenval in (sorted(jordan_block_structures, key=default_sort_key)):
            l_jordan_chains = jordan_block_structures[eigenval]
            for s in sorted(l_jordan_chains, reverse=True):  # Start with the biggest block
                s_chains = l_jordan_chains[s]
                block = self.jordan_cell(eigenval, s)
                number_of_s_chains = len(s_chains)
                for i in range(number_of_s_chains):
                    Jcells.append(type(self)(block))
                    chain_vectors = s_chains[i]
                    lc = len(chain_vectors)
                    assert lc == s
                    for j in range(lc):
                        generalized_eigen_vector = chain_vectors[j]
                        Pcols_new.append(generalized_eigen_vector)

        if calc_transformation:
            P = MutableMatrix.zeros(n)
            for j in range(n):
                P[:, j] = Pcols_new[j]
            return Jcells, type(self)(P)
        else:
            return Jcells

    def has(self, *patterns):
        """Test whether any subexpression matches any of the patterns.

        Examples
        ========

        >>> A = Matrix(((1, x), (0.2, 3)))
        >>> A.has(x)
        True
        >>> A.has(y)
        False
        >>> A.has(Float)
        True

        """
        return any(a.has(*patterns) for a in self._mat)

    def dual(self):
        """Returns the dual of a matrix, which is:

        `(1/2)*levicivita(i, j, k, l)*M(k, l)` summed over indices `k` and `l`

        Since the levicivita method is anti_symmetric for any pairwise
        exchange of indices, the dual of a symmetric matrix is the zero
        matrix. Strictly speaking the dual defined here assumes that the
        'matrix' `M` is a contravariant anti_symmetric second rank tensor,
        so that the dual is a covariant second rank tensor.

        """
        from ..functions import LeviCivita
        from . import zeros

        M, n = self[:, :], self.rows
        work = zeros(n)
        if self.is_symmetric():
            return work

        for i in range(1, n):
            for j in range(1, n):
                acum = 0
                for k in range(1, n):
                    acum += LeviCivita(i, j, 0, k)*M[0, k]
                work[i, j] = acum
                work[j, i] = -acum

        for l in range(1, n):
            acum = 0
            for a in range(1, n):
                for b in range(1, n):
                    acum += LeviCivita(0, l, a, b)*M[a, b]
            acum /= 2
            work[0, l] = -acum
            work[l, 0] = acum

        return work

    @classmethod
    def hstack(cls, *args):
        """Return a matrix formed by joining args horizontally (i.e.
        by repeated application of row_join).

        Examples
        ========

        >>> Matrix.hstack(eye(2), 2*eye(2))
        Matrix([
        [1, 0, 2, 0],
        [0, 1, 0, 2]])

        """
        _cls = type(args[0])
        return functools.reduce(_cls.row_join, args)

    @classmethod
    def vstack(cls, *args):
        """Return a matrix formed by joining args vertically (i.e.
        by repeated application of col_join).

        Examples
        ========

        >>> Matrix.vstack(eye(2), 2*eye(2))
        Matrix([
        [1, 0],
        [0, 1],
        [2, 0],
        [0, 2]])

        """
        _cls = type(args[0])
        return functools.reduce(_cls.col_join, args)

    def row_join(self, rhs):
        """Concatenates two matrices along self's last and rhs's first column

        Examples
        ========

        >>> M = zeros(3)
        >>> V = ones(3, 1)
        >>> M.row_join(V)
        Matrix([
        [0, 0, 0, 1],
        [0, 0, 0, 1],
        [0, 0, 0, 1]])

        See Also
        ========

        col_join

        """
        from . import MutableMatrix

        if not self:
            return type(self)(rhs)
        if self.rows != rhs.rows:
            raise ShapeError('`self` and `rhs` must have the same'
                             ' number of rows.')

        newmat = MutableMatrix.zeros(self.rows, self.cols + rhs.cols)
        newmat[:, :self.cols] = self
        newmat[:, self.cols:] = rhs
        return type(self)(newmat)

    def col_join(self, bott):
        """Concatenates two matrices along self's last and bott's first row

        Examples
        ========

        >>> M = zeros(3)
        >>> V = ones(1, 3)
        >>> M.col_join(V)
        Matrix([
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 0],
        [1, 1, 1]])

        See Also
        ========

        row_join

        """
        from . import MutableMatrix

        if not self:
            return type(self)(bott)
        if self.cols != bott.cols:
            raise ShapeError('`self` and `bott` must have the same'
                             ' number of columns.')

        newmat = MutableMatrix.zeros(self.rows + bott.rows, self.cols)
        newmat[:self.rows, :] = self
        newmat[self.rows:, :] = bott
        return type(self)(newmat)

    def row_insert(self, pos, mti):
        """Insert one or more rows at the given row position.

        Examples
        ========

        >>> M = zeros(3)
        >>> V = ones(1, 3)
        >>> M.row_insert(1, V)
        Matrix([
        [0, 0, 0],
        [1, 1, 1],
        [0, 0, 0],
        [0, 0, 0]])

        See Also
        ========

        col_insert

        """
        if not self:
            return type(self)(mti)

        if pos == 0:
            return mti.col_join(self)
        elif pos < 0:
            pos = self.rows + pos
        if pos < 0:
            pos = 0
        elif pos > self.rows:
            pos = self.rows

        if self.cols != mti.cols:
            raise ShapeError(
                '`self` and `mti` must have the same number of columns.')

        newmat = self.zeros(self.rows + mti.rows, self.cols)
        i, j = pos, pos + mti.rows
        newmat[:i, :] = self[:i, :]
        newmat[i: j, :] = mti
        newmat[j:, :] = self[i:, :]
        return newmat

    def col_insert(self, pos, mti):
        """Insert one or more columns at the given column position.

        Examples
        ========

        >>> M = zeros(3)
        >>> V = ones(3, 1)
        >>> M.col_insert(1, V)
        Matrix([
        [0, 1, 0, 0],
        [0, 1, 0, 0],
        [0, 1, 0, 0]])

        See Also
        ========

        row_insert

        """
        from . import MutableMatrix

        if not self:
            return type(self)(mti)

        if pos == 0:
            return mti.row_join(self)
        elif pos < 0:
            pos = self.cols + pos
        if pos < 0:
            pos = 0
        elif pos > self.cols:
            pos = self.cols

        if self.rows != mti.rows:
            raise ShapeError('self and mti must have the same number of rows.')

        newmat = MutableMatrix.zeros(self.rows, self.cols + mti.cols)
        i, j = pos, pos + mti.cols
        newmat[:, :i] = self[:, :i]
        newmat[:, i:j] = mti
        newmat[:, j:] = self[:, i:]
        return type(self)(newmat)

    def replace(self, F, G, exact=False):
        """Replaces Function F in Matrix entries with Function G.

        Examples
        ========

        >>> F, G = symbols('F, G', cls=Function)
        >>> M = Matrix(2, 2, lambda i, j: F(i+j))
        >>> M
        Matrix([
        [F(0), F(1)],
        [F(1), F(2)]])
        >>> N = M.replace(F, G)
        >>> N
        Matrix([
        [G(0), G(1)],
        [G(1), G(2)]])

        """
        M = self[:, :]

        return M.applyfunc(lambda x: x.replace(F, G))

    def pinv(self):
        """Calculate the Moore-Penrose pseudoinverse of the matrix.

        The Moore-Penrose pseudoinverse exists and is unique for any matrix.
        If the matrix is invertible, the pseudoinverse is the same as the
        inverse.

        Examples
        ========

        >>> Matrix([[1, 2, 3], [4, 5, 6]]).pinv()
        Matrix([
        [-17/18,  4/9],
        [  -1/9,  1/9],
        [ 13/18, -2/9]])

        See Also
        ========

        diofant.matrices.matrices.MatrixBase.inv
        pinv_solve

        References
        ==========

        * https://en.wikipedia.org/wiki/Moore-Penrose_pseudoinverse

        """
        A = self
        AH = self.H
        # Trivial case: pseudoinverse of all-zero matrix is its transpose.
        if A.is_zero:
            return AH
        if self.rows >= self.cols:
            return (AH * A).inv() * AH
        else:
            return AH * (A * AH).inv()

    def pinv_solve(self, B, arbitrary_matrix=None):
        """Solve Ax = B using the Moore-Penrose pseudoinverse.

        There may be zero, one, or infinite solutions.  If one solution
        exists, it will be returned.  If infinite solutions exist, one will
        be returned based on the value of arbitrary_matrix.  If no solutions
        exist, the least-squares solution is returned.

        Parameters
        ==========

        B : Matrix
            The right hand side of the equation to be solved for.  Must have
            the same number of rows as matrix A.
        arbitrary_matrix : Matrix
            If the system is underdetermined (e.g. A has more columns than
            rows), infinite solutions are possible, in terms of an arbitrary
            matrix.  This parameter may be set to a specific matrix to use
            for that purpose; if so, it must be the same shape as x, with as
            many rows as matrix A has columns, and as many columns as matrix
            B.  If left as None, an appropriate matrix containing dummy
            symbols in the form of ``wn_m`` will be used, with n and m being
            row and column position of each symbol.

        Returns
        =======

        x : Matrix
            The matrix that will satisfy Ax = B.  Will have as many rows as
            matrix A has columns, and as many columns as matrix B.

        Examples
        ========

        >>> A = Matrix([[1, 2, 3], [4, 5, 6]])
        >>> B = Matrix([7, 8])
        >>> A.pinv_solve(B)
        Matrix([
        [ _w0_0/6 - _w1_0/3 + _w2_0/6 - 55/18],
        [-_w0_0/3 + 2*_w1_0/3 - _w2_0/3 + 1/9],
        [ _w0_0/6 - _w1_0/3 + _w2_0/6 + 59/18]])
        >>> A.pinv_solve(B, arbitrary_matrix=Matrix([0, 0, 0]))
        Matrix([
        [-55/18],
        [   1/9],
        [ 59/18]])

        See Also
        ========

        lower_triangular_solve
        upper_triangular_solve
        cholesky_solve
        diagonal_solve
        LDLsolve
        LUsolve
        QRsolve
        pinv

        Notes
        =====

        This may return either exact solutions or least squares solutions.
        To determine which, check ``A * A.pinv() * B == B``.  It will be
        True if exact solutions exist, and False if only a least-squares
        solution exists.  Be aware that the left hand side of that equation
        may need to be simplified to correctly compare to the right hand
        side.

        References
        ==========

        * https://en.wikipedia.org/wiki/Moore-Penrose_pseudoinverse#Obtaining_all_solutions_of_a_linear_system

        """
        from . import eye
        A = self
        A_pinv = self.pinv()
        if arbitrary_matrix is None:
            rows, cols = A.cols, B.cols
            w = symbols(f'w:{rows}_:{cols}', cls=Dummy)
            arbitrary_matrix = self.__class__(cols, rows, w).T
        return A_pinv * B + (eye(A.cols) - A_pinv*A) * arbitrary_matrix


def classof(A, B):
    """
    Get the type of the result when combining matrices of different types.

    Currently the strategy is that immutability is contagious.

    Examples
    ========

    >>> M = Matrix([[1, 2], [3, 4]])  # a Mutable Matrix
    >>> IM = ImmutableMatrix([[1, 2], [3, 4]])
    >>> classof(M, IM)
    <class 'diofant.matrices.immutable.ImmutableMatrix'>

    """
    try:
        if A._class_priority > B._class_priority:
            return A.__class__
        else:
            return B.__class__
    except AttributeError:
        pass
    raise TypeError(f'Incompatible classes {A.__class__}, {B.__class__}')


def a2idx(j, n=None):
    """Return integer after making positive and validating against n."""
    if type(j) is not int:
        try:
            j = getattr(j, '__index__')
            j = j()
        except AttributeError as exc:
            raise IndexError(f'Invalid index a[{j!r}]') from exc
    if n is not None:
        if j < 0:
            j += n
        if not 0 <= j < n:
            raise IndexError(f'Index out of range: a[{j}]')
    return int(j)


def mgamma(mu, lower=False):
    r"""Returns a Dirac gamma matrix `\gamma^\mu` in the standard
    (Dirac) representation.

    If you want `\gamma_\mu`, use ``gamma(mu, True)``.

    We use a convention:

    `\gamma^5 = i \cdot \gamma^0 \cdot \gamma^1 \cdot \gamma^2 \cdot \gamma^3`

    `\gamma_5 = i \cdot \gamma_0 \cdot \gamma_1 \cdot \gamma_2 \cdot \gamma_3 = - \gamma^5`

    References
    ==========

    * https://en.wikipedia.org/wiki/Gamma_matrices

    """
    from . import Matrix
    if mu not in [0, 1, 2, 3, 5]:
        raise IndexError('Invalid Dirac index')
    if mu == 0:
        mat = (
            (1, 0, 0, 0),
            (0, 1, 0, 0),
            (0, 0, -1, 0),
            (0, 0, 0, -1)
        )
    elif mu == 1:
        mat = (
            (0, 0, 0, 1),
            (0, 0, 1, 0),
            (0, -1, 0, 0),
            (-1, 0, 0, 0)
        )
    elif mu == 2:
        mat = (
            (0, 0, 0, -I),
            (0, 0, I, 0),
            (0, I, 0, 0),
            (-I, 0, 0, 0)
        )
    elif mu == 3:
        mat = (
            (0, 0, 1, 0),
            (0, 0, 0, -1),
            (-1, 0, 0, 0),
            (0, 1, 0, 0)
        )
    else:  # mu == 5
        mat = (
            (0, 0, 1, 0),
            (0, 0, 0, 1),
            (1, 0, 0, 0),
            (0, 1, 0, 0)
        )
    m = Matrix(mat)
    if lower:
        if mu in [1, 2, 3, 5]:
            m = -m
    return m
