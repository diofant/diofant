import typing

from ...core import AtomicExpr, Expr, Integer, Symbol, Tuple
from ...core.assumptions import StdFactKB
from ...core.decorators import _sympifyit, call_highest_priority
from ...core.logic import fuzzy_bool
from ...core.sympify import sympify
from ...functions import adjoint, conjugate
from ...logic import false
from ...simplify import simplify
from ..matrices import ShapeError


class MatrixExpr(Expr):
    """Superclass for Matrix Expressions

    MatrixExprs represent abstract matrices, linear transformations represented
    within a particular basis.

    Examples
    ========

    >>> A = MatrixSymbol('A', 3, 3)
    >>> y = MatrixSymbol('y', 3, 1)
    >>> x = (A.T*A).inverse() * A * y

    See Also
    ========
    MatrixSymbol
    MatAdd
    MatMul
    Transpose
    Inverse

    """

    _op_priority = 11.0

    is_Matrix = True
    is_MatrixExpr = True
    is_Identity: typing.Optional[bool] = None
    is_Inverse = False
    is_Transpose = False
    is_ZeroMatrix = False
    is_MatAdd = False
    is_MatMul = False

    def __new__(cls, *args, **kwargs):
        args = map(sympify, args)
        return Expr.__new__(cls, *args, **kwargs)

    # The following is adapted from the core Expr object
    def __neg__(self):
        from .matmul import MatMul
        return MatMul(-1, self).doit()

    def __abs__(self):
        raise NotImplementedError

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__radd__')
    def __add__(self, other):
        from .matadd import MatAdd
        return MatAdd(self, other).doit()

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__add__')
    def __radd__(self, other):
        from .matadd import MatAdd
        return MatAdd(other, self).doit()

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rsub__')
    def __sub__(self, other):
        from .matadd import MatAdd
        return MatAdd(self, -other).doit()

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__sub__')
    def __rsub__(self, other):
        from .matadd import MatAdd
        return MatAdd(other, -self).doit()

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rmul__')
    def __mul__(self, other):
        from .matmul import MatMul
        return MatMul(self, other).doit()

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__mul__')
    def __rmul__(self, other):
        from .matmul import MatMul
        return MatMul(other, self).doit()

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rpow__')
    def __pow__(self, other):
        from .inverse import Inverse
        from .matpow import MatPow
        if not self.is_square:
            raise ShapeError(f'Power of non-square matrix {self}')
        elif self.is_Identity:
            return self
        elif other == -1:
            return Inverse(self)
        elif other == 0:
            return Identity(self.rows)
        elif other == 1:
            return self
        return MatPow(self, other)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__pow__')
    def __rpow__(self, other):  # pragma: no cover
        raise NotImplementedError('Matrix Power not defined')

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rtruediv__')
    def __truediv__(self, other):
        return self * other**Integer(-1)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__truediv__')
    def __rtruediv__(self, other):
        raise NotImplementedError()
        # return MatMul(other, Pow(self, -1))

    @property
    def rows(self):
        return self.shape[0]

    @property
    def cols(self):
        return self.shape[1]

    @property
    def is_square(self):
        return self.rows == self.cols

    def _eval_conjugate(self):
        from .adjoint import Adjoint
        from .transpose import Transpose
        return Adjoint(Transpose(self))

    def _eval_inverse(self):
        from .inverse import Inverse
        return Inverse(self)

    def _eval_transpose(self):
        from .transpose import Transpose
        return Transpose(self)

    def _eval_power(self, exp):
        from .matpow import MatPow
        return MatPow(self, exp)

    def _eval_simplify(self, ratio, measure):
        if self.is_Atom:
            return self
        else:
            return self.__class__(*[simplify(x, ratio=ratio, measure=measure)
                                    for x in self.args])

    def _eval_adjoint(self):
        from .adjoint import Adjoint
        return Adjoint(self)

    def _entry(self, i, j):  # pragma: no cover
        raise NotImplementedError('Indexing not implemented '
                                  f'for {self.__class__.__name__}')

    def adjoint(self):
        return adjoint(self)

    def conjugate(self):
        return conjugate(self)

    def transpose(self):
        from .transpose import transpose
        return transpose(self)

    T = property(transpose, None, None, 'Matrix transposition.')

    def inverse(self):
        return self._eval_inverse()

    def valid_index(self, i, j):
        def is_valid(idx):
            return isinstance(idx, (int, Integer, Symbol, Expr))
        return (is_valid(i) and is_valid(j) and
                (0 <= i) != false and (i < self.rows) != false and
                (0 <= j) != false and (j < self.cols) != false)

    def __getitem__(self, key):
        if not isinstance(key, tuple) and isinstance(key, slice):
            from .slice import MatrixSlice
            return MatrixSlice(self, key, (0, None, 1))
        if isinstance(key, tuple) and len(key) == 2:
            i, j = key
            if isinstance(i, slice) or isinstance(j, slice):
                from .slice import MatrixSlice
                return MatrixSlice(self, i, j)
            i, j = sympify(i), sympify(j)
            if self.valid_index(i, j) is not False:
                return self._entry(i, j)
            else:
                raise IndexError(f'Invalid indices ({i}, {j})')
        elif isinstance(key, (int, Integer)):
            # row-wise decomposition of matrix
            rows, cols = self.shape
            if not (isinstance(rows, Integer) and isinstance(cols, Integer)):
                raise IndexError('Single index only supported for '
                                 'non-symbolic matrix shapes.')
            key = sympify(key)
            i = key // cols
            j = key % cols
            if self.valid_index(i, j) is not False:
                return self._entry(i, j)
            else:
                raise IndexError(f'Invalid index {key}')
        elif isinstance(key, (Symbol, Expr)):
            raise IndexError('Single index only supported for '
                             'non-symbolic indices.')
        raise IndexError(f'Invalid index, wanted {self}[i,j]')

    def as_explicit(self):
        """
        Returns a dense Matrix with elements represented explicitly

        Returns an object of type ImmutableMatrix.

        Examples
        ========

        >>> I = Identity(3)
        >>> I
        I
        >>> I.as_explicit()
        Matrix([
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1]])

        See Also
        ========
        as_mutable: returns mutable Matrix type

        """
        from ..immutable import ImmutableMatrix
        return ImmutableMatrix([[    self[i, j]
                                     for j in range(self.cols)]
                                for i in range(self.rows)])

    def as_mutable(self):
        """
        Returns a dense, mutable matrix with elements represented explicitly

        Examples
        ========

        >>> I = Identity(3)
        >>> I
        I
        >>> I.shape
        (3, 3)
        >>> I.as_mutable()
        Matrix([
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1]])

        See Also
        ========
        as_explicit: returns ImmutableMatrix

        """
        return self.as_explicit().as_mutable()

    def __array__(self):
        from numpy import empty
        a = empty(self.shape, dtype=object)
        for i in range(self.rows):
            for j in range(self.cols):
                a[i, j] = self[i, j]
        return a

    def equals(self, other):
        """
        Test elementwise equality between matrices, potentially of different
        types

        >>> Identity(3).equals(eye(3))
        True

        """
        if all(x.is_Integer for x in self.shape):
            return self.as_explicit().equals(other)

    def canonicalize(self):
        return self

    def as_coeff_mmul(self):
        from .matmul import MatMul
        return 1, MatMul(self)


class MatrixElement(Expr):
    """Element of the matrix expression."""

    parent = property(lambda self: self.args[0])
    i = property(lambda self: self.args[1])
    j = property(lambda self: self.args[2])
    _diff_wrt = True

    def __new__(cls, name, n, m):
        n, m = map(sympify, (n, m))
        from .. import MatrixBase
        if isinstance(name, MatrixBase):
            if n.is_Integer and m.is_Integer:
                return name[n, m]
        name = sympify(name)
        return Expr.__new__(cls, name, n, m)

    def xreplace(self, rule):
        if self in rule:
            return rule[self]
        else:
            return self


class MatrixSymbol(MatrixExpr, AtomicExpr):
    """Symbolic representation of a Matrix object

    Creates a Diofant Symbol to represent a Matrix. This matrix has a shape and
    can be included in Matrix Expressions

    >>> A = MatrixSymbol('A', 3, 4)  # A 3 by 4 Matrix
    >>> B = MatrixSymbol('B', 4, 3)  # A 4 by 3 Matrix
    >>> A.shape
    (3, 4)
    >>> 2*A*B + Identity(3)
    I + 2*A*B

    """

    is_Atom = True

    is_number = False

    def __new__(cls, name, n, m, **assumptions):
        n, m = sympify(n), sympify(m)
        is_commutative = fuzzy_bool(assumptions.get('commutative', False))
        assumptions['commutative'] = is_commutative
        obj = Expr.__new__(cls)
        obj._name = name
        obj._shape = (n, m)
        obj._assumptions = StdFactKB(assumptions)
        return obj

    def _hashable_content(self):
        return ((self.name, self.shape) +
                tuple(sorted((k, v) for k, v in self._assumptions.items()
                             if v is not None)))

    @property
    def shape(self):
        return self._shape

    @property
    def name(self):
        return self._name

    def _eval_subs(self, old, new):
        # only do substitutions in shape
        shape = Tuple(*self.shape)._subs(old, new)
        return MatrixSymbol(self.name, *shape)

    def __call__(self, *args):
        raise TypeError( f'{self.__class__} object is not callable' )

    def _entry(self, i, j):
        return MatrixElement(self, i, j)

    @property
    def free_symbols(self):
        return {self}

    def doit(self, **hints):
        if hints.get('deep', True):
            return type(self)(self.name,
                              *(_.doit(**hints) for _ in self.shape),
                              **self._assumptions._generator)
        else:
            return self


class Identity(MatrixExpr):
    """The Matrix Identity I - multiplicative identity

    >>> A = MatrixSymbol('A', 3, 5)
    >>> I = Identity(3)
    >>> I*A
    A

    """

    is_Identity = True

    def __new__(cls, n):
        return super().__new__(cls, sympify(n))

    @property
    def rows(self):
        return self.args[0]

    @property
    def cols(self):
        return self.args[0]

    @property
    def shape(self):
        return self.args[0], self.args[0]

    def _eval_transpose(self):
        return self

    def _eval_trace(self):
        return self.rows

    def _eval_inverse(self):
        return self

    def conjugate(self):
        return self

    def _entry(self, i, j):
        if i == j:
            return Integer(1)
        else:
            return Integer(0)

    def _eval_determinant(self):
        return Integer(1)


class ZeroMatrix(MatrixExpr):
    """The Matrix Zero 0 - additive identity

    >>> A = MatrixSymbol('A', 3, 5)
    >>> Z = ZeroMatrix(3, 5)
    >>> A+Z
    A
    >>> Z*A.T
    0

    """

    is_ZeroMatrix = True

    def __new__(cls, m, n):
        return super().__new__(cls, m, n)

    @property
    def shape(self):
        return self.args[0], self.args[1]

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rpow__')
    def __pow__(self, other):
        if other != 1 and not self.is_square:
            raise ShapeError(f'Power of non-square matrix {self}')
        if other == 0:
            return Identity(self.rows)
        if other < 1:
            raise ValueError('Matrix det == 0; not invertible.')
        return self

    def _eval_transpose(self):
        return ZeroMatrix(self.cols, self.rows)

    def _eval_trace(self):
        return Integer(0)

    def _eval_determinant(self):
        return Integer(0)

    def conjugate(self):
        return self

    def _entry(self, i, j):
        return Integer(0)

    def __bool__(self):
        return False
