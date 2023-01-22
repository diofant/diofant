from ...core.logic import fuzzy_and
from ...core.sympify import sympify
from ..matrices import MatrixBase, ShapeError
from .matexpr import Identity, MatrixExpr, ZeroMatrix
from .matmul import MatMul


class MatPow(MatrixExpr):
    """Power of matrix expression."""

    def __new__(cls, base, exp):
        base = sympify(base, strict=True)
        if not base.is_Matrix:
            raise TypeError('Function parameter should be a matrix')
        exp = sympify(exp, strict=True)
        return super().__new__(cls, base, exp)

    @property
    def base(self):
        return self.args[0]

    @property
    def exp(self):
        return self.args[1]

    @property
    def shape(self):
        return self.base.shape

    def _entry(self, i, j):
        A = self.doit()
        if isinstance(A, MatPow):
            # We still have a MatPow, make an explicit MatMul out of it.
            if not A.base.is_square:
                raise ShapeError(f'Power of non-square matrix {A.base}')
            if A.exp.is_Integer and A.exp.is_positive:
                A = MatMul(*[A.base for k in range(A.exp)])
            # elif A.exp.is_Integer and self.exp.is_negative:
            # Note: possible future improvement: in principle we can take
            # positive powers of the inverse, but carefully avoid recursion,
            # perhaps by adding `_entry` to Inverse (as it is our subclass).
            # T = A.base.as_explicit().inverse()
            # A = MatMul(*[T for k in range(-A.exp)])
            else:
                raise NotImplementedError((f'({int(i):d}, {int(j):d}) entry')
                                          + ' of matrix power either not defined or not implemented')
        return A._entry(i, j)

    def _eval_is_commutative(self):
        return fuzzy_and([self.base.is_commutative, self.exp.is_commutative])

    def doit(self, **kwargs):
        deep = kwargs.get('deep', True)
        if deep:
            args = [arg.doit(**kwargs) for arg in self.args]
        else:
            args = self.args
        base = args[0]
        exp = args[1]
        if exp.is_zero and base.is_square:
            if isinstance(base, MatrixBase):
                return base.func(Identity(base.shape[0]))
            return Identity(base.shape[0])
        if isinstance(base, ZeroMatrix) and exp.is_negative:
            raise ValueError('Matrix det == 0; not invertible.')
        if isinstance(base, (Identity, ZeroMatrix)):
            return base
        if isinstance(base, MatrixBase) and exp.is_number:
            if exp == 1:
                return base
            return base**exp
        # Note: just evaluate cases we know, return unevaluated on others.
        # E.g., MatrixSymbol('x', n, m) to power 0 is not an error.
        if exp == 1:
            return base
        return MatPow(base, exp)
