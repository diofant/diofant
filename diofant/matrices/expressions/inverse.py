from ...core import Expr, Integer
from ...core.sympify import sympify
from ..matrices import ShapeError
from .matpow import MatPow


class Inverse(MatPow):
    """
    The multiplicative inverse of a matrix expression

    This is a symbolic object that simply stores its argument without
    evaluating it. To actually compute the inverse, use the ``.inverse()``
    method of matrices.

    Examples
    ========

    >>> A = MatrixSymbol('A', 3, 3)
    >>> B = MatrixSymbol('B', 3, 3)
    >>> Inverse(A)
    A^-1
    >>> A.inverse() == Inverse(A)
    True
    >>> (A*B).inverse()
    B^-1*A^-1
    >>> Inverse(A*B)
    (A*B)^-1

    """

    is_Inverse = True
    exp = Integer(-1)

    def __new__(cls, mat):
        mat = sympify(mat, strict=True)
        if not mat.is_Matrix:
            raise TypeError('mat should be a matrix')
        if not mat.is_square:
            raise ShapeError(f'Inverse of non-square matrix {mat}')
        return Expr.__new__(cls, mat)

    @property
    def arg(self):
        return self.args[0]

    @property
    def shape(self):
        return self.arg.shape

    def _eval_inverse(self):
        return self.arg

    def _eval_determinant(self):
        from .determinant import det
        return 1/det(self.arg)

    def doit(self, **hints):
        if hints.get('deep', True):
            return self.arg.doit(**hints).inverse()
        else:
            return self.arg.inverse()
