from ...core import Expr
from ...core.sympify import sympify
from ..matrices import ShapeError


class Determinant(Expr):
    """Matrix Determinant

    Represents the determinant of a matrix expression.

    Examples
    ========

    >>> A = MatrixSymbol('A', 3, 3)
    >>> Determinant(A)
    Determinant(A)

    >>> Determinant(eye(3)).doit()
    1

    Determinant of the empty matrix:

    >>> Determinant(Matrix()).doit()
    1

    """

    def __new__(cls, mat):
        mat = sympify(mat)
        if not mat.is_Matrix:
            raise TypeError(f'Input to Determinant, {mat!s}, not a matrix')

        if not mat.is_square:
            raise ShapeError('Det of a non-square matrix')

        return Expr.__new__(cls, mat)

    @property
    def arg(self):
        return self.args[0]

    def doit(self, **hints):
        try:
            return self.arg._eval_determinant()
        except (AttributeError, NotImplementedError):
            return self


def det(matexpr):
    """Compute the Matrix Determinant

    >>> A = MatrixSymbol('A', 3, 3)
    >>> det(A)
    Determinant(A)

    >>> det(eye(3))
    1

    See Also
    ========

    Determinant

    """
    return Determinant(matexpr).doit()
