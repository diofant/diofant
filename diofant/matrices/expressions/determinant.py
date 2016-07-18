from diofant import Basic, Expr, sympify
from .matexpr import ShapeError


class Determinant(Expr):
    """Matrix Determinant

    Represents the determinant of a matrix expression.

    Examples
    ========

    >>> from diofant import MatrixSymbol, Determinant, eye, Matrix
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
            raise TypeError("Input to Determinant, %s, not a matrix" % str(mat))

        if not mat.is_square:
            raise ShapeError("Det of a non-square matrix")

        return Basic.__new__(cls, mat)

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

    >>> from diofant import MatrixSymbol, det, eye
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
