from ...core import Expr, sympify


class Rank(Expr):
    """Matrix Rank

    Represents the rank of a matrix expression.

    Examples
    ========

    >>> A = MatrixSymbol('A', 3, 3)
    >>> Rank(A)
    Rank(A)

    >>> Rank(eye(3)).doit()
    3

    """

    def __new__(cls, mat):
        mat = sympify(mat)
        if not mat.is_Matrix:
            raise TypeError("Input to Rank, %s, not a matrix" % str(mat))

        return Expr.__new__(cls, mat)

    @property
    def arg(self):
        return self.args[0]

    def doit(self, **hints):
        try:
            return self.arg._eval_rank()
        except (AttributeError, NotImplementedError):
            return self


def rank(matexpr):
    """Compute the Matrix Rank

    >>> A = MatrixSymbol('A', 3, 3)
    >>> rank(A)
    Rank(A)

    >>> rank(eye(3))
    3

    See Also
    ========

    Rank

    """

    return Rank(matexpr).doit()
