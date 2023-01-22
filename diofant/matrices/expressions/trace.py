from ...core import Dummy, Expr
from ...core.sympify import sympify
from ..matrices import MatrixBase, ShapeError


class Trace(Expr):
    """Matrix Trace

    Represents the trace of a matrix expression.

    >>> A = MatrixSymbol('A', 3, 3)
    >>> Trace(A)
    Trace(A)

    """

    is_Trace = True

    def __new__(cls, mat):
        mat = sympify(mat)
        if not mat.is_Matrix:
            raise TypeError(f'input to Trace, {mat!s}, is not a matrix')

        if not mat.is_square:
            raise ShapeError('Trace of a non-square matrix')

        return Expr.__new__(cls, mat)

    def _eval_transpose(self):
        return self

    @property
    def arg(self):
        return self.args[0]

    def doit(self, **kwargs):
        if kwargs.get('deep', True):
            arg = self.arg.doit(**kwargs)
            try:
                return arg._eval_trace()
            except (AttributeError, NotImplementedError):
                return Trace(arg)
        else:
            # _eval_trace would go too deep here
            if isinstance(self.arg, MatrixBase):
                return trace(self.arg)
            return Trace(self.arg)

    def _eval_rewrite_as_Sum(self, arg):
        from ...concrete import Sum
        i = Dummy('i')
        return Sum(self.arg[i, i], (i, 0, self.arg.rows-1)).doit()


def trace(expr):
    """Trace of a Matrix.  Sum of the diagonal elements

    >>> X = MatrixSymbol('X', n, n)  # A square matrix
    >>> trace(2*X)
    2*Trace(X)

    >>> trace(eye(3))
    3

    See Also
    ========

    Trace

    """
    return Trace(expr).doit()
