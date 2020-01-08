from ...core import Integer
from .matexpr import MatrixExpr


class DiagonalMatrix(MatrixExpr):
    """Diagonal matrix class."""

    arg = property(lambda self: self.args[0])
    shape = property(lambda self: (self.arg.shape[0], self.arg.shape[0]))

    def _entry(self, i, j):
        return Integer(0) if i != j else self.arg[i, 0]


class DiagonalOf(MatrixExpr):
    """Represents the matrix diagonal."""

    arg = property(lambda self: self.args[0])
    shape = property(lambda self: (self.arg.shape[0], Integer(1)))

    def _entry(self, i, j):
        return self.arg[i, i]
