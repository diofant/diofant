from ...core import I, pi
from ...functions import exp, sqrt
from .matexpr import MatrixExpr


class DFT(MatrixExpr):
    """Discrete Fourier Transform."""

    n = property(lambda self: self.args[0])
    shape = property(lambda self: (self.n, self.n))

    def _entry(self, i, j):
        w = exp(-2*pi*I/self.n)
        return w**(i*j) / sqrt(self.n)

    def _eval_inverse(self):
        return IDFT(self.n)


class IDFT(DFT):
    """Inverse Discrete Fourier Transform."""

    def _entry(self, i, j):
        w = exp(-2*pi*I/self.n)
        return w**(-i*j) / sqrt(self.n)

    def _eval_inverse(self):
        return DFT(self.n)
