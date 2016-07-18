from diofant import S, I, Abs, simplify, exp, sqrt
from diofant.matrices.expressions.fourier import DFT, IDFT
from diofant.matrices import det, Matrix, Identity
from diofant.abc import n, i, j


def test_dft():
    assert DFT(4).shape == (4, 4)
    assert Abs(simplify(det(Matrix(DFT(4))))) == 1
    assert DFT(n)*IDFT(n) == Identity(n)
    assert DFT(n)[i, j] == exp(-2*S.Pi*I/n)**(i*j) / sqrt(n)
