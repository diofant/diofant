from diofant import pi, I, Abs, simplify, exp, sqrt
from diofant.matrices.expressions.fourier import DFT, IDFT
from diofant.matrices import det, Matrix, Identity
from diofant.abc import n, i, j


def test_dft():
    assert DFT(4).shape == (4, 4)
    assert Abs(simplify(det(Matrix(DFT(4))))) == 1
    assert DFT(n)*IDFT(n) == Identity(n)
    assert DFT(n)[i, j] == exp(-2*pi*I/n)**(i*j) / sqrt(n)
    assert DFT(n).inverse() == IDFT(n)


def test_idft():
    assert IDFT(3).shape == (3, 3)
    assert IDFT(n)[i, j] == (exp(-2*I*pi/n))**(-i*j)/sqrt(n)
    assert IDFT(n).inverse() == DFT(n)
