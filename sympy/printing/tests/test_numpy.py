from sympy import Piecewise, MatrixSymbol
from sympy.printing.lambdarepr import NumPyPrinter

from sympy.abc import x


def test_numpy_piecewise_regression():
    """
    NumPyPrinter needs to print Piecewise()'s choicelist as a list to avoid
    breaking compatibility with numpy 1.8. This is not necessary in numpy 1.9+.
    See sympy/sympy#9747 and sympy/sympy#9749 for details.
    """
    p = Piecewise((1, x < 0), (0, True))
    assert NumPyPrinter().doprint(p) == 'select([x < 0,True], [1,0], default=nan)'


def test_numpyprinter():
    A, B = MatrixSymbol('A', 2, 2), MatrixSymbol('B', 2, 2)
    assert NumPyPrinter().doprint(A*B) == '(A).dot(B)'
