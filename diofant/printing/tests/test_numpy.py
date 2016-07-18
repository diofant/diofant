from diofant import Piecewise, MatrixSymbol, Equality, Unequality, Mod
from diofant.printing.lambdarepr import NumPyPrinter

from diofant.abc import x, a, b


def test_numpy_piecewise_regression():
    """
    NumPyPrinter needs to print Piecewise()'s choicelist as a list to avoid
    breaking compatibility with numpy 1.8. This is not necessary in numpy 1.9+.
    See sympy/sympy#9747 and sympy/sympy#9749 for details.
    """
    p = Piecewise((1, x < 0), (0, True))
    assert NumPyPrinter().doprint(p) == 'select([less(x, 0),True], [1,0], default=nan)'


def test_numpyprinter():
    A, B = MatrixSymbol('A', 2, 2), MatrixSymbol('B', 2, 2)
    assert NumPyPrinter().doprint(A*B) == '(A).dot(B)'


def test_relational():
    p = NumPyPrinter()

    e = Equality(x, 1)
    assert p.doprint(e) == 'equal(x, 1)'

    e = Unequality(x, 1)
    assert p.doprint(e) == 'not_equal(x, 1)'

    e = (x < 1)
    assert p.doprint(e) == 'less(x, 1)'

    e = (x <= 1)
    assert p.doprint(e) == 'less_equal(x, 1)'

    e = (x > 1)
    assert p.doprint(e) == 'greater(x, 1)'

    e = (x >= 1)
    assert p.doprint(e) == 'greater_equal(x, 1)'


def test_mod():
    p = NumPyPrinter()

    e = Mod(a, b)
    assert p.doprint(e) == 'a%b'
