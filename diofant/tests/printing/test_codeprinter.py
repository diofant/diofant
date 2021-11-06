"""Tests for generic code printing."""

import pytest

from diofant import Dummy, Idx, IndexedBase, Matrix, MatrixSymbol, symbols
from diofant.printing.codeprinter import Assignment, CodePrinter


__all__ = ()


def setup_test_printer(**kwargs):
    p = CodePrinter(settings=kwargs)
    p._not_supported = set()
    p._number_symbols = set()
    return p


def test_print_Dummy():
    d = Dummy('d')
    p = setup_test_printer()
    assert p._print_Dummy(d) == f'd_{d.dummy_index:d}'


def test_Assignment():
    x, y = symbols('x, y')
    A = MatrixSymbol('A', 3, 1)
    mat = Matrix([1, 2, 3])
    B = IndexedBase('B')
    n = symbols('n', integer=True)
    i = Idx('i', n)
    # Here we just do things to show they don't error
    Assignment(x, y)
    Assignment(x, 0)
    Assignment(A, mat)
    Assignment(A[1, 0], 0)
    Assignment(A[1, 0], x)
    Assignment(B[i], x)
    Assignment(B[i], 0)
    # Here we test things to show that they error
    # Matrix to scalar
    pytest.raises(ValueError, lambda: Assignment(B[i], A))
    pytest.raises(ValueError, lambda: Assignment(B[i], mat))
    pytest.raises(ValueError, lambda: Assignment(x, mat))
    pytest.raises(ValueError, lambda: Assignment(x, A))
    pytest.raises(ValueError, lambda: Assignment(A[1, 0], mat))
    # Scalar to matrix
    pytest.raises(ValueError, lambda: Assignment(A, x))
    pytest.raises(ValueError, lambda: Assignment(A, 0))
    # Non-atomic lhs
    pytest.raises(TypeError, lambda: Assignment(mat, A))
    pytest.raises(TypeError, lambda: Assignment(0, x))
    pytest.raises(TypeError, lambda: Assignment(x*x, 1))
    pytest.raises(TypeError, lambda: Assignment(A + A, mat))
    pytest.raises(TypeError, lambda: Assignment(B, 0))


def test_print_Symbol():
    x, y = symbols('x, if')

    p = setup_test_printer()
    assert p._print(x) == 'x'
    assert p._print(y) == 'if'

    p.reserved_words.update(['if'])
    assert p._print(y) == 'if_'

    p = setup_test_printer(error_on_reserved=True)
    p.reserved_words.update(['if'])
    with pytest.raises(ValueError):
        p._print(y)

    p = setup_test_printer(reserved_word_suffix='_He_Man')
    p.reserved_words.update(['if'])
    assert p._print(y) == 'if_He_Man'
