import pytest

from diofant import Function, I, MatrixSymbol, asinh, sin
from diofant.abc import x, y
from diofant.printing.lambdarepr import NumExprPrinter
from diofant.utilities.lambdify import implemented_function


__all__ = ()


def test_numexprprinter():
    p = NumExprPrinter()
    M = MatrixSymbol('M', 1, 2)

    pytest.raises(TypeError, lambda: p.doprint(M))
    pytest.raises(TypeError, lambda: p.doprint([x, y]))

    assert p.doprint(I) == "evaluate('1j')"
    assert p.doprint(sin(x)) == "evaluate('sin(x)')"
    assert p.doprint(asinh(x)) == "evaluate('arcsinh(x)')"

    f = implemented_function('f', lambda x: 2*x)
    assert p.doprint(f(x)) == "evaluate('(2*x)')"

    g = Function('g')
    pytest.raises(TypeError, lambda: p.doprint(g(x)))
