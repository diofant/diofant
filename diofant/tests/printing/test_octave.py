import pytest

from diofant import octave_code as mcode
from diofant.abc import n, x, y, z
from diofant.core import (Catalan, E, EulerGamma, Function, GoldenRatio, I,
                          Integer, Lambda, Rational, Symbol, Tuple, nan, oo,
                          pi, zoo)
from diofant.functions import (Chi, Ci, LambertW, Piecewise, Shi, Si, airyai,
                               airyaiprime, airybi, airybiprime, besseli,
                               besselj, besselk, bessely, ceiling, cos, exp,
                               hankel1, hankel2, jn, laguerre, li, loggamma,
                               lowergamma, polygamma, sin, sqrt, uppergamma,
                               yn, zeta)
from diofant.logic import false, true
from diofant.matrices import (HadamardProduct, Identity, Matrix, MatrixSymbol,
                              SparseMatrix, eye)
from diofant.utilities.lambdify import implemented_function


__all__ = ()


def test_Integer():
    assert mcode(Integer(67)) == "67"
    assert mcode(Integer(-1)) == "-1"


def test_Rational():
    assert mcode(Rational(3, 7)) == "3/7"
    assert mcode(Rational(18, 9)) == "2"
    assert mcode(Rational(3, -7)) == "-3/7"
    assert mcode(Rational(-3, -7)) == "3/7"
    assert mcode(x + Rational(3, 7)) == "x + 3/7"
    assert mcode(Rational(3, 7)*x) == "3*x/7"


def test_Function():
    assert mcode(sin(x) ** cos(x)) == "sin(x).^cos(x)"
    assert mcode(abs(x)) == "abs(x)"
    assert mcode(ceiling(x)) == "ceil(x)"


def test_Pow():
    assert mcode(x**3) == "x.^3"
    assert mcode(x**(y**3)) == "x.^(y.^3)"
    assert mcode(x**Rational(2, 3)) == 'x.^(2/3)'
    g = implemented_function('g', Lambda(x, 2*x))
    assert mcode(1/(g(x)*3.5)**(x - y**x)/(x**2 + y)) == \
        "(3.5*2*x).^(-x + y.^x)./(x.^2 + y)"


def test_basic_ops():
    assert mcode(x*y) == "x.*y"
    assert mcode(x + y) == "x + y"
    assert mcode(x - y) == "x - y"
    assert mcode(-x) == "-x"


def test_1_over_x_and_sqrt():
    # 1.0 and 0.5 would do something different in regular StrPrinter,
    # but these are exact in IEEE floating point so no different here.
    assert mcode(1/x) == '1./x'
    assert mcode(x**-1) == mcode(x**-1.0) == '1./x'
    assert mcode(1/sqrt(x)) == '1./sqrt(x)'
    assert mcode(x**Rational(-1, 2)) == mcode(x**-0.5) == '1./sqrt(x)'
    assert mcode(sqrt(x)) == mcode(x**0.5) == 'sqrt(x)'
    assert mcode(1/pi) == '1/pi'
    assert mcode(pi**-1) == mcode(pi**-1.0) == '1/pi'
    assert mcode(pi**-0.5) == '1/sqrt(pi)'


def test_mix_number_mult_symbols():
    assert mcode(3*x) == "3*x"
    assert mcode(pi*x) == "pi*x"
    assert mcode(3/x) == "3./x"
    assert mcode(pi/x) == "pi./x"
    assert mcode(x/3) == "x/3"
    assert mcode(x/pi) == "x/pi"
    assert mcode(x*y) == "x.*y"
    assert mcode(x*y, order="none") == "x.*y"
    assert mcode(3*x*y) == "3*x.*y"
    assert mcode(3*pi*x*y) == "3*pi*x.*y"
    assert mcode(x/y) == "x./y"
    assert mcode(x*y**-2) == "x./y.^2"
    assert mcode(3*x/y) == "3*x./y"
    assert mcode(x*y/z) == "x.*y./z"
    assert mcode(x/y*z) == "x.*z./y"
    assert mcode(1/x/y) == "1./(x.*y)"
    assert mcode(2*pi*x/y/z) == "2*pi*x./(y.*z)"
    assert mcode(3*pi/x) == "3*pi./x"
    assert mcode(Rational(3, 5)) == "3/5"
    assert mcode(Rational(3, 5)*x) == "3*x/5"
    assert mcode(x/y/z) == "x./(y.*z)"
    assert mcode((x+y)/z) == "(x + y)./z"
    assert mcode((x+y)/(z+x)) == "(x + y)./(x + z)"
    assert mcode((x+y)/EulerGamma) == "(x + y)/0.5772156649015329"
    assert mcode(x/3/pi) == "x/(3*pi)"
    assert mcode(Rational(3, 5)*x*y/pi) == "3*x.*y/(5*pi)"


def test_mix_number_pow_symbols():
    assert mcode(pi**3) == 'pi^3'
    assert mcode(x**2) == 'x.^2'
    assert mcode(x**(pi**3)) == 'x.^(pi^3)'
    assert mcode(x**y) == 'x.^y'
    assert mcode(x**(y**z)) == 'x.^(y.^z)'
    assert mcode((x**y)**z) == '(x.^y).^z'


def test_imag():
    assert mcode(I) == "1i"
    assert mcode(5*I) == "5i"
    assert mcode(Rational(3, 2)*I) == "3*1i/2"
    assert mcode(3+4*I) == "3 + 4i"


def test_constants():
    assert mcode(pi) == "pi"
    assert mcode(oo) == "inf"
    assert mcode(-oo) == "-inf"
    assert mcode(nan) == "NaN"
    assert mcode(E) == "exp(1)"
    assert mcode(exp(1)) == "exp(1)"
    assert mcode(true) == "true"
    assert mcode(false) == "false"


def test_constants_other():
    assert mcode(2*GoldenRatio) == "2*(1+sqrt(5))/2"
    assert mcode(2*Catalan) == "2*0.915965594177219"
    assert mcode(2*EulerGamma) == "2*0.5772156649015329"


def test_boolean():
    assert mcode(True) == "true"
    assert mcode(False) == "false"
    assert mcode(x & y) == "x & y"
    assert mcode(x | y) == "x | y"
    assert mcode(~x) == "~x"
    assert mcode(x & y & z) == "x & y & z"
    assert mcode(x | y | z) == "x | y | z"
    assert mcode((x & y) | z) == "z | x & y"
    assert mcode((x | y) & z) == "z & (x | y)"


def test_Matrices():
    assert mcode(Matrix(1, 1, [10])) == "10"
    A = Matrix([[1, sin(x/2), abs(x)],
                [0, 1, pi],
                [0, exp(1), ceiling(x)]])
    expected = ("[1 sin(x/2)  abs(x);\n"
                "0        1      pi;\n"
                "0   exp(1) ceil(x)]")
    assert mcode(A) == expected
    # row and columns
    assert mcode(A[:, 0]) == "[1; 0; 0]"
    assert mcode(A[0, :]) == "[1 sin(x/2) abs(x)]"
    # empty matrices
    assert mcode(Matrix(0, 0, [])) == '[]'
    assert mcode(Matrix(0, 3, [])) == 'zeros(0, 3)'
    # annoying to read but correct
    assert mcode(Matrix([[x, x - y, -y]])) == "[x x - y -y]"


def test_vector_entries_hadamard():
    # For a row or column, user might to use the other dimension
    A = Matrix([[1, sin(2/x), 3*pi/x/5]])
    assert mcode(A) == "[1 sin(2./x) 3*pi./(5*x)]"
    assert mcode(A.T) == "[1; sin(2./x); 3*pi./(5*x)]"


def test_MatrixSymbol():
    n = Symbol('n', integer=True)
    A = MatrixSymbol('A', n, n)
    B = MatrixSymbol('B', n, n)
    assert mcode(A*B) == "A*B"
    assert mcode(B*A) == "B*A"
    assert mcode(2*A*B) == "2*A*B"
    assert mcode(B*2*A) == "2*B*A"
    assert mcode(A*(B + 3*Identity(n))) == "A*(B + 3*eye(n))"
    assert mcode(A**(x**2)) == "A^(x.^2)"
    assert mcode(A**3) == "A^3"
    assert mcode(A**Rational(1, 2)) == "A^(1/2)"


def test_special_matrices():
    assert mcode(6*Identity(3)) == "6*eye(3)"


def test_containers():
    assert mcode([1, 2, 3, [4, 5, [6, 7]], 8, [9, 10], 11]) == \
        "{1, 2, 3, {4, 5, {6, 7}}, 8, {9, 10}, 11}"
    assert mcode((1, 2, (3, 4))) == "{1, 2, {3, 4}}"
    assert mcode([1]) == "{1}"
    assert mcode((1,)) == "{1}"
    assert mcode(Tuple(*[1, 2, 3])) == "{1, 2, 3}"
    assert mcode((1, x*y, (3, x**2))) == "{1, x.*y, {3, x.^2}}"
    # scalar, matrix, empty matrix and empty list
    assert mcode((1, eye(3), Matrix(0, 0, []), [])) == "{1, [1 0 0;\n0 1 0;\n0 0 1], [], {}}"


def test_octave_noninline():
    source = mcode((x+y)/Catalan, assign_to='me', inline=False)
    expected = (
        "Catalan = 0.915965594177219;\n"
        "me = (x + y)/Catalan;"
    )
    assert source == expected


def test_octave_piecewise():
    expr = Piecewise((x, x < 1), (x**2, True))
    assert mcode(expr) == "((x < 1).*(x) + (~(x < 1)).*(x.^2))"
    assert mcode(expr, assign_to="r") == (
        "r = ((x < 1).*(x) + (~(x < 1)).*(x.^2));")
    assert mcode(expr, assign_to="r", inline=False) == (
        "if (x < 1)\n"
        "  r = x;\n"
        "else\n"
        "  r = x.^2;\n"
        "end")
    expr = Piecewise((x**2, x < 1), (x**3, x < 2), (x**4, x < 3), (x**5, True))
    expected = ("((x < 1).*(x.^2) + (~(x < 1)).*( ...\n"
                "(x < 2).*(x.^3) + (~(x < 2)).*( ...\n"
                "(x < 3).*(x.^4) + (~(x < 3)).*(x.^5))))")
    assert mcode(expr) == expected
    assert mcode(expr, assign_to="r") == "r = " + expected + ";"
    assert mcode(expr, assign_to="r", inline=False) == (
        "if (x < 1)\n"
        "  r = x.^2;\n"
        "elseif (x < 2)\n"
        "  r = x.^3;\n"
        "elseif (x < 3)\n"
        "  r = x.^4;\n"
        "else\n"
        "  r = x.^5;\n"
        "end")
    # Check that Piecewise without a True (default) condition error
    expr = Piecewise((x, x < 1), (x**2, x > 1), (sin(x), x > 0))
    pytest.raises(ValueError, lambda: mcode(expr))


def test_octave_piecewise_times_const():
    pw = Piecewise((x, x < 1), (x**2, True))
    assert mcode(2*pw) == "2*((x < 1).*(x) + (~(x < 1)).*(x.^2))"
    assert mcode(pw/x) == "((x < 1).*(x) + (~(x < 1)).*(x.^2))./x"
    assert mcode(pw/(x*y)) == "((x < 1).*(x) + (~(x < 1)).*(x.^2))./(x.*y)"
    assert mcode(pw/3) == "((x < 1).*(x) + (~(x < 1)).*(x.^2))/3"


def test_octave_matrix_assign_to():
    A = Matrix([[1, 2, 3]])
    assert mcode(A, assign_to='a') == "a = [1 2 3];"
    A = Matrix([[1, 2], [3, 4]])
    assert mcode(A, assign_to='A') == "A = [1 2;\n3 4];"


def test_octave_matrix_assign_to_more():
    # assigning to Symbol or MatrixSymbol requires lhs/rhs match
    A = Matrix([[1, 2, 3]])
    B = MatrixSymbol('B', 1, 3)
    C = MatrixSymbol('C', 2, 3)
    assert mcode(A, assign_to=B) == "B = [1 2 3];"
    pytest.raises(ValueError, lambda: mcode(A, assign_to=x))
    pytest.raises(ValueError, lambda: mcode(A, assign_to=C))


def test_octave_matrix_1x1():
    A = Matrix([[3]])
    B = MatrixSymbol('B', 1, 1)
    C = MatrixSymbol('C', 1, 2)
    assert mcode(A, assign_to=B) == "B = 3;"
    # FIXME?
    # assert mcode(A, assign_to=x) == "x = 3;"
    pytest.raises(ValueError, lambda: mcode(A, assign_to=C))


def test_octave_matrix_elements():
    A = Matrix([[x, 2, x*y]])
    assert mcode(A[0, 0]**2 + A[0, 1] + A[0, 2]) == "x.^2 + x.*y + 2"
    A = MatrixSymbol('AA', 1, 3)
    assert mcode(A) == "AA"
    assert mcode(A[0, 0]**2 + sin(A[0, 1]) + A[0, 2]) == \
        "sin(AA(1, 2)) + AA(1, 1).^2 + AA(1, 3)"
    assert mcode(sum(A)) == "AA(1, 1) + AA(1, 2) + AA(1, 3)"


def test_octave_boolean():
    assert mcode(True) == "true"
    assert mcode(true) == "true"
    assert mcode(False) == "false"
    assert mcode(false) == "false"


def test_octave_not_supported():
    assert mcode(zoo) == (
        "% Not supported in Octave:\n"
        "% ComplexInfinity\n"
        "zoo"
    )
    f = Function('f')
    assert mcode(f(x).diff(x)) == (
        "% Not supported in Octave:\n"
        "% Derivative\n"
        "Derivative(f(x), x)"
    )


def test_trick_indent_with_end_else_words():
    # words starting with "end" or "else" do not confuse the indenter
    t1 = Symbol('endless')
    t2 = Symbol('elsewhere')
    pw = Piecewise((t1, x < 0), (t2, x <= 1), (1, True))
    assert mcode(pw, inline=False) == (
        "if (x < 0)\n"
        "  endless\n"
        "elseif (x <= 1)\n"
        "  elsewhere\n"
        "else\n"
        "  1\n"
        "end")


def test_haramard():
    A = MatrixSymbol('A', 3, 3)
    B = MatrixSymbol('B', 3, 3)
    v = MatrixSymbol('v', 3, 1)
    h = MatrixSymbol('h', 1, 3)
    C = HadamardProduct(A, B)
    assert mcode(C) == "A.*B"
    assert mcode(C*v) == "(A.*B)*v"
    assert mcode(h*C*v) == "h*(A.*B)*v"
    assert mcode(C*A) == "(A.*B)*A"
    # mixing Hadamard and scalar strange b/c we vectorize scalars
    assert mcode(C*x*y) == "(x.*y)*(A.*B)"


def test_sparse():
    M = SparseMatrix(5, 6, {})
    M[2, 2] = 10
    M[1, 2] = 20
    M[1, 3] = 22
    M[0, 3] = 30
    M[3, 0] = x*y
    assert mcode(M) == (
        "sparse([4 2 3 1 2], [1 3 3 4 4], [x.*y 20 10 30 22], 5, 6)"
    )


def test_specfun():
    for f in [besselj, bessely, besseli, besselk]:
        assert mcode(f(n, x)) == f.__name__ + '(n, x)'
    assert mcode(hankel1(n, x)) == 'besselh(n, 1, x)'
    assert mcode(hankel2(n, x)) == 'besselh(n, 2, x)'
    assert mcode(airyai(x)) == 'airy(0, x)'
    assert mcode(airyaiprime(x)) == 'airy(1, x)'
    assert mcode(airybi(x)) == 'airy(2, x)'
    assert mcode(airybiprime(x)) == 'airy(3, x)'
    assert mcode(uppergamma(n, x)) == 'gammainc(x, n, \'upper\')'
    assert mcode(lowergamma(n, x)) == 'gammainc(x, n, \'lower\')'
    assert mcode(jn(n, x)) == 'sqrt(2)*sqrt(pi)*sqrt(1./x).*besselj(n + 1/2, x)/2'
    assert mcode(yn(n, x)) == 'sqrt(2)*sqrt(pi)*sqrt(1./x).*bessely(n + 1/2, x)/2'
    assert mcode(Chi(x)) == 'coshint(x)'
    assert mcode(Ci(x)) == 'cosint(x)'
    assert mcode(laguerre(n, x)) == 'laguerreL(n, x)'
    assert mcode(li(x)) == 'logint(x)'
    assert mcode(loggamma(x)) == 'gammaln(x)'
    assert mcode(polygamma(n, x)) == 'psi(n, x)'
    assert mcode(Shi(x)) == 'sinhint(x)'
    assert mcode(Si(x)) == 'sinint(x)'
    assert mcode(LambertW(x)) == 'lambertw(x)'
    assert mcode(LambertW(x, n)) == 'lambertw(n, x)'
    assert mcode(zeta(x)) == 'zeta(x)'
    assert mcode(zeta(x, y)) == '% Not supported in Octave:\n% zeta\nzeta(x, y)'
