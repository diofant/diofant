"""Octave printing tests."""

import pytest

from diofant import (Catalan, Chi, Ci, E, EulerGamma, Function, GoldenRatio,
                     HadamardProduct, I, Identity, Integer, Lambda, LambertW,
                     Matrix, MatrixSymbol, Piecewise, Rational, Shi, Si,
                     SparseMatrix, Symbol, Tuple, airyai, airyaiprime, airybi,
                     airybiprime, besseli, besselj, besselk, bessely, ceiling,
                     cos, exp, eye, false, hankel1, hankel2, jn, laguerre, li,
                     loggamma, lowergamma, nan, octave_code, oo, pi, polygamma,
                     sin, sqrt, true, uppergamma, yn, zeta, zoo)
from diofant.abc import n, x, y, z
from diofant.utilities.lambdify import implemented_function


__all__ = ()


def test_Integer():
    assert octave_code(Integer(67)) == '67'
    assert octave_code(Integer(-1)) == '-1'


def test_Rational():
    assert octave_code(Rational(3, 7)) == '3/7'
    assert octave_code(Rational(18, 9)) == '2'
    assert octave_code(Rational(3, -7)) == '-3/7'
    assert octave_code(Rational(-3, -7)) == '3/7'
    assert octave_code(x + Rational(3, 7)) == 'x + 3/7'
    assert octave_code(Rational(3, 7)*x) == '3*x/7'


def test_Function():
    assert octave_code(sin(x) ** cos(x)) == 'sin(x).^cos(x)'
    assert octave_code(abs(x)) == 'abs(x)'
    assert octave_code(ceiling(x)) == 'ceil(x)'


def test_Pow():
    assert octave_code(x**3) == 'x.^3'
    assert octave_code(x**(y**3)) == 'x.^(y.^3)'
    assert octave_code(x**Rational(2, 3)) == 'x.^(2/3)'
    g = implemented_function('g', Lambda(x, 2*x))
    assert octave_code(1/(g(x)*3.5)**(x - y**x)/(x**2 + y)) == \
        '(3.5*2*x).^(-x + y.^x)./(x.^2 + y)'


def test_basic_ops():
    assert octave_code(x*y) == 'x.*y'
    assert octave_code(x + y) == 'x + y'
    assert octave_code(x - y) == 'x - y'
    assert octave_code(-x) == '-x'


def test_1_over_x_and_sqrt():
    # 1.0 and 0.5 would do something different in regular StrPrinter,
    # but these are exact in IEEE floating point so no different here.
    assert octave_code(1/x) == '1./x'
    assert octave_code(x**-1) == octave_code(x**-1.0) == '1./x'
    assert octave_code(1/sqrt(x)) == '1./sqrt(x)'
    assert octave_code(x**Rational(-1, 2)) == octave_code(x**-0.5) == '1./sqrt(x)'
    assert octave_code(sqrt(x)) == octave_code(x**0.5) == 'sqrt(x)'
    assert octave_code(1/pi) == '1/pi'
    assert octave_code(pi**-1) == octave_code(pi**-1.0) == '1/pi'
    assert octave_code(pi**-0.5) == '1/sqrt(pi)'


def test_mix_number_mult_symbols():
    assert octave_code(3*x) == '3*x'
    assert octave_code(pi*x) == 'pi*x'
    assert octave_code(3/x) == '3./x'
    assert octave_code(pi/x) == 'pi./x'
    assert octave_code(x/3) == 'x/3'
    assert octave_code(x/pi) == 'x/pi'
    assert octave_code(x*y) == 'x.*y'
    assert octave_code(x*y, order='none') == 'x.*y'
    assert octave_code(3*x*y) == '3*x.*y'
    assert octave_code(3*pi*x*y) == '3*pi*x.*y'
    assert octave_code(x/y) == 'x./y'
    assert octave_code(x*y**-2) == 'x./y.^2'
    assert octave_code(3*x/y) == '3*x./y'
    assert octave_code(x*y/z) == 'x.*y./z'
    assert octave_code(x/y*z) == 'x.*z./y'
    assert octave_code(1/x/y) == '1./(x.*y)'
    assert octave_code(2*pi*x/y/z) == '2*pi*x./(y.*z)'
    assert octave_code(3*pi/x) == '3*pi./x'
    assert octave_code(Rational(3, 5)) == '3/5'
    assert octave_code(Rational(3, 5)*x) == '3*x/5'
    assert octave_code(x/y/z) == 'x./(y.*z)'
    assert octave_code((x+y)/z) == '(x + y)./z'
    assert octave_code((x+y)/(z+x)) == '(x + y)./(x + z)'
    assert octave_code((x+y)/EulerGamma) == '(x + y)/0.5772156649015329'
    assert octave_code(x/3/pi) == 'x/(3*pi)'
    assert octave_code(Rational(3, 5)*x*y/pi) == '3*x.*y/(5*pi)'


def test_mix_number_pow_symbols():
    assert octave_code(pi**3) == 'pi^3'
    assert octave_code(x**2) == 'x.^2'
    assert octave_code(x**(pi**3)) == 'x.^(pi^3)'
    assert octave_code(x**y) == 'x.^y'
    assert octave_code(x**(y**z)) == 'x.^(y.^z)'
    assert octave_code((x**y)**z) == '(x.^y).^z'


def test_imag():
    assert octave_code(I) == '1i'
    assert octave_code(5*I) == '5i'
    assert octave_code(Rational(3, 2)*I) == '3*1i/2'
    assert octave_code(3+4*I) == '3 + 4i'


def test_constants():
    assert octave_code(pi) == 'pi'
    assert octave_code(oo) == 'inf'
    assert octave_code(-oo) == '-inf'
    assert octave_code(nan) == 'NaN'
    assert octave_code(E) == 'exp(1)'
    assert octave_code(exp(1)) == 'exp(1)'
    assert octave_code(true) == 'true'
    assert octave_code(false) == 'false'


def test_constants_other():
    assert octave_code(2*GoldenRatio) == '2*(1+sqrt(5))/2'
    assert octave_code(2*Catalan) == '2*0.915965594177219'
    assert octave_code(2*EulerGamma) == '2*0.5772156649015329'


def test_boolean():
    assert octave_code(True) == 'true'
    assert octave_code(False) == 'false'
    assert octave_code(x & y) == 'x & y'
    assert octave_code(x | y) == 'x | y'
    assert octave_code(~x) == '~x'
    assert octave_code(x & y & z) == 'x & y & z'
    assert octave_code(x | y | z) == 'x | y | z'
    assert octave_code((x & y) | z) == 'z | x & y'
    assert octave_code((x | y) & z) == 'z & (x | y)'


def test_Matrices():
    assert octave_code(Matrix(1, 1, [10])) == '10'
    A = Matrix([[1, sin(x/2), abs(x)],
                [0, 1, pi],
                [0, exp(1), ceiling(x)]])
    expected = ('[1 sin(x/2)  abs(x);\n'
                '0        1      pi;\n'
                '0   exp(1) ceil(x)]')
    assert octave_code(A) == expected
    # row and columns
    assert octave_code(A[:, 0]) == '[1; 0; 0]'
    assert octave_code(A[0, :]) == '[1 sin(x/2) abs(x)]'
    # empty matrices
    assert octave_code(Matrix(0, 0, [])) == '[]'
    assert octave_code(Matrix(0, 3, [])) == 'zeros(0, 3)'
    # annoying to read but correct
    assert octave_code(Matrix([[x, x - y, -y]])) == '[x x - y -y]'


def test_vector_entries_hadamard():
    # For a row or column, user might to use the other dimension
    A = Matrix([[1, sin(2/x), 3*pi/x/5]])
    assert octave_code(A) == '[1 sin(2./x) 3*pi./(5*x)]'
    assert octave_code(A.T) == '[1; sin(2./x); 3*pi./(5*x)]'


def test_MatrixSymbol():
    n = Symbol('n', integer=True)
    A = MatrixSymbol('A', n, n)
    B = MatrixSymbol('B', n, n)
    assert octave_code(A*B) == 'A*B'
    assert octave_code(B*A) == 'B*A'
    assert octave_code(2*A*B) == '2*A*B'
    assert octave_code(B*2*A) == '2*B*A'
    assert octave_code(A*(B + 3*Identity(n))) == 'A*(B + 3*eye(n))'
    assert octave_code(A**(x**2)) == 'A^(x.^2)'
    assert octave_code(A**3) == 'A^3'
    assert octave_code(A**Rational(1, 2)) == 'A^(1/2)'


def test_special_matrices():
    assert octave_code(6*Identity(3)) == '6*eye(3)'


def test_containers():
    assert octave_code([1, 2, 3, [4, 5, [6, 7]], 8, [9, 10], 11]) == \
        '{1, 2, 3, {4, 5, {6, 7}}, 8, {9, 10}, 11}'
    assert octave_code((1, 2, (3, 4))) == '{1, 2, {3, 4}}'
    assert octave_code([1]) == '{1}'
    assert octave_code((1,)) == '{1}'
    assert octave_code(Tuple(*[1, 2, 3])) == '{1, 2, 3}'
    assert octave_code((1, x*y, (3, x**2))) == '{1, x.*y, {3, x.^2}}'
    # scalar, matrix, empty matrix and empty list
    assert octave_code((1, eye(3), Matrix(0, 0, []), [])) == '{1, [1 0 0;\n0 1 0;\n0 0 1], [], {}}'


def test_octave_noninline():
    source = octave_code((x+y)/Catalan, assign_to='me', inline=False)
    expected = (
        'Catalan = 0.915965594177219;\n'
        'me = (x + y)/Catalan;'
    )
    assert source == expected


def test_octave_piecewise():
    expr = Piecewise((x, x < 1), (x**2, True))
    assert octave_code(expr) == '((x < 1).*(x) + (~(x < 1)).*(x.^2))'
    assert octave_code(expr, assign_to='r') == (
        'r = ((x < 1).*(x) + (~(x < 1)).*(x.^2));')
    assert octave_code(expr, assign_to='r', inline=False) == (
        'if (x < 1)\n'
        '  r = x;\n'
        'else\n'
        '  r = x.^2;\n'
        'end')
    expr = Piecewise((x**2, x < 1), (x**3, x < 2), (x**4, x < 3), (x**5, True))
    expected = ('((x < 1).*(x.^2) + (~(x < 1)).*( ...\n'
                '(x < 2).*(x.^3) + (~(x < 2)).*( ...\n'
                '(x < 3).*(x.^4) + (~(x < 3)).*(x.^5))))')
    assert octave_code(expr) == expected
    assert octave_code(expr, assign_to='r') == 'r = ' + expected + ';'
    assert octave_code(expr, assign_to='r', inline=False) == (
        'if (x < 1)\n'
        '  r = x.^2;\n'
        'elseif (x < 2)\n'
        '  r = x.^3;\n'
        'elseif (x < 3)\n'
        '  r = x.^4;\n'
        'else\n'
        '  r = x.^5;\n'
        'end')


def test_octave_piecewise_times_const():
    pw = Piecewise((x, x < 1), (x**2, True))
    assert octave_code(2*pw) == '2*((x < 1).*(x) + (~(x < 1)).*(x.^2))'
    assert octave_code(pw/x) == '((x < 1).*(x) + (~(x < 1)).*(x.^2))./x'
    assert octave_code(pw/(x*y)) == '((x < 1).*(x) + (~(x < 1)).*(x.^2))./(x.*y)'
    assert octave_code(pw/3) == '((x < 1).*(x) + (~(x < 1)).*(x.^2))/3'


def test_octave_matrix_assign_to():
    A = Matrix([[1, 2, 3]])
    assert octave_code(A, assign_to='a') == 'a = [1 2 3];'
    A = Matrix([[1, 2], [3, 4]])
    assert octave_code(A, assign_to='A') == 'A = [1 2;\n3 4];'


def test_octave_matrix_assign_to_more():
    # assigning to Symbol or MatrixSymbol requires lhs/rhs match
    A = Matrix([[1, 2, 3]])
    B = MatrixSymbol('B', 1, 3)
    C = MatrixSymbol('C', 2, 3)
    assert octave_code(A, assign_to=B) == 'B = [1 2 3];'
    pytest.raises(ValueError, lambda: octave_code(A, assign_to=x))
    pytest.raises(ValueError, lambda: octave_code(A, assign_to=C))


def test_octave_matrix_1x1():
    A = Matrix([[3]])
    B = MatrixSymbol('B', 1, 1)
    C = MatrixSymbol('C', 1, 2)
    assert octave_code(A, assign_to=B) == 'B = 3;'
    pytest.raises(ValueError, lambda: octave_code(A, assign_to=C))


def test_octave_matrix_elements():
    A = Matrix([[x, 2, x*y]])
    assert octave_code(A[0, 0]**2 + A[0, 1] + A[0, 2]) == 'x.^2 + x.*y + 2'
    A = MatrixSymbol('AA', 1, 3)
    assert octave_code(A) == 'AA'
    assert octave_code(A[0, 0]**2 + sin(A[0, 1]) + A[0, 2]) == \
        'sin(AA(1, 2)) + AA(1, 1).^2 + AA(1, 3)'
    assert octave_code(sum(A)) == 'AA(1, 1) + AA(1, 2) + AA(1, 3)'


def test_octave_boolean():
    assert octave_code(True) == 'true'
    assert octave_code(true) == 'true'
    assert octave_code(False) == 'false'
    assert octave_code(false) == 'false'


def test_octave_not_supported():
    assert octave_code(zoo) == (
        '% Not supported in Octave:\n'
        '% ComplexInfinity\n'
        'zoo'
    )
    f = Function('f')
    assert octave_code(f(x).diff(x)) == (
        '% Not supported in Octave:\n'
        '% Derivative\n'
        'Derivative(f(x), x)'
    )


def test_trick_indent_with_end_else_words():
    # words starting with 'end' or 'else' do not confuse the indenter
    t1 = Symbol('endless')
    t2 = Symbol('elsewhere')
    pw = Piecewise((t1, x < 0), (t2, x <= 1), (1, True))
    assert octave_code(pw, inline=False) == (
        'if (x < 0)\n'
        '  endless\n'
        'elseif (x <= 1)\n'
        '  elsewhere\n'
        'else\n'
        '  1\n'
        'end')


def test_haramard():
    A = MatrixSymbol('A', 3, 3)
    B = MatrixSymbol('B', 3, 3)
    v = MatrixSymbol('v', 3, 1)
    h = MatrixSymbol('h', 1, 3)
    C = HadamardProduct(A, B)
    assert octave_code(C) == 'A.*B'
    assert octave_code(C*v) == '(A.*B)*v'
    assert octave_code(h*C*v) == 'h*(A.*B)*v'
    assert octave_code(C*A) == '(A.*B)*A'
    # mixing Hadamard and scalar strange b/c we vectorize scalars
    assert octave_code(C*x*y) == '(x.*y)*(A.*B)'


def test_sparse():
    M = SparseMatrix(5, 6, {})
    M[2, 2] = 10
    M[1, 2] = 20
    M[1, 3] = 22
    M[0, 3] = 30
    M[3, 0] = x*y
    assert octave_code(M) == (
        'sparse([4 2 3 1 2], [1 3 3 4 4], [x.*y 20 10 30 22], 5, 6)'
    )


def test_specfun():
    for f in [besselj, bessely, besseli, besselk]:
        assert octave_code(f(n, x)) == f.__name__ + '(n, x)'
    assert octave_code(hankel1(n, x)) == 'besselh(n, 1, x)'
    assert octave_code(hankel2(n, x)) == 'besselh(n, 2, x)'
    assert octave_code(airyai(x)) == 'airy(0, x)'
    assert octave_code(airyaiprime(x)) == 'airy(1, x)'
    assert octave_code(airybi(x)) == 'airy(2, x)'
    assert octave_code(airybiprime(x)) == 'airy(3, x)'
    assert octave_code(uppergamma(n, x)) == "gammainc(x, n, 'upper')"
    assert octave_code(lowergamma(n, x)) == "gammainc(x, n, 'lower')"
    assert octave_code(jn(n, x)) == 'sqrt(2)*sqrt(pi)*sqrt(1./x).*besselj(n + 1/2, x)/2'
    assert octave_code(yn(n, x)) == 'sqrt(2)*sqrt(pi)*sqrt(1./x).*bessely(n + 1/2, x)/2'
    assert octave_code(Chi(x)) == 'coshint(x)'
    assert octave_code(Ci(x)) == 'cosint(x)'
    assert octave_code(laguerre(n, x)) == 'laguerreL(n, x)'
    assert octave_code(li(x)) == 'logint(x)'
    assert octave_code(loggamma(x)) == 'gammaln(x)'
    assert octave_code(polygamma(n, x)) == 'psi(n, x)'
    assert octave_code(Shi(x)) == 'sinhint(x)'
    assert octave_code(Si(x)) == 'sinint(x)'
    assert octave_code(LambertW(x)) == 'lambertw(x)'
    assert octave_code(LambertW(x, n)) == 'lambertw(n, x)'
    assert octave_code(zeta(x)) == 'zeta(x)'
    assert octave_code(zeta(x, y)) == '% Not supported in Octave:\n% zeta\nzeta(x, y)'
