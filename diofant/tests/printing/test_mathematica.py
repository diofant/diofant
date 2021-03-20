"""Mathematica code printing tests."""

import pytest

from diofant import (QQ, Catalan, Derivative, Dummy, E, Eq, EulerGamma,
                     Function, Gt, Heaviside, Integer, Integral, Lambda, Le,
                     Limit, Matrix, Max, Min, Ne, Or, Piecewise, Poly,
                     Rational, RootOf, RootSum, SparseMatrix, Sum, Tuple, acos,
                     acosh, acot, acoth, asin, asinh, atan, atanh, binomial,
                     conjugate, cos, cosh, cot, coth, csch, erfc, exp,
                     factorial, factorial2, false, fibonacci, gamma, hyper, im,
                     log, loggamma, mathematica_code, meijerg, oo, pi,
                     polygamma, polylog, re, rf, sech, sign, sin, sinh,
                     symbols, tan, tanh, true, zeta)
from diofant.abc import x, y, z


__all__ = ()

f = Function('f')


def test_Integer():
    assert mathematica_code(Integer(67)) == '67'
    assert mathematica_code(Integer(-1)) == '-1'


def test_Rational():
    assert mathematica_code(Rational(3, 7)) == '3/7'
    assert mathematica_code(Rational(18, 9)) == '2'
    assert mathematica_code(Rational(3, -7)) == '-3/7'
    assert mathematica_code(Rational(-3, -7)) == '3/7'
    assert mathematica_code(x + Rational(3, 7)) == 'x + 3/7'
    assert mathematica_code(Rational(3, 7)*x) == '(3/7)*x'


def test_symbols():
    assert mathematica_code(x) == 'x'
    d = Dummy('d')
    assert mathematica_code(d) == f'd{d.dummy_index}'


def test_Function():
    assert mathematica_code(f(x, y, z)) == 'f[x, y, z]'
    assert mathematica_code(sin(x) ** cos(x)) == 'Sin[x]^Cos[x]'
    assert mathematica_code(sign(x)) == 'Sign[x]'

    assert mathematica_code(atanh(x), user_functions={'atanh': 'ArcTanh'}) == 'ArcTanh[x]'

    assert (mathematica_code(meijerg(((1, 1), (3, 4)), ((1,), ()), x)) ==
            'MeijerG[{{1, 1}, {3, 4}}, {{1}, {}}, x]')
    assert (mathematica_code(hyper((1, 2, 3), (3, 4), x)) ==
            'HypergeometricPFQ[{1, 2, 3}, {3, 4}, x]')

    assert mathematica_code(Min(x, y)) == 'Min[x, y]'
    assert mathematica_code(Max(x, y)) == 'Max[x, y]'
    assert mathematica_code(Max(x, 2)) == 'Max[2, x]'  # issue sympy/sympy#15344

    assert mathematica_code(binomial(x, y)) == 'Binomial[x, y]'

    assert mathematica_code(log(x)) == 'Log[x]'
    assert mathematica_code(tan(x)) == 'Tan[x]'
    assert mathematica_code(cot(x)) == 'Cot[x]'
    assert mathematica_code(asin(x)) == 'ArcSin[x]'
    assert mathematica_code(acos(x)) == 'ArcCos[x]'
    assert mathematica_code(atan(x)) == 'ArcTan[x]'
    assert mathematica_code(acot(x)) == 'ArcCot[x]'
    assert mathematica_code(sinh(x)) == 'Sinh[x]'
    assert mathematica_code(cosh(x)) == 'Cosh[x]'
    assert mathematica_code(tanh(x)) == 'Tanh[x]'
    assert mathematica_code(coth(x)) == 'Coth[x]'
    assert mathematica_code(asinh(x)) == 'ArcSinh[x]'
    assert mathematica_code(acosh(x)) == 'ArcCosh[x]'
    assert mathematica_code(atanh(x)) == 'ArcTanh[x]'
    assert mathematica_code(acoth(x)) == 'ArcCoth[x]'
    assert mathematica_code(sech(x)) == 'Sech[x]'
    assert mathematica_code(csch(x)) == 'Csch[x]'
    assert mathematica_code(erfc(x)) == 'Erfc[x]'
    assert mathematica_code(conjugate(x)) == 'Conjugate[x]'
    assert mathematica_code(re(x)) == 'Re[x]'
    assert mathematica_code(im(x)) == 'Im[x]'
    assert mathematica_code(polygamma(x, y)) == 'PolyGamma[x, y]'
    assert mathematica_code(factorial(x)) == 'Factorial[x]'
    assert mathematica_code(factorial2(x)) == 'Factorial2[x]'
    assert mathematica_code(rf(x, y)) == 'Pochhammer[x, y]'
    assert mathematica_code(gamma(x)) == 'Gamma[x]'
    assert mathematica_code(zeta(x)) == 'Zeta[x]'
    assert mathematica_code(Heaviside(x)) == 'UnitStep[x]'
    assert mathematica_code(fibonacci(x)) == 'Fibonacci[x]'
    assert mathematica_code(polylog(x, y)) == 'PolyLog[x, y]'
    assert mathematica_code(loggamma(x)) == 'LogGamma[x]'

    class MyFunc1(Function):
        @classmethod
        def eval(cls, x):
            pass

    class MyFunc2(Function):
        @classmethod
        def eval(cls, x, y):
            pass

    pytest.raises(ValueError,
                  lambda: mathematica_code(MyFunc1(x),
                                           user_functions={'MyFunc1':
                                                           ['Myfunc1']}))
    assert mathematica_code(MyFunc1(x),
                            user_functions={'MyFunc1':
                                            'Myfunc1'}) == 'Myfunc1[x]'
    assert mathematica_code(MyFunc2(x, y),
                            user_functions={'MyFunc2':
                                            [(lambda *x: False,
                                              'Myfunc2')]}) == 'MyFunc2[x, y]'


def test_Lambda():
    f1 = Lambda(x, x**2)
    assert mathematica_code(f1) == 'Function[{x}, x^2]'
    f2 = Lambda((x, y), x + 2*y)
    assert mathematica_code(f2) == 'Function[{x, y}, x + 2*y]'


def test_Derivative():
    assert mathematica_code(Derivative(f(x), x, x)) == 'Hold[D[f[x], x, x]]'
    assert mathematica_code(Derivative(sin(x), x)) == 'Hold[D[Sin[x], x]]'
    assert mathematica_code(Derivative(x, x)) == 'Hold[D[x, x]]'
    assert mathematica_code(Derivative(sin(x)*y**4, (x, 2))) == 'Hold[D[y^4*Sin[x], x, x]]'
    assert mathematica_code(Derivative(sin(x)*y**4, x, y, x)) == 'Hold[D[y^4*Sin[x], x, y, x]]'
    assert mathematica_code(Derivative(sin(x)*y**4, x, (y, 3), x)) == 'Hold[D[y^4*Sin[x], x, y, y, y, x]]'


def test_Pow():
    assert mathematica_code(x**3) == 'x^3'
    assert mathematica_code(x**(y**3)) == 'x^(y^3)'
    assert mathematica_code(1/(f(x)*3.5)**(x - y**x)/(x**2 + y)) == \
        '(3.5*f[x])^(-x + y^x)/(x^2 + y)'
    assert mathematica_code(x**-1.0) == 'x^(-1.0)'
    assert mathematica_code(x**Rational(2, 3)) == 'x^(2/3)'


def test_Mul():
    A, B, C, D = symbols('A B C D', commutative=False)
    assert mathematica_code(x*y*z) == 'x*y*z'
    assert mathematica_code(x*y*A) == 'x*y*A'
    assert mathematica_code(x*y*A*B) == 'x*y*A**B'
    assert mathematica_code(x*y*A*B*C) == 'x*y*A**B**C'
    assert mathematica_code(x*A*B*(C + D)*A*y) == 'x*y*A**B**(C + D)**A'


def test_constants():
    assert mathematica_code(pi) == 'Pi'
    assert mathematica_code(oo) == 'Infinity'
    assert mathematica_code(-oo) == '-Infinity'
    assert mathematica_code(EulerGamma) == 'EulerGamma'
    assert mathematica_code(Catalan) == 'Catalan'
    assert mathematica_code(E) == 'E'


def test_containers():
    assert mathematica_code([1, 2, 3, [4, 5, [6, 7]], 8, [9, 10], 11]) == \
        '{1, 2, 3, {4, 5, {6, 7}}, 8, {9, 10}, 11}'
    assert mathematica_code((1, 2, (3, 4))) == '{1, 2, {3, 4}}'
    assert mathematica_code([1]) == '{1}'
    assert mathematica_code((1,)) == '{1}'
    assert mathematica_code(Tuple(*[1, 2, 3])) == '{1, 2, 3}'


def test_Integral():
    assert mathematica_code(Integral(sin(sin(x)), x)) == 'Hold[Integrate[Sin[Sin[x]], x]]'
    assert mathematica_code(Integral(exp(-x**2 - y**2),
                                     (x, -oo, oo),
                                     (y, -oo, oo))) == \
        'Hold[Integrate[E^(-x^2 - y^2), {x, -Infinity, Infinity}, ' \
        '{y, -Infinity, Infinity}]]'


def test_Sum():
    assert mathematica_code(Sum(sin(x), (x, 0, 10))) == 'Hold[Sum[Sin[x], {x, 0, 10}]]'
    assert mathematica_code(Sum(exp(-x**2 - y**2),
                                (x, -oo, oo),
                                (y, -oo, oo))) == \
        'Hold[Sum[E^(-x^2 - y^2), {x, -Infinity, Infinity}, ' \
        '{y, -Infinity, Infinity}]]'


def test_Matrix():
    assert mathematica_code(Matrix()) == '{}'

    m = Matrix([[1, 2], [3, 4444]])
    assert mathematica_code(m) == mathematica_code(m.as_immutable()) == '{{1, 2}, {3, 4444}}'

    m = SparseMatrix(m)
    assert mathematica_code(m) == mathematica_code(m.as_immutable()) == '{{1, 2}, {3, 4444}}'


def test_Relational():
    assert mathematica_code(Eq(x, y)) == 'x == y'
    assert mathematica_code(Ne(x, y/(1 + y**2))) == 'x != (y/(y^2 + 1))'
    assert mathematica_code(Le(0, x**2)) == '0 <= x^2'
    assert mathematica_code(Gt(pi, 3, evaluate=False)) == 'Pi > 3'


def test_Booleans():
    assert mathematica_code(true) == 'True'
    assert mathematica_code(false) == 'False'


def test_Piecewise():
    g = Piecewise((0, Or(x <= -1, x >= 1)), (1 - x, x > 0), (1 + x, True))

    assert (mathematica_code(g) ==
            'Piecewise[{{0, x >= 1 || x <= -1}, '
            '{-x + 1, x > 0}, {x + 1, True}}]')


def test_RootOf():
    p = Poly(x**3 + y*x + 1, x)
    assert mathematica_code(RootOf(p, 0)) == 'Root[#^3 + #*y + 1 &, 1]'


def test_RootSum():
    r = RootSum(x**3 + x + 3, Lambda(y, log(y*z)))
    assert mathematica_code(r) == ('RootSum[Function[{x}, x^3 + x + 3], '
                                   'Function[{y}, Log[y*z]]]')


def test_AlgebraicElement():
    r = RootOf(x**7 + 3*x - 1, 3)
    K = QQ.algebraic_field(r)
    a = K([1, 0, 3, 2, 1])
    assert mathematica_code(a) == ('AlgebraicNumber[Root[#^7 + 3*# - 1 &, 4],'
                                   ' {1, 0, 3, 2, 1}]')


def test_Limit():
    e = Limit(sin(x)/x, x, 0)
    assert mathematica_code(e) == 'Hold[Limit[Sin[x]/x, x -> 0, Direction -> -1]]'
    e = Limit(sin(x)/x, x, 0, '-')
    assert mathematica_code(e) == 'Hold[Limit[Sin[x]/x, x -> 0, Direction -> 1]]'
    e = Limit(sin(x)/x, x, 0, 'real')
    assert mathematica_code(e) == 'Hold[Limit[Sin[x]/x, x -> 0, Direction -> Reals]]'
