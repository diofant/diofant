import pytest

from diofant import QQ
from diofant import mathematica_code as mcode
from diofant.abc import x, y, z
from diofant.concrete import Sum
from diofant.core import (Catalan, Derivative, Dummy, E, Eq, EulerGamma,
                          Function, Gt, Integer, Lambda, Le, Ne, Rational,
                          Tuple, oo, pi, symbols)
from diofant.functions import (Heaviside, Max, Min, Piecewise, acos, asin,
                               asinh, atan, atanh, binomial, conjugate, cos,
                               cosh, cot, coth, csch, erfc, exp, factorial,
                               factorial2, fibonacci, gamma, hyper, im, log,
                               meijerg, polygamma, polylog, re, rf, sech, sign,
                               sin, sinh, tan, tanh, zeta)
from diofant.integrals import Integral
from diofant.logic import Or, false, true
from diofant.matrices import Matrix, SparseMatrix
from diofant.polys import Poly, RootOf, RootSum
from diofant.series import Limit


__all__ = ()

f = Function('f')


def test_Integer():
    assert mcode(Integer(67)) == "67"
    assert mcode(Integer(-1)) == "-1"


def test_Rational():
    assert mcode(Rational(3, 7)) == "3/7"
    assert mcode(Rational(18, 9)) == "2"
    assert mcode(Rational(3, -7)) == "-3/7"
    assert mcode(Rational(-3, -7)) == "3/7"
    assert mcode(x + Rational(3, 7)) == "x + 3/7"
    assert mcode(Rational(3, 7)*x) == "(3/7)*x"


def test_symbols():
    assert mcode(x) == "x"
    d = Dummy("d")
    assert mcode(d) == "d%s" % d.dummy_index


def test_Function():
    assert mcode(f(x, y, z)) == "f[x, y, z]"
    assert mcode(sin(x) ** cos(x)) == "Sin[x]^Cos[x]"
    assert mcode(sign(x)) == "Sign[x]"

    assert mcode(atanh(x), user_functions={"atanh": "ArcTanh"}) == "ArcTanh[x]"

    assert (mcode(meijerg(((1, 1), (3, 4)), ((1,), ()), x)) ==
            "MeijerG[{{1, 1}, {3, 4}}, {{1}, {}}, x]")
    assert (mcode(hyper((1, 2, 3), (3, 4), x)) ==
            "HypergeometricPFQ[{1, 2, 3}, {3, 4}, x]")

    assert mcode(Min(x, y)) == "Min[x, y]"
    assert mcode(Max(x, y)) == "Max[x, y]"
    assert mcode(Max(x, 2)) == "Max[2, x]"  # issue sympy/sympy#15344

    assert mcode(binomial(x, y)) == "Binomial[x, y]"

    assert mcode(log(x)) == "Log[x]"
    assert mcode(tan(x)) == "Tan[x]"
    assert mcode(cot(x)) == "Cot[x]"
    assert mcode(asin(x)) == "ArcSin[x]"
    assert mcode(acos(x)) == "ArcCos[x]"
    assert mcode(atan(x)) == "ArcTan[x]"
    assert mcode(sinh(x)) == "Sinh[x]"
    assert mcode(cosh(x)) == "Cosh[x]"
    assert mcode(tanh(x)) == "Tanh[x]"
    assert mcode(coth(x)) == "Coth[x]"
    assert mcode(sech(x)) == "Sech[x]"
    assert mcode(csch(x)) == "Csch[x]"
    assert mcode(erfc(x)) == "Erfc[x]"
    assert mcode(conjugate(x)) == "Conjugate[x]"
    assert mcode(re(x)) == "Re[x]"
    assert mcode(im(x)) == "Im[x]"
    assert mcode(polygamma(x, y)) == "PolyGamma[x, y]"
    assert mcode(factorial(x)) == "Factorial[x]"
    assert mcode(factorial2(x)) == "Factorial2[x]"
    assert mcode(rf(x, y)) == "Pochhammer[x, y]"
    assert mcode(gamma(x)) == "Gamma[x]"
    assert mcode(zeta(x)) == "Zeta[x]"
    assert mcode(asinh(x)) == "ArcSinh[x]"
    assert mcode(Heaviside(x)) == "UnitStep[x]"
    assert mcode(fibonacci(x)) == "Fibonacci[x]"
    assert mcode(polylog(x, y)) == "PolyLog[x, y]"
    assert mcode(atanh(x)) == "ArcTanh[x]"

    class myfunc1(Function):
        @classmethod
        def eval(cls, x):
            pass

    class myfunc2(Function):
        @classmethod
        def eval(cls, x, y):
            pass

    pytest.raises(ValueError,
                  lambda: mcode(myfunc1(x),
                                user_functions={"myfunc1": ["Myfunc1"]}))
    assert mcode(myfunc1(x),
                 user_functions={"myfunc1": "Myfunc1"}) == "Myfunc1[x]"
    assert mcode(myfunc2(x, y),
                 user_functions={"myfunc2": [(lambda *x: False,
                                              "Myfunc2")]}) == "myfunc2[x, y]"


def test_Lambda():
    f1 = Lambda(x, x**2)
    assert mcode(f1) == "Function[{x}, x^2]"
    f2 = Lambda((x, y), x + 2*y)
    assert mcode(f2) == "Function[{x, y}, x + 2*y]"


def test_Derivative():
    assert mcode(Derivative(f(x), x, x)) == 'Hold[D[f[x], x, x]]'
    assert mcode(Derivative(sin(x), x)) == "Hold[D[Sin[x], x]]"
    assert mcode(Derivative(x, x)) == "Hold[D[x, x]]"
    assert mcode(Derivative(sin(x)*y**4, x, 2)) == "Hold[D[y^4*Sin[x], x, x]]"
    assert mcode(Derivative(sin(x)*y**4, x, y, x)) == "Hold[D[y^4*Sin[x], x, y, x]]"
    assert mcode(Derivative(sin(x)*y**4, x, y, 3, x)) == "Hold[D[y^4*Sin[x], x, y, y, y, x]]"


def test_Pow():
    assert mcode(x**3) == "x^3"
    assert mcode(x**(y**3)) == "x^(y^3)"
    assert mcode(1/(f(x)*3.5)**(x - y**x)/(x**2 + y)) == \
        "(3.5*f[x])^(-x + y^x)/(x^2 + y)"
    assert mcode(x**-1.0) == 'x^(-1.0)'
    assert mcode(x**Rational(2, 3)) == 'x^(2/3)'


def test_Mul():
    A, B, C, D = symbols('A B C D', commutative=False)
    assert mcode(x*y*z) == "x*y*z"
    assert mcode(x*y*A) == "x*y*A"
    assert mcode(x*y*A*B) == "x*y*A**B"
    assert mcode(x*y*A*B*C) == "x*y*A**B**C"
    assert mcode(x*A*B*(C + D)*A*y) == "x*y*A**B**(C + D)**A"


def test_constants():
    assert mcode(pi) == "Pi"
    assert mcode(oo) == "Infinity"
    assert mcode(-oo) == "-Infinity"
    assert mcode(EulerGamma) == "EulerGamma"
    assert mcode(Catalan) == "Catalan"
    assert mcode(E) == "E"


def test_containers():
    assert mcode([1, 2, 3, [4, 5, [6, 7]], 8, [9, 10], 11]) == \
        "{1, 2, 3, {4, 5, {6, 7}}, 8, {9, 10}, 11}"
    assert mcode((1, 2, (3, 4))) == "{1, 2, {3, 4}}"
    assert mcode([1]) == "{1}"
    assert mcode((1,)) == "{1}"
    assert mcode(Tuple(*[1, 2, 3])) == "{1, 2, 3}"


def test_Integral():
    assert mcode(Integral(sin(sin(x)), x)) == "Hold[Integrate[Sin[Sin[x]], x]]"
    assert mcode(Integral(exp(-x**2 - y**2),
                          (x, -oo, oo),
                          (y, -oo, oo))) == \
        "Hold[Integrate[E^(-x^2 - y^2), {x, -Infinity, Infinity}, " \
        "{y, -Infinity, Infinity}]]"


def test_Sum():
    assert mcode(Sum(sin(x), (x, 0, 10))) == "Hold[Sum[Sin[x], {x, 0, 10}]]"
    assert mcode(Sum(exp(-x**2 - y**2),
                     (x, -oo, oo),
                     (y, -oo, oo))) == \
        "Hold[Sum[E^(-x^2 - y^2), {x, -Infinity, Infinity}, " \
        "{y, -Infinity, Infinity}]]"


def test_Matrix():
    assert mcode(Matrix()) == '{}'

    m = Matrix([[1, 2], [3, 4444]])
    assert mcode(m) == mcode(m.as_immutable()) == '{{1, 2}, {3, 4444}}'

    m = SparseMatrix(m)
    assert mcode(m) == mcode(m.as_immutable()) == '{{1, 2}, {3, 4444}}'


def test_Relational():
    assert mcode(Eq(x, y)) == 'x == y'
    assert mcode(Ne(x, y/(1 + y**2))) == 'x != y/(y^2 + 1)'
    assert mcode(Le(0, x**2)) == '0 <= x^2'
    assert mcode(Gt(pi, 3, evaluate=False)) == 'Pi > 3'


def test_Booleans():
    assert mcode(true) == "True"
    assert mcode(false) == "False"


def test_Piecewise():
    g = Piecewise((0, Or(x <= -1, x >= 1)), (1 - x, x > 0), (1 + x, True))

    assert (mcode(g) ==
            'Piecewise[{{0, x >= 1 || x <= -1}, '
            '{-x + 1, x > 0}, {x + 1, True}}]')


def test_RootOf():
    p = Poly(x**3 + y*x + 1, x)
    assert mcode(RootOf(p, 0)) == 'Root[#^3 + #*y + 1 &, 1]'


def test_RootSum():
    r = RootSum(x**3 + x + 3, Lambda(y, log(y*z)))
    assert mcode(r) == ("RootSum[Function[{x}, x^3 + x + 3], "
                        "Function[{y}, Log[y*z]]]")


def test_AlgebraicElement():
    r = RootOf(x**7 + 3*x - 1, 3)
    K = QQ.algebraic_field(r)
    a = K([1, 2, 3, 0, 1])
    assert mcode(a) == ('AlgebraicNumber[Root[#^7 + 3*# - 1 &, 4],'
                        ' {1, 0, 3, 2, 1}]')


def test_Limit():
    e = Limit(sin(x)/x, x, 0)
    assert mcode(e) == "Hold[Limit[Sin[x]/x, x -> 0, Direction -> -1]]"
    e = Limit(sin(x)/x, x, 0, "-")
    assert mcode(e) == "Hold[Limit[Sin[x]/x, x -> 0, Direction -> 1]]"
    e = Limit(sin(x)/x, x, 0, "real")
    assert mcode(e) == "Hold[Limit[Sin[x]/x, x -> 0, Direction -> Reals]]"
