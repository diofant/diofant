from diofant.core import (S, pi, oo, symbols, Function,
                          Rational, Integer, Tuple, Derivative,
                          Eq, Ne, Le, Gt)
from diofant.integrals import Integral
from diofant.concrete import Sum
from diofant.functions import (exp, sin, cos, sign, atanh, meijerg, hyper,
                               Min, Max, Piecewise)
from diofant.matrices import Matrix
from diofant.logic import Or, true, false

from diofant import mathematica_code as mcode

x, y, z = symbols('x,y,z')
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


def test_Derivative():
    assert mcode(Derivative(f(x), x, x)) == 'D[f[x], x, x]'


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
    assert mcode(S.NegativeInfinity) == "-Infinity"
    assert mcode(S.EulerGamma) == "EulerGamma"
    assert mcode(S.Catalan) == "Catalan"
    assert mcode(S.Exp1) == "E"


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
    assert mcode(Matrix([[1, 2], [3, 4444]])) == '{{1, 2}, {3, 4444}}'


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
