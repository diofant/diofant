import pytest

from diofant import (Add, And, Ci, Derivative, DiracDelta, E, Eq, EulerGamma,
                     Expr, Function, I, Integral, Interval, Lambda, LambertW,
                     Matrix, Max, Min, Mul, Ne, O, Piecewise, Poly, Rational,
                     Si, Sum, Symbol, Tuple, acos, acosh, asin, asinh, atan,
                     cbrt, cos, cosh, diff, erf, erfi, exp, expand_func,
                     expand_mul, factor, fresnels, gamma, im, integrate, log,
                     lowergamma, meijerg, nan, oo, pi, polar_lift, polygamma,
                     re, sign, simplify, sin, sinh, sqrt, sstr, symbols,
                     sympify, tan, tanh, trigsimp)
from diofant.abc import A, L, R, a, b, c, h, i, k, m, s, t, w, x, y, z
from diofant.functions.elementary.complexes import periodic_argument
from diofant.integrals.heurisch import heurisch
from diofant.integrals.risch import NonElementaryIntegral
from diofant.utilities.randtest import verify_numerically


__all__ = ()

x_1, x_2 = symbols('x_1 x_2')
n = Symbol('n', integer=True)
f = Function('f')


def diff_test(i):
    """Return the set of symbols, s, which were used in testing that
    i.diff(s) agrees with i.doit().diff(s). If there is an error then
    the assertion will fail, causing the test to fail.
    """
    syms = i.free_symbols
    for s in syms:
        assert (i.diff(s).doit() - i.doit().diff(s)).expand() == 0
    return syms


def test_improper_integral():
    assert integrate(log(x), (x, 0, 1)) == -1
    assert integrate(x**(-2), (x, 1, oo)) == 1

    # issue sympy/sympy#10445:
    assert integrate(1/(1 + exp(x)), (x, 0, oo)) == log(2)

    # issue sympy/sympy#8945
    assert integrate(sin(x)**3/x, (x, 0, 1)) == -Si(3)/4 + 3*Si(1)/4
    assert integrate(sin(x)**3/x, (x, 0, oo)) == pi/4

    # issue sympy/sympy#4527
    k, m = symbols('k m', integer=True)
    assert integrate(sin(k*x)*sin(m*x), (x, 0, pi)) == Piecewise(
        (0, And(Eq(k, 0), Eq(m, 0))),
        (-pi/2, Eq(k, -m)),
        (pi/2, Eq(k, m)),
        (0, True))

    # issue sympy/sympy#10211
    assert integrate((1/sqrt(((y - x)**2 + h**2))**3),
                     (x, 0, w), (y, 0, w)) == 2*sqrt(1 + w**2/h**2)/h - 2/h


def test_constructor():
    # this is shared by Sum, so testing Integral's constructor
    # is equivalent to testing Sum's
    s1 = Integral(n, n)
    assert s1.limits == (Tuple(n),)
    s2 = Integral(n, (n,))
    assert s2.limits == (Tuple(n),)
    s3 = Integral(Sum(x, (x, 1, y)))
    assert s3.limits == (Tuple(y),)
    s4 = Integral(n, Tuple(n,))
    assert s4.limits == (Tuple(n),)

    s5 = Integral(n, (n, Interval(1, 2)))
    assert s5.limits == (Tuple(n, 1, 2),)


def test_basics():
    assert Integral(0, x) != 0
    assert Integral(x, (x, 1, 1)) != 0
    assert Integral(oo, x) != oo
    assert Integral(nan, x) == nan

    assert diff(Integral(y, y), x) == 0
    assert diff(Integral(x, (x, 0, 1)), x) == 0
    assert diff(Integral(x, x), x) == x
    assert diff(Integral(t, (t, 0, x)), x) == x + Integral(0, (t, 0, x))

    e = (t + 1)**2
    assert diff(integrate(e, (t, 0, x)), x) == \
        diff(Integral(e, (t, 0, x)), x).doit().expand() == \
        ((1 + x)**2).expand()
    assert diff(integrate(e, (t, 0, x)), t) == \
        diff(Integral(e, (t, 0, x)), t) == 0
    assert diff(integrate(e, (t, 0, x)), a) == \
        diff(Integral(e, (t, 0, x)), a) == 0
    assert diff(integrate(e, t), a) == diff(Integral(e, t), a) == 0

    assert integrate(e, (t, a, x)).diff(x) == \
        Integral(e, (t, a, x)).diff(x).doit().expand()
    assert Integral(e, (t, a, x)).diff(x).doit() == ((1 + x)**2)
    assert integrate(e, (t, x, a)).diff(x).doit() == (-(1 + x)**2).expand()

    assert integrate(t**2, (t, x, 2*x)).diff(x) == 7*x**2

    assert Integral(x, x).atoms() == {x}
    assert Integral(f(x), (x, 0, 1)).atoms() == {0, 1, x}

    assert diff_test(Integral(x, (x, 3*y))) == {y}
    assert diff_test(Integral(x, (a, 3*y))) == {x, y}

    assert integrate(x, (x, oo, oo)) == 0  # issue sympy/sympy#8171
    assert integrate(x, (x, -oo, -oo)) == 0

    # sum integral of terms
    assert integrate(y + x + exp(x), x) == x*y + x**2/2 + exp(x)

    assert Integral(x).is_commutative
    n = Symbol('n', commutative=False)
    assert Integral(n + x, x).is_commutative is False


def test_diff_wrt():
    class Test(Expr):
        _diff_wrt = True
        is_commutative = True

    t = Test()
    assert integrate(t + 1, t) == t**2/2 + t
    assert integrate(t + 1, (t, 0, 1)) == Rational(3, 2)

    pytest.raises(ValueError, lambda: integrate(x + 1, x + 1))
    pytest.raises(ValueError, lambda: integrate(x + 1, (x + 1, 0, 1)))


def test_basics_multiple():
    assert diff_test(Integral(x, (x, 3*x, 5*y), (y, x, 2*x))) == {x}
    assert diff_test(Integral(x, (x, 5*y), (y, x, 2*x))) == {x}
    assert diff_test(Integral(x, (x, 5*y), (y, y, 2*x))) == {x, y}
    assert diff_test(Integral(y, y, x)) == {x, y}
    assert diff_test(Integral(y*x, x, y)) == {x, y}
    assert diff_test(Integral(x + y, y, (y, 1, x))) == {x}
    assert diff_test(Integral(x + y, (x, x, y), (y, y, x))) == {x, y}


def test_conjugate_transpose():
    A, B = symbols('A B', commutative=False)

    x = Symbol('x', complex=True)
    p = Integral(A*B, (x,))
    assert p.adjoint().doit() == p.doit().adjoint()
    assert p.conjugate().doit() == p.doit().conjugate()
    assert p.transpose().doit() == p.doit().transpose()

    x = Symbol('x', extended_real=True)
    p = Integral(A*B, (x,))
    assert p.adjoint().doit() == p.doit().adjoint()
    assert p.conjugate().doit() == p.doit().conjugate()
    assert p.transpose().doit() == p.doit().transpose()


def test_integration():
    assert integrate(0, (t, 0, x)) == 0
    assert integrate(3, (t, 0, x)) == 3*x
    assert integrate(t, (t, 0, x)) == x**2/2
    assert integrate(3*t, (t, 0, x)) == 3*x**2/2
    assert integrate(3*t**2, (t, 0, x)) == x**3
    assert integrate(1/t, (t, 1, x)) == log(x)
    assert integrate(-1/t**2, (t, 1, x)) == 1/x - 1
    assert integrate(t**2 + 5*t - 8, (t, 0, x)) == x**3/3 + 5*x**2/2 - 8*x
    assert integrate(x**2, x) == x**3/3
    assert integrate((3*t*x)**5, x) == (3*t)**5 * x**6 / 6

    b = Symbol('b')
    c = Symbol('c')
    assert integrate(a*t, (t, 0, x)) == a*x**2/2
    assert integrate(a*t**4, (t, 0, x)) == a*x**5/5
    assert integrate(a*t**2 + b*t + c, (t, 0, x)) == a*x**3/3 + b*x**2/2 + c*x

    # issue sympy/sympy#6253
    # Note: this used to raise NotImplementedError
    # Note: psi in _check_antecedents becomes NaN.
    assert integrate((sqrt(1 - x) + sqrt(1 + x))**2/x, x, meijerg=True) == \
        Integral((sqrt(-x + 1) + sqrt(x + 1))**2/x, x)

    # issue sympy/sympy#8945
    assert integrate(cos(x)**2/x**2, x) == -Si(2*x) - cos(2*x)/(2*x) - 1/(2*x)

    # issue sympy/sympy#10680
    integrate(x**log(x**log(x**log(x))), x)  # not raises

    # issue sympy/sympy#4890
    assert integrate(exp(-log(x)**2), x) == \
        sqrt(pi)*exp(Rational(1, 4))*erf(log(x) - Rational(1, 2))/2
    assert integrate(exp(log(x)**2), x) == \
        sqrt(pi)*exp(-Rational(1, 4))*erfi(log(x) + Rational(1, 2))/2

    # issue sympy/sympy#8901
    assert integrate(tanh(x)) == x - log(tanh(x) + 1)

    # issue sympy/sympy#4403
    z = Symbol('z', positive=True)
    assert integrate(sqrt(x**2 + z**2), x) == \
        z**2*asinh(x/z)/2 + x*sqrt(x**2 + z**2)/2
    assert integrate(sqrt(x**2 - z**2), x) == \
        -z**2*acosh(x/z)/2 + x*sqrt(x**2 - z**2)/2
    assert integrate(sqrt(-x**2 - 4), x) == \
        -2*atan(x/sqrt(-4 - x**2)) + x*sqrt(-4 - x**2)/2


def test_multiple_integration():
    assert integrate((x**2)*(y**2), (x, 0, 1), (y, -1, 2)) == 1
    assert integrate((y**2)*(x**2), x, y) == Rational(1, 9)*(x**3)*(y**3)
    assert integrate(1/(x + 3)/(1 + x)**3, x) == \
        -Rational(1, 8)*log(3 + x) + Rational(1, 8)*log(1 + x) + x/(4 + 8*x + 4*x**2)

    # issue sympy/sympy#5178
    assert (integrate(sin(x)*f(y, z), (x, 0, pi), (y, 0, pi), (z, 0, pi)) ==
            2*Integral(f(y, z), (y, 0, pi), (z, 0, pi)))

    # issue sympy/sympy#5167
    assert Integral(f(x), (x, 1, 2), (w, 1, x), (z, 1, y)).doit() == \
        y*(x - 1)*Integral(f(x), (x, 1, 2)) - (x - 1)*Integral(f(x), (x, 1, 2))

    # issue sympy/sympy#14782
    assert integrate(sqrt(-x**2 + 1)*(-x**2 + x), (x, -1, 1)) != 0


def test_sympyissue_3532():
    assert integrate(exp(-x), (x, 0, oo)) == 1


def test_sympyissue_3560():
    assert integrate(sqrt(x)**3, x) == 2*sqrt(x)**5/5
    assert integrate(sqrt(x), x) == 2*sqrt(x)**3/3
    assert integrate(1/sqrt(x)**3, x) == -2/sqrt(x)


def test_integrate_poly():
    p = (x + x**2*y + y**3).as_poly()

    qx = integrate(p, x)
    qy = integrate(p, y)

    assert isinstance(qx, Poly) is True
    assert isinstance(qy, Poly) is True

    assert qx.gens == (x, y)
    assert qy.gens == (x, y)

    assert qx.as_expr() == x**2/2 + x**3*y/3 + x*y**3
    assert qy.as_expr() == x*y + x**2*y**2/2 + y**4/4


def test_integrate_poly_defined():
    p = (x + x**2*y + y**3).as_poly()

    Qx = integrate(p, (x, 0, 1))
    Qy = integrate(p, (y, 0, pi))

    assert isinstance(Qx, Poly) is True
    assert isinstance(Qy, Poly) is True

    assert Qx.gens == (y,)
    assert Qy.gens == (x,)

    assert Qx.as_expr() == Rational(1, 2) + y/3 + y**3
    assert Qy.as_expr() == pi**4/4 + pi*x + pi**2*x**2/2


def test_integrate_omit_var():
    assert integrate(x) == x**2/2

    pytest.raises(ValueError, lambda: integrate(2))
    pytest.raises(ValueError, lambda: integrate(x*y))


def test_integrate_poly_accurately():
    assert integrate(x*sin(y), x) == x**2*sin(y)/2

    # when passed to risch_norman, this will be a CPU hog, so this really
    # checks, that integrated function is recognized as polynomial
    assert integrate(x**1000*sin(y), x) == x**1001*sin(y)/1001


def test_sympyissue_3635():
    assert integrate(x**2, y) == x**2*y
    assert integrate(x**2, (y, -1, 1)) == 2*x**2


def test_integrate_linearterm_pow():
    # check integrate((a*x+b)^c, x)  --  issue sympy/sympy#3499
    y = Symbol('y', positive=True)
    # TODO: Remove conds='none' below, let the assumption take care of it.
    assert integrate(x**y, x, conds='none') == x**(y + 1)/(y + 1)
    assert integrate((exp(y)*x + 1/y)**(1 + sin(y)), x, conds='none') == \
        exp(-y)*(exp(y)*x + 1/y)**(2 + sin(y)) / (2 + sin(y))


def test_sympyissue_3618():
    assert integrate(pi*sqrt(x), x) == 2*pi*sqrt(x)**3/3
    assert integrate(pi*sqrt(x) + E*sqrt(x)**3, x) == \
        2*pi*sqrt(x)**3/3 + 2*E * sqrt(x)**5/5


def test_sympyissue_3623():
    assert integrate(cos((n + 1)*x), x) == Piecewise(
        (x, Eq(n + 1, 0)), (sin((n + 1)*x)/(n + 1), True))
    assert integrate(cos((n - 1)*x), x) == Piecewise(
        (x, Eq(n - 1, 0)), (sin((n - 1)*x)/(n - 1), True))
    assert integrate(cos((n + 1)*x) + cos((n - 1)*x), x) == \
        Piecewise((x, Eq(n + 1, 0)), (sin((n + 1)*x)/(n + 1), True)) + \
        Piecewise((x, Eq(n - 1, 0)), (sin((n - 1)*x)/(n - 1), True))


def test_sympyissue_3664():
    n = Symbol('n', integer=True, nonzero=True)
    assert integrate(-1./2 * x * sin(n * pi * x/2), [x, -2, 0]) == \
        2*cos(pi*n)/(pi*n)
    assert integrate(-x*sin(n*pi*x/2)/2, [x, -2, 0]) == \
        2*cos(pi*n)/(pi*n)


def test_sympyissue_3679():
    # definite integration of rational functions gives wrong answers
    assert NS(Integral(1/(x**2 - 8*x + 17), (x, 2, 4))) == '1.10714871779409'


def test_sympyissue_3686():  # remove this when fresnel itegrals are implemented
    assert expand_func(integrate(sin(x**2), x)) == \
        sqrt(2)*sqrt(pi)*fresnels(sqrt(2)*x/sqrt(pi))/2


def test_transcendental_functions():
    assert integrate(LambertW(2*x), x) == \
        -x + x*LambertW(2*x) + x/LambertW(2*x)


def test_sympyissue_3740():
    f = 4*log(x) - 2*log(x)**2
    fid = diff(integrate(f, x), x)
    assert abs(f.subs({x: 42}).evalf() - fid.subs({x: 42}).evalf()) < 1e-10


def test_sympyissue_3788():
    assert integrate(1/(1 + x**2), x) == atan(x)


def test_sympyissue_3952():
    f = sin(x)
    assert integrate(f, x) == -cos(x)
    pytest.raises(ValueError, lambda: integrate(f, 2*x))


def test_sympyissue_4516():
    assert integrate(2**x - 2*x, x) == 2**x/log(2) - x**2


def test_sympyissue_7450():
    ans = integrate(exp(-(1 + I)*x), (x, 0, oo))
    assert re(ans) == Rational(1, 2) and im(ans) == Rational(-1, 2)


def test_matrices():
    M = Matrix(2, 2, lambda i, j: (i + j + 1)*sin((i + j + 1)*x))

    assert integrate(M, x) == Matrix([
        [-cos(x), -cos(2*x)],
        [-cos(2*x), -cos(3*x)],
    ])


def test_integrate_functions():
    # issue sympy/sympy#4111
    assert integrate(f(x), x) == Integral(f(x), x)
    assert integrate(f(x), (x, 0, 1)) == Integral(f(x), (x, 0, 1))

    assert integrate(Derivative(f(y), y), x) == x*Derivative(f(y), y)


@pytest.mark.xfail
def test_integrate_functions_1():
    assert integrate(f(x)*diff(f(x), x), x) == f(x)**2/2
    assert integrate(diff(f(x), x)/f(x), x) == log(f(x))


@pytest.mark.xfail
def test_integrate_derivatives():
    assert integrate(Derivative(f(x), x), x) == f(x)


def test_transform():
    a = Integral(x**2 + 1, (x, -1, 2))
    fx = x
    fy = 3*y + 1
    assert a.doit() == a.transform(fx, fy).doit()
    assert a.transform(fx, fy).transform(fy, fx) == a
    fx = 3*x + 1
    fy = y
    assert a.transform(fx, fy).transform(fy, fx) == a
    a = Integral(sin(1/x), (x, 0, 1))
    assert a.transform(x, 1/y) == Integral(sin(y)/y**2, (y, 1, oo))
    assert a.transform(x, 1/y).transform(y, 1/x) == a
    a = Integral(exp(-x**2), (x, -oo, oo))
    assert a.transform(x, 2*y) == Integral(2*exp(-4*y**2), (y, -oo, oo))
    # < 3 arg limit handled properly
    assert Integral(x, x).transform(x, a*y).doit() == \
        Integral(y*a**2, y).doit()
    assert Integral(x, (x, 0, -3)).transform(x, 1/y).doit() == \
        Integral(-1/x**3, (x, -oo, Rational(-1, 3))).doit()
    assert Integral(x, (x, 0, 3)).transform(x, 1/y) == \
        Integral(y**(-3), (y, Rational(1, 3), oo))
    # issue sympy/sympy#8400
    i = Integral(x + y, (x, 1, 2), (y, 1, 2))
    assert i.transform(x, (x + 2*y, x)).doit() == \
        i.transform(x, (x + 2*z, x)).doit() == 3

    assert a.transform(b, x) == a
    assert a.transform(x, y) == Integral(exp(-y**2), (y, -oo, oo))

    i2 = Integral(cos(x**2 - 1), (x, 0, y))
    i = i2.subs({y: 1})

    pytest.raises(ValueError, lambda: i.transform(x**2 - 1, y))
    pytest.raises(ValueError, lambda: i.transform(x, y*z))
    pytest.raises(ValueError, lambda: i.transform(x, (y, y + z)))
    pytest.raises(ValueError, lambda: i2.transform(x, (z*y, y)))
    pytest.raises(ValueError, lambda: i.transform(x, (sin(y), y)))


def test_sympyissue_4052():
    f = Rational(1, 2)*asin(x) + x*sqrt(1 - x**2)/2

    assert integrate(cos(asin(x)), x) == f
    assert integrate(sin(acos(x)), x) == f


def NS(e, n=15, **options):
    return sstr(sympify(e).evalf(n, **options), full_prec=True)


def test_evalf_integrals():
    assert NS(Integral(cos(x)/x, (x, 1, oo)), quad='osc') == '-0.337403922900968'

    pytest.raises(ValueError,
                  lambda: NS(Integral(acos(x)/x, (x, 1, oo)), quad='osc'))

    assert NS(Integral(sin(x + I), (x, 0, pi/2))) == '1.54308063481524 + 1.17520119364380*I'

    assert Integral(pi, (x, y, z)).evalf() == Integral(pi, (x, y, z))
    assert Integral(pi, (x, y, y + z)).evalf() == Integral(pi, (x, y, y + z))
    #
    # Endpoints causing trouble (rounding error in integration points -> complex log)
    assert NS(
        2 + Integral(log(2*cos(x/2)), (x, -pi, pi)), 17, chop=True) == NS(2, 17)
    assert NS(
        2 + Integral(log(2*cos(x/2)), (x, -pi, pi)), 20, chop=True) == NS(2, 20)
    assert NS(
        2 + Integral(log(2*cos(x/2)), (x, -pi, pi)), 22, chop=True) == NS(2, 22)

    # issue sympy/sympy#4038

    # The output form of an integral may differ by a step function between
    # revisions, making this test a bit useless. This can't be said about
    # other two tests. For now, all values of this evaluation are used here,
    # but in future this should be reconsidered.
    assert NS(integrate(1/(x**5 + 1), x).subs({x: 4}), chop=True) in \
        ['-0.000976138910649103', '0.965906660135753', '1.93278945918216']

    assert NS(Integral(1/(x**5 + 1), (x, 2, 4))) == '0.0144361088886740'
    assert NS(
        integrate(1/(x**5 + 1), (x, 2, 4)), chop=True) == '0.0144361088886740'


@pytest.mark.slow
def test_evalf_integrals_slow():
    assert NS(Integral(x, (x, 2, 5)), 15) == '10.5000000000000'
    gauss = Integral(exp(-x**2), (x, -oo, oo))
    assert NS(gauss, 15) == '1.77245385090552'
    assert NS(gauss**2 - pi + E*Rational(
        1, 10**20), 15) in ('2.71828182845904e-20', '2.71828182845905e-20')
    # A monster of an integral from https://mathworld.wolfram.com/DefiniteIntegral.html
    a = 8*sqrt(3)/(1 + 3*t**2)
    b = 16*sqrt(2)*(3*t + 1)*sqrt(4*t**2 + t + 1)**3
    c = (3*t**2 + 1)*(11*t**2 + 2*t + 3)**2
    d = sqrt(2)*(249*t**2 + 54*t + 65)/(11*t**2 + 2*t + 3)**2
    f = a - b/c - d
    assert NS(Integral(f, (t, 0, 1)), 50) == \
        NS((3*sqrt(2) - 49*pi + 162*atan(sqrt(2)))/12, 50)
    # https://mathworld.wolfram.com/VardisIntegral.html
    assert NS(Integral(log(log(1/x))/(1 + x + x**2), (x, 0, 1)), 15) == \
        NS('pi/sqrt(3) * log(2*pi**(5/6) / gamma(1/6))', 15)
    # https://mathworld.wolfram.com/AhmedsIntegral.html
    assert NS(Integral(atan(sqrt(x**2 + 2))/(sqrt(x**2 + 2)*(x**2 + 1)), (x,
                                                                          0, 1)), 15) == NS(5*pi**2/96, 15)
    # https://mathworld.wolfram.com/AbelsIntegral.html
    assert NS(Integral(x/((exp(pi*x) - exp(
        -pi*x))*(x**2 + 1)), (x, 0, oo)), 15) == NS('log(2)/2-1/4', 15)
    # Complex part trimming
    # https://mathworld.wolfram.com/VardisIntegral.html
    assert NS(Integral(log(log(sin(x)/cos(x))), (x, pi/4, pi/2)), 15, chop=True) == \
        NS('pi/4*log(4*pi**3/gamma(1/4)**4)', 15)
    # Needs zero handling
    assert NS(pi - 4*Integral(sqrt(1 - x**2), (x, 0, 1)),
              15, maxn=30, chop=True, strict=False) in ('0.0', '0')
    # Oscillatory quadrature
    a = Integral(sin(x)/x**2, (x, 1, oo)).evalf(maxn=15, strict=False)
    assert 0.49 < a < 0.51
    assert NS(
        Integral(sin(x)/x**2, (x, 1, oo)), quad='osc') == '0.504067061906928'
    assert NS(Integral(
        cos(pi*x + 1)/x, (x, -oo, -1)), quad='osc') == '0.276374705640365'
    # indefinite integrals aren't evaluated
    assert NS(Integral(x, x)) == 'Integral(x, x)'
    assert NS(Integral(x, (x, y))) == 'Integral(x, (x, y))'


@pytest.mark.xfail
def test_failing_integrals():
    assert isinstance(NS(Integral(sqrt(x) + x*y, (x, 1, 2), (y, -1, 1))),
                      Float)  # == '2.43790283299492'
    assert isinstance(NS(Integral(sin(x + x*y), (x, -1, 1), (y, -1, 1))),
                      Float)  # == '0.0'


def test_integrate_DiracDelta():
    # This is here to check that deltaintegrate is being called, but also
    # to test definite integrals. More tests are in test_deltafunctions.py
    assert integrate(DiracDelta(x) * f(x), (x, -oo, oo)) == f(0)
    assert integrate(DiracDelta(x) * f(x), (x, 0, oo)) == f(0)/2
    assert integrate(DiracDelta(x)**2, (x, -oo, oo)) == DiracDelta(0)
    # issue sympy/sympy#4522
    assert integrate(integrate((4 - 4*x + x*y - 4*y) *
                               DiracDelta(x)*DiracDelta(y - 1), (x, 0, 1)), (y, 0, 1)) == 0
    # issue sympy/sympy#5729
    p = exp(-(x**2 + y**2))/pi
    assert integrate(p*DiracDelta(x - 10*y), (x, -oo, oo), (y, -oo, oo)) == \
        integrate(p*DiracDelta(x - 10*y), (y, -oo, oo), (x, -oo, oo)) == \
        integrate(p*DiracDelta(10*x - y), (x, -oo, oo), (y, -oo, oo)) == \
        integrate(p*DiracDelta(10*x - y), (y, -oo, oo), (x, -oo, oo)) == \
        1/sqrt(101*pi)


@pytest.mark.xfail
@pytest.mark.slow
def test_integrate_DiracDelta_fails():
    # issue sympy/sympy#6427
    integrate(DiracDelta(x - y - z), (x, 0, 1),
              (y, 0, 1), (z, 0, oo))  # = Rational(1, 2)


def test_integrate_returns_piecewise():
    assert integrate(x**y, x) == Piecewise(
        (log(x), Eq(y, -1)), (x**(y + 1)/(y + 1), True))
    assert integrate(x**y, y) == Piecewise(
        (y, Eq(log(x), 0)), (x**y/log(x), True))
    assert integrate(exp(n*x), x) == Piecewise(
        (x, Eq(n, 0)), (exp(n*x)/n, True))
    assert integrate(x*exp(n*x), x) == Piecewise(
        (x**2/2, Eq(n**3, 0)), ((x*n**2 - n)*exp(n*x)/n**3, True))
    assert integrate(x**(n*y), x) == Piecewise(
        (log(x), Eq(n*y, -1)), (x**(n*y + 1)/(n*y + 1), True))
    assert integrate(x**(n*y), y) == Piecewise(
        (y, Eq(n*log(x), 0)), (x**(n*y)/(n*log(x)), True))
    assert integrate(cos(n*x), x) == Piecewise(
        (x, Eq(n, 0)), (sin(n*x)/n, True))
    assert integrate(cos(n*x)**2, x) == Piecewise(
        (x, Eq(n, 0)), ((n*x/2 + sin(n*x)*cos(n*x)/2)/n, True))
    assert integrate(x*cos(n*x), x) == Piecewise(
        (x**2/2, Eq(n, 0)), (x*sin(n*x)/n + cos(n*x)/n**2, True))
    assert integrate(sin(n*x), x) == Piecewise(
        (0, Eq(n, 0)), (-cos(n*x)/n, True))
    assert integrate(sin(n*x)**2, x) == Piecewise(
        (0, Eq(n, 0)), ((n*x/2 - sin(n*x)*cos(n*x)/2)/n, True))
    assert integrate(x*sin(n*x), x) == Piecewise(
        (0, Eq(n, 0)), (-x*cos(n*x)/n + sin(n*x)/n**2, True))
    assert integrate(exp(x*y), (x, 0, z)) == Piecewise(
        (z, Eq(y, 0)), (exp(y*z)/y - 1/y, True))


def test_subs1():
    e = Integral(exp(x - y), x)
    assert e.subs({y: 3}) == Integral(exp(x - 3), x)
    e = Integral(exp(x - y), (x, 0, 1))
    assert e.subs({y: 3}) == Integral(exp(x - 3), (x, 0, 1))
    f = Lambda(x, exp(-x**2))
    conv = Integral(f(x - y)*f(y), (y, -oo, oo))
    assert conv.subs({x: 0}) == Integral(exp(-2*y**2), (y, -oo, oo))


def test_subs2():
    e = Integral(exp(x - y), x, t)
    assert e.subs({y: 3}) == Integral(exp(x - 3), x, t)
    e = Integral(exp(x - y), (x, 0, 1), (t, 0, 1))
    assert e.subs({y: 3}) == Integral(exp(x - 3), (x, 0, 1), (t, 0, 1))
    f = Lambda(x, exp(-x**2))
    conv = Integral(f(x - y)*f(y), (y, -oo, oo), (t, 0, 1))
    assert conv.subs({x: 0}) == Integral(exp(-2*y**2), (y, -oo, oo), (t, 0, 1))


def test_subs3():
    e = Integral(exp(x - y), (x, 0, y), (t, y, 1))
    assert e.subs({y: 3}) == Integral(exp(x - 3), (x, 0, 3), (t, 3, 1))
    f = Lambda(x, exp(-x**2))
    conv = Integral(f(x - y)*f(y), (y, -oo, oo), (t, x, 1))
    assert conv.subs({x: 0}) == Integral(exp(-2*y**2), (y, -oo, oo), (t, 0, 1))


def test_subs4():
    e = Integral(exp(x), (x, 0, y), (t, y, 1))
    assert e.subs({y: 3}) == Integral(exp(x), (x, 0, 3), (t, 3, 1))
    f = Lambda(x, exp(-x**2))
    conv = Integral(f(y)*f(y), (y, -oo, oo), (t, x, 1))
    assert conv.subs({x: 0}) == Integral(exp(-2*y**2), (y, -oo, oo), (t, 0, 1))


def test_subs5():
    e = Integral(exp(-x**2), (x, -oo, oo))
    assert e.subs({x: 5}) == e
    e = Integral(exp(-x**2 + y), x)
    assert e.subs({y: 5}) == Integral(exp(-x**2 + 5), x)
    e = Integral(exp(-x**2 + y), (x, x))
    assert e.subs({x: 5}) == Integral(exp(y - x**2), (x, 5))
    assert e.subs({y: 5}) == Integral(exp(-x**2 + 5), x)
    e = Integral(exp(-x**2 + y), (y, -oo, oo), (x, -oo, oo))
    assert e.subs({x: 5}) == e
    assert e.subs({y: 5}) == e
    # Test evaluation of antiderivatives
    e = Integral(exp(-x**2), (x, x))
    assert e.subs({x: 5}) == Integral(exp(-x**2), (x, 5))
    e = Integral(exp(x), x)
    assert (e.subs({x: 1}) - e.subs({x: 0}) -
            Integral(exp(x), (x, 0, 1))).doit().is_zero


def test_subs6():
    e = Integral(x*y, (x, f(x), f(y)))
    assert e.subs({x: 1}) == Integral(x*y, (x, f(1), f(y)))
    assert e.subs({y: 1}) == Integral(x, (x, f(x), f(1)))
    e = Integral(x*y, (x, f(x), f(y)), (y, f(x), f(y)))
    assert e.subs({x: 1}) == Integral(x*y, (x, f(1), f(y)), (y, f(1), f(y)))
    assert e.subs({y: 1}) == Integral(x*y, (x, f(x), f(y)), (y, f(x), f(1)))
    e = Integral(x*y, (x, f(x), f(a)), (y, f(x), f(a)))
    assert e.subs({a: 1}) == Integral(x*y, (x, f(x), f(1)), (y, f(x), f(1)))


def test_subs7():
    e = Integral(x, (x, 1, y), (y, 1, 2))
    assert e.subs({x: 1, y: 2}) == e
    e = Integral(sin(x) + sin(y), (x, sin(x), sin(y)),
                 (y, 1, 2))
    assert e.subs({sin(y): 1}) == e
    assert e.subs({sin(x): 1}) == Integral(sin(x) + sin(y), (x, 1, sin(y)),
                                           (y, 1, 2))


def test_expand():
    e = Integral(f(x)+f(x**2), (x, 1, y))
    assert e.expand() == Integral(f(x), (x, 1, y)) + Integral(f(x**2), (x, 1, y))


def test_integration_variable():
    pytest.raises(ValueError, lambda: Integral(exp(-x**2), 3))
    pytest.raises(ValueError, lambda: Integral(exp(-x**2), (3, -oo, oo)))


def test_expand_integral():
    assert Integral(cos(x**2)*(sin(x**2) + 1), (x, 0, 1)).expand() == \
        Integral(cos(x**2)*sin(x**2), (x, 0, 1)) + \
        Integral(cos(x**2), (x, 0, 1))
    assert Integral(cos(x**2)*(sin(x**2) + 1), x).expand() == \
        Integral(cos(x**2)*sin(x**2), x) + \
        Integral(cos(x**2), x)


def test_as_sum_midpoint1():
    e = Integral(sqrt(x**3 + 1), (x, 2, 10))
    assert e.as_sum(1, method='midpoint') == 8*sqrt(217)
    assert e.as_sum(2, method='midpoint') == 4*sqrt(65) + 12*sqrt(57)
    assert e.as_sum(3, method='midpoint') == 8*sqrt(217)/3 + \
        8*sqrt(3081)/27 + 8*sqrt(52809)/27
    assert e.as_sum(4, method='midpoint') == 2*sqrt(730) + \
        4*sqrt(7) + 4*sqrt(86) + 6*sqrt(14)
    assert abs(e.as_sum(4, method='midpoint').evalf() - e.evalf()) < 0.5

    e = Integral(sqrt(x**3 + y**3), (x, 2, 10), (y, 0, 10))
    pytest.raises(NotImplementedError, lambda: e.as_sum(4))


def test_as_sum_midpoint2():
    e = Integral((x + y)**2, (x, 0, 1))
    assert e.as_sum(1, method='midpoint').expand() == Rational(1, 4) + y + y**2
    assert e.as_sum(2, method='midpoint').expand() == Rational(5, 16) + y + y**2
    assert e.as_sum(3, method='midpoint').expand() == Rational(35, 108) + y + y**2
    assert e.as_sum(4, method='midpoint').expand() == Rational(21, 64) + y + y**2


def test_as_sum_trapezoid():
    e = Integral(sin(x), (x, 3, 7))
    assert e.as_sum(2, 'trapezoid') == 2*sin(5) + sin(3) + sin(7)


def test_as_sum_left():
    e = Integral((x + y)**2, (x, 0, 1))
    assert e.as_sum(1, method='left').expand() == y**2
    assert e.as_sum(2, method='left').expand() == Rational(1, 8) + y/2 + y**2
    assert e.as_sum(3, method='left').expand() == Rational(5, 27) + 2*y/3 + y**2
    assert e.as_sum(4, method='left').expand() == Rational(7, 32) + 3*y/4 + y**2


def test_as_sum_right():
    e = Integral((x + y)**2, (x, 0, 1))
    assert e.as_sum(1, method='right').expand() == 1 + 2*y + y**2
    assert e.as_sum(2, method='right').expand() == Rational(5, 8) + 3*y/2 + y**2
    assert e.as_sum(3, method='right').expand() == Rational(14, 27) + 4*y/3 + y**2
    assert e.as_sum(4, method='right').expand() == Rational(15, 32) + 5*y/4 + y**2


def test_as_sum_raises():
    e = Integral((x + y)**2, (x, 0, 1))
    pytest.raises(ValueError, lambda: e.as_sum(-1))
    pytest.raises(ValueError, lambda: e.as_sum(0))
    pytest.raises(ValueError, lambda: Integral(x).as_sum(3))
    pytest.raises(NotImplementedError, lambda: e.as_sum(oo))
    pytest.raises(NotImplementedError, lambda: e.as_sum(3, method='xxxx2'))


def test_integrate_conds():
    assert integrate(x**a*exp(-x), (x, 0, oo),
                     conds='separate') == (gamma(a + 1), -re(a) < 1)

    # issue sympy/sympy#4199
    ypos = Symbol('y', positive=True)
    # TODO: Remove conds='none' below, let the assumption take care of it.
    assert (integrate(exp(-I*2*pi*ypos*x)*x, (x, -oo, oo), conds='none') ==
            Integral(exp(-I*2*pi*ypos*x)*x, (x, -oo, oo)))
    assert (integrate(exp(-I*2*pi*ypos*x)*x, (x, 0, oo), conds='none') ==
            Integral(exp(-2*I*pi*x*ypos)*x, (x, 0, oo)))


def test_nested_doit():
    e = Integral(Integral(x, x), x)
    f = Integral(x, x, x)
    assert e.doit() == f.doit()


def test_sympyissue_4665():
    # Allow only upper or lower limit evaluation
    e = Integral(x**2, (x, None, 1))
    f = Integral(x**2, (x, 1, None))
    assert e.doit() == Rational(1, 3)
    assert f.doit() == Rational(-1, 3)
    assert Integral(x*y, (x, None, y)).subs({y: t}) == Integral(x*t, (x, None, t))
    assert Integral(x*y, (x, y, None)).subs({y: t}) == Integral(x*t, (x, t, None))
    assert integrate(x**2, (x, None, 1)) == Rational(1, 3)
    assert integrate(x**2, (x, 1, None)) == Rational(-1, 3)
    assert integrate('x**2', ('x', '1', None)) == Rational(-1, 3)


def test_integral_reconstruct():
    e = Integral(x**2, (x, -1, 1))
    assert e == Integral(*e.args)


def test_doit_integrals():
    e = Integral(Integral(2*x), (x, 0, 1))
    assert e.doit() == Rational(1, 3)
    assert e.doit(deep=False) == Rational(1, 3)
    f = Function('f')
    # doesn't matter if the integral can't be performed
    assert Integral(f(x), (x, 1, 1)).doit() == 0
    # doesn't matter if the limits can't be evaluated
    assert Integral(0, (x, 1, Integral(f(x), x))).doit() == 0
    assert Integral(x, (a, 0)).doit() == 0
    limits = ((a, 1, exp(x)), (x, 0))
    assert Integral(a, *limits).doit() == Rational(1, 4)
    assert Integral(a, *list(reversed(limits))).doit() == 0


def test_as_dummy():
    assert Integral(x, x).as_dummy() == Integral(x, x)


def test_sympyissue_4884():
    assert integrate(sqrt(x)*(1 + x)).simplify() == \
        Piecewise((2*sqrt(x)**3*(3*x + 5)/15, abs(x + 1) > 1),
                  (2*I*x*sqrt(-x)*(3*x + 5)/15, True))
    assert integrate(x**x*(1 + log(x))) == x**x


def test_is_number():
    assert Integral(x).is_number is False
    assert Integral(1, x).is_number is False
    assert Integral(1, (x, 1)).is_number is True
    assert Integral(1, (x, 1, 2)).is_number is True
    assert Integral(1, (x, 1, y)).is_number is False
    assert Integral(1, (x, y)).is_number is False
    assert Integral(x, y).is_number is False
    assert Integral(x, (y, 1, x)).is_number is False
    assert Integral(x, (y, 1, 2)).is_number is False
    assert Integral(x, (x, 1, 2)).is_number is True
    # `foo.is_number` should always be eqivalent to `not foo.free_symbols`
    # in each of these cases, there are pseudo-free symbols
    i = Integral(x, (y, 1, 1))
    assert i.is_number is False and i.evalf() == 0
    i = Integral(x, (y, z, z))
    assert i.is_number is False and i.evalf() == 0
    i = Integral(1, (y, z, z + 2))
    assert i.is_number is False and i.evalf() == 2

    assert Integral(x*y, (x, 1, 2), (y, 1, 3)).is_number is True
    assert Integral(x*y, (x, 1, 2), (y, 1, z)).is_number is False
    assert Integral(x, (x, 1)).is_number is True
    assert Integral(x, (x, 1, Integral(y, (y, 1, 2)))).is_number is True
    assert Integral(Sum(z, (z, 1, 2)), (x, 1, 2)).is_number is True
    # it is possible to get a false negative if the integrand is
    # actually an unsimplified zero, but this is true of is_number in general.
    assert Integral(sin(x)**2 + cos(x)**2 - 1, x).is_number is False
    assert Integral(f(x), (x, 0, 1)).is_number is True


def test_symbols():
    assert Integral(0, x).free_symbols == {x}
    assert Integral(x).free_symbols == {x}
    assert Integral(x, (x, None, y)).free_symbols == {y}
    assert Integral(x, (x, y, None)).free_symbols == {y}
    assert Integral(x, (x, 1, y)).free_symbols == {y}
    assert Integral(x, (x, y, 1)).free_symbols == {y}
    assert Integral(x, (x, x, y)).free_symbols == {x, y}
    assert Integral(x, x, y).free_symbols == {x, y}
    assert Integral(x, (x, 1, 2)).free_symbols == set()
    assert Integral(x, (y, 1, 2)).free_symbols == {x}
    # pseudo-free in this case
    assert Integral(x, (y, z, z)).free_symbols == {x, z}
    assert Integral(x, (y, 1, 2), (y, None, None)).free_symbols == {x, y}
    assert Integral(x, (y, 1, 2), (x, 1, y)).free_symbols == {y}
    assert Integral(2, (y, 1, 2), (y, 1, x), (x, 1, 2)).free_symbols == set()
    assert Integral(2, (y, x, 2), (y, 1, x), (x, 1, 2)).free_symbols == set()
    assert Integral(2, (x, 1, 2), (y, x, 2), (y, 1, 2)).free_symbols == \
        {x}


def test_is_zero():
    assert Integral(0, (x, 1, x)).is_zero
    assert Integral(1, (x, 1, 1)).is_zero
    assert Integral(1, (x, m)).is_zero is None
    assert Integral(1, (x, 1, 2), (y, 2)).is_nonzero
    assert Integral(x, (m, 0)).is_zero
    assert Integral(x + m, (m, 0)).is_zero is None
    i = Integral(m, (m, 1, exp(x)), (x, 0))
    assert i.is_zero is None
    assert Integral(m, (x, 0), (m, 1, exp(x))).is_zero

    assert Integral(x, (x, oo, oo)).is_zero  # issue sympy/sympy#8171
    assert Integral(x, (x, -oo, -oo)).is_zero

    # this is zero but is beyond the scope of what is_zero
    # should be doing
    assert Integral(sin(x), (x, 0, 2*pi)).is_zero is None


def test_is_real():
    assert Integral(x**3, (x, 1, 3)).is_real
    assert Integral(1/(x - 1), (x, -1, 1)).is_real is not True


def test_series():
    i = Integral(cos(x), (x, x))
    e = i.lseries(x)
    s1 = i.nseries(x, n=8).removeO().doit()
    s2 = Add(*[next(e) for j in range(4)])
    assert s1 == s2


def test_sympyissue_4403():
    x = Symbol('x', extended_real=True)
    y = Symbol('y', positive=True)
    assert integrate(1/(x**2 + y**2)**Rational(3, 2), x) == \
        x/(y**2*sqrt(x**2 + y**2))


def test_sympyissue_4100():
    R = Symbol('R', positive=True)
    assert integrate(sqrt(R**2 - x**2), (x, 0, R)) == pi*R**2/4


def test_sympyissue_5167():
    assert Integral(Integral(f(x), x), x) == Integral(f(x), x, x)
    assert Integral(f(x)).args == (f(x), Tuple(x))
    assert Integral(Integral(f(x))).args == (f(x), Tuple(x), Tuple(x))
    assert Integral(Integral(f(x)), y).args == (f(x), Tuple(x), Tuple(y))
    assert Integral(Integral(f(x), z), y).args == (f(x), Tuple(z), Tuple(y))
    assert Integral(Integral(Integral(f(x), x), y), z).args == \
        (f(x), Tuple(x), Tuple(y), Tuple(z))
    assert integrate(Integral(f(x), x), x) == Integral(f(x), x, x)
    assert integrate(Integral(f(x), y), x) == y*Integral(f(x), x)
    assert integrate(Integral(f(x), x), y) in [Integral(y*f(x), x), y*Integral(f(x), x)]
    assert integrate(Integral(2, x), x) == x**2
    assert integrate(Integral(2, x), y) == 2*x*y
    # don't re-order given limits
    assert Integral(1, x, y).args != Integral(1, y, x).args
    # do as many as possibble
    assert Integral(f(x), y, x, y, x).doit() == y**2*Integral(f(x), x, x)/2


def test_sympyissue_4890():
    z = Symbol('z', positive=True)
    assert integrate(exp(-z*log(x)**2), x) == \
        sqrt(pi)*exp(1/(4*z))*erf(sqrt(z)*log(x) - 1/(2*sqrt(z)))/(2*sqrt(z))


def test_sympyissue_4376():
    n = Symbol('n', integer=True, positive=True)
    assert simplify(integrate(n*(x**(1/n) - 1), (x, 0, Rational(1, 2))) -
                    (n**2 - 2**(1/n)*n**2 - n*2**(1/n))/(2**(1 + 1/n) + n*2**(1 + 1/n))) == 0


@pytest.mark.slow
def test_sympyissue_4517():
    assert integrate((sqrt(x) - x**3)/cbrt(x), x) == \
        6*x**Rational(7, 6)/7 - 3*x**Rational(11, 3)/11


def test_sympyissue_4527():
    k, m = symbols('k m', integer=True)
    assert integrate(sin(k*x)*sin(m*x), (x,)) == Piecewise(
        (0, And(Eq(k, 0), Eq(m, 0))),
        (-x*sin(m*x)**2/2 - x*cos(m*x)**2/2 + sin(m*x)*cos(m*x)/(2*m), Eq(k, -m)),
        (x*sin(m*x)**2/2 + x*cos(m*x)**2/2 - sin(m*x)*cos(m*x)/(2*m), Eq(k, m)),
        (m*sin(k*x)*cos(m*x)/(k**2 - m**2) -
         k*sin(m*x)*cos(k*x)/(k**2 - m**2), True))


@pytest.mark.slow
def test_sympyissue_3940():
    a, b, c, d = symbols('a:d', positive=True, finite=True)
    assert integrate(exp(-x**2 + I*c*x), x) == \
        -sqrt(pi)*exp(-c**2/4)*erf(I*c/2 - x)/2
    assert integrate(exp(a*x**2 + b*x + c), x) == \
        sqrt(pi)*exp(c)*exp(-b**2/(4*a))*erfi(sqrt(a)*x + b/(2*sqrt(a)))/(2*sqrt(a))

    assert expand_mul(integrate(exp(-x**2)*exp(I*k*x), (x, -oo, oo))) == \
        sqrt(pi)*exp(-k**2/4)
    a, d = symbols('a d', positive=True)
    assert expand_mul(integrate(exp(-a*x**2 + 2*d*x), (x, -oo, oo))) == \
        sqrt(pi)*exp(d**2/a)/sqrt(a)


def test_sympyissue_5413():
    # Note that this is not the same as testing ratint() becuase integrate()
    # pulls out the coefficient.
    assert integrate(-a/(a**2 + x**2), x) == I*log(-I*a + x)/2 - I*log(I*a + x)/2


def test_sympyissue_5907():
    a = Symbol('a', real=True)
    assert (integrate(1/(x**2 + a**2)**2, x) ==
            x/(2*a**4 + 2*a**2*x**2) + atan(x/a)/(2*a**3))


def test_sympyissue_4892a():
    c = Symbol('c', nonzero=True)
    P1 = -A*exp(-z)
    P2 = -A/(c*t)*(sin(x)**2 + cos(y)**2)

    h1 = -sin(x)**2 - cos(y)**2
    h2 = -sin(x)**2 + sin(y)**2 - 1

    # there is still some non-deterministic behavior in integrate
    # or trigsimp which permits one of the following
    assert integrate(c*(P2 - P1), t) in [
        c*(-A*(-h1)*log(c*t)/c + A*t*exp(-z)),
        c*(-A*(-h2)*log(c*t)/c + A*t*exp(-z)),
        c*(A*h1*log(c*t)/c + A*t*exp(-z)),
        c*(A*h2*log(c*t)/c + A*t*exp(-z)),
        (A*c*t - A*(-h1)*log(t)*exp(z))*exp(-z),
        (A*c*t - A*(-h2)*log(t)*exp(z))*exp(-z),
    ]


def test_sympyissue_4892b():
    # Issues relating to issue sympy/sympy#4596 are making the actual result of this hard
    # to test.  The answer should be something like
    #
    # (-sin(y) + sqrt(-72 + 48*cos(y) - 8*cos(y)**2)/2)*log(x + sqrt(-72 +
    # 48*cos(y) - 8*cos(y)**2)/(2*(3 - cos(y)))) + (-sin(y) - sqrt(-72 +
    # 48*cos(y) - 8*cos(y)**2)/2)*log(x - sqrt(-72 + 48*cos(y) -
    # 8*cos(y)**2)/(2*(3 - cos(y)))) + x**2*sin(y)/2 + 2*x*cos(y)

    expr = (sin(y)*x**3 + 2*cos(y)*x**2 + 12)/(x**2 + 2)
    assert trigsimp(factor(integrate(expr, x).diff(x) - expr)) == 0


def test_integrate_series():
    f = sin(x).series(x, 0, 10)
    g = x**2/2 - x**4/24 + x**6/720 - x**8/40320 + \
        x**10/3628800 + Integral(O(x**10), x)

    assert integrate(f, x) == g

    assert integrate(O(x**5), x) == Integral(O(x**5), x)


def test_atom_bug():
    assert heurisch(meijerg([], [], [1], [], x), x) is None


def test_limit_bug():
    z = Symbol('z', nonzero=True)
    assert integrate(sin(x*y*z), (x, 0, pi), (y, 0, pi)) == \
        (log(z**2) + 2*EulerGamma + 2*log(pi))/(2*z) - \
        (-log(pi*z) + log(pi**2*z**2)/2 + Ci(pi**2*z))/z + log(pi)/z


def test_sympyissue_4703():
    g = Function('g')
    assert integrate(exp(x)*g(x), x).has(Integral)


def test_sympyissue_1888():
    f = Function('f')
    assert integrate(f(x).diff(x)**2, x).has(Integral)

# The following tests work using meijerint.


def test_sympyissue_3558():
    assert integrate(cos(x*y), (x, -pi/2, pi/2), (y, 0, pi)) == 2*Si(pi**2/2)


def test_sympyissue_4422():
    assert integrate(1/sqrt(16 + 4*x**2), x) == asinh(x/2) / 2


def test_sympyissue_4493():
    assert simplify(integrate(x*sqrt(1 + 2*x), x)) == \
        sqrt(2*x + 1)*(6*x**2 + x - 1)/15


def test_sympyissue_4737():
    assert integrate(sin(x)/x, (x, -oo, oo)) == pi
    assert integrate(sin(x)/x, (x, 0, oo)) == pi/2
    assert integrate(sin(x)/x, x) == Si(x)


def test_sympyissue_4992():
    # Note: psi in _check_antecedents becomes NaN.
    a = Symbol('a', positive=True)
    assert simplify(expand_func(integrate(exp(-x)*log(x)*x**a, (x, 0, oo)))) == \
        (a*polygamma(0, a) + 1)*gamma(a)


def test_sympyissue_4487():
    assert simplify(integrate(exp(-x)*x**y, x)) == lowergamma(y + 1, x)


@pytest.mark.xfail
def test_sympyissue_4215():
    assert integrate(1/x**2, (x, -1, 1)) == oo


def test_sympyissue_4400():
    n = Symbol('n', integer=True, positive=True)
    assert integrate((x**n)*log(x), x) == \
        n*x*x**n*log(x)/(n**2 + 2*n + 1) + x*x**n*log(x)/(n**2 + 2*n + 1) - \
        x*x**n/(n**2 + 2*n + 1)


def test_sympyissue_4153():
    assert integrate(1/(1 + x + y + z), (x, 0, 1), (y, 0, 1), (z, 0, 1)) in [
        -12*log(3) - 3*log(6)/2 + 3*log(8)/2 + 5*log(2) + 7*log(4),
        6*log(2) + 8*log(4) - 27*log(3)/2, 22*log(2) - 27*log(3)/2,
        -12*log(3) - 3*log(6)/2 + 47*log(2)/2]


def test_sympyissue_4326():
    # It doesn't matter if we can do the integral.  Just make sure the result
    # doesn't contain nan.  This is really a test against _eval_interval.
    assert not integrate(((h*(x - R + b))/b)*sqrt(R**2 - x**2), (x, R - b, R)).has(nan)


def test_powers():
    assert integrate(2**x + 3**x, x) == 2**x/log(2) + 3**x/log(3)


def test_risch_option():
    # risch=True only allowed on indefinite integrals
    pytest.raises(ValueError, lambda: integrate(1/log(x), (x, 0, oo), risch=True))
    assert integrate(exp(-x**2), x, risch=True) == NonElementaryIntegral(exp(-x**2), x)
    assert integrate(log(1/x)*y, x, y, risch=True) == y**2*(x*log(1/x)/2 + x/2)
    assert integrate(erf(x), x, risch=True) == Integral(erf(x), x)
    # TODO: How to test risch=False?

    # issue sympy/sympy#2708
    # This test needs to use an integration function that can
    # not be evaluated in closed form.  Update as needed.

    f = 1/(a + z + log(z))
    integral_f = NonElementaryIntegral(f, (z, 2, 3))
    assert Integral(f, (z, 2, 3)).doit() == integral_f
    assert integrate(f + exp(z), (z, 2, 3)) == integral_f - exp(2) + exp(3)

    assert integrate(2*f + exp(z), (z, 2, 3)) == 2*integral_f - exp(2) + exp(3)
    assert (integrate(exp(1.2*n*s*z*(-t + z)/t), (z, 0, x)) ==
            1.0*NonElementaryIntegral(exp(-1.2*n*s*z)*exp(1.2*n*s*z**2/t),
                                      (z, 0, x)))


def test_sympyissue_6828():
    f = 1/(1.08*x**2 - 4.3)
    g = integrate(f, x).diff(x)
    assert verify_numerically(f, g, tol=1e-12)


@pytest.mark.xfail
def test_integrate_Piecewise_rational_over_reals():
    f = Piecewise(
        (0,                                              t - 478.515625*pi < 0),
        (13.2075145209219*pi/(0.000871222*t + 0.995)**2, t - 478.515625*pi >= 0))

    assert integrate(f, (t, 0, oo)) != 0  # ~20664.5


def test_sympyissue_4803():
    x_max = Symbol('x_max')
    assert integrate(y/pi*exp(-(x_max - x)/cos(a)), x) == \
        y*exp((x - x_max)/cos(a))*cos(a)/pi


def test_sympyissue_4234():
    assert integrate(1/sqrt(1 + tan(x)**2)) == tan(x) / sqrt(1 + tan(x)**2)


def test_sympyissue_4492():
    assert simplify(integrate(x**2 * sqrt(5 - x**2), x)) == Piecewise(
        (I*(2*x**5 - 15*x**3 + 25*x - 25*sqrt(x**2 - 5)*acosh(sqrt(5)*x/5)) /
            (8*sqrt(x**2 - 5)), 1 < abs(x**2)/5),
        ((-2*x**5 + 15*x**3 - 25*x + 25*sqrt(-x**2 + 5)*asin(sqrt(5)*x/5)) /
            (8*sqrt(-x**2 + 5)), True))


def test_sympyissue_8368():
    assert integrate(exp(-s*x)*cosh(x), (x, 0, oo)) == \
        Piecewise((pi*Piecewise((-s/(pi*(-s**2 + 1)), abs(s**2) < 1),
                                (1/(pi*s*(1 - 1/s**2)), abs(s**-2) < 1), (meijerg(((Rational(1, 2),), (0, 0)),
                                                                                  ((0, Rational(1, 2)), (0,)), polar_lift(s)**2), True)),
                   And(abs(periodic_argument(polar_lift(s)**2, oo)) < pi, Ne(s**2, 1),
                       cos(abs(periodic_argument(polar_lift(s)**2, oo))/2)*sqrt(abs(s**2)) -
                       1 > 0)), (Integral(exp(-s*x)*cosh(x), (x, 0, oo)), True))
    assert integrate(exp(-s*x)*sinh(x), (x, 0, oo)) == \
        Piecewise((pi*Piecewise((2/(pi*(2*s**2 - 2)), abs(s**2) < 1),
                                (-2/(pi*s**2*(-2 + 2/s**2)), abs(s**-2) < 1),
                                (meijerg(((0,), (Rational(-1, 2), Rational(1, 2))),
                                         ((0, Rational(1, 2)), (Rational(-1, 2),)),
                                         polar_lift(s)**2), True)),
                   And(abs(periodic_argument(polar_lift(s)**2, oo)) < pi, Ne(s**2, 1),
                       cos(abs(periodic_argument(polar_lift(s)**2, oo))/2)*sqrt(abs(s**2)) - 1 > 0)),
                  (Integral(E**(-s*x)*sinh(x), (x, 0, oo)), True))


def test_sympyissue_8901():
    assert integrate(sinh(1.0*x)) == 1.0*cosh(1.0*x)
    assert integrate(tanh(1.0*x)) == 1.0*x - 1.0*log(tanh(1.0*x) + 1)


@pytest.mark.slow
def test_sympyissue_7130():
    integrand = (cos(pi*i*x/L)**2 / (a + b*x)).rewrite(exp)
    assert x not in integrate(integrand, (x, 0, L)).free_symbols


def test_sympyissue_4950():
    assert integrate((-60*exp(x) - 19.2*exp(4*x))*exp(4*x), x) ==\
        -2.4*exp(8*x) - 12.0*exp(5*x)


def test_sympyissue_4968():
    assert integrate(sin(log(x**2))) == x*sin(2*log(x))/5 - 2*x*cos(2*log(x))/5


def test_sympyissue_7098():
    assert integrate(1/sqrt(x) * 1/sqrt(1 - x), (x, 0, 1)) == pi


def test_sympyissue_4187():
    assert integrate(log(x)*exp(-x), (x, 0, oo)) == -EulerGamma


def test_sympyissue_10567():
    vt = Matrix([a*t, b, c])
    assert integrate(vt, t) == Integral(vt, t).doit()
    assert integrate(vt, t) == Matrix([[a*t**2/2], [b*t], [c*t]])


def test_sympyissue_4231():
    e = (1 + 2*x + sqrt(x + log(x))*(1 + 3*x) +
         x**2)/(x*(x + sqrt(x + log(x)))*sqrt(x + log(x)))
    assert integrate(e, x) == 2*log(x + sqrt(x + log(x))) + 2*sqrt(x + log(x))


def test_sympyissue_11045():
    e = 1/(x*sqrt(x**2 - 1))
    assert integrate(e, (x, 1, 2)) == pi/3


def test_definite_integrals_abs():
    # issue sympy/sympy#8430
    assert integrate(abs(x), (x, 0, 1)) == Rational(1, 2)
    # issue sympy/sympy#7165
    r = Symbol('r', real=True)
    assert (integrate(abs(x - r**2), (x, 0, 2)) ==
            r**2*Max(0, Min(2, r**2)) + r**2*Min(2, r**2) -
            2*r**2 - Max(0, Min(2, r**2))**2/2 -
            Min(2, r**2)**2/2 + 2)
    # issue sympy/sympy#8733
    assert integrate(abs(x + 1), (x, 0, 1)) == Rational(3, 2)

    e = x*abs(x**2 - 9)
    assert integrate(e, (x, -2, 2)) == 0
    assert integrate(e, (x, -1, 2)) == Rational(39, 4)
    assert integrate(e, (x, -2, 7)) == Rational(1625, 4)
    assert integrate(e, (x, -3, 11)) == 3136
    assert integrate(e, (x, -17, -2)) == Rational(-78425, 4)
    assert integrate(e, (x, -17, 20)) == Rational(74481, 4)


def test_definite_integrals_other():
    # issue diofant/diofant#303
    assert integrate((cos(x)/x)**2, (x, pi, 2*pi)) == Si(2*pi) - Si(4*pi) + 1/pi/2


def test_sympyissue_12081():
    assert integrate(x**(-Rational(3, 2))*exp(-x), (x, 0, oo)) == oo


def test_sympyissue_7163():
    integrate((sign(x - 1) - sign(x - 2))*cos(x), x)  # not raises


def test_sympyissue_12221():
    e = sqrt(1 - x)/x
    r = 2*I*(-sqrt(2) - asin(sqrt(3)/3) + asin(sqrt(5)/5) + 2)
    assert integrate(e, (x, 3, 5)).simplify() == r


def test_sympyissue_12582():
    assert integrate(abs(x**2 - 3*x), (x, -15, 15)) == 2259


def test_sympyissue_13312():
    assert integrate(exp(-a*t), (t, b, oo)) == Piecewise((-b + oo, Eq(a, 0)),
                                                         (-exp(-oo*sign(a))/a + exp(-a*b)/a, True))


def test_sympyissue_13501():
    a = Symbol('a', real=True)
    assert integrate(1/(1 + a**2*x**2), x) == atan(a*x)/a


def test_diofantissue_447():
    assert integrate(1/(2*sin(x) + cos(x)),
                     x) == (sqrt(5)*log(tan(x/2) - 2 + sqrt(5))/5 -
                            sqrt(5)*log(tan(x/2) - sqrt(5) - 2)/5)


def test_sympyissue_4511():
    assert simplify(integrate(cos(x)**2/(1 - sin(x)), x)) in [x - cos(x),
                                                              +1 - cos(x) + x,
                                                              -1 - cos(x) + x,
                                                              -2/(tan(x/2)**2 + 1) + x]


def test_sympyissue_4551():
    assert (integrate(1/(x*sqrt(1 - x**2)), x) ==
            Piecewise((-acosh(1/x), abs(x**-2) > 1), (I*asin(1/x), True)))


@pytest.mark.xfail
def test_sympyissue_4212():
    assert not isinstance(integrate(sign(x), x), Integral)


@pytest.mark.xfail
def test_sympyissue_4491():
    assert not isinstance(integrate(x*sqrt(x**2 + 2*x + 4), x), Integral)


@pytest.mark.xfail
def test_sympyissue_4514():
    assert not isinstance(integrate(sin(2*x)/sin(x)), Integral)


def test_sympyissue_15810():
    assert integrate(1/(2**(2*x/3) + 1), (x, 0, oo)) == Rational(3, 2)


def test_diofantissue_149():
    a = Symbol('a', positive=True)
    expr = (2 - x)**a*sin(a/x)
    res = sqrt(pi)*a*meijerg(((), (a/2 + 1/2, a/2 + 1)),
                             ((0, 0, Rational(1, 2)), (Rational(-1, 2),)),
                             a**2/16)*gamma(a + 1)/4
    assert integrate(expr, (x, 0, 2)) == res


def test_sympyissue_7337():
    assert integrate(x*sqrt(2*x + 3), (x, -1, 1)) == Rational(2, 5)


def test_sympyissue_11877():
    assert integrate(log(Rational(1, 2) - x),
                     (x, 0, Rational(1, 2))) == -Rational(1, 2) - log(2)/2


def test_sympyissue_17841():
    e = 1/(x**2 + x + I)
    assert integrate(e.diff(x), x) == e


def test_sympyissue_18384():
    e = abs(sin(x)*cos(x))
    assert integrate(e, (x, pi, 2*pi)) == 1
    assert integrate(e, (x, 0, pi/2)) == Rational(1, 2)
    assert integrate(e, (x, pi/2, pi)) == Rational(1, 2)
    assert integrate(e, (x, pi, 3*pi/2)) == Rational(1, 2)
    assert integrate(e, (x, 3*pi/2, 2*pi)) == Rational(1, 2)


def test_sympyissue_20360():
    e = exp(pi*x*(n - Rational(1, 2)))
    r = Piecewise((y, Eq(2*pi*n - pi, 0)),
                  (2*exp(pi*y*(n - Rational(1, 2)))/(2*pi*n - pi) -
                   2/(2*pi*n - pi), True))
    assert integrate(e, (x, 0, y)) == r


def test_sympyissue_20941():
    assert integrate(x**2*sqrt(1 - x**2), (x, 0, 1)) == pi/16


@pytest.mark.slow
def test_sympyissue_21034():
    f1 = x*(-x**4/asin(5)**4 - x*sinh(x + log(asin(5))) + 5)
    f2 = (x + cosh(cos(4)))/(x*(x + 1/(12*x)))

    assert (f1.integrate(x).diff(x) - f1).simplify() == 0
    assert (f2.integrate(x).diff(x) -
            f2).simplify().rewrite(exp).simplify() == 0


def test_sympyissue_21041():
    eq = sin(k*x)*exp(-x**2)

    assert integrate(eq.subs({k: 2}),
                     (x, 0, oo)) == -I*sqrt(pi)*erf(I)/(2*E)


def test_sympyissue_21063():
    assert integrate(exp(-x**2),
                     (x, Mul(-1, oo, evaluate=False), oo)) == sqrt(pi)


def test_sympyissue_21091():
    assert integrate(exp(-x**2)*sin(x), (x, -oo, oo)) == 0


def test_sympyissue_21132():
    f = exp(I*sqrt(k)*t)*cos(a*t)
    r = [Piecewise((t, Eq(a, 0) & Eq(k, 0)),
                   (exp(I*sqrt(k)*t)*(-I*cos(sqrt(k)*t) +
                                      exp(-I*sqrt(k)*t)*sqrt(k)*t)/(2*sqrt(k)),
                    Eq(a, sqrt(k)) | Eq(a, -sqrt(k))),
                   (exp(I*sqrt(k)*t)*(a*sin(a*t) +
                                      I*sqrt(k)*cos(a*t))/(a**2 - k), True)),
         Piecewise((t, Eq(a, 0) & Eq(k, 0)),
                   (exp(I*sqrt(k)*t)*(sin(sqrt(k)*t) +
                                      exp(-I*sqrt(k)*t)*sqrt(k)*t)/(2*sqrt(k)),
                    Eq(a, sqrt(k)) | Eq(a, -sqrt(k))),
                   (exp(I*sqrt(k)*t)*(a*sin(a*t) +
                                      I*sqrt(k)*cos(a*t))/(a**2 - k), True))]
    ans = f.integrate(t)
    assert ans.simplify() in r
    assert ans.subs({k: 0}).subs({a: 0}) == t


def test_sympyissue_21342():
    assert (1/(exp(I*x) - 2)).integrate((x, 0, 2*pi)) == -pi


def test_sympyissue_21024():
    assert ((x + exp(3))/x**2).integrate(x) == log(x) - exp(3)/x
    assert ((x**2 + exp(5))/x).integrate(x) == x**2/2 + exp(5)*log(x)
    assert ((x/(2*x + tanh(1))).integrate(x) ==
            x/2 - (-1 + E)*(1 + E)*log(2*x + tanh(1))/Mul(4, 1 + E**2,
                                                          evaluate=False))
    assert ((log(x)*log(4*x) + log(3*x + exp(2))).integrate(x) ==
            x*log(x)**2 + x*log(3*x + E**2) - x + x*(-2*log(2) + 2) +
            (-2*x + 2*x*log(2))*log(x) + E**2*log(3*x + E**2)/3)
