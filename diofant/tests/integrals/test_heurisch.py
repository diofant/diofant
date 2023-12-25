import pytest

from diofant import (Derivative, E, Ei, Eq, Function, I, Integral, LambertW,
                     Piecewise, Rational, Sum, Symbol, acos, asin, asinh,
                     besselj, cos, cosh, diff, erf, erfi, exp, li, log, pi,
                     ratsimp, root, simplify, sin, sinh, sqrt, tan)
from diofant.abc import nu, x, y, z
from diofant.integrals.heurisch import components, heurisch, heurisch_wrapper


__all__ = ()

f = Function('f')


def test_components():
    assert components(x*y, x) == {x}
    assert components(1/(x + y), x) == {x}
    assert components(sin(x), x) == {sin(x), x}
    assert components(sin(x)*sqrt(log(x)), x) == \
        {log(x), sin(x), sqrt(log(x)), x}
    assert components(x*sin(exp(x)*y), x) == \
        {sin(y*exp(x)), x, exp(x)}
    assert components(x**Rational(17, 54)/sqrt(sin(x)), x) == \
        {sin(x), root(x, 54), sqrt(sin(x)), x}

    assert components(f(x), x) == \
        {x, f(x)}
    assert components(Derivative(f(x), x), x) == \
        {x, f(x), Derivative(f(x), x)}
    assert components(f(x)*diff(f(x), x), x) == \
        {x, f(x), Derivative(f(x), x), Derivative(f(x), x)}


def test_heurisch_polynomials():
    assert heurisch(1, x) == x
    assert heurisch(x, x) == x**2/2
    assert heurisch(x**17, x) == x**18/18


def test_heurisch_fractions():
    assert heurisch(1/x, x) == log(x)
    assert heurisch(1/(2 + x), x) == log(x + 2)
    assert heurisch(1/(x + sin(y)), x) == log(x + sin(y))

    # Up to a constant, where C = 5*pi*I/12, Mathematica gives identical
    # result in the first case. The difference is because diofant changes
    # signs of expressions without any care.
    # XXX ^ ^ ^ is this still correct?
    assert heurisch(5*x**5/(
        2*x**6 - 5), x) in [5*log(2*x**6 - 5) / 12, 5*log(-2*x**6 + 5) / 12]
    assert heurisch(5*x**5/(2*x**6 + 5), x) == 5*log(2*x**6 + 5) / 12

    assert heurisch(1/x**2, x) == -1/x
    assert heurisch(-1/x**5, x) == 1/(4*x**4)


def test_heurisch_log():
    assert heurisch(log(x), x) == x*log(x) - x
    assert heurisch(log(3*x), x) == -x + x*log(3) + x*log(x)
    assert heurisch(log(x**2), x) in [x*log(x**2) - 2*x, 2*x*log(x) - 2*x]


def test_heurisch_exp():
    assert heurisch(exp(x), x) == exp(x)
    assert heurisch(exp(-x), x) == -exp(-x)
    assert heurisch(exp(17*x), x) == exp(17*x) / 17
    assert heurisch(x*exp(x), x) == x*exp(x) - exp(x)
    assert heurisch(x*exp(x**2), x) == exp(x**2) / 2

    assert heurisch(exp(-x**2), x) is None

    assert heurisch(2**x, x) == 2**x/log(2)
    assert heurisch(x*2**x, x) == x*2**x/log(2) - 2**x*log(2)**(-2)

    assert heurisch(Integral(x**z*y, (y, 1, 2), (z, 2, 3)).function, x) == (x*x**z*y)/(z+1)
    assert heurisch(Sum(x**z, (z, 1, 2)).function, z) == x**z/log(x)


def test_heurisch_trigonometric():
    assert heurisch(sin(x), x) == -cos(x)
    assert heurisch(pi*sin(x) + 1, x) == x - pi*cos(x)

    assert heurisch(cos(x), x) == sin(x)
    assert heurisch(tan(x), x) in [
        log(1 + tan(x)**2)/2,
        log(tan(x) + I) + I*x,
        log(tan(x) - I) - I*x,
    ]

    assert heurisch(sin(x)*sin(y), x) == -cos(x)*sin(y)
    assert heurisch(sin(x)*sin(y), y) == -cos(y)*sin(x)

    assert heurisch(sin(x)*cos(x), x) in [sin(x)**2/2, -cos(x)**2/2]
    assert heurisch(cos(x)/sin(x), x) == log(sin(x))

    assert heurisch(x*sin(7*x), x) == sin(7*x) / 49 - x*cos(7*x) / 7
    assert heurisch(1/pi/4 * x**2*cos(x), x) == 1/pi/4*(x**2*sin(x) -
                                                        2*sin(x) + 2*x*cos(x))

    assert heurisch(acos(x/4) * asin(x/4), x) == 2*x - (sqrt(16 - x**2))*asin(x/4) \
        + (sqrt(16 - x**2))*acos(x/4) + x*asin(x/4)*acos(x/4)

    assert heurisch(1/sin(1/x)/x**2, x) == -log(tan(1/x/2))


def test_heurisch_hyperbolic():
    assert heurisch(sinh(x), x) == cosh(x)
    assert heurisch(cosh(x), x) == sinh(x)

    assert heurisch(x*sinh(x), x) == x*cosh(x) - sinh(x)
    assert heurisch(x*cosh(x), x) == x*sinh(x) - cosh(x)

    assert heurisch(
        x*asinh(x/2), x) == x**2*asinh(x/2)/2 + asinh(x/2) - x*sqrt(4 + x**2)/4


def test_heurisch_mixed():
    assert heurisch(sin(x)*exp(x), x) == exp(x)*sin(x)/2 - exp(x)*cos(x)/2


def test_heurisch_radicals():
    assert heurisch(1/sqrt(x), x) == 2*sqrt(x)
    assert heurisch(1/sqrt(x)**3, x) == -2/sqrt(x)
    assert heurisch(sqrt(x)**3, x) == 2*sqrt(x)**5/5

    assert heurisch(sin(x)*sqrt(cos(x)), x) == -2*sqrt(cos(x))**3/3
    y = Symbol('y')
    assert heurisch(sin(y*sqrt(x)), x) == 2/y**2*sin(y*sqrt(x)) - \
        2*sqrt(x)*cos(y*sqrt(x))/y
    assert heurisch_wrapper(sin(y*sqrt(x)), x) == Piecewise(
        (0, Eq(y, 0)),
        (-2*sqrt(x)*cos(sqrt(x)*y)/y + 2*sin(sqrt(x)*y)/y**2, True))
    y = Symbol('y', positive=True)
    assert heurisch_wrapper(sin(y*sqrt(x)), x) == 2/y**2*sin(y*sqrt(x)) - \
        2*sqrt(x)*cos(y*sqrt(x))/y


def test_heurisch_special():
    assert heurisch(erf(x), x) == x*erf(x) + exp(-x**2)/sqrt(pi)
    assert heurisch(exp(-x**2)*erf(x), x) == sqrt(pi)*erf(x)**2 / 4


def test_heurisch_symbolic_coeffs():
    assert heurisch(1/(x + y), x) == log(x + y)
    assert heurisch(1/(x + sqrt(2)), x) == log(x + sqrt(2))
    assert simplify(diff(heurisch(log(x + y + z), y), y)) == log(x + y + z)


def test_heurisch_symbolic_coeffs_1130():
    y = Symbol('y')
    assert (heurisch_wrapper(1/(x**2 + y), x) ==
            Piecewise((-1/x, Eq(y, 0)),
                      (log(x - sqrt(-y))/(2*sqrt(-y)) -
                       log(x + sqrt(-y))/(2*sqrt(-y)), True)))
    y = Symbol('y', positive=True)
    assert heurisch_wrapper(1/(x**2 + y), x) in [I/sqrt(y)*log(x + sqrt(-y))/2 -
                                                 I/sqrt(y)*log(x - sqrt(-y))/2, I*log(x + I*sqrt(y)) /
                                                 (2*sqrt(y)) - I*log(x - I*sqrt(y))/(2*sqrt(y))]


def test_heurisch_hacking():
    assert (heurisch(sqrt(1 + 7*x**2), x, hints=[]) ==
            x*sqrt(1 + 7*x**2)/2 + sqrt(7)*asinh(sqrt(7)*x)/14)
    assert (heurisch(sqrt(1 - 7*x**2), x, hints=[]) ==
            x*sqrt(1 - 7*x**2)/2 + sqrt(7)*asin(sqrt(7)*x)/14)

    assert heurisch(sqrt(y*x**2 - 1), x, hints=[]) is None

    assert (heurisch(1/sqrt(1 + 7*x**2), x, hints=[]) ==
            sqrt(7)*asinh(sqrt(7)*x)/7)
    assert (heurisch(1/sqrt(1 - 7*x**2), x, hints=[]) ==
            sqrt(7)*asin(sqrt(7)*x)/7)

    assert heurisch(exp(-7*x**2), x, hints=[]) == sqrt(7*pi)*erf(sqrt(7)*x)/14
    assert heurisch(exp(2*x**2), x,
                    hints=[]) == sqrt(2)*sqrt(pi)*erfi(sqrt(2)*x)/4

    assert (heurisch(exp(2*x**2 - 3*x), x, hints=[]) ==
            sqrt(2)*sqrt(pi)*erfi(sqrt(2)*x - 3*sqrt(2)/4)/(4*E**Rational(9, 8)))

    assert heurisch(1/sqrt(9 - 4*x**2), x, hints=[]) == asin(2*x/3)/2
    assert heurisch(1/sqrt(9 + 4*x**2), x, hints=[]) == asinh(2*x/3)/2

    assert heurisch(li(x), x, hints=[]) == x*li(x) - Ei(2*log(x))
    assert heurisch(li(log(x)), x, hints=[]) is None

    assert (heurisch(sqrt(1 + x), x, hints=[x, sqrt(1 + x)]) ==
            2*x*sqrt(x + 1)/3 + 2*sqrt(x + 1)/3)


def test_heurisch_function():
    assert heurisch(f(x), x) is None


def test_heurisch_wrapper():
    assert heurisch_wrapper(1, x) == x
    f = 1/(y + x)
    assert heurisch_wrapper(f, x) == log(x + y)
    f = 1/(y - x)
    assert heurisch_wrapper(f, x) == -log(x - y)
    f = 1/((y - x)*(y + x))
    assert heurisch_wrapper(f, x) == \
        Piecewise((1/x, Eq(y, 0)), (log(x + y)/2/y - log(x - y)/2/y, True))
    # issue sympy/sympy#6926
    f = sqrt(x**2/((y - x)*(y + x)))
    assert (heurisch_wrapper(f, x) ==
            x*sqrt(x**2/(-x**2 + y**2)) - y**2*sqrt(x**2/(-x**2 + y**2))/x)


def test_sympyissue_3609():
    assert heurisch(1/(x * (1 + log(x)**2)), x) == I*log(log(x) + I)/2 - \
        I*log(log(x) - I)/2

# These are examples from the Poor Man's Integrator
# http://www-sop.inria.fr/cafe/Manuel.Bronstein/pmint/examples/


def test_pmint_rat():
    f = (x**7 - 24*x**4 - 4*x**2 + 8*x - 8)/(x**8 + 6*x**6 + 12*x**4 + 8*x**2)
    g = (4 + 8*x**2 + 6*x + 3*x**3)/(x**5 + 4*x**3 + 4*x) + log(x)

    assert ratsimp(heurisch(f, x)) == g


def test_pmint_trig():
    f = (x - tan(x)) / tan(x)**2 + tan(x)
    g = -x**2/2 - x/tan(x) + log(tan(x)**2 + 1)/2

    assert heurisch(f, x) == g


@pytest.mark.slow  # 8 seconds on 3.4 GHz
def test_pmint_logexp():
    f = (1 + x + x*exp(x))*(x + log(x) + exp(x) - 1)/(x + log(x) + exp(x))**2/x
    g = log(x**2 + 2*x*exp(x) + 2*x*log(x) + exp(2*x) + 2*exp(x)*log(x) + log(x)**2)/2 + 1/(x + exp(x) + log(x))

    # TODO: Optimal solution is g = 1/(x + log(x) + exp(x)) + log(x + log(x) + exp(x)),
    # but Diofant requires a lot of guidance to properly simplify heurisch() output.

    assert ratsimp(heurisch(f, x)) == g


@pytest.mark.slow  # 8 seconds on 3.4 GHz
def test_pmint_erf():
    f = exp(-x**2)*erf(x)/(erf(x)**3 - erf(x)**2 - erf(x) + 1)
    g = sqrt(pi)*log(erf(x) - 1)/8 - sqrt(pi)*log(erf(x) + 1)/8 - sqrt(pi)/(4*erf(x) - 4)

    assert ratsimp(heurisch(f, x)) == g


def test_pmint_LambertW():
    f = LambertW(x)
    g = x*LambertW(x) - x + x/LambertW(x)

    assert heurisch(f, x) == g


@pytest.mark.xfail
def test_pmint_besselj():
    # TODO: in both cases heurisch() gives None. Wrong besselj() derivative?

    f = besselj(nu + 1, x)/besselj(nu, x)
    g = nu*log(x) - log(besselj(nu, x))

    assert simplify(heurisch(f, x) - g) == 0

    f = (nu*besselj(nu, x) - x*besselj(nu + 1, x))/x
    g = besselj(nu, x)

    assert simplify(heurisch(f, x) - g) == 0


@pytest.mark.slow
def test_pmint_WrightOmega():
    def omega(x):
        return LambertW(exp(x))

    f = (1 + omega(x) * (2 + cos(omega(x)) * (x + omega(x))))/(1 + omega(x))/(x + omega(x))
    g = log(x + omega(x)) + sin(omega(x))

    assert heurisch(f, x) == g


def test_RR():
    # Make sure the algorithm does the right thing if the ring is RR. See
    # issue sympy/sympy#8685.
    assert heurisch(sqrt(1 + 0.25*x**2), x, hints=[]) == \
        0.5*x*sqrt(0.25*x**2 + 1) + 1.0*asinh(0.5*x)

# TODO: convert the rest of PMINT tests:
# Airy functions
# f = (x - AiryAi(x)*AiryAi(1, x)) / (x**2 - AiryAi(x)**2)
# g = Rational(1,2)*ln(x + AiryAi(x)) + Rational(1,2)*ln(x - AiryAi(x))
# f = x**2 * AiryAi(x)
# g = -AiryAi(x) + AiryAi(1, x)*x
# Whittaker functions
# f = WhittakerW(mu + 1, nu, x) / (WhittakerW(mu, nu, x) * x)
# g = x/2 - mu*ln(x) - ln(WhittakerW(mu, nu, x))
