"""
This test suite is testing the Gruntz algorithm implementation.

See the documentation in the gruntz module.  The
algorithm itself is highly recursive by nature, so ``compare()`` is
logically the lowest part of the algorithm, yet in some sense it's the most
complex part, because it needs to calculate a limit to return the result.
"""

import pytest

from diofant import (Add, E, Ei, EulerGamma, GoldenRatio, I, Integer, Li,
                     Limit, Max, Min, Mul, Pow, Rational, Symbol, acosh, acot,
                     airyai, airybi, atan, binomial, cbrt, cos, cosh, coth,
                     digamma, erf, exp, factorial, fibonacci, gamma, li, limit,
                     log, loggamma, oo, pi, root, sign, sin, sinh, sqrt, tan,
                     tanh, zeta)
from diofant.abc import a, n, y
from diofant.calculus.gruntz import compare, leadterm, mrv, rewrite, signinf


__all__ = ()

x = Symbol('x', real=True, positive=True)
m = Symbol('m', real=True, positive=True)


@pytest.mark.slow
def test_gruntz_evaluation():
    # Gruntz' thesis, problems 8.1-20, pp. 122 to 123
    assert limit(exp(x)*(exp(1/x - exp(-x)) - exp(1/x)), x, oo) == -1
    assert limit(exp(x)*(exp(1/x + exp(-x) + exp(-x**2))
                         - exp(1/x - exp(-exp(x)))), x, oo) == 1
    assert limit(exp(exp(x - exp(-x))/(1 - 1/x)) - exp(exp(x)), x, oo) == oo
    assert limit(exp(exp(exp(x)/(1 - 1/x)))
                 - exp(exp(exp(x)/(1 - 1/x - log(x)**(-log(x))))), x, oo) == -oo
    assert limit(exp(exp(exp(x + exp(-x)))) / exp(exp(exp(x))), x, oo) == oo
    assert limit(exp(exp(exp(x))) / exp(exp(exp(x - exp(-exp(x))))),
                 x, oo) == oo
    assert limit(exp(exp(exp(x))) / exp(exp(exp(x - exp(-exp(exp(x)))))),
                 x, oo) == 1
    assert limit(exp(exp(x)) / exp(exp(x - exp(-exp(exp(x))))), x, oo) == 1
    assert limit(log(x)**2*exp(sqrt(log(x))*(log(log(x)))**2 *
                 exp(sqrt(log(log(x)))*(log(log(log(x))))**3))/sqrt(x),
                 x, oo) == 0
    assert limit((x*log(x)*(log(x*exp(x) - x**2))**2) /
                 (log(log(x**2 + 2*exp(exp(3*x**3*log(x)))))),
                 x, oo) == Rational(1, 3)
    assert limit((exp(x*exp(-x)/(exp(-x) + exp(-2*x**2/(x + 1)))) -
                 exp(x))/x, x, oo) == -exp(2)
    assert limit((3**x + 5**x)**(1/x), x, oo) == 5
    assert limit(x/log(x**(log(x**(log(2)/log(x))))), x, oo) == oo
    assert limit(exp(exp(2*log(x**5 + x)*log(log(x)))) /
                 exp(exp(10*log(x)*log(log(x)))), x, oo) == oo
    assert limit(exp(exp(Rational(5, 2)*x**(-Rational(5, 7)) +
                         Rational(21, 8)*x**Rational(6, 11) + 2*x**(-8) +
                         Rational(54, 17)*x**Rational(49, 45)))**8 /
                 log(log(-log(Rational(4, 3) *
                     x**(-Rational(5, 14)))))**Rational(7, 6), x, oo) == oo
    assert limit((exp(4*x*exp(-x)/(1/exp(x) + 1/exp(2*x**2/(x + 1)))) -
                  exp(x))/exp(x)**4, x, oo) == 1
    assert limit(exp(x*exp(-x)/(exp(-x) +
                                exp(-2*x**2/(x + 1))))/exp(x), x, oo) == 1
    assert limit((exp(exp(-x/(1 + exp(-x)))) *
                  exp(-x/(1 + exp(-x/(1 + exp(-x))))) *
                  exp(exp(-x + exp(-x/(1 + exp(-x)))))) /
                 (exp(-x/(1 + exp(-x))))**2 - exp(x) + x, x, oo) == 2
    assert limit(log(x)*(log(log(x) + log(log(x))) - log(log(x))) /
                 (log(log(x) + log(log(log(x))))), x, oo) == 1
    assert limit(exp((log(log(x + exp(log(x)*log(log(x)))))) /
                     log(log(log(exp(x) + x + log(x))))), x, oo) == E
    # Another
    assert limit(exp(exp(exp(x + exp(-x))))/exp(exp(x)), x, oo) == oo
    assert limit(log(log(x*exp(x*exp(x)) + 1)) -
                 exp(exp(log(log(x)) + 1/x)), x, oo) == 0


def test_gruntz_eval_special():
    # Gruntz, p. 126
    assert limit(exp(x)*(sin(1/x + exp(-x)) - sin(1/x + exp(-x**2))),
                 x, oo) == 1
    assert limit((erf(x - exp(-exp(x))) - erf(x))*exp(exp(x))*exp(x**2),
                 x, oo) == -2/sqrt(pi)
    assert limit(exp(exp(x))*(exp(sin(1/x + exp(-exp(x)))) -
                              exp(sin(1/x))), x, oo) == 1
    assert limit(exp(x)*(gamma(x + exp(-x)) - gamma(x)), x, oo) == oo
    assert limit(exp(exp(digamma(digamma(x))))/x, x, oo) == exp(-Rational(1, 2))
    assert limit(exp(exp(digamma(log(x))))/x, x, oo) == exp(-Rational(1, 2))
    assert limit(digamma(digamma(digamma(x))), x, oo) == oo
    assert limit(loggamma(loggamma(x)), x, oo) == oo
    assert limit(((gamma(x + 1/gamma(x)) - gamma(x))/log(x) - cos(1/x)) *
                 x*log(x), x, oo) == -Rational(1, 2)
    assert limit(x*(gamma(x - 1/gamma(x)) - gamma(x) + log(x)),
                 x, oo) == Rational(1, 2)
    assert limit((gamma(x + 1/gamma(x)) - gamma(x))/log(x), x, oo) == 1
    assert limit(erf(log(1/x)), x, oo) == -1


@pytest.mark.slow
def test_gruntz_eval_special_slow():
    assert limit(gamma(x + 1)/sqrt(2*pi) - exp(-x)*(x**(x + Rational(1, 2)) +
                 x**(x - Rational(1, 2))/12), x, oo) == oo
    assert limit(exp(exp(exp(digamma(digamma(digamma(x))))))/x, x, oo) == 0
    assert limit(exp(gamma(x - exp(-x))*exp(1/x)) - exp(gamma(x)), x, oo) == oo
    assert limit((Ei(x - exp(-exp(x))) - Ei(x)) *
                 exp(-x)*exp(exp(x))*x, x, oo) == -1
    assert limit(exp((log(2) + 1)*x)*(zeta(x + exp(-x)) - zeta(x)),
                 x, oo) == -log(2)
    # TODO 8.36 (bessel)
    assert limit(Max(x, exp(x))/log(Min(exp(-x), exp(-exp(x)))), x, oo) == -1


def test_gruntz_other():
    assert limit(sqrt(log(x + 1)) - sqrt(log(x)), x, oo) == 0  # p12, 2.5
    assert limit(((1 + 1/x)**y - 1)*x, x, oo) == y  # p12, 2.6
    assert limit(x**n/exp(x), x, oo) == 0  # p14, 2.9
    assert limit((1 + 1/x)*x - 1/log(1 + 1/x),
                 x, oo) == Rational(1, 2)  # p15, 2.10
    assert limit((root(1 + 1/x, n) - 1)/(root(1 + 1/x, m) - 1),
                 x, oo) == m/n  # p13, 2.7


def test_gruntz_hyperbolic():
    assert limit(cosh(x), x, oo) == oo
    assert limit(cosh(-x), x, oo) == oo
    assert limit(sinh(x), x, oo) == oo
    assert limit(sinh(-x), x, oo) == -oo
    assert limit(2*cosh(x)*exp(x), x, oo) == oo
    assert limit(2*cosh(-x)*exp(-x), x, oo) == 1
    assert limit(2*sinh(x)*exp(x), x, oo) == oo
    assert limit(2*sinh(-x)*exp(-x), x, oo) == -1
    assert limit(tanh(x), x, oo) == 1
    assert limit(tanh(-x), x, oo) == -1
    assert limit(coth(x), x, oo) == 1
    assert limit(coth(-x), x, oo) == -1


def test_compare():
    assert compare(Integer(2), x, x) < 0
    assert compare(x, exp(x), x) < 0
    assert compare(exp(x), exp(x**2), x) < 0
    assert compare(exp(x**2), exp(exp(x)), x) < 0
    assert compare(Integer(1), exp(exp(x)), x) < 0

    assert compare(x, Integer(2), x) > 0
    assert compare(exp(x), x, x) > 0
    assert compare(exp(x**2), exp(x), x) > 0
    assert compare(exp(exp(x)), exp(x**2), x) > 0
    assert compare(exp(exp(x)), Integer(1), x) > 0

    assert compare(Integer(2), Integer(3), x) == 0
    assert compare(Integer(3), Integer(-5), x) == 0
    assert compare(Integer(2), Integer(-5), x) == 0

    assert compare(x, x**2, x) == 0
    assert compare(x**2, x**3, x) == 0
    assert compare(x**3, 1/x, x) == 0
    assert compare(1/x, x**m, x) == 0
    assert compare(x**m, -x, x) == 0

    assert compare(exp(x), exp(-x), x) == 0
    assert compare(exp(-x), exp(2*x), x) == 0
    assert compare(exp(2*x), exp(x)**2, x) == 0
    assert compare(exp(x)**2, exp(x + exp(-x)), x) == 0
    assert compare(exp(x), exp(x + exp(-x)), x) == 0

    assert compare(exp(x**2), 1/exp(x**2), x) == 0

    assert compare(exp(x), x**5, x) > 0
    assert compare(exp(x**2), exp(x)**2, x) > 0
    assert compare(exp(x), exp(x + exp(-x)), x) == 0
    assert compare(exp(x + exp(-x)), exp(x), x) == 0
    assert compare(exp(x + exp(-x)), exp(-x), x) == 0
    assert compare(exp(-x), x, x) > 0
    assert compare(x, exp(-x), x) < 0
    assert compare(exp(x + 1/x), x, x) > 0
    assert compare(exp(-exp(x)), exp(x), x) > 0
    assert compare(exp(exp(-exp(x)) + x), exp(-exp(x)), x) < 0

    assert compare(exp(exp(x)), exp(x + exp(-exp(x))), x) > 0


def test_signinf():
    assert signinf(Integer(0), x) == 0
    assert signinf(Integer(3), x) == 1
    assert signinf(Integer(-5), x) == -1
    assert signinf(log(x), x) == 1
    assert signinf(exp(-x), x) == 1
    assert signinf(exp(x), x) == 1
    assert signinf(-exp(x), x) == -1
    assert signinf(3 - 1/x, x) == 1
    assert signinf(-3 - 1/x, x) == -1
    assert signinf(sin(1/x), x) == 1
    assert signinf(x**2, x) == 1
    assert signinf(x**5, x) == 1

    assert signinf(x, x) == 1
    assert signinf(-x, x) == -1
    y = Symbol('y', positive=True)
    assert signinf(y, x) == 1
    assert signinf(-y, x) == -1
    assert signinf(y*x, x) == 1
    assert signinf(-y*x, x) == -1


def test_mrv():
    assert mrv(x, x) == {x}
    assert mrv(x + 1/x, x) == {x}
    assert mrv(x**2, x) == {x}
    assert mrv(log(x), x) == {x}
    assert mrv(exp(x), x) == {exp(x)}
    assert mrv(exp(-x), x) == {exp(-x)}
    assert mrv(exp(x**2), x) == {exp(x**2)}
    assert mrv(-exp(1/x), x) == {x}
    assert mrv(exp(x + 1/x), x) == {exp(x + 1/x)}

    assert mrv(exp(x + exp(-exp(x))), x) == {exp(-exp(x))}
    assert mrv(exp(x + exp(-x)), x) == {exp(x + exp(-x)), exp(-x)}
    assert mrv(exp(1/x + exp(-x)), x) == {exp(-x)}

    assert mrv(exp(x + exp(-x**2)), x) == {exp(-x**2)}

    assert mrv(exp(-x + 1/x**2) - exp(x + 1/x), x) == {exp(x + 1/x),
                                                       exp(1/x**2 - x)}

    assert mrv(exp(x**2) + x*exp(x) + exp(x*log(log(x)))/x, x) == {exp(x**2)}
    assert mrv(exp(x)*(exp(1/x + exp(-x)) - exp(1/x)), x) == {exp(x), exp(-x)}
    assert mrv(log(x**2 + 2*exp(exp(3*x**3*log(x)))),
               x) == {exp(exp(3*x**3*log(x)))}
    assert mrv(log(x - log(x))/log(x), x) == {x}
    assert mrv((exp(1/x - exp(-x)) - exp(1/x))*exp(x), x) == {exp(x), exp(-x)}
    assert mrv(1/exp(-x + exp(-x)) - exp(x), x) == {exp(x), exp(-x),
                                                    exp(x - exp(-x))}
    assert mrv(log(log(x*exp(x*exp(x)) + 1)), x) == {exp(x*exp(x))}
    assert mrv(exp(exp(log(log(x) + 1/x))), x) == {x}

    assert mrv((log(log(x) + log(log(x))) -
                log(log(x)))/log(log(x) + log(log(log(x))))*log(x), x) == {x}
    assert mrv(log(log(x*exp(x*exp(x)) + 1)) - exp(exp(log(log(x) + 1/x))),
               x) == {exp(x*exp(x))}

    # Gruntz: p47, 3.21
    h = exp(-x/(1 + exp(-x)))
    e = Mul(exp(h), exp(-x/(1 + h)), exp(exp(-x + h)),
            Pow(h, -2, evaluate=False), evaluate=False) - exp(x) + x
    assert mrv(e, x) == {exp(-x + h), exp(-x/(1 + h)), h, exp(x), exp(-x)}


def test_rewrite():
    assert rewrite(Integer(1), x, m) == (1, None)

    e = exp(x)
    assert rewrite(e, x, m) == (1/m, -x)
    e = exp(x**2)
    assert rewrite(e, x, m) == (1/m, -x**2)
    e = exp(x + 1/x)
    assert rewrite(e, x, m) == (1/m, -x - 1/x)
    e = 1/exp(-x + exp(-x)) - exp(x)
    assert rewrite(e, x, m) == (Add(-1/m, 1/(m*exp(m)), evaluate=False), -x)

    e = exp(x)*log(log(exp(x)))
    assert mrv(e, x) == {exp(x)}
    assert rewrite(e, x, m) == (1/m*log(x), -x)

    e = exp(-x + 1/x**2) - exp(x + 1/x)
    assert rewrite(e, x, m) == (Add(m, Mul(-1, exp(1/x + x**(-2))/m,
                                           evaluate=False),
                                    evaluate=False), -x + 1/x**2)


def test_leadterm():
    assert leadterm(Integer(1), x) == (1, 0)

    assert leadterm(-exp(1/x), x) == (-1, 0)
    assert leadterm(1/exp(-x + exp(-x)) - exp(x), x) == (-1, 0)
    assert leadterm((exp(1/x - exp(-x)) - exp(1/x))*exp(x), x) == (-exp(1/x), 0)

    # Gruntz: p51, 3.25
    assert leadterm((log(exp(x) + x) - x)/log(exp(x) + log(x))*exp(x),
                    x) == (1, 0)

    # Gruntz: p56, 3.27
    assert mrv(exp(-x + exp(-x)*exp(-x*log(x))), x) == {exp(-x*log(x))}
    assert leadterm(exp(-x + exp(-x)*exp(-x*log(x))), x) == (exp(-x), 0)


def test_limit():
    assert limit(x, x, oo) == oo
    assert limit(-x, x, oo) == -oo
    assert limit(-x, x, oo) == -oo
    assert limit((-x)**2, x, oo) == oo
    assert limit(-x**2, x, oo) == -oo
    assert limit((1/x)*log(1/x), x, oo) == 0  # Gruntz: p15, 2.11
    assert limit(1/x, x, oo) == 0  # issue sympy/sympy#11667
    assert limit(exp(x), x, oo) == oo
    assert limit(-exp(x), x, oo) == -oo
    assert limit(exp(x)/x, x, oo) == oo
    assert limit(1/x - exp(-x), x, oo) == 0
    assert limit(x + 1/x, x, oo) == oo
    assert limit(exp(1/log(1 - 1/x)), x, oo) == 0

    assert limit((1/x)**(1/x), x, oo) == 1  # Gruntz: p15, 2.11
    assert limit((exp(1/x) - 1)*x, x, oo) == 1
    assert limit(1 + 1/x, x, oo) == 1
    assert limit(-exp(1/x), x, oo) == -1
    assert limit(x + exp(-x), x, oo) == oo
    assert limit(x + exp(-x**2), x, oo) == oo
    assert limit(x + exp(-exp(x)), x, oo) == oo
    assert limit(13 + 1/x - exp(-x), x, oo) == 13

    assert limit(x - log(1 + exp(x)), x, oo) == 0
    assert limit(x - log(a + exp(x)), x, oo) == 0
    assert limit(exp(x)/(1 + exp(x)), x, oo) == 1
    assert limit(exp(x)/(a + exp(x)), x, oo) == 1

    assert limit((3**x + 5**x)**(1/x), x, oo) == 5  # issue sympy/sympy#3463

    assert limit(Ei(x + exp(-x))*exp(-x)*x, x, oo) == 1

    assert limit(1/li(x), x, oo) == 0
    assert limit(1/Li(x), x, oo) == 0

    # issue diofant/diofant#56
    assert limit((log(E + 1/x) - 1)**(1 - sqrt(E + 1/x)), x, oo) == oo

    # issue sympy/sympy#9471
    assert limit(27**log(x, 3)/x**3, x, oo) == 1
    assert limit(27**(log(x, 3) + 1)/x**3, x, oo) == 27

    # issue sympy/sympy#9449
    assert limit(x*(abs(1/x + y) - abs(y - 1/x))/2, x, oo) == sign(y)**-1

    # issue sympy/sympy#8481
    assert limit(m**x * exp(-m) / factorial(x), x, oo) == 0

    # issue sympy/sympy#4187
    assert limit(exp(1/x)*log(1/x) - Ei(1/x), x, oo) == -EulerGamma
    assert limit(exp(x)*log(x) - Ei(x), x, oo) == oo

    # issue sympy/sympy#10382
    assert limit(fibonacci(x + 1)/fibonacci(x), x, oo) == GoldenRatio

    assert limit(zeta(x), x, oo) == 1
    assert limit(zeta(m)*zeta(x), x, oo) == zeta(m)


def test_I():
    assert limit(I*x, x, oo) == I*oo
    assert limit(y*I*x, x, oo) == sign(y)*I*oo
    assert limit(y*3*I*x, x, oo) == sign(y)*I*oo
    assert limit(y*3*sin(I)*x, x, oo).simplify() == sign(y)*I*oo


def test_sympyissue_4814():
    assert limit((x + 1)**(1/log(x + 1)), x, oo) == E


def test_intractable():
    assert limit(1/gamma(x), x, oo) == 0
    assert limit(1/loggamma(x), x, oo) == 0
    assert limit(gamma(x)/loggamma(x), x, oo) == oo
    assert limit(exp(gamma(x))/gamma(x), x, oo) == oo
    assert limit(gamma(3 + 1/x), x, oo) == 2
    assert limit(gamma(Rational(1, 7) + 1/x), x, oo) == gamma(Rational(1, 7))
    assert limit(log(x**x)/log(gamma(x)), x, oo) == 1
    assert limit(log(gamma(gamma(x)))/exp(x), x, oo) == oo
    assert limit(acosh(1 + 1/x)*sqrt(x), x, oo) == sqrt(2)

    # issue sympy/sympy#10804
    assert limit(2*airyai(x)*root(x, 4) *
                 exp(2*x**Rational(3, 2)/3), x, oo) == 1/sqrt(pi)
    assert limit(airybi(x)*root(x, 4) *
                 exp(-2*x**Rational(3, 2)/3), x, oo) == 1/sqrt(pi)
    assert limit(airyai(1/x), x, oo) == (3**Rational(5, 6) *
                                         gamma(Rational(1, 3))/(6*pi))
    assert limit(airybi(1/x), x, oo) == cbrt(3)*gamma(Rational(1, 3))/(2*pi)
    assert limit(airyai(2 + 1/x), x, oo) == airyai(2)
    assert limit(airybi(2 + 1/x), x, oo) == airybi(2)

    # issue sympy/sympy#10976
    assert limit(erf(m/x)/erf(1/x), x, oo) == m

    assert limit(Max(x**2, x, exp(x))/x, x, oo) == oo


def test_branch_cuts():
    assert limit(sqrt(-1 + I/x), x, oo) == +I
    assert limit(sqrt(-1 - I/x), x, oo) == -I
    assert limit(log(-1 + I/x), x, oo) == +I*pi
    assert limit(log(-1 - I/x), x, oo) == -I*pi
    # Gruntz: p23
    assert limit(atan(2*I - x), x, 0, -1) == I*log(3)/2 - pi/2
    assert limit(atan(2*I - x), x, 0, +1) == I*log(3)/2 + pi/2


def test_aseries_trig():
    assert limit(1/log(atan(x)), x, oo) == -1/(-log(pi) + log(2))
    assert limit(1/acot(-x), x, oo) == -oo


def test_exp_log_series():
    assert limit(x/log(log(x*exp(x))), x, oo) == oo


def test_sympyissue_3644():
    assert limit(((x**7 + x + 1)/(2**x + x**2))**(-1/x), x, oo) == 2


def test_sympyissue_6843():
    n = Symbol('n', integer=True, positive=True)
    r = (n + 1)*(1 + 1/x)**(n + 1)/((1 + 1/x)**(n + 1) - 1) - (1 + 1/x)*x
    assert limit(r, x, oo) == n/2


def test_sympyissue_4190():
    assert limit(x - gamma(1/x), x, oo) == EulerGamma


def test_sympyissue_4109():
    assert limit(1/gamma(1/x), x, oo) == 0
    assert limit(gamma(1/x)/x, x, oo) == 1


def test_sympyissue_6682():
    assert limit(exp(2*Ei(-1/x))*x**2, x, oo) == exp(2*EulerGamma)


def test_sympyissue_7096():
    assert limit((-1/x)**-pi, x, oo) == oo*sign((-1)**(-pi))


def test_sympyissue_8462():
    assert limit(binomial(x, x/2), x, oo) == oo
    # issue sympy/sympy#10801
    assert limit(16**x/(x*binomial(2*x, x)**2), x, oo) == pi


def test_issue_74():
    assert limit(sign(log(1 + 1/x)), x, oo) == +1
    assert limit(sign(log(1 - 1/x)), x, oo) == -1
    assert limit(sign(sin(+1/x)), x, oo) == +1
    assert limit(sign(sin(-1/x)), x, oo) == -1
    assert limit(sign(tan(+1/x)), x, oo) == +1
    assert limit(sign(tan(-1/x)), x, oo) == -1
    assert limit(sign(cos(pi/2 + 1/x)), x, oo) == -1
    assert limit(sign(cos(pi/2 - 1/x)), x, oo) == +1


def test_issue_75():
    assert limit(abs(log(x)), x, oo) == oo
    assert limit(tan(abs(pi/2 + 1/x))/acosh(pi/2 + 1/x), x, oo) == -oo
    assert limit(tan(abs(pi/2 - 1/x))/acosh(pi/2 - 1/x), x, oo) == +oo

    assert limit(abs(log(2 + 1/x)) - log(2 + 1/x), x, oo) == 0
    assert limit(abs(log(2 - 1/x)) - log(2 - 1/x), x, oo) == 0


def test_sympyissue_8241():
    e = x/log(x)**(log(x)/(m*log(log(x))))
    assert isinstance(limit(e, x, oo), Limit)


def test_sympyissue_23845():
    e = ((-sqrt(5)*(-sqrt(5)/2 + Rational(3, 2))**x +
          5*(-sqrt(5)/2 + Rational(3, 2))**x +
          (sqrt(5)/2 + Rational(3, 2))**x*(sqrt(5) + 5)) /
         (-3*sqrt(5)*(-sqrt(5)/2 + Rational(3, 2))**x +
          5*(-sqrt(5)/2 + Rational(3, 2))**x +
          (5 + 3*sqrt(5))*(sqrt(5)/2 + Rational(3, 2))**x))
    assert limit(e, x, oo) == (sqrt(5) + 5)/(5 + 3*sqrt(5))
