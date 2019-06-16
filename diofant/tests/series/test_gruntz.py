"""
This test suite is testing the Gruntz algorithm implementation using the
bottom up approach.  See the documentation in the gruntz module.  The
algorithm itself is highly recursive by nature, so ``compare()`` is
logically the lowest part of the algorithm, yet in some sense it's the most
complex part, because it needs to calculate a limit to return the result.
"""

import pytest

from diofant import (Add, E, Ei, EulerGamma, GoldenRatio, I, Integer, Li,
                     Limit, Mul, Pow, Rational, Symbol, acosh, acot, airyai,
                     airybi, atan, binomial, cbrt, cos, cosh, coth, digamma,
                     erf, exp, factorial, fibonacci, gamma, li, log, loggamma,
                     oo, pi, root, sign, sin, sinh, sqrt, tan, tanh, zeta)
from diofant.series.gruntz import compare
from diofant.series.gruntz import limitinf as gruntz
from diofant.series.gruntz import mrv, mrv_leadterm, rewrite
from diofant.series.gruntz import sign as mrv_sign


__all__ = ()

x = Symbol('x', real=True, positive=True)
m = Symbol('m', real=True, positive=True)


@pytest.mark.slow
def test_gruntz_evaluation():
    # Gruntz' thesis pp. 122 to 123
    # 8.1
    assert gruntz(exp(x)*(exp(1/x - exp(-x)) - exp(1/x)), x) == -1
    # 8.2
    assert gruntz(exp(x)*(exp(1/x + exp(-x) + exp(-x**2))
                          - exp(1/x - exp(-exp(x)))), x) == 1
    # 8.3
    assert gruntz(exp(exp(x - exp(-x))/(1 - 1/x)) - exp(exp(x)), x) == oo
    # 8.4
    assert gruntz(exp(exp(exp(x)/(1 - 1/x)))
                  - exp(exp(exp(x)/(1 - 1/x - log(x)**(-log(x))))), x) == -oo
    # 8.5
    assert gruntz(exp(exp(exp(x + exp(-x)))) / exp(exp(exp(x))), x) == oo
    # 8.6
    assert gruntz(exp(exp(exp(x))) / exp(exp(exp(x - exp(-exp(x))))),
                  x) == oo
    # 8.7
    assert gruntz(exp(exp(exp(x))) / exp(exp(exp(x - exp(-exp(exp(x)))))),
                  x) == 1
    # 8.8
    assert gruntz(exp(exp(x)) / exp(exp(x - exp(-exp(exp(x))))), x) == 1
    # 8.9
    assert gruntz(log(x)**2 * exp(sqrt(log(x))*(log(log(x)))**2
                                  * exp(sqrt(log(log(x))) * (log(log(log(x))))**3)) / sqrt(x),
                  x) == 0
    # 8.10
    assert gruntz((x*log(x)*(log(x*exp(x) - x**2))**2)
                  / (log(log(x**2 + 2*exp(exp(3*x**3*log(x)))))), x) == Rational(1, 3)
    # 8.11
    assert gruntz((exp(x*exp(-x)/(exp(-x) + exp(-2*x**2/(x + 1)))) - exp(x))/x,
                  x) == -exp(2)
    # 8.12
    assert gruntz((3**x + 5**x)**(1/x), x) == 5
    # 8.13
    assert gruntz(x/log(x**(log(x**(log(2)/log(x))))), x) == oo
    # 8.14
    assert gruntz(exp(exp(2*log(x**5 + x)*log(log(x))))
                  / exp(exp(10*log(x)*log(log(x)))), x) == oo
    # 8.15
    assert gruntz(exp(exp(Rational(5, 2)*x**(-Rational(5, 7)) + Rational(21, 8)*x**Rational(6, 11)
                          + 2*x**(-8) + Rational(54, 17)*x**Rational(49, 45)))**8
                  / log(log(-log(Rational(4, 3)*x**(-Rational(5, 14)))))**Rational(7, 6), x) == oo
    # 8.16
    assert gruntz((exp(4*x*exp(-x)/(1/exp(x) + 1/exp(2*x**2/(x + 1)))) - exp(x))
                  / exp(x)**4, x) == 1
    # 8.17
    assert gruntz(exp(x*exp(-x)/(exp(-x) + exp(-2*x**2/(x + 1))))/exp(x), x) \
        == 1
    # 8.18
    assert gruntz((exp(exp(-x/(1 + exp(-x))))*exp(-x/(1 + exp(-x/(1 + exp(-x)))))
                   * exp(exp(-x + exp(-x/(1 + exp(-x))))))
                  / (exp(-x/(1 + exp(-x))))**2 - exp(x) + x, x) == 2
    # 8.19
    assert gruntz(log(x)*(log(log(x) + log(log(x))) - log(log(x)))
                  / (log(log(x) + log(log(log(x))))), x) == 1
    # 8.20
    assert gruntz(exp((log(log(x + exp(log(x)*log(log(x))))))
                      / (log(log(log(exp(x) + x + log(x)))))), x) == E
    # Another
    assert gruntz(exp(exp(exp(x + exp(-x)))) / exp(exp(x)), x) == oo


def test_gruntz_eval_special():
    # Gruntz, p. 126
    assert gruntz(exp(x)*(sin(1/x + exp(-x)) - sin(1/x + exp(-x**2))), x) == 1
    assert gruntz((erf(x - exp(-exp(x))) - erf(x)) * exp(exp(x)) * exp(x**2),
                  x) == -2/sqrt(pi)
    assert gruntz(exp(exp(x)) * (exp(sin(1/x + exp(-exp(x)))) - exp(sin(1/x))),
                  x) == 1
    assert gruntz(exp(x)*(gamma(x + exp(-x)) - gamma(x)), x) == oo
    assert gruntz(exp(exp(digamma(digamma(x))))/x, x) == exp(-Rational(1, 2))
    assert gruntz(exp(exp(digamma(log(x))))/x, x) == exp(-Rational(1, 2))
    assert gruntz(digamma(digamma(digamma(x))), x) == oo
    assert gruntz(loggamma(loggamma(x)), x) == oo
    assert gruntz(((gamma(x + 1/gamma(x)) - gamma(x))/log(x) - cos(1/x))
                  * x*log(x), x) == -Rational(1, 2)
    assert gruntz(x * (gamma(x - 1/gamma(x)) - gamma(x) + log(x)), x) \
        == Rational(1, 2)
    assert gruntz((gamma(x + 1/gamma(x)) - gamma(x)) / log(x), x) == 1


@pytest.mark.slow
def test_gruntz_eval_special_slow():
    assert gruntz(gamma(x + 1)/sqrt(2*pi)
                  - exp(-x)*(x**(x + Rational(1, 2)) + x**(x - Rational(1, 2))/12), x) == oo
    assert gruntz(exp(exp(exp(digamma(digamma(digamma(x))))))/x, x) == 0
    assert gruntz(exp(gamma(x - exp(-x))*exp(1/x)) - exp(gamma(x)), x) == oo
    assert gruntz(
        (Ei(x - exp(-exp(x))) - Ei(x))*exp(-x)*exp(exp(x))*x, x) == -1
    assert gruntz(
        exp((log(2) + 1)*x) * (zeta(x + exp(-x)) - zeta(x)), x) == -log(2)

    # TODO 8.36 - 8.37 (bessel, max-min)


def test_gruntz_other():
    assert gruntz(sqrt(log(x + 1)) - sqrt(log(x)), x) == 0  # p12, 2.5
    y = Symbol('y')
    assert gruntz(((1 + 1/x)**y - 1)*x, x) == y  # p12, 2.6
    n = Symbol('n', integer=True)
    assert gruntz(x**n/exp(x), x) == 0  # p14, 2.9
    assert gruntz((1 + 1/x)*x - 1/log(1 + 1/x), x) == Rational(1, 2)  # p15, 2.10
    m = Symbol('m', integer=True)
    assert gruntz((root(1 + 1/x, n) - 1)/(root(1 + 1/x, m) - 1),
                  x) == m/n  # p13, 2.7


def test_gruntz_hyperbolic():
    assert gruntz(cosh(x), x) == oo
    assert gruntz(cosh(-x), x) == oo
    assert gruntz(sinh(x), x) == oo
    assert gruntz(sinh(-x), x) == -oo
    assert gruntz(2*cosh(x)*exp(x), x) == oo
    assert gruntz(2*cosh(-x)*exp(-x), x) == 1
    assert gruntz(2*sinh(x)*exp(x), x) == oo
    assert gruntz(2*sinh(-x)*exp(-x), x) == -1
    assert gruntz(tanh(x), x) == 1
    assert gruntz(tanh(-x), x) == -1
    assert gruntz(coth(x), x) == 1
    assert gruntz(coth(-x), x) == -1


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


def test_sign():
    assert mrv_sign(Integer(0), x) == 0
    assert mrv_sign(Integer(3), x) == 1
    assert mrv_sign(Integer(-5), x) == -1
    assert mrv_sign(log(x), x) == 1
    assert mrv_sign(exp(-x), x) == 1
    assert mrv_sign(exp(x), x) == 1
    assert mrv_sign(-exp(x), x) == -1
    assert mrv_sign(3 - 1/x, x) == 1
    assert mrv_sign(-3 - 1/x, x) == -1
    assert mrv_sign(sin(1/x), x) == 1
    assert mrv_sign(x**2, x) == 1
    assert mrv_sign(x**5, x) == 1

    assert mrv_sign(x, x) == 1
    assert mrv_sign(-x, x) == -1
    y = Symbol("y", positive=True)
    assert mrv_sign(y, x) == 1
    assert mrv_sign(-y, x) == -1
    assert mrv_sign(y*x, x) == 1
    assert mrv_sign(-y*x, x) == -1


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

    assert mrv(
        exp(-x + 1/x**2) - exp(x + 1/x), x) == {exp(x + 1/x), exp(1/x**2 - x)}

    assert mrv(exp(x**2) + x*exp(x) + exp(x*log(log(x)))/x, x) == {exp(x**2)}
    assert mrv(
        exp(x)*(exp(1/x + exp(-x)) - exp(1/x)), x) == {exp(x), exp(-x)}
    assert mrv(log(
        x**2 + 2*exp(exp(3*x**3*log(x)))), x) == {exp(exp(3*x**3*log(x)))}
    assert mrv(log(x - log(x))/log(x), x) == {x}
    assert mrv(
        (exp(1/x - exp(-x)) - exp(1/x))*exp(x), x) == {exp(x), exp(-x)}
    assert mrv(
        1/exp(-x + exp(-x)) - exp(x), x) == {exp(x), exp(-x), exp(x - exp(-x))}
    assert mrv(log(log(x*exp(x*exp(x)) + 1)), x) == {exp(x*exp(x))}
    assert mrv(exp(exp(log(log(x) + 1/x))), x) == {x}

    assert mrv((log(log(x) + log(log(x))) - log(log(x)))
               / log(log(x) + log(log(log(x))))*log(x), x) == {x}
    assert mrv(log(log(x*exp(x*exp(x)) + 1)) - exp(exp(log(log(x) + 1/x))), x) == \
        {exp(x*exp(x))}

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


def test_mrv_leadterm():
    assert mrv_leadterm(Integer(1), x) == (1, 0)

    assert mrv_leadterm(-exp(1/x), x) == (-1, 0)
    assert mrv_leadterm(1/exp(-x + exp(-x)) - exp(x), x) == (-1, 0)
    assert mrv_leadterm(
        (exp(1/x - exp(-x)) - exp(1/x))*exp(x), x) == (-exp(1/x), 0)

    # Gruntz: p51, 3.25
    assert mrv_leadterm((log(exp(x) + x) - x)/log(exp(x) + log(x))*exp(x), x) == \
        (1, 0)

    # Gruntz: p56, 3.27
    assert mrv(exp(-x + exp(-x)*exp(-x*log(x))), x) == {exp(-x*log(x))}
    assert mrv_leadterm(exp(-x + exp(-x)*exp(-x*log(x))), x) == (exp(-x), 0)


def test_limit():
    assert gruntz(x, x) == oo
    assert gruntz(-x, x) == -oo
    assert gruntz(-x, x) == -oo
    assert gruntz((-x)**2, x) == oo
    assert gruntz(-x**2, x) == -oo
    assert gruntz((1/x)*log(1/x), x) == 0  # Gruntz: p15, 2.11
    assert gruntz(1/x, x) == 0
    assert gruntz(exp(x), x) == oo
    assert gruntz(-exp(x), x) == -oo
    assert gruntz(exp(x)/x, x) == oo
    assert gruntz(1/x - exp(-x), x) == 0
    assert gruntz(x + 1/x, x) == oo

    assert gruntz((1/x)**(1/x), x) == 1  # Gruntz: p15, 2.11
    assert gruntz((exp(1/x) - 1)*x, x) == 1
    assert gruntz(1 + 1/x, x) == 1
    assert gruntz(-exp(1/x), x) == -1
    assert gruntz(x + exp(-x), x) == oo
    assert gruntz(x + exp(-x**2), x) == oo
    assert gruntz(x + exp(-exp(x)), x) == oo
    assert gruntz(13 + 1/x - exp(-x), x) == 13

    a = Symbol('a')
    assert gruntz(x - log(1 + exp(x)), x) == 0
    assert gruntz(x - log(a + exp(x)), x) == 0
    assert gruntz(exp(x)/(1 + exp(x)), x) == 1
    assert gruntz(exp(x)/(a + exp(x)), x) == 1

    assert gruntz((3**x + 5**x)**(1/x), x) == 5  # issue sympy/sympy#3463

    assert gruntz(Ei(x + exp(-x))*exp(-x)*x, x) == 1

    assert gruntz(1/li(x), x) == 0
    assert gruntz(1/Li(x), x) == 0

    # issue diofant/diofant#56
    assert gruntz((log(E + 1/x) - 1)**(1 - sqrt(E + 1/x)), x) == oo

    # issue sympy/sympy#9471
    assert gruntz((((27**(log(x, 3))))/x**3), x) == 1
    assert gruntz((((27**(log(x, 3) + 1)))/x**3), x) == 27

    # issue sympy/sympy#9449
    y = Symbol('y')
    assert gruntz(x*(abs(1/x + y) - abs(y - 1/x))/2, x) == sign(y)

    # issue sympy/sympy#8481
    assert gruntz(m**x * exp(-m) / factorial(x), x) == 0

    # issue sympy/sympy#4187
    assert gruntz(exp(1/x)*log(1/x) - Ei(1/x), x) == -EulerGamma
    assert gruntz(exp(x)*log(x) - Ei(x), x) == oo

    # issue sympy/sympy#10382
    assert gruntz(fibonacci(x + 1)/fibonacci(x), x) == GoldenRatio

    assert gruntz(zeta(x), x) == 1
    assert gruntz(zeta(m)*zeta(x), x) == zeta(m)


def test_I():
    y = Symbol("y")
    assert gruntz(I*x, x) == I*oo
    assert gruntz(y*I*x, x) == sign(y)*I*oo
    assert gruntz(y*3*I*x, x) == sign(y)*I*oo
    assert gruntz(y*3*sin(I)*x, x).simplify() == sign(y)*I*oo


def test_sympyissue_4814():
    assert gruntz((x + 1)**(1/log(x + 1)), x) == E


def test_intractable():
    assert gruntz(1/gamma(x), x) == 0
    assert gruntz(1/loggamma(x), x) == 0
    assert gruntz(gamma(x)/loggamma(x), x) == oo
    assert gruntz(exp(gamma(x))/gamma(x), x) == oo
    assert gruntz(gamma(3 + 1/x), x) == 2
    assert gruntz(gamma(Rational(1, 7) + 1/x), x) == gamma(Rational(1, 7))
    assert gruntz(log(x**x)/log(gamma(x)), x) == 1
    assert gruntz(log(gamma(gamma(x)))/exp(x), x) == oo

    # issue sympy/sympy#10804
    assert gruntz(2*airyai(x)*root(x, 4) *
                  exp(2*x**Rational(3, 2)/3), x) == 1/sqrt(pi)
    assert gruntz(airybi(x)*root(x, 4) *
                  exp(-2*x**Rational(3, 2)/3), x) == 1/sqrt(pi)
    assert gruntz(airyai(1/x), x) == (3**Rational(5, 6) *
                                      gamma(Rational(1, 3))/(6*pi))
    assert gruntz(airybi(1/x), x) == cbrt(3)*gamma(Rational(1, 3))/(2*pi)
    assert gruntz(airyai(2 + 1/x), x) == airyai(2)
    assert gruntz(airybi(2 + 1/x), x) == airybi(2)


def test_branch_cuts():
    assert gruntz(sqrt(-1 + I/x), x) == +I
    assert gruntz(sqrt(-1 - I/x), x) == -I
    assert gruntz(log(-1 + I/x), x) == +I*pi
    assert gruntz(log(-1 - I/x), x) == -I*pi


def test_aseries_trig():
    assert gruntz(1/log(atan(x)), x) == -1/(-log(pi) + log(2))
    assert gruntz(1/acot(-x), x) == -oo


def test_exp_log_series():
    assert gruntz(x/log(log(x*exp(x))), x) == oo


def test_sympyissue_3644():
    assert gruntz(((x**7 + x + 1)/(2**x + x**2))**(-1/x), x) == 2


def test_sympyissue_6843():
    n = Symbol('n', integer=True, positive=True)
    r = (n + 1)*(1 + 1/x)**(n + 1)/((1 + 1/x)**(n + 1) - 1) - (1 + 1/x)*x
    assert gruntz(r, x).simplify() == n/2


def test_sympyissue_4190():
    assert gruntz(x - gamma(1/x), x) == EulerGamma


def test_sympyissue_5172():
    n = Symbol('n', real=True, positive=True)
    r = Symbol('r', positive=True)
    c = Symbol('c')
    p = Symbol('p', positive=True)
    m = Symbol('m', negative=True)
    expr = ((2*n*(n - r + 1)/(n + r*(n - r + 1)))**c +
            (r - 1)*(n*(n - r + 2)/(n + r*(n - r + 1)))**c - n)/(n**c - n)
    expr = expr.subs({c: c + 1})
    assert gruntz(expr.subs({c: m}), n) == 1
    assert gruntz(expr.subs({c: p}), n).simplify() == \
        (2**(p + 1) + r - 1)/(r + 1)**(p + 1)


def test_sympyissue_4109():
    assert gruntz(1/gamma(1/x), x) == 0
    assert gruntz(gamma(1/x)/x, x) == 1


def test_sympyissue_6682():
    assert gruntz(exp(2*Ei(-1/x))*x**2, x) == exp(2*EulerGamma)


def test_sympyissue_7096():
    assert gruntz((-1/x)**-pi, x) == oo*sign((-1)**(-pi))


def test_sympyissue_8462():
    assert gruntz(binomial(x, x/2), x) == oo
    # issue sympy/sympy#10801
    assert gruntz(16**x/(x*binomial(2*x, x)**2), x) == pi


def test_diofantissue_74():
    assert gruntz(sign(log(1 + 1/x)), x) == +1
    assert gruntz(sign(log(1 - 1/x)), x) == -1
    assert gruntz(sign(sin(+1/x)), x) == +1
    assert gruntz(sign(sin(-1/x)), x) == -1
    assert gruntz(sign(tan(+1/x)), x) == +1
    assert gruntz(sign(tan(-1/x)), x) == -1
    assert gruntz(sign(cos(pi/2 + 1/x)), x) == -1
    assert gruntz(sign(cos(pi/2 - 1/x)), x) == +1


def test_diofantissue_75():
    assert gruntz(abs(log(x)), x) == oo
    assert gruntz(tan(abs(pi/2 + 1/x))/acosh(pi/2 + 1/x), x) == -oo
    assert gruntz(tan(abs(pi/2 - 1/x))/acosh(pi/2 - 1/x), x) == +oo

    assert gruntz(abs(log(2 + 1/x)) - log(2 + 1/x), x) == 0
    assert gruntz(abs(log(2 - 1/x)) - log(2 - 1/x), x) == 0


def test_sympyissue_8241():
    e = x/log(x)**(log(x)/(m*log(log(x))))
    pytest.raises(NotImplementedError, lambda: gruntz(e, x))
    assert isinstance(e.limit(x, oo), Limit)


def test_sympyissue_10976():
    assert gruntz(erf(m/x)/erf(1/x), x) == m
