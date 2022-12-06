from collections import defaultdict
from random import randrange, uniform

import pytest

from diofant import (E1, RR, Abs, Add, And, Chi, Ci, Ei, Heaviside, I, Integer,
                     Integral, Mul, Piecewise, Rational, Shi, Si, Symbol,
                     acosh, acoth, arg, asin, atan, besseli, besselj, cbrt,
                     combsimp, cos, cosh, default_sort_key, erf, exp,
                     exp_polar, expand, expand_func, expand_mul, expint,
                     fourier_transform, fresnelc, fresnels, gamma, hyper,
                     hyperexpand, integrate, laplace_transform, log,
                     lowergamma, meijerg, nan, oo, pi, piecewise_fold,
                     polygamma, powdenest, powsimp, re, simplify, sin, sinh,
                     sqrt, symbols, unpolarify)
from diofant.abc import R, a, b, c, d, r, s, t, x, y, z
from diofant.integrals.meijerint import (_create_lookup_table, _inflate_g,
                                         _rewrite1, _rewrite_single,
                                         meijerint_definite,
                                         meijerint_indefinite,
                                         meijerint_inversion)
from diofant.integrals.meijerint import z as z_dummy
from diofant.utilities.randtest import random_complex_number as randcplx
from diofant.utilities.randtest import verify_numerically


__all__ = ()


def test_rewrite_single():
    def t(expr, c, m):
        e = _rewrite_single(meijerg([a], [b], [c], [d], expr), x)
        assert e is not None
        assert isinstance(e[0][0][2], meijerg)
        assert e[0][0][2].argument.as_coeff_mul(x) == (c, (m,))

    def tn(expr):
        assert _rewrite_single(meijerg([a], [b], [c], [d], expr), x) is None

    t(x, 1, x)
    t(x**2, 1, x**2)
    t(x**2 + y*x**2, y + 1, x**2)
    tn(x**2 + x)
    tn(x**y)

    def u(expr, x):
        r = _rewrite_single(expr, x)
        e = Add(*[res[0]*res[2] for res in r[0]]).replace(
            exp_polar, exp)  # XXX Hack?
        assert verify_numerically(e, expr, x)

    u(exp(-x)*sin(x), x)

    # The following has stopped working because hyperexpand changed slightly.
    # It is probably not worth fixing
    # u(exp(-x)*sin(x)*cos(x), x)

    # This one cannot be done numerically, since it comes out as a g-function
    # of argument 4*pi
    # NOTE This also tests a bug in inverse mellin transform (which used to
    #      turn exp(4*pi*I*t) into a factor of exp(4*pi*I)**t instead of
    #      exp_polar).
    # u(exp(x)*sin(x), x)
    assert _rewrite_single(exp(x)*sin(x), x) == \
        ([(-sqrt(2)/(2*sqrt(pi)), 0,
           meijerg(((-Rational(1, 2), 0, Rational(1, 4), Rational(1, 2), Rational(3, 4)), (1,)),
                   ((), (-Rational(1, 2), 0)), 64*exp_polar(-4*I*pi)/x**4))], True)


def test_rewrite1():
    assert _rewrite1(x**3*meijerg([a], [b], [c], [d], x**2 + y*x**2)*5, x) == \
        (5, x**3, [(1, 0, meijerg([a], [b], [c], [d], x**2*(y + 1)))], True)


def test_meijerint_indefinite_numerically():
    def t(fac, arg):
        g = meijerg([a], [b], [c], [d], arg)*fac
        subs = {a: randcplx()/10, b: randcplx()/10 + I,
                c: randcplx(), d: randcplx()}
        integral = meijerint_indefinite(g, x)
        assert integral is not None
        assert verify_numerically(g.subs(subs), integral.diff(x).subs(subs), x)
    t(1, x)
    t(2, x)
    t(1, 2*x)
    t(1, x**2)
    t(5, x**Rational(3, 2))
    t(x**3, x)
    t(3*x**Rational(3, 2), 4*x**Rational(7, 3))


def test_meijerint_definite():
    v, b = meijerint_definite(x, x, 0, 0)
    assert v.is_zero
    assert b is True
    v, b = meijerint_definite(x, x, oo, oo)
    assert v.is_zero
    assert b is True


def test_inflate():
    subs = {a: randcplx()/10, b: randcplx()/10 + I, c: randcplx(),
            d: randcplx(), y: randcplx()/10}

    def t(a, b, arg, n):
        m1 = meijerg(a, b, arg)
        m2 = Mul(*_inflate_g(m1, n))
        # NOTE: (the random number)**9 must still be on the principal sheet.
        # Thus make b&d small to create random numbers of small imaginary part.
        return verify_numerically(m1.subs(subs), m2.subs(subs), x, b=0.1, d=-0.1)
    assert t([[a], [b]], [[c], [d]], x, 3)
    assert t([[a, y], [b]], [[c], [d]], x, 3)
    assert t([[a], [b]], [[c, y], [d]], 2*x**3, 3)


def test_recursive():
    a, b, c = symbols('a b c', positive=True)
    r = exp(-(x - a)**2)*exp(-(x - b)**2)
    e = integrate(r, (x, 0, oo), meijerg=True)
    assert simplify(e.expand()) == (
        sqrt(2)*sqrt(pi)*(
            (erf(sqrt(2)*(a + b)/2) + 1)*exp(-a**2/2 + a*b - b**2/2))/4)
    e = integrate(exp(-(x - a)**2)*exp(-(x - b)**2)*exp(c*x), (x, 0, oo), meijerg=True)
    assert simplify(e) == (
        sqrt(2)*sqrt(pi)*(erf(sqrt(2)*(2*a + 2*b + c)/4) + 1)*exp(-a**2 - b**2
                                                                  + (2*a + 2*b + c)**2/8)/4)
    assert simplify(integrate(exp(-(x - a - b - c)**2), (x, 0, oo), meijerg=True)) == \
        sqrt(pi)/2*(1 + erf(a + b + c))
    assert simplify(integrate(exp(-(x + a + b + c)**2), (x, 0, oo), meijerg=True)) == \
        sqrt(pi)/2*(1 - erf(a + b + c))


def test_meijerint():
    s, t, mu = symbols('s t mu', extended_real=True)
    assert integrate(meijerg([], [], [0], [], s*t)
                     * meijerg([], [], [mu/2], [-mu/2], t**2/4),
                     (t, 0, oo)).is_Piecewise
    s = symbols('s', positive=True)
    assert integrate(x**s*meijerg([[], []], [[0], []], x), (x, 0, oo)) == \
        gamma(s + 1)
    assert integrate(x**s*meijerg([[], []], [[0], []], x), (x, 0, oo),
                     meijerg=True) == gamma(s + 1)
    assert isinstance(integrate(x**s*meijerg([[], []], [[0], []], x),
                                (x, 0, oo), meijerg=False),
                      Integral)

    assert meijerint_indefinite(exp(x), x) == exp(x)

    # issue sympy/sympy#8368
    assert meijerint_indefinite(cosh(x)*exp(-x*t), x) == (
        (exp(x*(-t + 1))*(-t - 1) + exp(-x*(t + 1))*(-t + 1))/(t**2 - 1)/2)

    # TODO what simplifications should be done automatically?
    # This tests "extra case" for antecedents_1.
    a, b = symbols('a b', positive=True)
    assert simplify(meijerint_definite(x**a, x, 0, b)[0]) == \
        b**(a + 1)/(a + 1)

    # This tests various conditions and expansions:
    assert meijerint_definite((x + 1)**3*exp(-x), x, 0, oo) == (16, True)

    # Again, how about simplifications?
    sigma, mu = symbols('sigma mu', positive=True)
    i, c = meijerint_definite(exp(-((x - mu)/(2*sigma))**2), x, 0, oo)
    assert simplify(i) == sqrt(pi)*sigma*(erf(mu/(2*sigma)) + 1)
    assert c

    i, _ = meijerint_definite(exp(-mu*x)*exp(sigma*x), x, 0, oo)
    # TODO it would be nice to test the condition
    assert simplify(i) == 1/(mu - sigma)

    # Test substitutions to change limits
    assert meijerint_definite(exp(x), x, -oo, 2) == (exp(2), True)
    # Note: causes a NaN in _check_antecedents
    assert expand(meijerint_definite(exp(x), x, 0, I)[0]) == exp(I) - 1
    assert expand(meijerint_definite(exp(-x), x, 0, x)[0]) == \
        1 - exp(-exp(I*arg(x))*abs(x))

    # Test -oo to oo
    assert meijerint_definite(exp(-x**2), x, -oo, oo) == (sqrt(pi), True)
    assert meijerint_definite(exp(-abs(x)), x, -oo, oo) == (2, True)
    assert meijerint_definite(exp(-(2*x - 3)**2), x, -oo, oo) == \
        (sqrt(pi)/2, True)
    assert meijerint_definite(exp(-abs(2*x - 3)), x, -oo, oo) == (1, True)
    assert meijerint_definite(exp(-((x - mu)/sigma)**2/2)/sqrt(2*pi*sigma**2),
                              x, -oo, oo) == (1, True)

    # Test one of the extra conditions for 2 g-functinos
    assert meijerint_definite(exp(-x)*sin(x), x, 0, oo) == (Rational(1, 2), True)

    # Test a bug
    def res(n):
        return (1/(1 + x**2)).diff((x, n)).subs({x: 1})*(-1)**n
    for n in range(6):
        assert integrate(exp(-x)*sin(x)*x**n, (x, 0, oo), meijerg=True) == \
            res(n)

    # This used to test trigexpand... now it is done by linear substitution
    assert simplify(integrate(exp(-x)*sin(x + a), (x, 0, oo), meijerg=True)
                    ) == sqrt(2)*sin(a + pi/4)/2

    # Test the condition 14 from prudnikov.
    # (This is besselj*besselj in disguise, to stop the product from being
    #  recognised in the tables.)
    a, b, s = symbols('a b s')
    assert meijerint_definite(meijerg([], [], [a/2], [-a/2], x/4)
                              * meijerg([], [], [b/2], [-b/2], x/4)*x**(s - 1), x, 0, oo) == \
        (4*2**(2*s - 2)*gamma(-2*s + 1)*gamma(a/2 + b/2 + s)
         / (gamma(-a/2 + b/2 - s + 1)*gamma(a/2 - b/2 - s + 1)
            * gamma(a/2 + b/2 - s + 1)),
            And(0 < -2*re(4*s) + 8, 0 < re(a/2 + b/2 + s), re(2*s) < 1))

    # test a bug
    assert integrate(sin(x**a)*sin(x**b), (x, 0, oo), meijerg=True) == \
        Integral(sin(x**a)*sin(x**b), (x, 0, oo))

    # test better hyperexpand
    assert integrate(exp(-x**2)*log(x), (x, 0, oo), meijerg=True) == \
        (sqrt(pi)*polygamma(0, Rational(1, 2))/4).expand()

    # Test hyperexpand bug.
    n = symbols('n', integer=True)
    assert simplify(integrate(exp(-x)*x**n, x, meijerg=True)) == \
        lowergamma(n + 1, x)

    # Test a bug with argument 1/x
    alpha = symbols('alpha', positive=True)
    assert meijerint_definite((2 - x)**alpha*sin(alpha/x), x, 0, 2) == \
        (sqrt(pi)*alpha*gamma(alpha + 1)*meijerg(((), (alpha/2 + Rational(1, 2),
                                                       alpha/2 + 1)), ((0, 0, Rational(1, 2)), (-Rational(1, 2),)), alpha**2/16)/4, True)

    # test a bug related to 3016
    a, s = symbols('a s', positive=True)
    assert simplify(integrate(x**s*exp(-a*x**2), (x, -oo, oo))) == \
        a**(-s/2 - Rational(1, 2))*((-1)**s + 1)*gamma(s/2 + Rational(1, 2))/2

    # issue sympy/sympy#6348
    assert integrate(exp(I*x)/(1 + x**2),
                     (x, -oo, oo)).simplify().rewrite(exp) == pi*exp(-1)

    # issue sympy/sympy#13536
    a = Symbol('a', real=True, positive=True)
    assert integrate(1/x**2, (x, oo, a)) == -1/a


def test_bessel():
    assert simplify(integrate(besselj(a, z)*besselj(b, z)/z, (z, 0, oo),
                              meijerg=True, conds='none')) == \
        2*sin(pi*(a/2 - b/2))/(pi*(a - b)*(a + b))
    assert simplify(integrate(besselj(a, z)*besselj(a, z)/z, (z, 0, oo),
                              meijerg=True, conds='none')) == 1/(2*a)

    # TODO more orthogonality integrals

    assert simplify(integrate(sin(z*x)*(x**2 - 1)**(-(y + Rational(1, 2))),
                              (x, 1, oo), meijerg=True, conds='none')
                    * 2/((z/2)**y*sqrt(pi)*gamma(Rational(1, 2) - y))) == \
        besselj(y, z)

    # Werner Rosenheinrich
    # SOME INDEFINITE INTEGRALS OF BESSEL FUNCTIONS

    assert integrate(x*besselj(0, x), x, meijerg=True) == x*besselj(1, x)
    assert integrate(x*besseli(0, x), x, meijerg=True) == x*besseli(1, x)
    # TODO can do higher powers, but come out as high order ... should they be
    #      reduced to order 0, 1?
    assert integrate(besselj(1, x), x, meijerg=True) == -besselj(0, x)
    assert integrate(besselj(1, x)**2/x, x, meijerg=True) == \
        -(besselj(0, x)**2 + besselj(1, x)**2)/2
    # TODO more besseli when tables are extended or recursive mellin works
    assert integrate(besselj(0, x)**2/x**2, x, meijerg=True) == \
        -2*x*besselj(0, x)**2 - 2*x*besselj(1, x)**2 \
        + 2*besselj(0, x)*besselj(1, x) - besselj(0, x)**2/x
    assert integrate(besselj(0, x)*besselj(1, x), x, meijerg=True) == \
        -besselj(0, x)**2/2
    assert integrate(x**2*besselj(0, x)*besselj(1, x), x, meijerg=True) == \
        x**2*besselj(1, x)**2/2
    assert integrate(besselj(0, x)*besselj(1, x)/x, x, meijerg=True) == \
        (x*besselj(0, x)**2 + x*besselj(1, x)**2 -
            besselj(0, x)*besselj(1, x))
    # TODO how does besselj(0, a*x)*besselj(0, b*x) work?
    # TODO how does besselj(0, x)**2*besselj(1, x)**2 work?
    # TODO sin(x)*besselj(0, x) etc come out a mess
    # TODO can x*log(x)*besselj(0, x) be done?
    # TODO how does besselj(1, x)*besselj(0, x+a) work?
    # TODO more indefinite integrals when struve functions etc are implemented

    # test a substitution
    assert integrate(besselj(1, x**2)*x, x, meijerg=True) == \
        -besselj(0, x**2)/2


def test_inversion():
    def inv(f):
        return piecewise_fold(meijerint_inversion(f, s, t))
    assert inv(1/(s**2 + 1)) == sin(t)*Heaviside(t)
    assert inv(s/(s**2 + 1)) == cos(t)*Heaviside(t)
    assert inv(exp(-s)/s) == Heaviside(t - 1)
    assert inv(1/sqrt(1 + s**2)) == besselj(0, t)*Heaviside(t)

    # Test some antcedents checking.
    assert meijerint_inversion(sqrt(s)/sqrt(1 + s**2), s, t) is None
    assert inv(exp(s**2)) is None
    assert meijerint_inversion(exp(-s**2), s, t) is None


@pytest.mark.slow
def test_lookup_table():
    table = defaultdict(list)
    _create_lookup_table(table)
    for _, l in sorted(table.items(), key=default_sort_key):
        for formula, terms, *_ in sorted(l, key=default_sort_key):
            subs = {}
            for a in list(formula.free_symbols) + [z_dummy]:
                if hasattr(a, 'properties') and a.properties:
                    # these Wilds match positive integers
                    subs[a] = randrange(1, 10)
                else:
                    subs[a] = uniform(1.5, 2.0)
            if not isinstance(terms, list):
                terms = terms(subs)

            # First test that hyperexpand can do this.
            expanded = [hyperexpand(g) for (_, g) in terms]
            assert all(x.is_Piecewise or not x.has(meijerg) for x in expanded)

            # Now test that the meijer g-function is indeed as advertised.
            expanded = Add(*[f*x for (f, x) in terms])
            a, b = formula.evalf(subs=subs, strict=False), expanded.evalf(subs=subs, strict=False)
            r = min(abs(a), abs(b))
            if r < 1:
                assert abs(a - b).evalf(strict=False) <= 1e-10
            else:
                assert (abs(a - b)/r).evalf(strict=False) <= 1e-10


def test_branch_bug():
    # TODO combsimp cannot prove that the factor is unity
    assert powdenest(integrate(erf(x**3), x, meijerg=True).diff(x),
                     polar=True) == 2*erf(x**3)*gamma(Rational(2, 3))/3/gamma(Rational(5, 3))
    assert integrate(erf(x**3), x, meijerg=True) == \
        2*x*erf(x**3)*gamma(Rational(2, 3))/(3*gamma(Rational(5, 3))) \
        - 2*gamma(Rational(2, 3))*lowergamma(Rational(2, 3), x**6)/(3*sqrt(pi)*gamma(Rational(5, 3)))


def test_linear_subs():
    assert integrate(sin(x - 1), x, meijerg=True) == -cos(1 - x)
    assert integrate(besselj(1, x - 1), x, meijerg=True) == -besselj(0, 1 - x)


@pytest.mark.slow
def test_probability():
    # various integrals from probability theory
    mu1, mu2 = symbols('mu1 mu2', real=True, nonzero=True)
    sigma1, sigma2 = symbols('sigma1 sigma2', real=True,
                             nonzero=True, positive=True)
    rate = Symbol('lambda', real=True, positive=True)

    def normal(x, mu, sigma):
        return 1/sqrt(2*pi*sigma**2)*exp(-(x - mu)**2/2/sigma**2)

    def exponential(x, rate):
        return rate*exp(-rate*x)

    assert integrate(normal(x, mu1, sigma1), (x, -oo, oo), meijerg=True) == 1
    assert integrate(x*normal(x, mu1, sigma1), (x, -oo, oo), meijerg=True) == \
        mu1
    assert integrate(x**2*normal(x, mu1, sigma1), (x, -oo, oo), meijerg=True) \
        == mu1**2 + sigma1**2
    assert integrate(x**3*normal(x, mu1, sigma1), (x, -oo, oo), meijerg=True) \
        == mu1**3 + 3*mu1*sigma1**2
    assert integrate(normal(x, mu1, sigma1)*normal(y, mu2, sigma2),
                     (x, -oo, oo), (y, -oo, oo), meijerg=True) == 1
    assert integrate(x*normal(x, mu1, sigma1)*normal(y, mu2, sigma2),
                     (x, -oo, oo), (y, -oo, oo), meijerg=True) == mu1
    assert integrate(y*normal(x, mu1, sigma1)*normal(y, mu2, sigma2),
                     (x, -oo, oo), (y, -oo, oo), meijerg=True) == mu2
    assert integrate(x*y*normal(x, mu1, sigma1)*normal(y, mu2, sigma2),
                     (x, -oo, oo), (y, -oo, oo), meijerg=True) == mu1*mu2
    assert integrate((x + y + 1)*normal(x, mu1, sigma1)*normal(y, mu2, sigma2),
                     (x, -oo, oo), (y, -oo, oo), meijerg=True) == 1 + mu1 + mu2
    assert integrate((x + y - 1)*normal(x, mu1, sigma1)*normal(y, mu2, sigma2),
                     (x, -oo, oo), (y, -oo, oo), meijerg=True) == \
        -1 + mu1 + mu2

    i = integrate(x**2*normal(x, mu1, sigma1)*normal(y, mu2, sigma2),
                  (x, -oo, oo), (y, -oo, oo), meijerg=True)
    assert not i.has(Abs)
    assert simplify(i) == mu1**2 + sigma1**2
    assert integrate(y**2*normal(x, mu1, sigma1)*normal(y, mu2, sigma2),
                     (x, -oo, oo), (y, -oo, oo), meijerg=True) == \
        sigma2**2 + mu2**2

    assert integrate(exponential(x, rate), (x, 0, oo), meijerg=True) == 1
    assert integrate(x*exponential(x, rate), (x, 0, oo), meijerg=True) == \
        1/rate
    assert integrate(x**2*exponential(x, rate), (x, 0, oo), meijerg=True) == \
        2/rate**2

    def E(expr):
        res1 = integrate(expr*exponential(x, rate)*normal(y, mu1, sigma1),
                         (x, 0, oo), (y, -oo, oo), meijerg=True)
        res2 = integrate(expr*exponential(x, rate)*normal(y, mu1, sigma1),
                         (y, -oo, oo), (x, 0, oo), meijerg=True)
        assert expand_mul(res1) == expand_mul(res2)
        return res1

    assert E(1) == 1
    assert E(x*y) == mu1/rate
    assert E(x*y**2) == mu1**2/rate + sigma1**2/rate
    ans = sigma1**2 + 1/rate**2
    assert simplify(E((x + y + 1)**2) - E(x + y + 1)**2) == ans
    assert simplify(E((x + y - 1)**2) - E(x + y - 1)**2) == ans
    assert simplify(E((x + y)**2) - E(x + y)**2) == ans

    # Beta' distribution
    alpha, beta = symbols('alpha beta', positive=True, real=True)
    betadist = x**(alpha - 1)*(1 + x)**(-alpha - beta)*gamma(alpha + beta) \
        / gamma(alpha)/gamma(beta)
    assert integrate(betadist, (x, 0, oo), meijerg=True) == 1
    i = integrate(x*betadist, (x, 0, oo), meijerg=True, conds='separate')
    assert (combsimp(i[0]), i[1]) == (alpha/(beta - 1), 1 < beta)
    j = integrate(x**2*betadist, (x, 0, oo), meijerg=True, conds='separate')
    assert j[1] == (1 < beta - 1)
    assert combsimp(j[0] - i[0]**2) == (alpha + beta - 1)*alpha \
        / (beta - 2)/(beta - 1)**2

    # Beta distribution
    # NOTE: this is evaluated using antiderivatives. It also tests that
    #       meijerint_indefinite returns the simplest possible answer.
    a, b = symbols('a b', positive=True)
    betadist = x**(a - 1)*(-x + 1)**(b - 1)*gamma(a + b)/(gamma(a)*gamma(b))
    assert simplify(integrate(betadist, (x, 0, 1), meijerg=True)) == 1
    assert simplify(integrate(x*betadist, (x, 0, 1), meijerg=True)) == \
        a/(a + b)
    assert simplify(integrate(x**2*betadist, (x, 0, 1), meijerg=True)) == \
        a*(a + 1)/(a + b)/(a + b + 1)
    assert simplify(integrate(x**y*betadist, (x, 0, 1), meijerg=True)) == \
        Piecewise((gamma(a + b)*gamma(a + y)/(gamma(a)*gamma(a + b + y)),
                   -a - re(y) + 1 < 1),
                  (Integral(x**(a + y - 1)*(-x + 1)**(b - 1)*gamma(a + b)/(gamma(a)*gamma(b)),
                            (x, 0, 1)), True))

    # Chi distribution
    k = Symbol('k', integer=True, positive=True)
    chi = 2**(1 - k/2)*x**(k - 1)*exp(-x**2/2)/gamma(k/2)
    assert powsimp(integrate(chi, (x, 0, oo), meijerg=True)) == 1
    assert simplify(integrate(x*chi, (x, 0, oo), meijerg=True)) == \
        sqrt(2)*gamma((k + 1)/2)/gamma(k/2)
    assert simplify(integrate(x**2*chi, (x, 0, oo), meijerg=True)) == k

    # Chi^2 distribution
    chisquared = 2**(-k/2)/gamma(k/2)*x**(k/2 - 1)*exp(-x/2)
    assert powsimp(integrate(chisquared, (x, 0, oo), meijerg=True)) == 1
    assert simplify(integrate(x*chisquared, (x, 0, oo), meijerg=True)) == k
    assert simplify(integrate(x**2*chisquared, (x, 0, oo), meijerg=True)) == \
        k*(k + 2)
    assert combsimp(integrate(((x - k)/sqrt(2*k))**3*chisquared, (x, 0, oo),
                              meijerg=True)) == 2*sqrt(2)/sqrt(k)

    # Dagum distribution
    a, b, p = symbols('a b p', positive=True, real=True)
    # XXX (x/b)**a does not work
    dagum = a*p/x*(x/b)**(a*p)/(1 + x**a/b**a)**(p + 1)
    assert simplify(integrate(dagum, (x, 0, oo), meijerg=True)) == 1
    # XXX conditions are a mess
    arg = x*dagum
    assert simplify(integrate(arg, (x, 0, oo), meijerg=True, conds='none')
                    ) == a*b*gamma(1 - 1/a)*gamma(p + 1 + 1/a)/(
        (a*p + 1)*gamma(p))
    assert simplify(integrate(x*arg, (x, 0, oo), meijerg=True, conds='none')
                    ) == a*b**2*gamma(1 - 2/a)*gamma(p + 1 + 2/a)/(
        (a*p + 2)*gamma(p))

    # F-distribution
    d1, d2 = symbols('d1 d2', positive=True)
    f = sqrt(((d1*x)**d1 * d2**d2)/(d1*x + d2)**(d1 + d2))/x \
        / gamma(d1/2)/gamma(d2/2)*gamma((d1 + d2)/2)
    assert simplify(integrate(f, (x, 0, oo), meijerg=True)) == 1
    # TODO conditions are a mess
    assert simplify(integrate(x*f, (x, 0, oo), meijerg=True, conds='none')
                    ) == d2/(d2 - 2)
    assert simplify(integrate(x**2*f, (x, 0, oo), meijerg=True, conds='none')
                    ) == d2**2*(d1 + 2)/d1/(d2 - 4)/(d2 - 2)

    # TODO gamma, rayleigh

    # inverse gaussian
    lamda, mu = symbols('lamda mu', positive=True)
    dist = sqrt(lamda/2/pi)*x**(-Rational(3, 2))*exp(-lamda*(x - mu)**2/x/2/mu**2)

    def mysimp(expr):
        return simplify(expr.rewrite(exp))

    assert mysimp(integrate(dist, (x, 0, oo))) == 1
    assert mysimp(integrate(x*dist, (x, 0, oo))) == mu
    assert mysimp(integrate((x - mu)**2*dist, (x, 0, oo))) == mu**3/lamda
    assert mysimp(integrate((x - mu)**3*dist, (x, 0, oo))) == 3*mu**5/lamda**2

    # Levi
    c = Symbol('c', positive=True)
    assert integrate(sqrt(c/2/pi)*exp(-c/2/(x - mu))/(x - mu)**Rational(3, 2),
                     (x, mu, oo)) == 1
    # higher moments oo

    # log-logistic
    distn = (beta/alpha)*x**(beta - 1)/alpha**(beta - 1) / \
        (1 + x**beta/alpha**beta)**2
    assert simplify(integrate(distn, (x, 0, oo))) == 1
    # NOTE the conditions are a mess, but correctly state beta > 1
    assert simplify(integrate(x*distn, (x, 0, oo), conds='none')) == \
        pi*alpha/beta/sin(pi/beta)
    # (similar comment for conditions applies)
    assert simplify(integrate(x**y*distn, (x, 0, oo), conds='none')) == \
        pi*alpha**y*y/beta/sin(pi*y/beta)

    # weibull
    k = Symbol('k', positive=True, real=True)
    n = Symbol('n', positive=True)
    distn = k/lamda*(x/lamda)**(k - 1)*exp(-(x/lamda)**k)
    assert simplify(integrate(distn, (x, 0, oo))) == 1
    assert simplify(integrate(x**n*distn, (x, 0, oo))) == \
        lamda**n*gamma(1 + n/k)

    # rice distribution
    nu, sigma = symbols('nu sigma', positive=True)
    rice = x/sigma**2*exp(-(x**2 + nu**2)/2/sigma**2)*besseli(0, x*nu/sigma**2)
    assert integrate(rice, (x, 0, oo), meijerg=True) == 1
    # can someone verify higher moments?

    # Laplace distribution
    mu = Symbol('mu', extended_real=True)
    b = Symbol('b', positive=True)
    laplace = exp(-abs(x - mu)/b)/2/b
    assert integrate(laplace, (x, -oo, oo), meijerg=True) == 1
    assert integrate(x*laplace, (x, -oo, oo), meijerg=True) == mu
    assert integrate(x**2*laplace, (x, -oo, oo), meijerg=True) == \
        2*b**2 + mu**2

    # TODO are there other distributions supported on (-oo, oo) that we can do?

    # misc tests
    k = Symbol('k', positive=True)
    assert combsimp(expand_mul(integrate(log(x)*x**(k - 1)*exp(-x)/gamma(k),
                                         (x, 0, oo)))) == polygamma(0, k)


@pytest.mark.slow
def test_expint():
    """Test various exponential integrals."""
    assert simplify(integrate(exp(-z*x)/x**y,
                              (x, 1, oo),
                              meijerg=True,
                              conds='none').rewrite(expint)) == expint(y, z)

    assert integrate(exp(-z*x)/x, (x, 1, oo), meijerg=True,
                     conds='none').rewrite(expint).expand() == \
        expint(1, z)
    assert integrate(exp(-z*x)/x**2, (x, 1, oo), meijerg=True,
                     conds='none').rewrite(expint).expand() == \
        expint(2, z).rewrite(Ei).rewrite(expint)
    assert integrate(exp(-z*x)/x**3, (x, 1, oo), meijerg=True,
                     conds='none').rewrite(expint).expand() == \
        expint(3, z).rewrite(Ei).rewrite(expint).expand()

    t = Symbol('t', positive=True)
    assert integrate(-cos(x)/x, (x, t, oo), meijerg=True).expand() == Ci(t)
    assert integrate(-sin(x)/x, (x, t, oo), meijerg=True).expand() == \
        Si(t) - pi/2
    assert integrate(sin(x)/x, (x, 0, z), meijerg=True) == Si(z)
    assert integrate(sinh(x)/x, (x, 0, z), meijerg=True) == Shi(z)
    assert integrate(exp(-x)/x, x, meijerg=True).expand().rewrite(expint) == \
        I*pi - expint(1, x)
    assert integrate(exp(-x)/x**2, x, meijerg=True).rewrite(expint).expand() \
        == expint(1, x) - exp(-x)/x - I*pi

    u = Symbol('u', polar=True)
    assert integrate(cos(u)/u, u, meijerg=True).expand().as_independent(u)[1] \
        == Ci(u)
    assert integrate(cosh(u)/u, u, meijerg=True).expand().as_independent(u)[1] \
        == Chi(u)

    assert (integrate(expint(1, x), x,
                      meijerg=True).rewrite(expint).expand() ==
            x*expint(1, x) - exp(-x))
    assert (integrate(expint(2, x), x,
                      meijerg=True).rewrite(expint).expand() ==
            -x**2*expint(1, x)/2 + x*exp(-x)/2 - exp(-x)/2)
    assert (simplify(unpolarify(integrate(expint(y, x), x,
                                          meijerg=True).rewrite(expint))) ==
            -expint(y + 1, x))

    assert integrate(Si(x), x, meijerg=True) == x*Si(x) + cos(x)
    assert integrate(Ci(u), u, meijerg=True).expand() == u*Ci(u) - sin(u)
    assert integrate(Shi(x), x, meijerg=True) == x*Shi(x) - cosh(x)
    assert integrate(Chi(u), u, meijerg=True).expand() == u*Chi(u) - sinh(u)

    assert integrate(Si(x)*exp(-x), (x, 0, oo), meijerg=True) == pi/4
    assert integrate(expint(1, x)*sin(x), (x, 0, oo), meijerg=True) == log(2)/2


def test_messy():
    assert laplace_transform(Si(x), x, s) == ((-atan(s) + pi/2)/s, 0, True)

    assert laplace_transform(Shi(x), x, s) == (acoth(s)/s, 1, True)

    # where should the logs be simplified?
    assert laplace_transform(Chi(x), x, s) == \
        ((log(s**(-2)) - log((s**2 - 1)/s**2))/(2*s), 1, True)

    # TODO maybe simplify the inequalities?
    assert laplace_transform(besselj(a, x), x, s)[1:] == \
        (0, And(Integer(0) < re(a/2) + Rational(1, 2), Integer(0) < re(a/2) + 1))

    # NOTE s < 0 can be done, but argument reduction is not good enough yet
    assert fourier_transform(besselj(1, x)/x, x, s, noconds=False) == \
        (Piecewise((0, 4*abs(pi**2*s**2) > 1),
                   (2*sqrt(-4*pi**2*s**2 + 1), True)), s > 0)
    # TODO FT(besselj(0,x)) - conditions are messy (but for acceptable reasons)
    #                       - folding could be better

    assert integrate(E1(x)*besselj(0, x), (x, 0, oo), meijerg=True) == \
        log(1 + sqrt(2))
    assert integrate(E1(x)*besselj(1, x), (x, 0, oo), meijerg=True) == \
        log(Rational(1, 2) + sqrt(2)/2)

    assert integrate(1/x/sqrt(1 - x**2), x, meijerg=True) == \
        Piecewise((-acosh(1/x), 1 < abs(x**(-2))), (I*asin(1/x), True))


def test_sympyissue_6122():
    assert integrate(exp(-I*x**2), (x, -oo, oo), meijerg=True) == \
        -I*sqrt(pi)*exp(I*pi/4)


def test_sympyissue_6252():
    expr = 1/x/cbrt(a + b*x)
    anti = integrate(expr, x, meijerg=True)
    assert not expr.has(hyper)
    # XXX the expression is a mess, but actually upon differentiation and
    # putting in numerical values seems to work...
    assert not anti.has(hyper)


def test_fresnel():
    assert expand_func(integrate(sin(pi*x**2/2), x)) == fresnels(x)
    assert expand_func(integrate(cos(pi*x**2/2), x)) == fresnelc(x)


def test_sympyissue_6860():
    assert meijerint_indefinite(x**x**x, x) is None


def test_meijerint_indefinite_abs():
    # issue sympy/sympy#4311
    assert meijerint_indefinite(x*abs(9 - x**2), x) is not nan
    # issue sympy/sympy#7165
    assert meijerint_indefinite(abs(y - x**2), y) is not nan
    # issue sympy/sympy#8733
    assert meijerint_indefinite(abs(x + 1), x) is not nan


def test_sympyissue_10681():
    f = integrate(r**2*(R**2 - r**2)**0.5, r, meijerg=True)
    g = (1.0/3)*R**1.0*r**3*hyper((-0.5, Rational(3, 2)), (Rational(5, 2),),
                                  r**2*exp_polar(2*I*pi)/R**2)
    assert RR.almosteq((f/g).evalf(), 1.0, 1e-12)


def test_sympyissue_11806():
    y, L = symbols('y L', positive=True)
    assert integrate(1/sqrt(x**2 + y**2)**3, (x, -L, L)) == \
        2*L/(y**2*sqrt(L**2 + y**2))


def test_meijerint_doc():
    from diofant.integrals.meijerint_doc import __doc__
    assert __doc__[:89] == r"""Elementary functions:

.. math::
  a = a {G_{1, 1}^{1, 0}\left(\begin{matrix}  & 1 \\0 & """
