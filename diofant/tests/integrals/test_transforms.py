import pytest

from diofant import (E1, And, Ci, CosineTransform, Ei, EulerGamma,
                     FourierTransform, Function, Heaviside, I, Integer,
                     Integral, InverseCosineTransform, InverseFourierTransform,
                     InverseLaplaceTransform, InverseSineTransform,
                     LaplaceTransform, Matrix, Max, MellinTransform, Min, Ne,
                     Or, Piecewise, Rational, Si, SineTransform, Symbol, atan,
                     atan2, besseli, besselj, besselk, bessely, combsimp, cos,
                     cosh, cosine_transform, cot, erf, exp, exp_polar, expand,
                     expand_complex, expand_mul, expand_trig, expint, eye,
                     factor_terms, factorial, fourier_transform, fresnelc,
                     fresnels, gamma, hankel_transform, hyperexpand,
                     inverse_cosine_transform, inverse_fourier_transform,
                     inverse_hankel_transform, inverse_laplace_transform,
                     inverse_mellin_transform, inverse_sine_transform,
                     laplace_transform, lerchphi, log, logcombine, meijerg,
                     mellin_transform, oo, periodic_argument, pi, polar_lift,
                     powsimp, re, simplify, sin, sine_transform, sinh, sqrt,
                     symbols, tan, trigsimp, unpolarify)
from diofant.abc import a, b, beta, c, d, nu, rho, s, t, w, x
from diofant.integrals.transforms import IntegralTransformError, _simplifyconds


__all__ = ()


def test_undefined_function():
    f = Function('f')
    assert mellin_transform(f(x), x, s) == MellinTransform(f(x), x, s)
    assert mellin_transform(f(x) + exp(-x), x, s) == \
        (MellinTransform(f(x), x, s) + gamma(s), (0, oo), True)

    assert laplace_transform(2*f(x), x, s) == 2*LaplaceTransform(f(x), x, s)
    # TODO test derivative and other rules when implemented


def test_free_symbols():
    f = Function('f')
    assert mellin_transform(f(x), x, s).free_symbols == {s}
    assert mellin_transform(f(x)*a, x, s).free_symbols == {s, a}


def test_as_integral():
    f = Function('f')
    assert mellin_transform(f(x), x, s).rewrite('Integral') == \
        Integral(x**(s - 1)*f(x), (x, 0, oo))
    assert fourier_transform(f(x), x, s).rewrite('Integral') == \
        Integral(f(x)*exp(-2*I*pi*s*x), (x, -oo, oo))
    assert laplace_transform(f(x), x, s).rewrite('Integral') == \
        Integral(f(x)*exp(-s*x), (x, 0, oo))
    assert str(inverse_mellin_transform(f(s), s, x, (a, b)).rewrite('Integral')) \
        == 'Integral(x**(-s)*f(s), (s, _c - oo*I, _c + oo*I))'
    assert str(inverse_laplace_transform(f(s), s, x).rewrite('Integral')) == \
        'Integral(E**(s*x)*f(s), (s, _c - oo*I, _c + oo*I))'
    assert inverse_fourier_transform(f(s), s, x).rewrite('Integral') == \
        Integral(f(s)*exp(2*I*pi*s*x), (s, -oo, oo))
    assert sine_transform(f(x), x, s).rewrite('Integral') == \
        Integral(sqrt(2)/sqrt(pi)*f(x)*sin(x*s), (x, 0, oo))


def test_mellin_transform():
    MT = mellin_transform

    bpos = symbols('b', positive=True)

    # 8.4.2
    assert MT(x**nu*Heaviside(x - 1), x, s) == \
        (-1/(nu + s), (-oo, -re(nu)), True)
    assert MT(x**nu*Heaviside(1 - x), x, s) == \
        (1/(nu + s), (-re(nu), oo), True)

    assert MT((1 - x)**(beta - 1)*Heaviside(1 - x), x, s) == \
        (gamma(beta)*gamma(s)/gamma(beta + s), (0, oo), re(-beta) < 0)
    assert MT((x - 1)**(beta - 1)*Heaviside(x - 1), x, s) == \
        (gamma(beta)*gamma(1 - beta - s)/gamma(1 - s),
            (-oo, -re(beta) + 1), re(-beta) < 0)

    assert MT((1 + x)**(-rho), x, s) == \
        (gamma(s)*gamma(rho - s)/gamma(rho), (0, re(rho)), True)

    # TODO also the conditions should be simplified
    assert MT(abs(1 - x)**(-rho), x, s) == (
        2*sin(pi*rho/2)*gamma(1 - rho)*cos(pi*(rho/2 - s))*gamma(s)*gamma(rho-s)/pi,
        (0, re(rho)), And(re(rho) - 1 < 0, re(rho) < 1))
    mt = MT((1 - x)**(beta - 1)*Heaviside(1 - x)
            + a*(x - 1)**(beta - 1)*Heaviside(x - 1), x, s)
    assert mt[1], mt[2] == ((0, -re(beta) + 1), True)

    assert MT((x**a - b**a)/(x - b), x, s)[0] == \
        pi*b**(a + s - 1)*sin(pi*a)/(sin(pi*s)*sin(pi*(a + s)))
    assert MT((x**a - bpos**a)/(x - bpos), x, s) == \
        (pi*bpos**(a + s - 1)*sin(pi*a)/(sin(pi*s)*sin(pi*(a + s))),
            (Max(-re(a), 0), Min(1 - re(a), 1)), True)

    expr = (sqrt(x + b**2) + b)**a
    assert MT(expr.subs({b: bpos}), x, s) == \
        (-a*(2*bpos)**(a + 2*s)*gamma(s)*gamma(-a - 2*s)/gamma(-a - s + 1),
         (0, -re(a)/2), True)
    expr = (sqrt(x + b**2) + b)**a/sqrt(x + b**2)
    assert MT(expr.subs({b: bpos}), x, s) == \
        (2**(a + 2*s)*bpos**(a + 2*s - 1)*gamma(s)
                                         * gamma(1 - a - 2*s)/gamma(1 - a - s),
            (0, -re(a)/2 + Rational(1, 2)), True)

    # 8.4.2
    assert MT(exp(-x), x, s) == (gamma(s), (0, oo), True)
    assert MT(exp(-1/x), x, s) == (gamma(-s), (-oo, 0), True)

    # 8.4.5
    assert MT(log(x)**4*Heaviside(1 - x), x, s) == (24/s**5, (0, oo), True)
    assert MT(log(x)**3*Heaviside(x - 1), x, s) == (6/s**4, (-oo, 0), True)
    assert MT(log(x + 1), x, s) == (pi/(s*sin(pi*s)), (-1, 0), True)
    assert MT(log(1/x + 1), x, s) == (pi/(s*sin(pi*s)), (0, 1), True)
    assert MT(log(abs(1 - x)), x, s) == (pi/(s*tan(pi*s)), (-1, 0), True)
    assert MT(log(abs(1 - 1/x)), x, s) == (pi/(s*tan(pi*s)), (0, 1), True)

    # TODO we cannot currently do these (needs summation of 3F2(-1))
    #      this also implies that they cannot be written as a single g-function
    #      (although this is possible)
    mt = MT(log(x)/(x + 1), x, s)
    assert mt[1:] == ((0, 1), True)
    assert not hyperexpand(mt[0], allow_hyper=True).has(meijerg)
    mt = MT(log(x)**2/(x + 1), x, s)
    assert mt[1:] == ((0, 1), True)
    assert not hyperexpand(mt[0], allow_hyper=True).has(meijerg)
    mt = MT(log(x)/(x + 1)**2, x, s)
    assert mt[1:] == ((0, 2), True)
    assert not hyperexpand(mt[0], allow_hyper=True).has(meijerg)

    # 8.4.14
    assert MT(erf(sqrt(x)), x, s) == \
        (-gamma(s + Rational(1, 2))/(sqrt(pi)*s), (-Rational(1, 2), 0), True)


@pytest.mark.slow
def test_mellin_transform_bessel():
    a, b = symbols('a, b', complex=True)
    MT = mellin_transform

    # 8.4.19
    assert MT(besselj(a, 2*sqrt(x)), x, s) == \
        (gamma(a/2 + s)/gamma(a/2 - s + 1), (-re(a)/2, Rational(3, 4)), True)
    assert MT(sin(sqrt(x))*besselj(a, sqrt(x)), x, s) == \
        (2**a*gamma(-2*s + Rational(1, 2))*gamma(a/2 + s + Rational(1, 2))/(
            gamma(-a/2 - s + 1)*gamma(a - 2*s + 1)), (
            -re(a)/2 - Rational(1, 2), Rational(1, 4)), True)
    assert MT(cos(sqrt(x))*besselj(a, sqrt(x)), x, s) == \
        (2**a*gamma(a/2 + s)*gamma(-2*s + Rational(1, 2))/(
            gamma(-a/2 - s + Rational(1, 2))*gamma(a - 2*s + 1)), (
            -re(a)/2, Rational(1, 4)), True)
    assert MT(besselj(a, sqrt(x))**2, x, s) == \
        (gamma(a + s)*gamma(Rational(1, 2) - s)
         / (sqrt(pi)*gamma(1 - s)*gamma(1 + a - s)),
            (-re(a), Rational(1, 2)), True)
    assert MT(besselj(a, sqrt(x))*besselj(-a, sqrt(x)), x, s) == \
        (gamma(s)*gamma(Rational(1, 2) - s)
         / (sqrt(pi)*gamma(1 - a - s)*gamma(1 + a - s)),
            (0, Rational(1, 2)), True)
    # NOTE: prudnikov gives the strip below as (1/2 - re(a), 1). As far as
    #       I can see this is wrong (since besselj(z) ~ 1/sqrt(z) for z large)
    assert MT(besselj(a - 1, sqrt(x))*besselj(a, sqrt(x)), x, s) == \
        (gamma(1 - s)*gamma(a + s - Rational(1, 2))
         / (sqrt(pi)*gamma(Rational(3, 2) - s)*gamma(a - s + Rational(1, 2))),
            (Rational(1, 2) - re(a), Rational(1, 2)), True)
    assert MT(besselj(a, sqrt(x))*besselj(b, sqrt(x)), x, s) == \
        (4**s*gamma(1 - 2*s)*gamma((a + b)/2 + s)
         / (gamma(1 - s + (b - a)/2)*gamma(1 - s + (a - b)/2)
            * gamma(1 - s + (a + b)/2)),
            (-(re(a) + re(b))/2, Rational(1, 2)), True)
    assert MT(besselj(a, sqrt(x))**2 + besselj(-a, sqrt(x))**2, x, s)[1:] == \
        ((Max(re(a), -re(a)), Rational(1, 2)), True)

    # Section 8.4.20
    assert MT(bessely(a, 2*sqrt(x)), x, s) == \
        (-cos(pi*(a/2 - s))*gamma(s - a/2)*gamma(s + a/2)/pi,
            (Max(-re(a)/2, re(a)/2), Rational(3, 4)), True)
    assert MT(sin(sqrt(x))*bessely(a, sqrt(x)), x, s) == \
        (-4**s*sin(pi*(a/2 - s))*gamma(Rational(1, 2) - 2*s)
         * gamma((1 - a)/2 + s)*gamma((1 + a)/2 + s)
         / (sqrt(pi)*gamma(1 - s - a/2)*gamma(1 - s + a/2)),
            (Max(-(re(a) + 1)/2, (re(a) - 1)/2), Rational(1, 4)), True)
    assert MT(cos(sqrt(x))*bessely(a, sqrt(x)), x, s) == \
        (-4**s*cos(pi*(a/2 - s))*gamma(s - a/2)*gamma(s + a/2)*gamma(Rational(1, 2) - 2*s)
         / (sqrt(pi)*gamma(Rational(1, 2) - s - a/2)*gamma(Rational(1, 2) - s + a/2)),
            (Max(-re(a)/2, re(a)/2), Rational(1, 4)), True)
    assert MT(besselj(a, sqrt(x))*bessely(a, sqrt(x)), x, s) == \
        (-cos(pi*s)*gamma(s)*gamma(a + s)*gamma(Rational(1, 2) - s)
         / (pi**Rational(3, 2)*gamma(1 + a - s)),
            (Max(-re(a), 0), Rational(1, 2)), True)
    assert MT(besselj(a, sqrt(x))*bessely(b, sqrt(x)), x, s) == \
        (-4**s*cos(pi*(a/2 - b/2 + s))*gamma(1 - 2*s)
         * gamma(a/2 - b/2 + s)*gamma(a/2 + b/2 + s)
         / (pi*gamma(a/2 - b/2 - s + 1)*gamma(a/2 + b/2 - s + 1)),
            (Max((-re(a) + re(b))/2, (-re(a) - re(b))/2), Rational(1, 2)), True)
    # NOTE bessely(a, sqrt(x))**2 and bessely(a, sqrt(x))*bessely(b, sqrt(x))
    # are a mess (no matter what way you look at it ...)
    assert MT(bessely(a, sqrt(x))**2, x, s)[1:] == \
             ((Max(-re(a), 0, re(a)), Rational(1, 2)), True)

    # Section 8.4.22
    # TODO we can't do any of these (delicate cancellation)

    # Section 8.4.23
    assert MT(besselk(a, 2*sqrt(x)), x, s) == \
        (gamma(
         s - a/2)*gamma(s + a/2)/2, (Max(-re(a)/2, re(a)/2), oo), True)
    assert MT(besselj(a, 2*sqrt(2*sqrt(x)))*besselk(
        a, 2*sqrt(2*sqrt(x))), x, s) == (4**(-s)*gamma(2*s) *
                                         gamma(a/2 + s)/(2*gamma(a/2 - s + 1)), (Max(0, -re(a)/2), oo), True)
    # TODO bessely(a, x)*besselk(a, x) is a mess
    assert MT(besseli(a, sqrt(x))*besselk(a, sqrt(x)), x, s) == \
        (gamma(s)*gamma(
            a + s)*gamma(-s + Rational(1, 2))/(2*sqrt(pi)*gamma(a - s + 1)),
         (Max(-re(a), 0), Rational(1, 2)), True)
    assert MT(besseli(b, sqrt(x))*besselk(a, sqrt(x)), x, s) == \
        (2**(2*s - 1)*gamma(-2*s + 1)*gamma(-a/2 + b/2 + s) *
         gamma(a/2 + b/2 + s)/(gamma(-a/2 + b/2 - s + 1) *
                               gamma(a/2 + b/2 - s + 1)), (Max(-re(a)/2 - re(b)/2,
                                                               re(a)/2 - re(b)/2), Rational(1, 2)), True)

    # TODO products of besselk are a mess

    mt = MT(exp(-x/2)*besselk(a, x/2), x, s)
    mt0 = combsimp(trigsimp(combsimp(mt[0].expand(func=True))))
    assert mt0 == 2*pi**Rational(3, 2)*cos(pi*s)*gamma(-s + Rational(1, 2))/(
        (cos(2*pi*a) - cos(2*pi*s))*gamma(-a - s + 1)*gamma(a - s + 1))
    assert mt[1:] == ((Max(-re(a), re(a)), oo), True)
    # TODO exp(x/2)*besselk(a, x/2) [etc] cannot currently be done
    # TODO various strange products of special orders


def test_expint():
    aneg = Symbol('a', negative=True)
    u = Symbol('u', polar=True)

    assert mellin_transform(E1(x), x, s) == (gamma(s)/s, (0, oo), True)
    assert inverse_mellin_transform(gamma(s)/s, s, x,
                                    (0, oo)).rewrite(expint).expand() == E1(x)
    assert mellin_transform(expint(a, x), x, s) == \
        (gamma(s)/(a + s - 1), (Max(1 - re(a), 0), oo), True)
    # XXX IMT has hickups with complicated strips ...
    assert simplify(unpolarify(
                    inverse_mellin_transform(gamma(s)/(aneg + s - 1), s, x,
                                             (1 - aneg, oo)).rewrite(expint).expand(func=True))) == \
        expint(aneg, x)

    assert mellin_transform(Si(x), x, s) == \
        (-2**s*sqrt(pi)*gamma(s/2 + Rational(1, 2))/(
            2*s*gamma(-s/2 + 1)), (-1, 0), True)
    assert inverse_mellin_transform(-2**s*sqrt(pi)*gamma((s + 1)/2)
                                    / (2*s*gamma(-s/2 + 1)), s, x, (-1, 0)) \
        == Si(x)

    assert mellin_transform(Ci(sqrt(x)), x, s) == \
        (-2**(2*s - 1)*sqrt(pi)*gamma(s)/(s*gamma(-s + Rational(1, 2))), (0, 1), True)
    assert inverse_mellin_transform(
        -4**s*sqrt(pi)*gamma(s)/(2*s*gamma(-s + Rational(1, 2))),
        s, u, (0, 1)).expand() == Ci(sqrt(u))

    # TODO LT of Si, Shi, Chi is a mess ...
    assert laplace_transform(Ci(x), x, s) == (-log(1 + s**2)/2/s, 0, True)
    assert laplace_transform(expint(a, x), x, s) == \
        (lerchphi(s*polar_lift(-1), 1, a), 0, Integer(0) < re(a))
    assert laplace_transform(expint(1, x), x, s) == (log(s + 1)/s, 0, True)
    assert laplace_transform(expint(2, x), x, s) == \
        ((s - log(s + 1))/s**2, 0, True)

    assert inverse_laplace_transform(-log(1 + s**2)/2/s, s, u).expand() == \
        Heaviside(u)*Ci(u)
    assert inverse_laplace_transform(log(s + 1)/s, s, x).rewrite(expint) == \
        Heaviside(x)*E1(x)
    assert inverse_laplace_transform((s - log(s + 1))/s**2, s,
                                     x).rewrite(expint).expand() == \
        (expint(2, x)*Heaviside(x)).rewrite(Ei).rewrite(expint).expand()


def test_inverse_mellin_transform():
    IMT = inverse_mellin_transform

    # test passing "None"
    assert IMT(1/(s**2 - 1), s, x, (-1, None)) == \
        -x*Heaviside(-x + 1)/2 - Heaviside(x - 1)/(2*x)
    assert IMT(1/(s**2 - 1), s, x, (None, 1)) == \
        (-x/2 + 1/(2*x))*Heaviside(-x + 1)

    def simp_pows(expr):
        return simplify(powsimp(expand_mul(expr, deep=False), force=True)).replace(exp_polar, exp)

    assert simp_pows(IMT(d**c*d**(s - 1)*sin(pi*c)
                         * gamma(s)*gamma(s + c)*gamma(1 - s)*gamma(1 - s - c)/pi,
                         s, x, (Max(-re(c), 0), Min(1 - re(c), 1)))) \
        == (x**c - d**c)/(x - d)

    # issue sympy/sympy#8882

    # This is the original test.
    # r = Symbol('r')
    # psi = 1/r*sin(r)*exp(-(a0*r))
    # h = -1/2*diff(psi, r, r) - 1/r*psi
    # f = 4*pi*psi*h*r**2
    # assert integrate(f, (r, -oo, 3), meijerg=True).has(Integral) == True

    # To save time, only the critical part is included.
    F = -a**(-s + 1)*(4 + 1/a**2)**(-s/2)*sqrt(1/a**2)*exp(-s*I*pi) * \
        sin(s*atan(sqrt(1/a**2)/2))*gamma(s)
    pytest.raises(IntegralTransformError, lambda:
                  inverse_mellin_transform(F, s, x, (-1, oo),
                                           **{'as_meijerg': True, 'needeval': True}))


@pytest.mark.slow
def test_inverse_mellin_transform2():
    IMT = inverse_mellin_transform

    assert IMT(gamma(s), s, x, (0, oo)) == exp(-x)
    assert IMT(gamma(-s), s, x, (-oo, 0)) == exp(-1/x)
    assert simplify(IMT(s/(2*s**2 - 2), s, x, (2, oo))) == \
        (x**2 + 1)*Heaviside(1 - x)/(4*x)

    # test expansion of sums
    assert IMT(gamma(s) + gamma(s - 1), s, x, (1, oo)) == (x + 1)*exp(-x)/x

    # test factorisation of polys
    r = symbols('r', extended_real=True)
    assert (IMT(1/(s**2 + 1), s, exp(-x),
            (None, oo)).subs({x: r}).rewrite(sin).simplify() ==
            sin(r)*Heaviside(1 - exp(-r)))

    # test multiplicative substitution
    _a, _b = symbols('a b', positive=True)
    assert IMT(_b**(-s/_a)*factorial(s/_a)/s, s, x, (0, oo)) == exp(-_b*x**_a)
    assert IMT(factorial(_a/_b + s/_b)/(_a + s), s, x, (-_a, oo)) == x**_a*exp(-x**_b)

    def simp_pows(expr):
        return simplify(powsimp(expand_mul(expr, deep=False), force=True)).replace(exp_polar, exp)

    # Now test the inverses of all direct transforms tested above

    # Section 8.4.2
    nu = symbols('nu', real=True)
    assert IMT(-1/(nu + s), s, x, (-oo, None)) == x**nu*Heaviside(x - 1)
    assert IMT(1/(nu + s), s, x, (None, oo)) == x**nu*Heaviside(1 - x)
    assert simp_pows(IMT(gamma(beta)*gamma(s)/gamma(s + beta), s, x, (0, oo))) \
        == (1 - x)**(beta - 1)*Heaviside(1 - x)
    assert simp_pows(IMT(gamma(beta)*gamma(1 - beta - s)/gamma(1 - s),
                         s, x, (-oo, None))) \
        == (x - 1)**(beta - 1)*Heaviside(x - 1)
    assert simp_pows(IMT(gamma(s)*gamma(rho - s)/gamma(rho), s, x, (0, None))) \
        == (1/(x + 1))**rho

    assert simplify(IMT(1/sqrt(pi)*(-c/2)*gamma(s)*gamma((1 - c)/2 - s)
                        * gamma(-c/2 - s)/gamma(1 - c - s),
                        s, x, (0, -re(c)/2))) == \
        (1 + sqrt(x + 1))**c
    assert simplify(IMT(2**(a + 2*s)*b**(a + 2*s - 1)*gamma(s)*gamma(1 - a - 2*s)
                        / gamma(1 - a - s), s, x, (0, (-re(a) + 1)/2))) == \
        b**(a - 1)*(sqrt(1 + x/b**2) + 1)**(a - 1)*(b**2*sqrt(1 + x/b**2) +
                                                    b**2 + x)/(b**2 + x)
    assert simplify(IMT(-2**(c + 2*s)*c*b**(c + 2*s)*gamma(s)*gamma(-c - 2*s)
                        / gamma(-c - s + 1), s, x, (0, -re(c)/2))) == \
        b**c*(sqrt(1 + x/b**2) + 1)**c

    # Section 8.4.5
    assert IMT(24/s**5, s, x, (0, oo)) == log(x)**4*Heaviside(1 - x)
    assert expand(IMT(6/s**4, s, x, (-oo, 0)), force=True) == \
        log(x)**3*Heaviside(x - 1)
    assert IMT(pi/(s*sin(pi*s)), s, x, (-1, 0)) == log(x + 1)
    assert IMT(pi/(s*sin(pi*s/2)), s, x, (-2, 0)) == log(x**2 + 1)
    assert IMT(pi/(s*sin(2*pi*s)), s, x, (-Rational(1, 2), 0)) == log(sqrt(x) + 1)
    assert IMT(pi/(s*sin(pi*s)), s, x, (0, 1)) == log(1 + 1/x)

    # TODO
    def mysimp(expr):
        return expand(
            powsimp(logcombine(expr, force=True), force=True, deep=True),
            force=True).replace(exp_polar, exp)

    assert mysimp(mysimp(IMT(pi/(s*tan(pi*s)), s, x, (-1, 0)))) in [
        log(1 - x)*Heaviside(1 - x) + log(x - 1)*Heaviside(x - 1),
        log(x)*Heaviside(x - 1) + log(1 - 1/x)*Heaviside(x - 1) + log(-x +
                                                                      1)*Heaviside(-x + 1)]
    # test passing cot
    assert mysimp(IMT(pi*cot(pi*s)/s, s, x, (0, 1))) in [
        log(1/x - 1)*Heaviside(1 - x) + log(1 - 1/x)*Heaviside(x - 1),
        -log(x)*Heaviside(-x + 1) + log(1 - 1/x)*Heaviside(x - 1) + log(-x +
                                                                        1)*Heaviside(-x + 1), ]

    # 8.4.14
    assert IMT(-gamma(s + Rational(1, 2))/(sqrt(pi)*s), s, x, (-Rational(1, 2), 0)) == \
        erf(sqrt(x))

    # 8.4.19
    assert simplify(IMT(gamma(a/2 + s)/gamma(a/2 - s + 1), s, x, (-re(a)/2, Rational(3, 4)))) \
        == besselj(a, 2*sqrt(x))
    assert simplify(IMT(2**a*gamma(Rational(1, 2) - 2*s)*gamma(s + (a + 1)/2)
                        / (gamma(1 - s - a/2)*gamma(1 - 2*s + a)),
                        s, x, (-(re(a) + 1)/2, Rational(1, 4)))) == \
        sin(sqrt(x))*besselj(a, sqrt(x))
    assert simplify(IMT(2**a*gamma(a/2 + s)*gamma(Rational(1, 2) - 2*s)
                        / (gamma(Rational(1, 2) - s - a/2)*gamma(1 - 2*s + a)),
                        s, x, (-re(a)/2, Rational(1, 4)))) == \
        cos(sqrt(x))*besselj(a, sqrt(x))
    # TODO this comes out as an amazing mess, but simplifies nicely
    assert simplify(IMT(gamma(a + s)*gamma(Rational(1, 2) - s)
                        / (sqrt(pi)*gamma(1 - s)*gamma(1 + a - s)),
                        s, x, (-re(a), Rational(1, 2)))) == \
        besselj(a, sqrt(x))**2
    assert simplify(IMT(gamma(s)*gamma(Rational(1, 2) - s)
                        / (sqrt(pi)*gamma(1 - s - a)*gamma(1 + a - s)),
                        s, x, (0, Rational(1, 2)))) == \
        besselj(-a, sqrt(x))*besselj(a, sqrt(x))
    assert simplify(IMT(4**s*gamma(-2*s + 1)*gamma(a/2 + b/2 + s)
                        / (gamma(-a/2 + b/2 - s + 1)*gamma(a/2 - b/2 - s + 1)
                           * gamma(a/2 + b/2 - s + 1)),
                        s, x, (-(re(a) + re(b))/2, Rational(1, 2)))) == \
        besselj(a, sqrt(x))*besselj(b, sqrt(x))

    # Section 8.4.20
    # TODO this can be further simplified!
    assert simplify(IMT(-2**(2*s)*cos(pi*a/2 - pi*b/2 + pi*s)*gamma(-2*s + 1) *
                        gamma(a/2 - b/2 + s)*gamma(a/2 + b/2 + s) /
                        (pi*gamma(a/2 - b/2 - s + 1)*gamma(a/2 + b/2 - s + 1)),
                        s, x,
                        (Max(-re(a)/2 - re(b)/2, -re(a)/2 + re(b)/2), Rational(1, 2)))) == \
        besselj(a, sqrt(x))*-(besselj(-b, sqrt(x)) -
                              besselj(b, sqrt(x))*cos(pi*b))/sin(pi*b)
    # TODO more

    # for coverage

    assert IMT(pi/cos(pi*s), s, x, (0, Rational(1, 2))) == sqrt(x)/(x + 1)


@pytest.mark.slow
def test_laplace_transform():
    LT = laplace_transform
    a, b, c, = symbols('a b c', positive=True)
    f = Function('f')

    # Test unevaluated form
    assert laplace_transform(f(t), t, w) == LaplaceTransform(f(t), t, w)
    assert inverse_laplace_transform(
        f(w), w, t, plane=0) == InverseLaplaceTransform(f(w), w, t, 0)

    # test a bug
    spos = symbols('s', positive=True)
    assert LT(exp(t), t, spos)[:2] == (1/(spos - 1), 1)

    # basic tests from wikipedia

    assert LT((t - a)**b*exp(-c*(t - a))*Heaviside(t - a), t, s) == \
        ((s + c)**(-b - 1)*exp(-a*s)*gamma(b + 1), -c, True)
    assert LT(t**a, t, s) == (s**(-a - 1)*gamma(a + 1), 0, True)
    assert LT(Heaviside(t), t, s) == (1/s, 0, True)
    assert LT(Heaviside(t - a), t, s) == (exp(-a*s)/s, 0, True)
    assert LT(1 - exp(-a*t), t, s) == (a/(s*(a + s)), 0, True)

    assert LT((exp(2*t) - 1)*exp(-b - t)*Heaviside(t)/2, t, s, noconds=True) \
        == exp(-b)/(s**2 - 1)

    assert LT(exp(t), t, s)[:2] == (1/(s - 1), 1)
    assert LT(exp(2*t), t, s)[:2] == (1/(s - 2), 2)
    assert LT(exp(a*t), t, s)[:2] == (1/(s - a), a)

    assert LT(log(t/a), t, s) == ((log(a*s) + EulerGamma)/s/-1, 0, True)

    assert LT(erf(t), t, s) == ((-erf(s/2) + 1)*exp(s**2/4)/s, 0, True)

    assert LT(sin(a*t), t, s) == (a/(a**2 + s**2), 0, True)
    assert LT(cos(a*t), t, s) == (s/(a**2 + s**2), 0, True)
    # TODO would be nice to have these come out better
    assert LT(
        exp(-a*t)*sin(b*t), t, s) == (b/(b**2 + (a + s)**2), -a, True)
    assert LT(exp(-a*t)*cos(b*t), t, s) == \
        ((a + s)/(b**2 + (a + s)**2), -a, True)

    assert LT(besselj(0, t), t, s) == (1/sqrt(1 + s**2), 0, True)
    assert LT(besselj(1, t), t, s) == (1 - 1/sqrt(1 + 1/s**2), 0, True)
    # TODO general order works, but is a *mess*
    # TODO besseli also works, but is an even greater mess

    # test a bug in conditions processing
    # TODO the auxiliary condition should be recognised/simplified
    assert LT(exp(t)*cos(t), t, s)[:-1] in [
        ((s - 1)/(s**2 - 2*s + 2), -oo),
        ((s - 1)/((s - 1)**2 + 1), -oo),
    ]

    # Fresnel functions
    assert laplace_transform(fresnels(t), t, s) == \
        ((-sin(s**2/(2*pi))*fresnels(s/pi) + sin(s**2/(2*pi))/2 -
            cos(s**2/(2*pi))*fresnelc(s/pi) + cos(s**2/(2*pi))/2)/s, 0, True)
    assert laplace_transform(fresnelc(t), t, s) == (
        (sin(s**2/(2*pi))*fresnelc(s/pi)/s - cos(s**2/(2*pi))*fresnels(s/pi)/s
         + sqrt(2)*cos(s**2/(2*pi) + pi/4)/(2*s), 0, True))

    assert LT(Matrix([[exp(t), t*exp(-t)], [t*exp(-t), exp(t)]]), t, s) ==\
        Matrix([
            [(1/(s - 1), 1, True), ((s + 1)**(-2), 0, True)],
            [((s + 1)**(-2), 0, True), (1/(s - 1), 1, True)]
        ])


def test_sympyissue_7173():
    LT = laplace_transform

    # hyperbolic
    assert LT(sinh(x + 3), x, s) == ((s*sinh(3) + cosh(3))/(s**2 - 1), 1, True)
    # trig (make sure they are not being rewritten in terms of exp)
    assert LT(cos(x + 3), x, s) == ((s*cos(3) - sin(3))/(s**2 + 1), 0, True)

    t = symbols('t', real=True, positive=True)
    w = symbols('w', real=True)

    assert LT(sinh(w*t)*cosh(w*t), t, s)[0] == w/(s**2 - 4*w**2)


def test_inverse_laplace_transform():
    ILT = inverse_laplace_transform
    a, b = symbols('a b', positive=True, finite=True)

    def simp_hyp(expr):
        return factor_terms(expand_mul(expr)).rewrite(sin)

    # just test inverses of all of the above
    assert ILT(1/s, s, t) == Heaviside(t)
    assert ILT(1/s**2, s, t) == t*Heaviside(t)
    assert ILT(1/s**5, s, t) == t**4*Heaviside(t)/24
    assert ILT(exp(-a*s)/s, s, t) == Heaviside(t - a)
    assert ILT(exp(-a*s)/(s + b), s, t) == exp(b*(a - t))*Heaviside(-a + t)
    assert ILT(a/(s**2 + a**2), s, t) == sin(a*t)*Heaviside(t)
    assert ILT(s/(s**2 + a**2), s, t) == cos(a*t)*Heaviside(t)
    # TODO is there a way around simp_hyp?
    assert simp_hyp(ILT(a/(s**2 - a**2), s, t)) == sinh(a*t)*Heaviside(t)
    assert simp_hyp(ILT(s/(s**2 - a**2), s, t)) == cosh(a*t)*Heaviside(t)
    assert ILT(a/((s + b)**2 + a**2), s, t) == exp(-b*t)*sin(a*t)*Heaviside(t)
    assert ILT(
        (s + b)/((s + b)**2 + a**2), s, t) == exp(-b*t)*cos(a*t)*Heaviside(t)
    # TODO sinh/cosh shifted come out a mess. also delayed trig is a mess
    # TODO should this simplify further?
    assert ILT(exp(-a*s)/s**b, s, t) == \
        (t - a)**(b - 1)*Heaviside(t - a)/gamma(b)

    assert ILT(exp(-a*s)/sqrt(1 + s**2), s, t) == \
        Heaviside(t - a)*besselj(0, a - t)  # note: besselj(0, x) is even

    # XXX ILT turns these branch factor into trig functions ...
    assert simplify(ILT(a**b*(s + sqrt(s**2 - a**2))**(-b)/sqrt(s**2 - a**2),
                        s, t).rewrite(exp)) == \
        Heaviside(t)*besseli(b, a*t)
    assert ILT(a**b*(s + sqrt(s**2 + a**2))**(-b)/sqrt(s**2 + a**2),
               s, t).rewrite(exp) == \
        Heaviside(t)*besselj(b, a*t)

    assert ILT(1/(s*sqrt(s + 1)), s, t) == Heaviside(t)*erf(sqrt(t))
    # TODO can we make erf(t) work?

    assert ILT(1/(s**2*(s**2 + 1)), s, t) == (t - sin(t))*Heaviside(t)

    assert ILT((s * eye(2) - Matrix([[1, 0], [0, 2]])).inv(), s, t) ==\
        Matrix([[exp(t)*Heaviside(t), 0], [0, exp(2*t)*Heaviside(t)]])


def test_fourier_transform():
    FT = fourier_transform
    IFT = inverse_fourier_transform

    def simp(x):
        return simplify(expand_trig(expand_complex(expand(x))))

    def sinc(x):
        return sin(pi*x)/(pi*x)
    k = symbols('k', extended_real=True)
    f = Function('f')

    # TODO for this to work with real a, need to expand abs(a*x) to abs(a)*abs(x)
    a = symbols('a', positive=True)
    b = symbols('b', positive=True)

    posk = symbols('posk', positive=True)

    # Test unevaluated form
    assert fourier_transform(f(x), x, k) == FourierTransform(f(x), x, k)
    assert inverse_fourier_transform(
        f(k), k, x) == InverseFourierTransform(f(k), k, x)

    # basic examples from wikipedia
    assert simp(FT(Heaviside(1 - abs(2*a*x)), x, k)) == sinc(k/a)/a
    # TODO IFT is a *mess*
    assert simp(FT(Heaviside(1 - abs(a*x))*(1 - abs(a*x)), x, k)) == sinc(k/a)**2/a
    # TODO IFT

    assert FT(exp(-a*x)*Heaviside(x), x, k) == 1/(a + 2*pi*I*k)
    # NOTE: the ift comes out in pieces
    assert IFT(1/(a + 2*pi*I*x), x, posk,
               noconds=False) == (exp(-a*posk), True)
    assert IFT(1/(a + 2*pi*I*x), x, -posk,
               noconds=False) == (0, True)
    assert IFT(1/(a + 2*pi*I*x), x, symbols('k', negative=True),
               noconds=False) == (0, True)
    # TODO IFT without factoring comes out as meijer g

    assert FT(x*exp(-a*x)*Heaviside(x), x, k) == 1/(a + 2*pi*I*k)**2
    assert FT(exp(-a*x)*sin(b*x)*Heaviside(x), x, k) == \
        b/(b**2 + (a + 2*I*pi*k)**2)

    assert FT(exp(-a*x**2), x, k) == sqrt(pi)*exp(-pi**2*k**2/a)/sqrt(a)
    assert IFT(sqrt(pi/a)*exp(-(pi*k)**2/a), k, x) == exp(-a*x**2)
    assert FT(exp(-a*abs(x)), x, k) == 2*a/(a**2 + 4*pi**2*k**2)
    # TODO IFT (comes out as meijer G)

    # TODO besselj(n, x), n an integer > 0 actually can be done...

    # TODO are there other common transforms (no distributions!)?


def test_sine_transform():
    f = Function('f')

    # Test unevaluated form
    assert sine_transform(f(t), t, w) == SineTransform(f(t), t, w)
    assert inverse_sine_transform(
        f(w), w, t) == InverseSineTransform(f(w), w, t)

    assert sine_transform(1/sqrt(t), t, w) == 1/sqrt(w)
    assert inverse_sine_transform(1/sqrt(w), w, t) == 1/sqrt(t)

    assert sine_transform(
        (1/sqrt(t))**3, t, w) == sqrt(w)*gamma(Rational(1, 4))/(2*gamma(Rational(5, 4)))

    assert sine_transform(t**(-a), t, w) == 2**(
        -a + Rational(1, 2))*w**(a - 1)*gamma(-a/2 + 1)/gamma((a + 1)/2)
    assert inverse_sine_transform(2**(-a +
                                      Rational(1, 2))*w**(a - 1)*gamma(-a/2 + 1)/gamma(a/2 + Rational(1, 2)), w, t) == t**(-a)

    assert sine_transform(
        exp(-a*t), t, w) == sqrt(2)*w/(sqrt(pi)*(a**2 + w**2))
    assert inverse_sine_transform(
        sqrt(2)*w/(sqrt(pi)*(a**2 + w**2)), w, t) == exp(-a*t)

    assert sine_transform(
        log(t)/t, t, w) == -sqrt(2)*sqrt(pi)*(log(w**2) + 2*EulerGamma)/4

    assert sine_transform(
        t*exp(-a*t**2), t, w) == sqrt(2)*w*exp(-w**2/(4*a))/(4*a**Rational(3, 2))
    assert inverse_sine_transform(
        sqrt(2)*w*exp(-w**2/(4*a))/(4*a**Rational(3, 2)), w, t) == t*exp(-a*t**2)


def test_cosine_transform():
    f = Function('f')

    # Test unevaluated form
    assert cosine_transform(f(t), t, w) == CosineTransform(f(t), t, w)
    assert inverse_cosine_transform(
        f(w), w, t) == InverseCosineTransform(f(w), w, t)

    assert cosine_transform(1/sqrt(t), t, w) == 1/sqrt(w)
    assert inverse_cosine_transform(1/sqrt(w), w, t) == 1/sqrt(t)

    assert cosine_transform(1/(
        a**2 + t**2), t, w) == sqrt(2)*sqrt(pi)*exp(-a*w)/(2*a)

    assert cosine_transform(t**(
        -a), t, w) == 2**(-a + Rational(1, 2))*w**(a - 1)*gamma((-a + 1)/2)/gamma(a/2)
    assert inverse_cosine_transform(2**(-a +
                                        Rational(1, 2))*w**(a - 1)*gamma(-a/2 + Rational(1, 2))/gamma(a/2), w, t) == t**(-a)

    assert cosine_transform(
        exp(-a*t), t, w) == sqrt(2)*a/(sqrt(pi)*(a**2 + w**2))
    assert inverse_cosine_transform(
        sqrt(2)*a/(sqrt(pi)*(a**2 + w**2)), w, t) == exp(-a*t)

    assert cosine_transform(exp(-a*sqrt(t))*cos(a*sqrt(
        t)), t, w) == a*exp(-a**2/(2*w))/(2*w**Rational(3, 2))

    assert cosine_transform(1/(a + t), t, w) == sqrt(2)*(
        (-2*Si(a*w) + pi)*sin(a*w)/2 - cos(a*w)*Ci(a*w))/sqrt(pi)
    assert inverse_cosine_transform(sqrt(2)*meijerg(((Rational(1, 2), 0), ()), (
        (Rational(1, 2), 0, 0), (Rational(1, 2),)), a**2*w**2/4)/(2*pi), w, t) == 1/(a + t)

    assert cosine_transform(1/sqrt(a**2 + t**2), t, w) == sqrt(2)*meijerg(
        ((Rational(1, 2),), ()), ((0, 0), (Rational(1, 2),)), a**2*w**2/4)/(2*sqrt(pi))
    assert inverse_cosine_transform(sqrt(2)*meijerg(((Rational(1, 2),), ()), ((0, 0), (Rational(1, 2),)), a**2*w**2/4)/(2*sqrt(pi)), w, t) == 1/(a*sqrt(1 + t**2/a**2))


def test_hankel_transform():
    r = Symbol('r')
    k = Symbol('k')
    nu = Symbol('nu')
    m = Symbol('m')

    assert hankel_transform(1/r, r, k, 0) == 1/k
    assert inverse_hankel_transform(1/k, k, r, 0) == 1/r

    assert hankel_transform(
        1/r**m, r, k, 0) == 2**(-m + 1)*k**(m - 2)*gamma(-m/2 + 1)/gamma(m/2)
    assert inverse_hankel_transform(
        2**(-m + 1)*k**(m - 2)*gamma(-m/2 + 1)/gamma(m/2), k, r, 0) == r**(-m)

    assert hankel_transform(1/r**m, r, k, nu) == (
        2*2**(-m)*k**(m - 2)*gamma(-m/2 + nu/2 + 1)/gamma(m/2 + nu/2))
    assert inverse_hankel_transform(2**(-m + 1)*k**(
        m - 2)*gamma(-m/2 + nu/2 + 1)/gamma(m/2 + nu/2), k, r, nu) == r**(-m)

    assert hankel_transform(r**nu*exp(-a*r), r, k, nu) == \
        2**(nu + 1)*a*k**(-nu - 3)*(a**2/k**2 + 1)**(-nu -
                                                     Rational(3, 2))*gamma(nu + Rational(3, 2))/sqrt(pi)
    assert inverse_hankel_transform(
        2**(nu + 1)*a*k**(-nu - 3)*(a**2/k**2 + 1)**(-nu - Rational(3, 2))*gamma(
            nu + Rational(3, 2))/sqrt(pi), k, r, nu) == r**nu*exp(-a*r)


def test_sympyissue_7181():
    assert mellin_transform(1/(1 - x), x, s) is not None


def test_laplace_transform_2():
    LT = laplace_transform
    # issue sympy/sympy#7173
    assert LT(sinh(a*x)*cosh(a*x), x, s) == \
        (a/(s**2 - 4*a**2), 0,
         And(Or(abs(periodic_argument(exp_polar(I*pi)*polar_lift(a), oo)) <
                pi/2, abs(periodic_argument(exp_polar(I*pi)*polar_lift(a), oo)) <=
                pi/2), Or(abs(periodic_argument(a, oo)) < pi/2,
                          abs(periodic_argument(a, oo)) <= pi/2)))

    # issues sympy/sympy#8368 and sympy/sympy#7173
    assert LT(sinh(x)*cosh(x), x, s) == (1/(s**2 - 4), 2, Ne(s/2, 1))


def test_sympyissue_8514():
    a, b, c, = symbols('a b c', positive=True, finite=True)
    t = symbols('t', positive=True)
    ft = simplify(inverse_laplace_transform(1/(a*s**2 + b*s + c), s, t))
    assert ft.rewrite(atan2) == ((exp(t*(exp(I*atan2(0, -4*a*c + b**2)/2) -
                                         exp(-I*atan2(0, -4*a*c + b**2)/2)) *
                                      sqrt(abs(4*a*c - b**2))/(4*a))*exp(t*cos(atan2(0, -4*a*c + b**2)/2)
                                                                         * sqrt(abs(4*a*c - b**2))/a) + I*sin(t*sin(atan2(0, -4*a*c + b**2)/2)
                                                                                                              * sqrt(abs(4*a*c - b**2))/(2*a)) - cos(t*sin(atan2(0, -4*a*c + b**2)/2)
                                                                                                                                                     * sqrt(abs(4*a*c - b**2))/(2*a)))*exp(-t*(b + cos(atan2(0, -4*a*c + b**2)/2)
                                                                                                                                                                                               * sqrt(abs(4*a*c - b**2)))/(2*a))/sqrt(-4*a*c + b**2))


def test__simplifyconds():
    assert _simplifyconds(1 < abs(x), x, 1) is True
    assert _simplifyconds(abs(x**2) < 1, x, 0) == (abs(x**2) < 1)


def test_sympyissue_21202():
    res = (Piecewise((s/(s**2 - 4), (4*abs(s**-2) < 1) | (abs(s**2)/4 < 1)),
                     (pi*meijerg(((Rational(1, 2),), (0, 0)), ((0, Rational(1, 2)),
                                 (0,)), s**2/4)/2, True)), 2, Ne(s**2/4, 1))
    assert laplace_transform(cosh(2*x), x, s) == res
