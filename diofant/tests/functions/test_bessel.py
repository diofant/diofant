from random import randint, uniform

import pytest

from diofant import (Float, I, Integer, Limit, O, Rational, Symbol, airyai,
                     airyaiprime, airybi, airybiprime, besseli, besselj,
                     besselk, besselsimp, bessely, cbrt, conjugate, cos, cosh,
                     diff, exp, exp_polar, expand_func, gamma, hankel1,
                     hankel2, hyper, im, jn, jn_zeros, log, nan, oo, pi,
                     polar_lift, re, root, simplify, sin, sinh, sqrt, yn, zoo)
from diofant.abc import k, n, x, y, z
from diofant.core.function import ArgumentIndexError
from diofant.functions.special.bessel import fn
from diofant.utilities.randtest import random_complex_number as randcplx
from diofant.utilities.randtest import verify_derivative_numerically as td
from diofant.utilities.randtest import verify_numerically as tn


__all__ = ()


def test_bessel_rand():
    for f in [besselj, bessely, besseli, besselk, hankel1, hankel2, jn, yn]:
        assert td(f(randcplx(), z), z)


def test_bessel_twoinputs():
    for f in [besselj, bessely, besseli, besselk, hankel1, hankel2, jn, yn]:
        pytest.raises(TypeError, lambda: f(1))
        pytest.raises(TypeError, lambda: f(1, 2, 3))


def test_diff():
    assert besselj(n, z).diff(z) == besselj(n - 1, z)/2 - besselj(n + 1, z)/2
    assert bessely(n, z).diff(z) == bessely(n - 1, z)/2 - bessely(n + 1, z)/2
    assert besseli(n, z).diff(z) == besseli(n - 1, z)/2 + besseli(n + 1, z)/2
    assert besselk(n, z).diff(z) == -besselk(n - 1, z)/2 - besselk(n + 1, z)/2
    assert hankel1(n, z).diff(z) == hankel1(n - 1, z)/2 - hankel1(n + 1, z)/2
    assert hankel2(n, z).diff(z) == hankel2(n - 1, z)/2 - hankel2(n + 1, z)/2

    pytest.raises(ArgumentIndexError, lambda: besselj(n, z).fdiff(3))
    pytest.raises(ArgumentIndexError, lambda: jn(n, z).fdiff(3))
    pytest.raises(ArgumentIndexError, lambda: airyai(z).fdiff(2))
    pytest.raises(ArgumentIndexError, lambda: airybi(z).fdiff(2))
    pytest.raises(ArgumentIndexError, lambda: airyaiprime(z).fdiff(2))
    pytest.raises(ArgumentIndexError, lambda: airybiprime(z).fdiff(2))


def test_rewrite():
    assert besselj(n, z).rewrite(jn) == sqrt(2*z/pi)*jn(n - Rational(1, 2), z)
    assert bessely(n, z).rewrite(yn) == sqrt(2*z/pi)*yn(n - Rational(1, 2), z)
    assert besseli(n, z).rewrite(besselj) == \
        exp(-I*n*pi/2)*besselj(n, polar_lift(I)*z)
    assert besselj(n, z).rewrite(besseli) == \
        exp(I*n*pi/2)*besseli(n, polar_lift(-I)*z)
    assert besselj(2, z).rewrite(bessely) == besselj(2, z)
    assert bessely(2, z).rewrite(besselj) == bessely(2, z)
    assert bessely(2, z).rewrite(besseli) == bessely(2, z)
    assert besselk(2, z).rewrite(besseli) == besselk(2, z)
    assert besselk(2, z).rewrite(besselj) == besselk(2, z)
    assert besselk(2, z).rewrite(bessely) == besselk(2, z)

    nu = randcplx()

    assert tn(besselj(nu, z), besselj(nu, z).rewrite(besseli), z)
    assert tn(besselj(nu, z), besselj(nu, z).rewrite(bessely), z)

    assert tn(besseli(nu, z), besseli(nu, z).rewrite(besselj), z)
    assert tn(besseli(nu, z), besseli(nu, z).rewrite(bessely), z)

    assert tn(bessely(nu, z), bessely(nu, z).rewrite(besselj), z)
    assert tn(bessely(nu, z), bessely(nu, z).rewrite(besseli), z)

    assert tn(besselk(nu, z), besselk(nu, z).rewrite(besselj), z)
    assert tn(besselk(nu, z), besselk(nu, z).rewrite(besseli), z)
    assert tn(besselk(nu, z), besselk(nu, z).rewrite(bessely), z)


def test_expand():
    assert expand_func(besselj(Rational(1, 2), z).rewrite(jn)) == \
        sqrt(2)*sin(z)/(sqrt(pi)*sqrt(z))
    assert expand_func(bessely(Rational(1, 2), z).rewrite(yn)) == \
        -sqrt(2)*cos(z)/(sqrt(pi)*sqrt(z))
    assert expand_func(besselj(I, z)) == besselj(I, z)

    # Test simplify helper
    assert simplify(besselj(Rational(1, 2), z)) == sqrt(2)*sin(z)/(sqrt(pi)*sqrt(z))

    # XXX: teach sin/cos to work around arguments like
    # x*exp_polar(I*pi*n/2).  Then change besselsimp -> expand_func
    assert besselsimp(besselj(Rational(1, 2), z)) == sqrt(2)*sin(z)/(sqrt(pi)*sqrt(z))
    assert besselsimp(besselj(Rational(-1, 2), z)) == sqrt(2)*cos(z)/(sqrt(pi)*sqrt(z))
    assert besselsimp(besselj(Rational(5, 2), z)) == \
        -sqrt(2)*(z**2*sin(z) + 3*z*cos(z) - 3*sin(z))/(sqrt(pi)*z**Rational(5, 2))
    assert besselsimp(besselj(-Rational(5, 2), z)) == \
        -sqrt(2)*(z**2*cos(z) - 3*z*sin(z) - 3*cos(z))/(sqrt(pi)*z**Rational(5, 2))

    assert besselsimp(bessely(Rational(1, 2), z)) == \
        -(sqrt(2)*cos(z))/(sqrt(pi)*sqrt(z))
    assert besselsimp(bessely(Rational(-1, 2), z)) == sqrt(2)*sin(z)/(sqrt(pi)*sqrt(z))
    assert besselsimp(bessely(Rational(5, 2), z)) == \
        sqrt(2)*(z**2*cos(z) - 3*z*sin(z) - 3*cos(z))/(sqrt(pi)*z**Rational(5, 2))
    assert besselsimp(bessely(Rational(-5, 2), z)) == \
        -sqrt(2)*(z**2*sin(z) + 3*z*cos(z) - 3*sin(z))/(sqrt(pi)*z**Rational(5, 2))

    assert besselsimp(besseli(Rational(1, 2), z)) == sqrt(2)*sinh(z)/(sqrt(pi)*sqrt(z))
    assert besselsimp(besseli(Rational(-1, 2), z)) == \
        sqrt(2)*cosh(z)/(sqrt(pi)*sqrt(z))
    assert besselsimp(besseli(Rational(5, 2), z)) == \
        sqrt(2)*(z**2*sinh(z) - 3*z*cosh(z) + 3*sinh(z))/(sqrt(pi)*z**Rational(5, 2))
    assert besselsimp(besseli(Rational(-5, 2), z)) == \
        sqrt(2)*(z**2*cosh(z) - 3*z*sinh(z) + 3*cosh(z))/(sqrt(pi)*z**Rational(5, 2))

    assert besselsimp(besselk(Rational(1, 2), z)) == \
        besselsimp(besselk(Rational(-1, 2), z)) == sqrt(pi)*exp(-z)/(sqrt(2)*sqrt(z))
    assert besselsimp(besselk(Rational(5, 2), z)) == \
        besselsimp(besselk(Rational(-5, 2), z)) == \
        sqrt(2)*sqrt(pi)*(z**2 + 3*z + 3)*exp(-z)/(2*z**Rational(5, 2))

    def check(eq, ans):
        return tn(eq, ans) and eq == ans

    rn = randcplx(a=1, b=0, d=0, c=2)

    for besselx in [besselj, bessely, besseli, besselk]:
        ri = Rational(2*randint(-11, 10) + 1, 2)  # half integer in [-21/2, 21/2]
        assert tn(besselsimp(besselx(ri, z)), besselx(ri, z))

    assert check(expand_func(besseli(rn, x)),
                 besseli(rn - 2, x) - 2*(rn - 1)*besseli(rn - 1, x)/x)
    assert check(expand_func(besseli(-rn, x)),
                 besseli(-rn + 2, x) + 2*(-rn + 1)*besseli(-rn + 1, x)/x)

    assert check(expand_func(besselj(rn, x)),
                 -besselj(rn - 2, x) + 2*(rn - 1)*besselj(rn - 1, x)/x)
    assert check(expand_func(besselj(-rn, x)),
                 -besselj(-rn + 2, x) + 2*(-rn + 1)*besselj(-rn + 1, x)/x)

    assert check(expand_func(besselk(rn, x)),
                 besselk(rn - 2, x) + 2*(rn - 1)*besselk(rn - 1, x)/x)
    assert check(expand_func(besselk(-rn, x)),
                 besselk(-rn + 2, x) - 2*(-rn + 1)*besselk(-rn + 1, x)/x)

    assert check(expand_func(bessely(rn, x)),
                 -bessely(rn - 2, x) + 2*(rn - 1)*bessely(rn - 1, x)/x)
    assert check(expand_func(bessely(-rn, x)),
                 -bessely(-rn + 2, x) + 2*(-rn + 1)*bessely(-rn + 1, x)/x)

    n = Symbol('n', integer=True, positive=True)

    assert expand_func(besseli(n + 2, z)) == \
        besseli(n, z) + (-2*n - 2)*(-2*n*besseli(n, z)/z + besseli(n - 1, z))/z
    assert expand_func(besselj(n + 2, z)) == \
        -besselj(n, z) + (2*n + 2)*(2*n*besselj(n, z)/z - besselj(n - 1, z))/z
    assert expand_func(besselk(n + 2, z)) == \
        besselk(n, z) + (2*n + 2)*(2*n*besselk(n, z)/z + besselk(n - 1, z))/z
    assert expand_func(bessely(n + 2, z)) == \
        -bessely(n, z) + (2*n + 2)*(2*n*bessely(n, z)/z - bessely(n - 1, z))/z

    assert expand_func(besseli(n + Rational(1, 2), z).rewrite(jn)) == \
        (sqrt(2)*sqrt(z)*exp(-I*pi*(n + Rational(1, 2))/2) *
         exp_polar(I*pi/4)*jn(n, z*exp_polar(I*pi/2))/sqrt(pi))
    assert expand_func(besselj(n + Rational(1, 2), z).rewrite(jn)) == \
        sqrt(2)*sqrt(z)*jn(n, z)/sqrt(pi)

    r = Symbol('r', extended_real=True)
    p = Symbol('p', positive=True)
    i = Symbol('i', integer=True)

    for besselx in [besselj, bessely, besseli, besselk]:
        assert besselx(i, p).is_extended_real
        assert besselx(i, x).is_extended_real is None
        assert besselx(x, z).is_extended_real is None

        # issue sympy/sympy#21486
        assert expand_func(besselx(oo, x)) == besselx(oo, x, evaluate=False)
        assert expand_func(besselx(-oo, x)) == besselx(-oo, x, evaluate=False)

    for besselx in [besselj, besseli]:
        assert besselx(i, r).is_extended_real
    for besselx in [bessely, besselk]:
        assert besselx(i, r).is_extended_real is None


def test_fn():
    assert fn(1, z) == 1/z**2
    assert fn(2, z) == -1/z + 3/z**3
    assert fn(3, z) == -6/z**2 + 15/z**4
    assert fn(4, z) == 1/z - 45/z**3 + 105/z**5


def mjn(n, z):
    return expand_func(jn(n, z))


def myn(n, z):
    return expand_func(yn(n, z))


def test_jn():
    assert mjn(0, z) == sin(z)/z
    assert mjn(1, z) == sin(z)/z**2 - cos(z)/z
    assert mjn(2, z) == (3/z**3 - 1/z)*sin(z) - (3/z**2) * cos(z)
    assert mjn(3, z) == (15/z**4 - 6/z**2)*sin(z) + (1/z - 15/z**3)*cos(z)
    assert mjn(4, z) == (1/z + 105/z**5 - 45/z**3)*sin(z) + \
        (-105/z**4 + 10/z**2)*cos(z)
    assert mjn(5, z) == (945/z**6 - 420/z**4 + 15/z**2)*sin(z) + \
        (-1/z - 945/z**5 + 105/z**3)*cos(z)
    assert mjn(6, z) == (-1/z + 10395/z**7 - 4725/z**5 + 210/z**3)*sin(z) + \
        (-10395/z**6 + 1260/z**4 - 21/z**2)*cos(z)

    assert expand_func(jn(n, z)) == jn(n, z)


def test_yn():
    assert myn(0, z) == -cos(z)/z
    assert myn(1, z) == -cos(z)/z**2 - sin(z)/z
    assert myn(2, z) == -((3/z**3 - 1/z)*cos(z) + (3/z**2)*sin(z))
    assert expand_func(yn(n, z)) == yn(n, z)


def test_sympify_yn():
    assert Integer(15) in myn(3, pi).atoms()
    assert myn(3, pi) == 15/pi**4 - 6/pi**2


def eq(a, b, tol=1e-6):
    for x, y in zip(a, b):
        if not abs(x - y) < tol:
            return False
    return True


def test_jn_zeros():
    assert eq(jn_zeros(0, 4), [3.141592, 6.283185, 9.424777, 12.566370])
    assert eq(jn_zeros(1, 4), [4.493409, 7.725251, 10.904121, 14.066193])
    assert eq(jn_zeros(2, 4), [5.763459, 9.095011, 12.322940, 15.514603])
    assert eq(jn_zeros(3, 4), [6.987932, 10.417118, 13.698023, 16.923621])
    assert eq(jn_zeros(4, 4), [8.182561, 11.704907, 15.039664, 18.301255])


def test_bessel_eval():
    n, m, k = Symbol('n', integer=True), Symbol('m'), Symbol('k', integer=True, zero=False)

    for f in [besselj, besseli]:
        assert f(0, 0) == 1
        assert f(2.1, 0) == 0
        assert f(-3, 0) == 0
        assert f(-10.2, 0) == zoo
        assert f(1 + 3*I, 0) == 0
        assert f(-3 + I, 0) == zoo
        assert f(-2*I, 0) == nan
        assert f(n, 0) not in (0, 1)
        assert f(m, 0) not in (0, 1)
        assert f(k, 0) == 0

    assert bessely(0, 0) == -oo
    assert besselk(0, 0) == +oo
    for f in [bessely, besselk]:
        assert f(1 + I, 0) == zoo
        assert f(I, 0) == nan

    for x in (oo, -oo, I*oo, -I*oo):
        assert besselk(m, x) == 0

    for f in [besselj, bessely]:
        assert f(m, +oo) == 0
        assert f(m, -oo) == 0

    for f in [besseli, besselk]:
        assert f(m, +I*oo) == 0
        assert f(m, -I*oo) == 0

    for f in [besseli, besselk]:
        assert f(-4, z) == f(4, z)
        assert f(-3, z) == f(3, z)
        assert f(-n, z) == f(n, z)
        assert f(-m, z) != f(m, z)

    for f in [besselj, bessely]:
        assert f(-4, z) == f(4, z)
        assert f(-3, z) == -f(3, z)
        assert f(-n, z) == (-1)**n*f(n, z)
        assert f(-m, z) != (-1)**m*f(m, z)

    for f in [besselj, besseli]:
        assert f(m, -z) == (-z)**m*z**(-m)*f(m, z)

    assert besseli(2, -z) == besseli(2, z)
    assert besseli(3, -z) == -besseli(3, z)

    assert besselj(0, -z) == besselj(0, z)
    assert besselj(1, -z) == -besselj(1, z)

    assert besseli(0, I*z) == besselj(0, z)
    assert besseli(1, I*z) == I*besselj(1, z)
    assert besselj(3, I*z) == -I*besseli(3, z)


def test_bessel_nan():
    for f in [besselj, bessely, besseli, besselk, hankel1, hankel2, yn, jn]:
        assert f(1, nan) == nan


def test_conjugate():
    n, z, x = Symbol('n'), Symbol('z', extended_real=False), Symbol('x', extended_real=True)
    y, t = Symbol('y', extended_real=True, positive=True), Symbol('t', negative=True)

    for f in [besseli, besselj, besselk, bessely, jn, yn, hankel1, hankel2]:
        assert f(n, -1).conjugate() != f(conjugate(n), -1)
        assert f(n, x).conjugate() != f(conjugate(n), x)
        assert f(n, t).conjugate() != f(conjugate(n), t)

    rz = randcplx(b=0.5)

    for f in [besseli, besselj, besselk, bessely, jn, yn]:
        assert f(n, 1 + I).conjugate() == f(conjugate(n), 1 - I)
        assert f(n, 0).conjugate() == f(conjugate(n), 0)
        assert f(n, 1).conjugate() == f(conjugate(n), 1)
        assert f(n, z).conjugate() == f(conjugate(n), conjugate(z))
        assert f(n, y).conjugate() == f(conjugate(n), y)
        assert tn(f(n, rz).conjugate(), f(conjugate(n), conjugate(rz)))

    assert hankel1(n, 1 + I).conjugate() == hankel2(conjugate(n), 1 - I)
    assert hankel1(n, 0).conjugate() == hankel2(conjugate(n), 0)
    assert hankel1(n, 1).conjugate() == hankel2(conjugate(n), 1)
    assert hankel1(n, y).conjugate() == hankel2(conjugate(n), y)
    assert hankel1(n, z).conjugate() == hankel2(conjugate(n), conjugate(z))
    assert tn(hankel1(n, rz).conjugate(), hankel2(conjugate(n), conjugate(rz)))

    assert hankel2(n, 1 + I).conjugate() == hankel1(conjugate(n), 1 - I)
    assert hankel2(n, 0).conjugate() == hankel1(conjugate(n), 0)
    assert hankel2(n, 1).conjugate() == hankel1(conjugate(n), 1)
    assert hankel2(n, y).conjugate() == hankel1(conjugate(n), y)
    assert hankel2(n, z).conjugate() == hankel1(conjugate(n), conjugate(z))
    assert tn(hankel2(n, rz).conjugate(), hankel1(conjugate(n), conjugate(rz)))


def test_branching():
    assert besselj(polar_lift(k), x) == besselj(k, x)
    assert besseli(polar_lift(k), x) == besseli(k, x)

    n = Symbol('n', integer=True)
    assert besselj(n, exp_polar(2*pi*I)*x) == besselj(n, x)
    assert besselj(n, polar_lift(x)) == besselj(n, x)
    assert besseli(n, exp_polar(2*pi*I)*x) == besseli(n, x)
    assert besseli(n, polar_lift(x)) == besseli(n, x)

    def tn(func, s):
        c = uniform(1, 5)
        expr = func(s, c*exp_polar(I*pi)) - func(s, c*exp_polar(-I*pi))
        eps = 1e-15
        expr2 = func(s + eps, -c + eps*I) - func(s + eps, -c - eps*I)
        return abs(expr - expr2).evalf(strict=False) < 1e-10

    nu = Symbol('nu')
    assert besselj(nu, exp_polar(2*pi*I)*x) == exp(2*pi*I*nu)*besselj(nu, x)
    assert besseli(nu, exp_polar(2*pi*I)*x) == exp(2*pi*I*nu)*besseli(nu, x)
    assert tn(besselj, 2)
    assert tn(besselj, pi)
    assert tn(besselj, I)
    assert tn(besseli, 2)
    assert tn(besseli, pi)
    assert tn(besseli, I)


def test_airy_base():
    z = Symbol('z')
    x = Symbol('x', extended_real=True)
    y = Symbol('y', extended_real=True)

    assert conjugate(airyai(z)) == airyai(conjugate(z))
    assert airyai(x).is_extended_real
    assert airyai(z).is_extended_real is None

    assert airyai(x+I*y).as_real_imag() == (
        airyai(x - I*x*abs(y)/abs(x))/2 + airyai(x + I*x*abs(y)/abs(x))/2,
        I*x*(airyai(x - I*x*abs(y)/abs(x)) -
             airyai(x + I*x*abs(y)/abs(x)))*abs(y)/(2*y*abs(x)))


def test_airyai():
    z = Symbol('z', extended_real=False)
    r = Symbol('r', extended_real=True)
    t = Symbol('t', negative=True)
    p = Symbol('p', positive=True)

    assert isinstance(airyai(z), airyai)

    assert airyai(0) == cbrt(3)/(3*gamma(Rational(2, 3)))
    assert airyai(oo) == 0
    assert airyai(-oo) == 0

    assert diff(airyai(z), z) == airyaiprime(z)

    assert airyai(z).series(z, 0, 3) == (
        3**Rational(5, 6)*gamma(Rational(1, 3))/(6*pi) - root(3, 6)*z*gamma(Rational(2, 3))/(2*pi) + O(z**3))

    l = Limit(airyai(I/x)/(exp(-Rational(2, 3)*(I/x)**Rational(3, 2))*sqrt(pi*sqrt(I/x))/2), x, 0)
    assert l.doit() == l  # cover _airyais._eval_aseries

    assert airyai(z).rewrite(hyper) == (
        -3**Rational(2, 3)*z*hyper((), (Rational(4, 3),), z**3/9)/(3*gamma(Rational(1, 3))) +
        cbrt(3)*hyper((), (Rational(2, 3),), z**3/9)/(3*gamma(Rational(2, 3))))

    assert isinstance(airyai(z).rewrite(besselj), airyai)
    assert airyai(t).rewrite(besselj) == (
        sqrt(-t)*(besselj(-Rational(1, 3), 2*(-t)**Rational(3, 2)/3) +
                  besselj(Rational(1, 3), 2*(-t)**Rational(3, 2)/3))/3)
    assert airyai(z).rewrite(besseli) == (
        -z*besseli(Rational(1, 3), 2*z**Rational(3, 2)/3)/(3*cbrt(z**Rational(3, 2))) +
        cbrt(z**Rational(3, 2))*besseli(-Rational(1, 3), 2*z**Rational(3, 2)/3)/3)
    assert airyai(p).rewrite(besseli) == (
        sqrt(p)*(besseli(-Rational(1, 3), 2*p**Rational(3, 2)/3) -
                 besseli(Rational(1, 3), 2*p**Rational(3, 2)/3))/3)

    assert expand_func(airyai(2*cbrt(3*z**5))) == (
        -sqrt(3)*(-1 + cbrt(z**5)/z**Rational(5, 3))*airybi(2*cbrt(3)*z**Rational(5, 3))/6 +
        (1 + cbrt(z**5)/z**Rational(5, 3))*airyai(2*cbrt(3)*z**Rational(5, 3))/2)
    assert expand_func(airyai(x*y)) == airyai(x*y)
    assert expand_func(airyai(log(x))) == airyai(log(x))
    assert expand_func(airyai(2*root(3*z**5, 5))) == airyai(2*root(3*z**5, 5))

    assert (airyai(r).as_real_imag() ==
            airyai(r).as_real_imag(deep=False) == (airyai(r), 0))
    assert airyai(x).as_real_imag() == airyai(x).as_real_imag(deep=False)
    assert (airyai(x).as_real_imag() ==
            (airyai(re(x) - I*re(x)*abs(im(x))/abs(re(x)))/2 +
             airyai(re(x) + I*re(x)*abs(im(x))/abs(re(x)))/2,
             I*(airyai(re(x) - I*re(x)*abs(im(x))/abs(re(x))) -
                airyai(re(x) + I*re(x)*abs(im(x))/abs(re(x)))) *
             re(x)*abs(im(x))/(2*im(x)*abs(re(x)))))

    assert airyai(x).taylor_term(-1, x) == 0


def test_airybi():
    z = Symbol('z', extended_real=False)
    t = Symbol('t', negative=True)
    p = Symbol('p', positive=True)

    assert isinstance(airybi(z), airybi)

    assert airybi(0) == 3**Rational(5, 6)/(3*gamma(Rational(2, 3)))
    assert airybi(oo) == oo
    assert airybi(-oo) == 0

    assert diff(airybi(z), z) == airybiprime(z)

    assert airybi(z).series(z, 0, 3) == (
        cbrt(3)*gamma(Rational(1, 3))/(2*pi) + 3**Rational(2, 3)*z*gamma(Rational(2, 3))/(2*pi) + O(z**3))
    l = Limit(airybi(I/x)/(exp(Rational(2, 3)*(I/x)**Rational(3, 2))*sqrt(pi*sqrt(I/x))), x, 0)
    assert l.doit() == l

    assert airybi(z).rewrite(hyper) == (
        root(3, 6)*z*hyper((), (Rational(4, 3),), z**3/9)/gamma(Rational(1, 3)) +
        3**Rational(5, 6)*hyper((), (Rational(2, 3),), z**3/9)/(3*gamma(Rational(2, 3))))

    assert isinstance(airybi(z).rewrite(besselj), airybi)
    assert (airybi(t).rewrite(besselj) ==
            sqrt(3)*sqrt(-t)*(besselj(-1/3, 2*(-t)**Rational(3, 2)/3) -
                              besselj(Rational(1, 3),
                                      2*(-t)**Rational(3, 2)/3))/3)
    assert airybi(z).rewrite(besseli) == (
        sqrt(3)*(z*besseli(Rational(1, 3), 2*z**Rational(3, 2)/3)/cbrt(z**Rational(3, 2)) +
                 cbrt(z**Rational(3, 2))*besseli(-Rational(1, 3), 2*z**Rational(3, 2)/3))/3)
    assert airybi(p).rewrite(besseli) == (
        sqrt(3)*sqrt(p)*(besseli(-Rational(1, 3), 2*p**Rational(3, 2)/3) +
                         besseli(Rational(1, 3), 2*p**Rational(3, 2)/3))/3)
    assert airybi(p).rewrite(besselj) == airybi(p)

    assert expand_func(airybi(2*cbrt(3*z**5))) == (
        sqrt(3)*(1 - cbrt(z**5)/z**Rational(5, 3))*airyai(2*cbrt(3)*z**Rational(5, 3))/2 +
        (1 + cbrt(z**5)/z**Rational(5, 3))*airybi(2*cbrt(3)*z**Rational(5, 3))/2)
    assert expand_func(airybi(x*y)) == airybi(x*y)
    assert expand_func(airybi(log(x))) == airybi(log(x))
    assert expand_func(airybi(2*root(3*z**5, 5))) == airybi(2*root(3*z**5, 5))

    assert airybi(x).taylor_term(-1, x) == 0


def test_airyaiprime():
    z = Symbol('z', extended_real=False)
    t = Symbol('t', negative=True)
    p = Symbol('p', positive=True)

    assert isinstance(airyaiprime(z), airyaiprime)

    assert airyaiprime(0) == -3**Rational(2, 3)/(3*gamma(Rational(1, 3)))
    assert airyaiprime(oo) == 0

    assert diff(airyaiprime(z), z) == z*airyai(z)

    assert airyaiprime(z).series(z, 0, 3) == (
        -3**Rational(2, 3)/(3*gamma(Rational(1, 3))) + cbrt(3)*z**2/(6*gamma(Rational(2, 3))) + O(z**3))

    assert airyaiprime(z).rewrite(hyper) == (
        cbrt(3)*z**2*hyper((), (Rational(5, 3),), z**3/9)/(6*gamma(Rational(2, 3))) -
        3**Rational(2, 3)*hyper((), (Rational(1, 3),), z**3/9)/(3*gamma(Rational(1, 3))))

    assert isinstance(airyaiprime(z).rewrite(besselj), airyaiprime)
    assert (airyaiprime(t).rewrite(besselj) ==
            t*(besselj(-Rational(2, 3), 2*(-t)**Rational(3, 2)/3) -
               besselj(Rational(2, 3), 2*(-t)**Rational(3, 2)/3))/3)
    assert airyaiprime(z).rewrite(besseli) == (
        z**2*besseli(Rational(2, 3), 2*z**Rational(3, 2)/3)/(3*(z**Rational(3, 2))**Rational(2, 3)) -
        (z**Rational(3, 2))**Rational(2, 3)*besseli(-Rational(1, 3), 2*z**Rational(3, 2)/3)/3)
    assert airyaiprime(p).rewrite(besseli) == (
        p*(-besseli(-Rational(2, 3), 2*p**Rational(3, 2)/3) + besseli(Rational(2, 3), 2*p**Rational(3, 2)/3))/3)
    assert airyaiprime(p).rewrite(besselj) == airyaiprime(p)

    assert expand_func(airyaiprime(2*cbrt(3*z**5))) == (
        sqrt(3)*(z**Rational(5, 3)/cbrt(z**5) - 1)*airybiprime(2*cbrt(3)*z**Rational(5, 3))/6 +
        (z**Rational(5, 3)/cbrt(z**5) + 1)*airyaiprime(2*cbrt(3)*z**Rational(5, 3))/2)
    assert expand_func(airyaiprime(x*y)) == airyaiprime(x*y)
    assert expand_func(airyaiprime(log(x))) == airyaiprime(log(x))
    assert expand_func(airyaiprime(2*root(3*z**5, 5))) == airyaiprime(2*root(3*z**5, 5))

    assert airyaiprime(-2).evalf(50) == Float('0.61825902074169104140626429133247528291577794512414753', dps=50)


def test_airybiprime():
    z = Symbol('z', extended_real=False)
    t = Symbol('t', negative=True)
    p = Symbol('p', positive=True)

    assert isinstance(airybiprime(z), airybiprime)

    assert airybiprime(0) == root(3, 6)/gamma(Rational(1, 3))
    assert airybiprime(oo) == oo
    assert airybiprime(-oo) == 0

    assert diff(airybiprime(z), z) == z*airybi(z)

    assert airybiprime(z).series(z, 0, 3) == (
        root(3, 6)/gamma(Rational(1, 3)) + 3**Rational(5, 6)*z**2/(6*gamma(Rational(2, 3))) + O(z**3))

    assert airybiprime(z).rewrite(hyper) == (
        3**Rational(5, 6)*z**2*hyper((), (Rational(5, 3),), z**3/9)/(6*gamma(Rational(2, 3))) +
        root(3, 6)*hyper((), (Rational(1, 3),), z**3/9)/gamma(Rational(1, 3)))

    assert isinstance(airybiprime(z).rewrite(besselj), airybiprime)
    assert (airybiprime(t).rewrite(besselj) ==
            -sqrt(3)*t*(besselj(-Rational(2, 3), 2*(-t)**Rational(3, 2)/3) +
                        besselj(Rational(2, 3), 2*(-t)**Rational(3, 2)/3))/3)
    assert airybiprime(z).rewrite(besseli) == (
        sqrt(3)*(z**2*besseli(Rational(2, 3), 2*z**Rational(3, 2)/3)/(z**Rational(3, 2))**Rational(2, 3) +
                 (z**Rational(3, 2))**Rational(2, 3)*besseli(-Rational(2, 3), 2*z**Rational(3, 2)/3))/3)
    assert airybiprime(p).rewrite(besseli) == (
        sqrt(3)*p*(besseli(-Rational(2, 3), 2*p**Rational(3, 2)/3) + besseli(Rational(2, 3), 2*p**Rational(3, 2)/3))/3)
    assert airybiprime(p).rewrite(besselj) == airybiprime(p)

    assert expand_func(airybiprime(2*cbrt(3*z**5))) == (
        sqrt(3)*(z**Rational(5, 3)/cbrt(z**5) - 1)*airyaiprime(2*cbrt(3)*z**Rational(5, 3))/2 +
        (z**Rational(5, 3)/cbrt(z**5) + 1)*airybiprime(2*cbrt(3)*z**Rational(5, 3))/2)
    assert expand_func(airybiprime(x*y)) == airybiprime(x*y)
    assert expand_func(airybiprime(log(x))) == airybiprime(log(x))
    assert expand_func(airybiprime(2*root(3*z**5, 5))) == airybiprime(2*root(3*z**5, 5))

    assert airybiprime(-2).evalf(50) == Float('0.27879516692116952268509756941098324140300059345163131', dps=50)
