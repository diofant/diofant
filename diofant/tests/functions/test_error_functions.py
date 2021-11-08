from random import uniform

import pytest

from diofant import (E1, Chi, Ci, E, Ei, EulerGamma, Float, I, Integer, Li,
                     Limit, O, Rational, Shi, Si, Symbol, conjugate, cos, cosh,
                     diff, erf, erf2, erf2inv, erfc, erfcinv, erfi, erfinv,
                     exp, exp_polar, expand, expand_func, expint, fresnelc,
                     fresnels, gamma, hyper, im, integrate, li, limit, log,
                     meijerg, nan, oo, pi, polar_lift, re, root, sign, sin,
                     sinh, sqrt, uppergamma)
from diofant.abc import x, y, z
from diofant.core.function import ArgumentIndexError
from diofant.functions.special.error_functions import _eis, _erfs
from diofant.utilities.randtest import (random_complex_number,
                                        verify_derivative_numerically,
                                        verify_numerically)


__all__ = ()

w = Symbol('w', extended_real=True)
n = Symbol('n', integer=True)


def test_erf():
    assert erf(nan) == nan

    assert erf(oo) == 1
    assert erf(-oo) == -1

    assert erf(0) == 0

    assert erf(I*oo) == oo*I
    assert erf(-I*oo) == -oo*I

    assert erf(-2) == -erf(2)
    assert erf(-x*y) == -erf(x*y)
    assert erf(-x - y) == -erf(x + y)

    assert erf(erfinv(x)) == x
    assert erf(erfcinv(x)) == 1 - x
    assert erf(erf2inv(0, x)) == x
    assert erf(erf2inv(0, erf(erfcinv(1 - erf(erfinv(x)))))) == x

    assert erf(I).is_extended_real is False
    assert erf(w).is_extended_real is True
    assert erf(z).is_extended_real is None

    assert conjugate(erf(z)) == erf(conjugate(z))

    assert erf(x).as_leading_term(x) == 2*x/sqrt(pi)
    assert erf(1/x).as_leading_term(x) == erf(1/x)

    assert erf(z).rewrite('uppergamma') == sqrt(z**2)*erf(sqrt(z**2))/z
    assert erf(z).rewrite('erfc') == 1 - erfc(z)
    assert erf(z).rewrite('erfi') == -I*erfi(I*z)
    assert erf(z).rewrite('fresnels') == (1 + I)*(fresnelc(z*(1 - I)/sqrt(pi)) -
                                                  I*fresnels(z*(1 - I)/sqrt(pi)))
    assert erf(z).rewrite('fresnelc') == (1 + I)*(fresnelc(z*(1 - I)/sqrt(pi)) -
                                                  I*fresnels(z*(1 - I)/sqrt(pi)))
    assert erf(z).rewrite('hyper') == 2*z*hyper([Rational(1, 2)], [Rational(3, 2)], -z**2)/sqrt(pi)
    assert erf(z).rewrite('meijerg') == z*meijerg([Rational(1, 2)], [], [0], [Rational(-1, 2)], z**2)/sqrt(pi)
    assert erf(z).rewrite('expint') == sqrt(z**2)/z - z*expint(Rational(1, 2), z**2)/sqrt(pi)

    assert limit(exp(x)*exp(x**2)*(erf(x + 1/exp(x)) - erf(x)), x, oo) == \
        2/sqrt(pi)
    assert limit((1 - erf(z))*exp(z**2)*z, z, oo) == 1/sqrt(pi)
    assert limit((1 - erf(x))*exp(x**2)*sqrt(pi)*x, x, oo) == 1
    assert limit(((1 - erf(x))*exp(x**2)*sqrt(pi)*x - 1)*2*x**2, x, oo) == -1

    l = Limit((1 - erf(y/x))*exp(y**2/x**2), x, 0)
    assert l.doit() == l  # cover _erfs._eval_aseries

    assert erf(x).as_real_imag() == \
        ((erf(re(x) - I*re(x)*abs(im(x))/abs(re(x)))/2 +
          erf(re(x) + I*re(x)*abs(im(x))/abs(re(x)))/2,
          I*(erf(re(x) - I*re(x)*abs(im(x))/abs(re(x))) -
             erf(re(x) + I*re(x)*abs(im(x))/abs(re(x)))) *
          re(x)*abs(im(x))/(2*im(x)*abs(re(x)))))
    assert erf(x).as_real_imag() == erf(x).as_real_imag(deep=False)
    assert erf(w).as_real_imag() == (erf(w), 0)
    assert erf(w).as_real_imag() == erf(w).as_real_imag(deep=False)
    assert erf(I).as_real_imag() == (0, erfi(1))

    pytest.raises(ArgumentIndexError, lambda: erf(x).fdiff(2))

    assert erf(x).taylor_term(3, x, *(2*x/sqrt(pi), 0)) == -2*x**3/3/sqrt(pi)


def test_erf_series():
    assert erf(x).series(x, 0, 7) == 2*x/sqrt(pi) - \
        2*x**3/3/sqrt(pi) + x**5/5/sqrt(pi) + O(x**7)


def test_erf_evalf():
    assert abs(erf(Float(2.0)) - 0.995322265) < 1E-8  # XXX


def test__erfs():
    assert _erfs(z).diff(z) == -2/sqrt(pi) + 2*z*_erfs(z)
    pytest.raises(ArgumentIndexError, lambda: _erfs(x).fdiff(2))

    assert _erfs(1/z).series(z) == \
        z/sqrt(pi) - z**3/(2*sqrt(pi)) + 3*z**5/(4*sqrt(pi)) + O(z**6)

    assert expand(erf(z).rewrite('tractable').diff(z).rewrite('intractable')) \
        == erf(z).diff(z)
    assert _erfs(z).rewrite('intractable') == (-erf(z) + 1)*exp(z**2)


def test_erfc():
    assert erfc(nan) == nan

    assert erfc(oo) == 0
    assert erfc(-oo) == 2

    assert erfc(0) == 1

    assert erfc(I*oo) == -oo*I
    assert erfc(-I*oo) == oo*I

    assert erfc(-x) == Integer(2) - erfc(x)
    assert erfc(erfcinv(x)) == x
    assert erfc(erfinv(x)) == 1 - x

    assert erfc(I).is_extended_real is False
    assert erfc(w).is_extended_real is True
    assert erfc(z).is_extended_real is None

    assert conjugate(erfc(z)) == erfc(conjugate(z))

    assert erfc(x).as_leading_term(x) == 1
    assert erfc(1/x).as_leading_term(x) == erfc(1/x)

    assert erfc(z).rewrite('erf') == 1 - erf(z)
    assert erfc(z).rewrite('erfi') == 1 + I*erfi(I*z)
    assert erfc(z).rewrite('fresnels') == 1 - (1 + I)*(fresnelc(z*(1 - I)/sqrt(pi)) -
                                                       I*fresnels(z*(1 - I)/sqrt(pi)))
    assert erfc(z).rewrite('fresnelc') == 1 - (1 + I)*(fresnelc(z*(1 - I)/sqrt(pi)) -
                                                       I*fresnels(z*(1 - I)/sqrt(pi)))
    assert erfc(z).rewrite('hyper') == 1 - 2*z*hyper([Rational(1, 2)], [Rational(3, 2)], -z**2)/sqrt(pi)
    assert erfc(z).rewrite('meijerg') == 1 - z*meijerg([Rational(1, 2)], [], [0], [Rational(-1, 2)], z**2)/sqrt(pi)
    assert erfc(z).rewrite('uppergamma') == 1 - sqrt(z**2)*erf(sqrt(z**2))/z
    assert erfc(z).rewrite('expint') == 1 - sqrt(z**2)/z + z*expint(Rational(1, 2), z**2)/sqrt(pi)

    assert erfc(x).as_real_imag() == \
        ((erfc(re(x) - I*re(x)*abs(im(x))/abs(re(x)))/2 +
          erfc(re(x) + I*re(x)*abs(im(x))/abs(re(x)))/2,
          I*(erfc(re(x) - I*re(x)*abs(im(x))/abs(re(x))) -
             erfc(re(x) + I*re(x)*abs(im(x))/abs(re(x)))) *
          re(x)*abs(im(x))/(2*im(x)*abs(re(x)))))
    assert erfc(x).as_real_imag(deep=False) == erfc(x).as_real_imag()
    assert erfc(w).as_real_imag() == (erfc(w), 0)
    assert erfc(w).as_real_imag(deep=False) == erfc(w).as_real_imag()
    assert erfc(I).as_real_imag() == (1, -erfi(1))

    pytest.raises(ArgumentIndexError, lambda: erfc(x).fdiff(2))

    assert erfc(x).taylor_term(3, x, *(-2*x/sqrt(pi), 0)) == 2*x**3/3/sqrt(pi)

    assert erfc(x).limit(x, oo) == 0

    assert erfc(x).diff(x) == -2*exp(-x**2)/sqrt(pi)


def test_erfc_series():
    assert erfc(x).series(x, 0, 7) == 1 - 2*x/sqrt(pi) + \
        2*x**3/3/sqrt(pi) - x**5/5/sqrt(pi) + O(x**7)


def test_erfc_evalf():
    assert abs(erfc(Float(2.0)) - 0.00467773) < 1E-8  # XXX


def test_erfi():
    assert erfi(nan) == nan

    assert erfi(+oo) == +oo
    assert erfi(-oo) == -oo

    assert erfi(0) == 0

    assert erfi(+I*oo) == I
    assert erfi(-I*oo) == -I

    assert erfi(-x) == -erfi(x)

    assert erfi(I*erfinv(x)) == I*x
    assert erfi(I*erfcinv(x)) == I*(1 - x)
    assert erfi(I*erf2inv(0, x)) == I*x

    assert erfi(I).is_extended_real is False
    assert erfi(w).is_extended_real is True
    assert erfi(z).is_extended_real is None

    assert conjugate(erfi(z)) == erfi(conjugate(z))

    assert erfi(z).rewrite('erf') == -I*erf(I*z)
    assert erfi(z).rewrite('erfc') == I*erfc(I*z) - I
    assert erfi(z).rewrite('fresnels') == (1 - I)*(fresnelc(z*(1 + I)/sqrt(pi)) -
                                                   I*fresnels(z*(1 + I)/sqrt(pi)))
    assert erfi(z).rewrite('fresnelc') == (1 - I)*(fresnelc(z*(1 + I)/sqrt(pi)) -
                                                   I*fresnels(z*(1 + I)/sqrt(pi)))
    assert erfi(z).rewrite('hyper') == 2*z*hyper([Rational(1, 2)], [Rational(3, 2)], z**2)/sqrt(pi)
    assert erfi(z).rewrite('meijerg') == z*meijerg([Rational(1, 2)], [], [0], [Rational(-1, 2)], -z**2)/sqrt(pi)
    assert erfi(z).rewrite('uppergamma') == (sqrt(-z**2)/z*(uppergamma(Rational(1, 2),
                                                                       -z**2)/sqrt(pi) - 1))
    assert erfi(z).rewrite('expint') == sqrt(-z**2)/z - z*expint(Rational(1, 2), -z**2)/sqrt(pi)

    assert erfi(x).as_real_imag() == \
        ((erfi(re(x) - I*re(x)*abs(im(x))/abs(re(x)))/2 +
          erfi(re(x) + I*re(x)*abs(im(x))/abs(re(x)))/2,
          I*(erfi(re(x) - I*re(x)*abs(im(x))/abs(re(x))) -
             erfi(re(x) + I*re(x)*abs(im(x))/abs(re(x)))) *
          re(x)*abs(im(x))/(2*im(x)*abs(re(x)))))
    assert erfi(x).as_real_imag(deep=False) == erfi(x).as_real_imag()
    assert erfi(w).as_real_imag() == (erfi(w), 0)
    assert erfi(w).as_real_imag(deep=False) == erfi(w).as_real_imag()
    assert erfi(I).as_real_imag() == (0, erf(1))

    pytest.raises(ArgumentIndexError, lambda: erfi(x).fdiff(2))

    assert erfi(x).taylor_term(3, x, *(2*x/sqrt(pi), 0)) == 2*x**3/3/sqrt(pi)

    assert erfi(x).limit(x, oo) == oo


def test_erfi_series():
    assert erfi(x).series(x, 0, 7) == 2*x/sqrt(pi) + \
        2*x**3/3/sqrt(pi) + x**5/5/sqrt(pi) + O(x**7)


def test_erfi_evalf():
    assert abs(erfi(Float(2.0)) - 18.5648024145756) < 1E-13  # XXX


def test_erf2():
    assert erf2(0, 0) == 0
    assert erf2(x, x) == 0
    assert erf2(nan, 0) == nan

    assert erf2(-oo, y) == erf(y) + 1
    assert erf2(+oo, y) == erf(y) - 1

    assert erf2(x, +oo) == +1 - erf(x)
    assert erf2(x, -oo) == -1 - erf(x)

    assert erf2(x, erf2inv(x, y)) == y

    assert erf2(-x, -y) == -erf2(x, y)
    assert erf2(-x, +y) == +erf(y) + erf(x)
    assert erf2(+x, -y) == -erf(y) - erf(x)
    assert erf2(x, y).rewrite('fresnels') == erf(y).rewrite(fresnels)-erf(x).rewrite(fresnels)
    assert erf2(x, y).rewrite('fresnelc') == erf(y).rewrite(fresnelc)-erf(x).rewrite(fresnelc)
    assert erf2(x, y).rewrite('hyper') == erf(y).rewrite(hyper)-erf(x).rewrite(hyper)
    assert erf2(x, y).rewrite('meijerg') == erf(y).rewrite(meijerg)-erf(x).rewrite(meijerg)
    assert erf2(x, y).rewrite('uppergamma') == erf(y).rewrite(uppergamma) - erf(x).rewrite(uppergamma)
    assert erf2(x, y).rewrite('expint') == erf(y).rewrite(expint)-erf(x).rewrite(expint)

    assert erf2(I, w).is_extended_real is False
    assert erf2(2*w, w).is_extended_real is True
    assert erf2(z, w).is_extended_real is None
    assert erf2(w, z).is_extended_real is None

    assert conjugate(erf2(x, y)) == erf2(conjugate(x), conjugate(y))

    assert erf2(x, y).rewrite('erf') == erf(y) - erf(x)
    assert erf2(x, y).rewrite('erfc') == erfc(x) - erfc(y)
    assert erf2(x, y).rewrite('erfi') == I*(erfi(I*x) - erfi(I*y))

    pytest.raises(ArgumentIndexError, lambda: erfi(x).fdiff(3))
    pytest.raises(ArgumentIndexError, lambda: erf2(x, y).fdiff(3))

    assert erf2(x, y).diff(x) == -2*exp(-x**2)/sqrt(pi)
    assert erf2(x, y).diff(y) == +2*exp(-y**2)/sqrt(pi)


def test_erfinv():
    assert erfinv(0) == 0
    assert erfinv(-1) == -oo
    assert erfinv(+1) == +oo
    assert erfinv(nan) == nan

    assert erfinv(erf(+w)) == w
    assert erfinv(erf(-w)) == -w

    assert erfinv(x).diff() == sqrt(pi)*exp(erfinv(x)**2)/2

    assert erfinv(z).rewrite('erfcinv') == erfcinv(1-z)

    pytest.raises(ArgumentIndexError, lambda: erfinv(x).fdiff(2))


def test_erfinv_evalf():
    assert abs(erfinv(Float(0.2)) - 0.179143454621292) < 1E-13


def test_erfcinv():
    assert erfcinv(1) == 0
    assert erfcinv(0) == oo
    assert erfcinv(nan) == nan

    assert erfcinv(x).diff() == -sqrt(pi)*exp(erfcinv(x)**2)/2

    assert erfcinv(z).rewrite('erfinv') == erfinv(1-z)

    pytest.raises(ArgumentIndexError, lambda: erfcinv(x).fdiff(2))


def test_erf2inv():
    assert erf2inv(0, 0) == 0
    assert erf2inv(0, 1) == oo
    assert erf2inv(1, 0) == 1
    assert erf2inv(0, y) == erfinv(y)
    assert erf2inv(oo, y) == erfcinv(-y)
    assert erf2inv(x, 0) == x
    assert erf2inv(x, oo) == erfinv(x)

    assert erf2inv(x, y).diff(x) == exp(-x**2 + erf2inv(x, y)**2)
    assert erf2inv(x, y).diff(y) == sqrt(pi)*exp(erf2inv(x, y)**2)/2

    pytest.raises(ArgumentIndexError, lambda: erf2inv(x, y).fdiff(3))


# NOTE we multiply by exp_polar(I*pi) and need this to be on the principal
# branch, hence take x in the lower half plane (d=0).


def mytn(expr1, expr2, expr3, x, d=0):
    subs = {}
    for a in expr1.free_symbols:
        if a != x:
            subs[a] = random_complex_number()
    return expr2 == expr3 and verify_numerically(expr1.subs(subs),
                                                 expr2.subs(subs), x, d=d)


def mytd(expr1, expr2, x):
    subs = {}
    for a in expr1.free_symbols:
        if a != x:
            subs[a] = random_complex_number()
    return expr1.diff(x) == expr2 and verify_derivative_numerically(expr1.subs(subs), x)


def tn_branch(func, s=None):
    def fn(x):
        if s is None:
            return func(x)
        return func(s, x)
    c = uniform(1, 5)
    expr = fn(c*exp_polar(I*pi)) - fn(c*exp_polar(-I*pi))
    eps = 1e-15
    expr2 = fn(-c + eps*I) - fn(-c - eps*I)
    return abs(expr - expr2).evalf(strict=False) < 1e-10


def test_ei():
    pos = Symbol('p', positive=True)
    neg = Symbol('n', negative=True)
    assert Ei(0) == -oo
    assert Ei(+oo) == oo
    assert Ei(-oo) == 0
    assert Ei(-pos) == Ei(polar_lift(-1)*pos) - I*pi
    assert Ei(neg) == Ei(polar_lift(neg)) - I*pi
    assert tn_branch(Ei)
    assert mytd(Ei(x), exp(x)/x, x)
    assert mytn(Ei(x), Ei(x).rewrite(uppergamma),
                -uppergamma(0, x*polar_lift(-1)) - I*pi, x)
    assert mytn(Ei(x), Ei(x).rewrite(expint),
                -expint(1, x*polar_lift(-1)) - I*pi, x)
    assert Ei(x).rewrite(expint).rewrite(Ei) == Ei(x)
    assert Ei(x*exp_polar(2*I*pi)) == Ei(x) + 2*I*pi
    assert Ei(x*exp_polar(-2*I*pi)) == Ei(x) - 2*I*pi

    assert mytn(Ei(x), Ei(x).rewrite(Shi), Chi(x) + Shi(x), x)
    assert mytn(Ei(x*polar_lift(I)), Ei(x*polar_lift(I)).rewrite(Si),
                Ci(x) + I*Si(x) + I*pi/2, x)

    assert Ei(log(x)).rewrite(li) == li(x)
    assert Ei(2*log(x)).rewrite(li) == li(x**2)

    assert Ei(x).series(x) == (EulerGamma + log(x) + x + x**2/4 +
                               x**3/18 + x**4/96 + x**5/600 + O(x**6))
    assert Ei(1 + x).series(x) == (Ei(1) + E*x + E*x**3/6 - E*x**4/12 +
                                   3*E*x**5/40 + O(x**6))

    pytest.raises(ArgumentIndexError, lambda: Ei(x).fdiff(2))

    assert (Ei(exp_polar(I*pi)).evalf() ==
            Float('-0.21938393439552029', dps=15) +
            I*Float('3.1415926535897931', dps=15))


def test_expint():
    assert mytn(expint(x, y), expint(x, y).rewrite(uppergamma),
                y**(x - 1)*uppergamma(1 - x, y), x)
    assert mytd(
        expint(x, y), -y**(x - 1)*meijerg([], [1, 1], [0, 0, 1 - x], [], y), x)
    assert mytd(expint(x, y), -expint(x - 1, y), y)
    assert mytn(expint(1, x), expint(1, x).rewrite(Ei),
                -Ei(x*polar_lift(-1)) + I*pi, x)

    assert expint(-4, x) == exp(-x)/x + 4*exp(-x)/x**2 + 12*exp(-x)/x**3 \
        + 24*exp(-x)/x**4 + 24*exp(-x)/x**5
    assert expint(-Rational(3, 2), x) == \
        exp(-x)/x + 3*exp(-x)/(2*x**2) - 3*sqrt(pi)*erf(sqrt(x))/(4*x**Rational(5, 2)) \
        + 3*sqrt(pi)/(4*x**Rational(5, 2))

    assert tn_branch(expint, 1)
    assert tn_branch(expint, 2)
    assert tn_branch(expint, 3)
    assert tn_branch(expint, 1.7)
    assert tn_branch(expint, pi)

    assert expint(y, x*exp_polar(2*I*pi)) == \
        x**(y - 1)*(exp(2*I*pi*y) - 1)*gamma(-y + 1) + expint(y, x)
    assert expint(y, x*exp_polar(-2*I*pi)) == \
        x**(y - 1)*(exp(-2*I*pi*y) - 1)*gamma(-y + 1) + expint(y, x)
    assert expint(2, x*exp_polar(2*I*pi)) == 2*I*pi*x + expint(2, x)
    assert expint(2, x*exp_polar(-2*I*pi)) == -2*I*pi*x + expint(2, x)
    assert (expint(n, x*exp_polar(2*I*pi)) ==
            expint(n, x*exp_polar(2*I*pi), evaluate=False))

    assert expint(1, x).rewrite(Ei).rewrite(expint) == expint(1, x)
    assert (expint(2, x, evaluate=False).rewrite(Shi) ==
            expint(2, x, evaluate=False))
    assert mytn(E1(x), E1(x).rewrite(Shi), Shi(x) - Chi(x), x)
    assert mytn(E1(polar_lift(I)*x), E1(polar_lift(I)*x).rewrite(Si),
                -Ci(x) + I*Si(x) - I*pi/2, x)

    assert mytn(expint(2, x), expint(2, x).rewrite(Ei).rewrite(expint),
                -x*E1(x) + exp(-x), x)
    assert mytn(expint(3, x), expint(3, x).rewrite(Ei).rewrite(expint),
                x**2*E1(x)/2 + (1 - x)*exp(-x)/2, x)

    assert expint(Rational(3, 2), z).series(z, n=6) == \
        2 + 2*z - z**2/3 + z**3/15 - z**4/84 + z**5/540 - \
        2*sqrt(pi)*sqrt(z) + O(z**6)

    assert E1(z).series(z) == -EulerGamma - log(z) + z - \
        z**2/4 + z**3/18 - z**4/96 + z**5/600 + O(z**6)

    assert expint(4, z).series(z) == Rational(1, 3) - z/2 + z**2/2 + \
        z**3*(log(z)/6 - Rational(11, 36) + EulerGamma/6) - z**4/24 + \
        z**5/240 + O(z**6)
    assert (expint(x, x).series(x, x0=1, n=2) ==
            expint(1, 1) + (x - 1)*(-meijerg(((), (1, 1)),
                                             ((0, 0, 0), ()), 1) - 1/E) +
            O((x - 1)**2, (x, 1)))

    pytest.raises(ArgumentIndexError, lambda: expint(x, y).fdiff(3))


def test__eis():
    assert _eis(z).diff(z) == -_eis(z) + 1/z
    pytest.raises(ArgumentIndexError, lambda: _eis(x).fdiff(2))

    assert _eis(1/z).series(z) == \
        z + z**2 + 2*z**3 + 6*z**4 + 24*z**5 + O(z**6)

    assert Ei(z).rewrite('tractable') == exp(z)*_eis(z)
    assert li(z).rewrite('tractable') == z*_eis(log(z))

    assert _eis(z).rewrite('intractable') == exp(-z)*Ei(z)

    assert expand(li(z).rewrite('tractable').diff(z).rewrite('intractable')) \
        == li(z).diff(z)

    assert expand(Ei(z).rewrite('tractable').diff(z).rewrite('intractable')) \
        == Ei(z).diff(z)

    assert _eis(z).series(z, n=2) == EulerGamma + log(z) + z*(-log(z) -
                                                              EulerGamma + 1) + z**2*(log(z)/2 - Rational(3, 4) + EulerGamma/2) + O(z**2)

    l = Limit(Ei(y/x)/exp(y/x), x, 0)
    assert l.doit() == l  # cover _eis._eval_aseries


def tn_arg(func):
    def test(arg, e1, e2):
        v = uniform(1, 5)
        v1 = func(arg*x).subs({x: v}).evalf(strict=False)
        v2 = func(e1*v + e2*1e-15).evalf(strict=False)
        return abs(v1 - v2).evalf(strict=False) < 1e-10
    return test(exp_polar(I*pi/2), I, 1) and \
        test(exp_polar(-I*pi/2), -I, 1) and \
        test(exp_polar(I*pi), -1, I) and \
        test(exp_polar(-I*pi), -1, -I)


def test_li():
    z = Symbol('z')
    zr = Symbol('z', extended_real=True)
    zp = Symbol('z', positive=True)
    zn = Symbol('z', negative=True)

    assert li(0) == 0
    assert li(1) == -oo
    assert li(oo) == oo

    assert isinstance(li(z), li)

    assert diff(li(z), z) == 1/log(z)
    pytest.raises(ArgumentIndexError, lambda: li(z).fdiff(2))

    assert conjugate(li(z)) == li(conjugate(z))
    assert conjugate(li(-zr)) == li(-zr)
    assert conjugate(li(-zp)) == conjugate(li(-zp))
    assert conjugate(li(zn)) == conjugate(li(zn))

    assert li(z).rewrite(Li) == Li(z) + li(2)
    assert li(z).rewrite(Ei) == Ei(log(z))
    assert li(z).rewrite(uppergamma) == (-log(1/log(z))/2 - log(-log(z)) +
                                         log(log(z))/2 - expint(1, -log(z)))
    assert li(z).rewrite(Si) == (-log(I*log(z)) - log(1/log(z))/2 +
                                 log(log(z))/2 + Ci(I*log(z)) + Shi(log(z)))
    assert li(z).rewrite(Ci) == (-log(I*log(z)) - log(1/log(z))/2 +
                                 log(log(z))/2 + Ci(I*log(z)) + Shi(log(z)))
    assert li(z).rewrite(Shi) == (-log(1/log(z))/2 + log(log(z))/2 +
                                  Chi(log(z)) - Shi(log(z)))
    assert li(z).rewrite(Chi) == (-log(1/log(z))/2 + log(log(z))/2 +
                                  Chi(log(z)) - Shi(log(z)))
    assert li(z).rewrite(hyper) == (log(z)*hyper((1, 1), (2, 2), log(z)) -
                                    log(1/log(z))/2 + log(log(z))/2 + EulerGamma)
    assert li(z).rewrite(meijerg) == (-log(1/log(z))/2 - log(-log(z)) + log(log(z))/2 -
                                      meijerg(((), (1,)), ((0, 0), ()), -log(z)))


def test_Li():
    assert Li(2) == 0
    assert Li(oo) == oo

    assert isinstance(Li(z), Li)

    assert diff(Li(z), z) == 1/log(z)
    pytest.raises(ArgumentIndexError, lambda: Li(z).fdiff(2))

    assert Li(z).rewrite(li) == li(z) - li(2)

    assert Li(4).evalf(30) == Float('1.92242131492155809316615998937961', dps=30)


def test_si():
    assert Si(I*x) == I*Shi(x)
    assert Shi(I*x) == I*Si(x)
    assert Si(-I*x) == -I*Shi(x)
    assert Shi(-I*x) == -I*Si(x)
    assert Si(-x) == -Si(x)
    assert Shi(-x) == -Shi(x)
    assert Si(exp_polar(2*pi*I)*x) == Si(x)
    assert Si(exp_polar(-2*pi*I)*x) == Si(x)
    assert Shi(exp_polar(2*pi*I)*x) == Shi(x)
    assert Shi(exp_polar(-2*pi*I)*x) == Shi(x)

    assert Si(oo) == pi/2
    assert Si(-oo) == -pi/2
    assert Shi(oo) == oo
    assert Shi(-oo) == -oo

    assert mytd(Si(x), sin(x)/x, x)
    assert mytd(Shi(x), sinh(x)/x, x)

    assert mytn(Si(x), Si(x).rewrite(Ei),
                -I*(-Ei(x*exp_polar(-I*pi/2))/2
                    + Ei(x*exp_polar(I*pi/2))/2 - I*pi) + pi/2, x)
    assert mytn(Si(x), Si(x).rewrite(expint),
                -I*(-expint(1, x*exp_polar(-I*pi/2))/2 +
                    expint(1, x*exp_polar(I*pi/2))/2) + pi/2, x)
    assert mytn(Shi(x), Shi(x).rewrite(Ei),
                Ei(x)/2 - Ei(x*exp_polar(I*pi))/2 + I*pi/2, x)
    assert mytn(Shi(x), Shi(x).rewrite(expint),
                expint(1, x)/2 - expint(1, x*exp_polar(I*pi))/2 - I*pi/2, x)

    assert tn_arg(Si)
    assert tn_arg(Shi)

    assert Si(x).series(x, n=9) == \
        x - x**3/18 + x**5/600 - x**7/35280 + O(x**9)
    assert Shi(x).series(x, n=9) == \
        x + x**3/18 + x**5/600 + x**7/35280 + O(x**9)
    assert Si(sin(x)).series(x, n=7) == x - 2*x**3/9 + 17*x**5/450 + O(x**7)
    assert Si(x).series(x, 1, n=3) == \
        Si(1) + (x - 1)*sin(1) + (x - 1)**2*(-sin(1)/2 + cos(1)/2) + O((x - 1)**3, (x, 1))

    pytest.raises(ArgumentIndexError, lambda: Si(z).fdiff(2))


def test_ci():
    m1 = exp_polar(I*pi)
    m1_ = exp_polar(-I*pi)
    pI = exp_polar(I*pi/2)
    mI = exp_polar(-I*pi/2)

    assert Ci(m1*x) == Ci(x) + I*pi
    assert Ci(m1_*x) == Ci(x) - I*pi
    assert Ci(pI*x) == Chi(x) + I*pi/2
    assert Ci(mI*x) == Chi(x) - I*pi/2
    assert Chi(m1*x) == Chi(x) + I*pi
    assert Chi(m1_*x) == Chi(x) - I*pi
    assert Chi(pI*x) == Ci(x) + I*pi/2
    assert Chi(mI*x) == Ci(x) - I*pi/2
    assert Ci(exp_polar(2*I*pi)*x) == Ci(x) + 2*I*pi
    assert Chi(exp_polar(-2*I*pi)*x) == Chi(x) - 2*I*pi
    assert Chi(exp_polar(2*I*pi)*x) == Chi(x) + 2*I*pi
    assert Ci(exp_polar(-2*I*pi)*x) == Ci(x) - 2*I*pi

    assert Ci(oo) == 0
    assert Ci(-oo) == I*pi
    assert Chi(oo) == oo
    assert Chi(-oo) == oo

    assert mytd(Ci(x), cos(x)/x, x)
    assert mytd(Chi(x), cosh(x)/x, x)

    assert mytn(Ci(x), Ci(x).rewrite(Ei),
                Ei(x*exp_polar(-I*pi/2))/2 + Ei(x*exp_polar(I*pi/2))/2, x)
    assert mytn(Chi(x), Chi(x).rewrite(Ei),
                Ei(x)/2 + Ei(x*exp_polar(I*pi))/2 - I*pi/2, x)

    assert tn_arg(Ci)
    assert tn_arg(Chi)

    assert Ci(x).series(x) == \
        EulerGamma + log(x) - x**2/4 + x**4/96 + O(x**6)
    assert Chi(x).series(x) == \
        EulerGamma + log(x) + x**2/4 + x**4/96 + O(x**6)
    assert limit(log(x) - Ci(2*x), x, 0) == -log(2) - EulerGamma


def test_fresnel():
    assert fresnels(0) == 0
    assert fresnels(+oo) == Rational(+1, 2)
    assert fresnels(-oo) == Rational(-1, 2)

    assert fresnels(z) == fresnels(z)
    assert fresnels(-z) == -fresnels(z)
    assert fresnels(I*z) == -I*fresnels(z)
    assert fresnels(-I*z) == I*fresnels(z)

    assert conjugate(fresnels(z)) == fresnels(conjugate(z))

    assert fresnels(z).diff(z) == sin(pi*z**2/2)

    assert fresnels(z).rewrite(erf) == (1 + I)/4 * (
        erf((1 + I)/2*sqrt(pi)*z) - I*erf((1 - I)/2*sqrt(pi)*z))

    assert fresnels(z).rewrite(hyper) == \
        pi*z**3/6 * hyper([Rational(3, 4)], [Rational(3, 2), Rational(7, 4)], -pi**2*z**4/16)

    assert fresnels(z).series(z, n=15) == \
        pi*z**3/6 - pi**3*z**7/336 + pi**5*z**11/42240 + O(z**15)

    assert fresnels(y/z).limit(z, 0) == fresnels(oo*sign(y))

    assert fresnels(x).taylor_term(-1, z) == 0
    assert fresnels(x).taylor_term(1, z, *(pi*z**3/6,)) == -pi**3*z**7/336
    assert fresnels(x).taylor_term(1, z) == -pi**3*z**7/336

    assert fresnels(w).is_extended_real is True
    assert fresnels(z).is_extended_real is None

    assert fresnels(z).as_real_imag() == \
        ((fresnels(re(z) - I*re(z)*abs(im(z))/abs(re(z)))/2 +
          fresnels(re(z) + I*re(z)*abs(im(z))/abs(re(z)))/2,
          I*(fresnels(re(z) - I*re(z)*abs(im(z))/abs(re(z))) -
             fresnels(re(z) + I*re(z)*abs(im(z))/abs(re(z)))) *
          re(z)*abs(im(z))/(2*im(z)*abs(re(z)))))
    assert fresnels(z).as_real_imag(deep=False) == fresnels(z).as_real_imag()
    assert fresnels(w).as_real_imag() == (fresnels(w), 0)
    assert fresnels(w).as_real_imag(deep=False) == fresnels(w).as_real_imag()
    assert (fresnels(I, evaluate=False).as_real_imag() ==
            (0, -erf(sqrt(pi)/2 + I*sqrt(pi)/2)/4 +
             I*(-erf(sqrt(pi)/2 + I*sqrt(pi)/2) + erf(sqrt(pi)/2 -
                I*sqrt(pi)/2))/4 - erf(sqrt(pi)/2 - I*sqrt(pi)/2)/4))

    assert fresnels(2 + 3*I).as_real_imag() == (
        fresnels(2 + 3*I)/2 + fresnels(2 - 3*I)/2,
        I*(fresnels(2 - 3*I) - fresnels(2 + 3*I))/2
    )

    assert expand_func(integrate(fresnels(z), z)) == \
        z*fresnels(z) + cos(pi*z**2/2)/pi

    assert fresnels(z).rewrite(meijerg) == sqrt(2)*pi*z**Rational(9, 4) * \
        meijerg(((), (1,)), ((Rational(3, 4),),
                             (Rational(1, 4), 0)), -pi**2*z**4/16)/(2*(-z)**Rational(3, 4)*(z**2)**Rational(3, 4))

    assert fresnelc(0) == 0
    assert fresnelc(+oo) == Rational(+1, 2)
    assert fresnelc(-oo) == Rational(-1, 2)

    assert fresnelc(z) == fresnelc(z)
    assert fresnelc(-z) == -fresnelc(z)
    assert fresnelc(I*z) == I*fresnelc(z)
    assert fresnelc(-I*z) == -I*fresnelc(z)

    assert conjugate(fresnelc(z)) == fresnelc(conjugate(z))

    assert fresnelc(z).diff(z) == cos(pi*z**2/2)
    pytest.raises(ArgumentIndexError, lambda: fresnels(z).fdiff(2))
    pytest.raises(ArgumentIndexError, lambda: fresnelc(z).fdiff(2))

    assert fresnelc(z).rewrite(erf) == (1 - I)/4 * (
        erf((1 + I)/2*sqrt(pi)*z) + I*erf((1 - I)/2*sqrt(pi)*z))

    assert fresnelc(z).rewrite(hyper) == \
        z * hyper([Rational(1, 4)], [Rational(1, 2), Rational(5, 4)], -pi**2*z**4/16)

    assert fresnelc(x).taylor_term(-1, z) == 0
    assert fresnelc(x).taylor_term(1, z, *(z,)) == -pi**2*z**5/40
    assert fresnelc(x).taylor_term(1, z) == -pi**2*z**5/40

    assert fresnelc(z).series(z, n=15) == \
        z - pi**2*z**5/40 + pi**4*z**9/3456 - pi**6*z**13/599040 + O(z**15)

    assert fresnelc(y/z).limit(z, 0) == fresnelc(oo*sign(y))

    # issue sympy/sympy#6510
    assert fresnels(z).series(z, oo) == \
        (-1/(pi**2*z**3) + O(z**(-6), (z, oo)))*sin(pi*z**2/2) + \
        (3/(pi**3*z**5) - 1/(pi*z) + O(z**(-6), (z, oo)))*cos(pi*z**2/2) + Rational(1, 2)
    assert fresnelc(z).series(z, oo) == \
        (-1/(pi**2*z**3) + O(z**(-6), (z, oo)))*cos(pi*z**2/2) + \
        (-3/(pi**3*z**5) + 1/(pi*z) + O(z**(-6), (z, oo)))*sin(pi*z**2/2) + Rational(1, 2)
    assert fresnels(1/z).series(z) == \
        (-z**3/pi**2 + O(z**6))*sin(pi/(2*z**2)) + (-z/pi + 3*z**5/pi**3 +
                                                    O(z**6))*cos(pi/(2*z**2)) + Rational(1, 2)
    assert fresnelc(1/z).series(z) == \
        (-z**3/pi**2 + O(z**6))*cos(pi/(2*z**2)) + (z/pi - 3*z**5/pi**3 +
                                                    O(z**6))*sin(pi/(2*z**2)) + Rational(1, 2)

    assert fresnelc(w).is_extended_real is True

    assert fresnelc(z).as_real_imag() == \
        ((fresnelc(re(z) - I*re(z)*abs(im(z))/abs(re(z)))/2 +
          fresnelc(re(z) + I*re(z)*abs(im(z))/abs(re(z)))/2,
          I*(fresnelc(re(z) - I*re(z)*abs(im(z))/abs(re(z))) -
             fresnelc(re(z) + I*re(z)*abs(im(z))/abs(re(z)))) *
          re(z)*abs(im(z))/(2*im(z)*abs(re(z)))))

    assert fresnelc(2 + 3*I).as_real_imag() == (
        fresnelc(2 - 3*I)/2 + fresnelc(2 + 3*I)/2,
        I*(fresnelc(2 - 3*I) - fresnelc(2 + 3*I))/2
    )

    assert expand_func(integrate(fresnelc(z), z)) == \
        z*fresnelc(z) - sin(pi*z**2/2)/pi

    assert fresnelc(z).rewrite(meijerg) == sqrt(2)*pi*z**Rational(3, 4) * \
        meijerg(((), (1,)), ((Rational(1, 4),),
                             (Rational(3, 4), 0)), -pi**2*z**4/16)/(2*root(-z, 4)*root(z**2, 4))

    verify_numerically(re(fresnels(z)), fresnels(z).as_real_imag()[0], z)
    verify_numerically(im(fresnels(z)), fresnels(z).as_real_imag()[1], z)
    verify_numerically(fresnels(z), fresnels(z).rewrite(hyper), z)
    verify_numerically(fresnels(z), fresnels(z).rewrite(meijerg), z)

    verify_numerically(re(fresnelc(z)), fresnelc(z).as_real_imag()[0], z)
    verify_numerically(im(fresnelc(z)), fresnelc(z).as_real_imag()[1], z)
    verify_numerically(fresnelc(z), fresnelc(z).rewrite(hyper), z)
    verify_numerically(fresnelc(z), fresnelc(z).rewrite(meijerg), z)
