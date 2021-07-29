import pytest

from diofant import (I, O, Rational, Symbol, atanh, conjugate, elliptic_e,
                     elliptic_f, elliptic_k, elliptic_pi, gamma, hyper,
                     meijerg, oo, pi, sin, sqrt, tan, zoo)
from diofant.abc import m, n, z
from diofant.core.function import ArgumentIndexError
from diofant.utilities.randtest import random_complex_number as randcplx
from diofant.utilities.randtest import verify_derivative_numerically as td
from diofant.utilities.randtest import verify_numerically as tn


__all__ = ()

i = Symbol('i', integer=True)
j = Symbol('k', integer=True, positive=True)


def test_elliptic_k():
    assert elliptic_k(0) == pi/2
    assert elliptic_k(Rational(1, 2)) == 8*pi**Rational(3, 2)/gamma(-Rational(1, 4))**2
    assert elliptic_k(1) == zoo
    assert elliptic_k(-1) == gamma(Rational(1, 4))**2/(4*sqrt(2*pi))
    assert elliptic_k(oo) == 0
    assert elliptic_k(-oo) == 0
    assert elliptic_k(I*oo) == 0
    assert elliptic_k(-I*oo) == 0
    assert elliptic_k(zoo) == 0

    assert elliptic_k(z).diff(z) == (elliptic_e(z) - (1 - z)*elliptic_k(z))/(2*z*(1 - z))
    assert td(elliptic_k(z), z)
    pytest.raises(ArgumentIndexError, lambda: elliptic_k(z).fdiff(2))

    zi = Symbol('z', extended_real=False)
    assert elliptic_k(zi).conjugate() == elliptic_k(zi.conjugate())
    zr = Symbol('z', extended_real=True, negative=True)
    assert elliptic_k(zr).conjugate() == elliptic_k(zr)
    assert elliptic_k(z).conjugate() == conjugate(elliptic_k(z), evaluate=False)

    assert elliptic_k(z).rewrite(hyper) == \
        (pi/2)*hyper((Rational(1, 2), Rational(1, 2)), (1,), z)
    assert tn(elliptic_k(z), (pi/2)*hyper((Rational(1, 2), Rational(1, 2)), (1,), z))
    assert elliptic_k(z).rewrite(meijerg) == \
        meijerg(((Rational(1, 2), Rational(1, 2)), []), ((0,), (0,)), -z)/2
    assert tn(elliptic_k(z), meijerg(((Rational(1, 2), Rational(1, 2)), []), ((0,), (0,)), -z)/2)

    assert elliptic_k(z).series(z) == pi/2 + pi*z/8 + 9*pi*z**2/128 + \
        25*pi*z**3/512 + 1225*pi*z**4/32768 + 3969*pi*z**5/131072 + O(z**6)


def test_elliptic_f():
    assert elliptic_f(z, 0) == z
    assert elliptic_f(0, m) == 0
    assert elliptic_f(pi*i/2, m) == i*elliptic_k(m)
    assert elliptic_f(z, oo) == 0
    assert elliptic_f(z, -oo) == 0

    assert elliptic_f(-z, m) == -elliptic_f(z, m)

    assert elliptic_f(z, m).diff(z) == 1/sqrt(1 - m*sin(z)**2)
    assert elliptic_f(z, m).diff(m) == elliptic_e(z, m)/(2*m*(1 - m)) - elliptic_f(z, m)/(2*m) - \
        sin(2*z)/(4*(1 - m)*sqrt(1 - m*sin(z)**2))
    r = randcplx()
    assert td(elliptic_f(z, r), z)
    assert td(elliptic_f(r, m), m)
    pytest.raises(ArgumentIndexError, lambda: elliptic_f(z, m).fdiff(3))

    mi = Symbol('m', extended_real=False)
    assert elliptic_f(z, mi).conjugate() == elliptic_f(z.conjugate(), mi.conjugate())
    mr = Symbol('m', extended_real=True, negative=True)
    assert elliptic_f(z, mr).conjugate() == elliptic_f(z.conjugate(), mr)
    assert elliptic_f(z, m).conjugate() == conjugate(elliptic_f(z, m), evaluate=False)

    assert elliptic_f(z, m).series(z) == \
        z + z**5*(3*m**2/40 - m/30) + m*z**3/6 + O(z**6)


def test_elliptic_e():
    assert elliptic_e(z, 0) == z
    assert elliptic_e(0, m) == 0
    assert elliptic_e(i*pi/2, m) == i*elliptic_e(m)
    assert elliptic_e(z, oo) == zoo
    assert elliptic_e(z, -oo) == zoo
    assert elliptic_e(0) == pi/2
    assert elliptic_e(1) == 1
    assert elliptic_e(oo) == I*oo
    assert elliptic_e(-oo) == oo
    assert elliptic_e(zoo) == zoo

    assert elliptic_e(-z, m) == -elliptic_e(z, m)

    assert elliptic_e(z, m).diff(z) == sqrt(1 - m*sin(z)**2)
    assert elliptic_e(z, m).diff(m) == (elliptic_e(z, m) - elliptic_f(z, m))/(2*m)
    assert elliptic_e(z).diff(z) == (elliptic_e(z) - elliptic_k(z))/(2*z)
    r = randcplx()
    assert td(elliptic_e(r, m), m)
    assert td(elliptic_e(z, r), z)
    assert td(elliptic_e(z), z)
    pytest.raises(ArgumentIndexError, lambda: elliptic_e(z, m).fdiff(3))
    pytest.raises(ArgumentIndexError, lambda: elliptic_e(z).fdiff(2))

    mi = Symbol('m', extended_real=False)
    assert elliptic_e(z, mi).conjugate() == elliptic_e(z.conjugate(), mi.conjugate())
    assert elliptic_e(mi).conjugate() == elliptic_e(mi.conjugate())
    mr = Symbol('m', extended_real=True, negative=True)
    assert elliptic_e(z, mr).conjugate() == elliptic_e(z.conjugate(), mr)
    assert elliptic_e(mr).conjugate() == elliptic_e(mr)
    assert elliptic_e(z, m).conjugate() == conjugate(elliptic_e(z, m))
    assert elliptic_e(z).conjugate() == conjugate(elliptic_e(z))

    assert elliptic_e(z).rewrite(hyper) == (pi/2)*hyper((Rational(-1, 2), Rational(1, 2)), (1,), z)
    assert elliptic_e(z, m).rewrite(hyper) == elliptic_e(z, m)
    assert tn(elliptic_e(z), (pi/2)*hyper((Rational(-1, 2), Rational(1, 2)), (1,), z))
    assert elliptic_e(z).rewrite(meijerg) == \
        -meijerg(((Rational(1, 2), Rational(3, 2)), []), ((0,), (0,)), -z)/4
    assert elliptic_e(z, m).rewrite(meijerg) == elliptic_e(z, m)
    assert tn(elliptic_e(z), -meijerg(((Rational(1, 2), Rational(3, 2)), []), ((0,), (0,)), -z)/4)

    assert elliptic_e(z, m).series(z) == \
        z + z**5*(-m**2/40 + m/30) - m*z**3/6 + O(z**6)
    assert elliptic_e(z).series(z) == pi/2 - pi*z/8 - 3*pi*z**2/128 - \
        5*pi*z**3/512 - 175*pi*z**4/32768 - 441*pi*z**5/131072 + O(z**6)


def test_elliptic_pi():
    assert elliptic_pi(0, z, m) == elliptic_f(z, m)
    assert elliptic_pi(1, z, m) == elliptic_f(z, m) + \
        (sqrt(1 - m*sin(z)**2)*tan(z) - elliptic_e(z, m))/(1 - m)
    assert elliptic_pi(n, i*pi/2, m) == i*elliptic_pi(n, m)
    assert elliptic_pi(n, z, 0) == atanh(sqrt(n - 1)*tan(z))/sqrt(n - 1)
    assert elliptic_pi(n, z, n) == elliptic_f(z, n) - elliptic_pi(1, z, n) + tan(z)/sqrt(1 - n*sin(z)**2)
    assert elliptic_pi(oo, z, m) == 0
    assert elliptic_pi(-oo, z, m) == 0
    assert elliptic_pi(n, z, oo) == 0
    assert elliptic_pi(n, z, -oo) == 0
    assert elliptic_pi(0, m) == elliptic_k(m)
    assert elliptic_pi(1, m) == zoo
    assert elliptic_pi(n, 0) == pi/(2*sqrt(1 - n))
    assert elliptic_pi(2, 1) == -oo
    assert elliptic_pi(-1, 1) == oo
    assert elliptic_pi(n, n) == elliptic_e(n)/(1 - n)
    assert elliptic_pi(oo, m) == 0
    assert elliptic_pi(n, oo) == 0

    assert elliptic_pi(n, -z, m) == -elliptic_pi(n, z, m)

    ni, mi = Symbol('n', extended_real=False), Symbol('m', extended_real=False)
    assert elliptic_pi(ni, z, mi).conjugate() == \
        elliptic_pi(ni.conjugate(), z.conjugate(), mi.conjugate())
    nr, mr = Symbol('n', extended_real=True, negative=True), \
        Symbol('m', extended_real=True, negative=True)
    assert elliptic_pi(nr, z, mr).conjugate() == elliptic_pi(nr, z.conjugate(), mr)
    assert elliptic_pi(n, m).conjugate() == elliptic_pi(n.conjugate(), m.conjugate())
    assert elliptic_pi(n, z, m).conjugate() == conjugate(elliptic_pi(n, z, m))

    assert elliptic_pi(n, z, m).diff(n) == (elliptic_e(z, m) + (m - n)*elliptic_f(z, m)/n +
                                            (n**2 - m)*elliptic_pi(n, z, m)/n - n*sqrt(1 -
                                                                                       m*sin(z)**2)*sin(2*z)/(2*(1 - n*sin(z)**2)))/(2*(m - n)*(n - 1))
    assert elliptic_pi(n, z, m).diff(z) == 1/(sqrt(1 - m*sin(z)**2)*(1 - n*sin(z)**2))
    assert elliptic_pi(n, z, m).diff(m) == (elliptic_e(z, m)/(m - 1) + elliptic_pi(n, z, m) -
                                            m*sin(2*z)/(2*(m - 1)*sqrt(1 - m*sin(z)**2)))/(2*(n - m))
    assert elliptic_pi(n, m).diff(n) == (elliptic_e(m) + (m - n)*elliptic_k(m)/n +
                                         (n**2 - m)*elliptic_pi(n, m)/n)/(2*(m - n)*(n - 1))
    assert elliptic_pi(n, m).diff(m) == (elliptic_e(m)/(m - 1) + elliptic_pi(n, m))/(2*(n - m))

    # workaround fredrik-johansson/mpmath#571, suggested by Kalevi Suominen
    # in https://github.com/sympy/sympy/issues/20933#issuecomment-779077562
    bounds = {'a': -0.9, 'b': -0.9, 'c': 0.9, 'd': 0.9}
    rx, ry = randcplx(**bounds), randcplx(**bounds)
    assert td(elliptic_pi(n, rx, ry), n, **bounds)
    assert td(elliptic_pi(rx, z, ry), z, **bounds)
    assert td(elliptic_pi(rx, ry, m), m, **bounds)

    pytest.raises(ArgumentIndexError, lambda: elliptic_pi(n, z, m).fdiff(4))
    pytest.raises(ArgumentIndexError, lambda: elliptic_pi(n, m).fdiff(3))

    assert elliptic_pi(n, z, m).series(z) == z + z**3*(m/6 + n/3) + \
        z**5*(3*m**2/40 + m*n/10 - m/30 + n**2/5 - n/15) + O(z**6)
