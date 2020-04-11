import pytest

from diofant import (I, Symbol, Ynm, Ynm_c, Znm, assoc_legendre, conjugate,
                     cos, cot, diff, exp, factorial, pi, sin, sqrt)
from diofant.abc import m, n
from diofant.core.function import ArgumentIndexError


__all__ = ()


def test_Ynm():
    # https://en.wikipedia.org/wiki/Spherical_harmonics
    th, ph = Symbol('theta', extended_real=True), Symbol('phi', extended_real=True)

    assert Ynm(0, 0, th, ph).expand(func=True) == 1/(2*sqrt(pi))
    assert Ynm(1, -1, th, ph) == -exp(-2*I*ph)*Ynm(1, 1, th, ph)
    assert Ynm(1, -1, th, ph).expand(func=True) == sqrt(6)*sin(th)*exp(-I*ph)/(4*sqrt(pi))
    assert Ynm(1, -1, th, ph).expand(func=True) == sqrt(6)*sin(th)*exp(-I*ph)/(4*sqrt(pi))
    assert Ynm(1, 0, th, ph).expand(func=True) == sqrt(3)*cos(th)/(2*sqrt(pi))
    assert Ynm(1, 1, th, ph).expand(func=True) == -sqrt(6)*sin(th)*exp(I*ph)/(4*sqrt(pi))
    assert Ynm(2, 0, th, ph).expand(func=True) == 3*sqrt(5)*cos(th)**2/(4*sqrt(pi)) - sqrt(5)/(4*sqrt(pi))
    assert Ynm(2, 1, th, ph).expand(func=True) == -sqrt(30)*sin(th)*exp(I*ph)*cos(th)/(4*sqrt(pi))
    assert Ynm(2, -2, th, ph).expand(func=True) == (-sqrt(30)*exp(-2*I*ph)*cos(th)**2/(8*sqrt(pi))
                                                    + sqrt(30)*exp(-2*I*ph)/(8*sqrt(pi)))
    assert Ynm(2, 2, th, ph).expand(func=True) == (-sqrt(30)*exp(2*I*ph)*cos(th)**2/(8*sqrt(pi))
                                                   + sqrt(30)*exp(2*I*ph)/(8*sqrt(pi)))

    assert diff(Ynm(n, m, th, ph), th) == (m*cot(th)*Ynm(n, m, th, ph)
                                           + sqrt((-m + n)*(m + n + 1))*exp(-I*ph)*Ynm(n, m + 1, th, ph))
    assert diff(Ynm(n, m, th, ph), ph) == I*m*Ynm(n, m, th, ph)
    pytest.raises(ArgumentIndexError, lambda: Ynm(n, m, th, ph).fdiff(1))

    assert conjugate(Ynm(n, m, th, ph)) == (-1)**(2*m)*exp(-2*I*m*ph)*Ynm(n, m, th, ph)

    assert Ynm(n, m, -th, ph) == Ynm(n, m, th, ph)
    assert Ynm(n, m, th, -ph) == exp(-2*I*m*ph)*Ynm(n, m, th, ph)
    assert Ynm(n, -m, th, ph) == (-1)**m*exp(-2*I*m*ph)*Ynm(n, m, th, ph)

    assert (Ynm(n, m, th, ph).rewrite(sin) ==
            Ynm(n, m, th, ph).rewrite(cos) ==
            exp(I*m*ph)*sqrt((2*n + 1)*factorial(-m + n)/factorial(m + n)) *
            assoc_legendre(n, m, cos(th))/(2*sqrt(pi)))
    assert (Ynm(n, m, th, ph).as_real_imag() ==
            (sqrt((2*n + 1)*factorial(-m + n)/factorial(m + n))*cos(m*ph) *
             assoc_legendre(n, m, cos(th))/(2*sqrt(pi)),
             sqrt((2*n + 1)*factorial(-m + n)/factorial(m + n))*sin(m*ph) *
             assoc_legendre(n, m, cos(th))/(2*sqrt(pi))))


def test_Ynm_c():
    th, ph = Symbol('theta', extended_real=True), Symbol('phi', extended_real=True)

    assert Ynm_c(n, m, th, ph) == (-1)**(2*m)*exp(-2*I*m*ph)*Ynm(n, m, th, ph)


def test_Znm():
    # https://en.wikipedia.org/wiki/Solid_harmonics#List_of_lowest_functions
    th, ph = Symbol('theta', extended_real=True), Symbol('phi', extended_real=True)

    assert Znm(0, 0, th, ph) == Ynm(0, 0, th, ph)
    assert Znm(1, -1, th, ph) == (-sqrt(2)*I*(Ynm(1, 1, th, ph)
                                              - exp(-2*I*ph)*Ynm(1, 1, th, ph))/2)
    assert Znm(1, 0, th, ph) == Ynm(1, 0, th, ph)
    assert Znm(1, 1, th, ph) == (sqrt(2)*(Ynm(1, 1, th, ph)
                                          + exp(-2*I*ph)*Ynm(1, 1, th, ph))/2)
    assert Znm(0, 0, th, ph).expand(func=True) == 1/(2*sqrt(pi))
    assert Znm(1, -1, th, ph).expand(func=True) == (sqrt(3)*I*sin(th)*exp(I*ph)/(4*sqrt(pi))
                                                    - sqrt(3)*I*sin(th)*exp(-I*ph)/(4*sqrt(pi)))
    assert Znm(1, 0, th, ph).expand(func=True) == sqrt(3)*cos(th)/(2*sqrt(pi))
    assert Znm(1, 1, th, ph).expand(func=True) == (-sqrt(3)*sin(th)*exp(I*ph)/(4*sqrt(pi))
                                                   - sqrt(3)*sin(th)*exp(-I*ph)/(4*sqrt(pi)))
    assert Znm(2, -1, th, ph).expand(func=True) == (sqrt(15)*I*sin(th)*exp(I*ph)*cos(th)/(4*sqrt(pi))
                                                    - sqrt(15)*I*sin(th)*exp(-I*ph)*cos(th)/(4*sqrt(pi)))
    assert Znm(2, 0, th, ph).expand(func=True) == 3*sqrt(5)*cos(th)**2/(4*sqrt(pi)) - sqrt(5)/(4*sqrt(pi))
    assert Znm(2, 1, th, ph).expand(func=True) == (-sqrt(15)*sin(th)*exp(I*ph)*cos(th)/(4*sqrt(pi))
                                                   - sqrt(15)*sin(th)*exp(-I*ph)*cos(th)/(4*sqrt(pi)))
