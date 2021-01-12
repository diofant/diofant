"""Tests for functions for generating interesting polynomials."""

import random

import pytest

from diofant import (ZZ, cyclotomic_poly, interpolating_poly, random_poly,
                     ring, swinnerton_dyer_poly, symbols, symmetric_poly)
from diofant.abc import x, y, z


__all__ = ()


def test_swinnerton_dyer_poly():
    pytest.raises(ValueError, lambda: swinnerton_dyer_poly(0, x))

    assert swinnerton_dyer_poly(1, x, polys=True) == (x**2 - 2).as_poly()
    assert swinnerton_dyer_poly(1, polys=True) == (x**2 - 2).as_poly()

    assert swinnerton_dyer_poly(1, x) == x**2 - 2
    assert swinnerton_dyer_poly(2, x) == x**4 - 10*x**2 + 1
    assert swinnerton_dyer_poly(3, x) == (x**8 - 40*x**6 +
                                          352*x**4 - 960*x**2 + 576)
    assert swinnerton_dyer_poly(4, x) == (x**16 - 136*x**14 + 6476*x**12 -
                                          141912*x**10 + 1513334*x**8 -
                                          7453176*x**6 + 13950764*x**4 -
                                          5596840*x**2 + 46225)


def test_cyclotomic_poly():
    pytest.raises(ValueError, lambda: cyclotomic_poly(0, x))

    assert cyclotomic_poly(1, x, polys=True) == (x - 1).as_poly()
    assert cyclotomic_poly(1, polys=True) == (x - 1).as_poly()

    assert cyclotomic_poly(1, x) == x - 1
    assert cyclotomic_poly(2, x) == x + 1
    assert cyclotomic_poly(3, x) == x**2 + x + 1
    assert cyclotomic_poly(4, x) == x**2 + 1
    assert cyclotomic_poly(5, x) == x**4 + x**3 + x**2 + x + 1
    assert cyclotomic_poly(6, x) == x**2 - x + 1


def test_symmetric_poly():
    pytest.raises(ValueError, lambda: symmetric_poly(-1, x, y, z))
    pytest.raises(ValueError, lambda: symmetric_poly(5, x, y, z))

    assert symmetric_poly(1, x, y, z, polys=True) == (x + y + z).as_poly()

    assert symmetric_poly(0, x, y, z) == 1
    assert symmetric_poly(1, x, y, z) == x + y + z
    assert symmetric_poly(2, x, y, z) == x*y + x*z + y*z
    assert symmetric_poly(3, x, y, z) == x*y*z


def test_random_poly():
    poly = random_poly(x, 10, -100, 100)

    assert poly.as_poly().degree() == 10
    assert all(-100 <= coeff <= 100 for coeff in poly.as_poly().coeffs()) is True

    poly = random_poly(x, 10, -100, 100, polys=True)

    assert poly.degree() == 10
    assert all(-100 <= coeff <= 100 for coeff in poly.coeffs()) is True

    poly = random_poly(x, 0, -10, 10, polys=True)

    assert poly.degree() == 0
    assert all(-10 <= c <= 10 for c in poly.all_coeffs())

    poly = random_poly(x, 1, -20, 20, polys=True)

    assert poly.degree() == 1
    assert all(-20 <= c <= 20 for c in poly.all_coeffs())

    poly = random_poly(x, 2, -30, 30, polys=True)

    assert poly.degree() == 2
    assert all(-30 <= c <= 30 for c in poly.all_coeffs())

    poly = random_poly(x, 3, -40, 40, polys=True)

    assert poly.degree() == 3
    assert all(-40 <= c <= 40 for c in poly.all_coeffs())

    poly = random_poly(x, 3, -400, 400, polys=True)

    assert poly.degree() == 3
    assert all(-400 <= c <= 400 for c in poly.all_coeffs())

    random.seed(11)
    assert random_poly(x, 10, -1, 1, polys=True).all_coeffs() == [1, 0, 1, 1,
                                                                  -1, 0, 0, -1,
                                                                  0, 0, 1]

    for i in range(10):
        poly = random_poly(x, 3, -10, 10, percent=50, polys=True)
        assert poly.all_coeffs()[-1]
        assert len([c for c in poly.all_coeffs() if c == 0]) == 2


def test_interpolating_poly():
    x0, x1, x2, y0, y1, y2 = symbols('x:3, y:3')

    assert interpolating_poly(0, x) == 0
    assert interpolating_poly(1, x) == y0

    assert interpolating_poly(2, x) == \
        y0*(x - x1)/(x0 - x1) + y1*(x - x0)/(x1 - x0)

    assert interpolating_poly(3, x) == \
        y0*(x - x1)*(x - x2)/((x0 - x1)*(x0 - x2)) + \
        y1*(x - x0)*(x - x2)/((x1 - x0)*(x1 - x2)) + \
        y2*(x - x0)*(x - x1)/((x2 - x0)*(x2 - x1))


def test_fateman_poly_F_1():
    R, x, y = ring('x y', ZZ)

    f, g, h = R.fateman_poly_F_1()

    assert f == (1 + sum(R.gens))*(2 + sum(R.gens))
    assert g == (1 + sum(_**2 for _ in R.gens))*(-3*y*x**2 + y**2 - 1)
    assert h == 1

    R, x, y, z, t = ring('x y z t', ZZ)

    f, g, h = R.fateman_poly_F_1()

    assert f == (1 + sum(R.gens))*(2 + sum(R.gens))
    assert g == (1 + sum(_**2 for _ in R.gens))*(-3*y*x**2 + y**2 - 1)


def test_fateman_poly_F_2():
    R, x, y = ring('x y', ZZ)

    f, g, h = R.fateman_poly_F_2()
    D = (1 + sum(R.gens))**2

    assert f == D*(-2 + x - sum(R.gens[1:]))**2
    assert g == D*(2 + sum(R.gens))**2
    assert h == D

    R, x, y, z, t = ring('x y z t', ZZ)

    f, g, h = R.fateman_poly_F_2()
    D = (1 + sum(R.gens))**2

    assert f == D*(-2 + x - sum(R.gens[1:]))**2
    assert g == D*(2 + sum(R.gens))**2
    assert h == D


def test_fateman_poly_F_3():
    R, x, y = ring('x y', ZZ)

    f, g, h = R.fateman_poly_F_3()
    D = (1 + sum(_**R.ngens for _ in R.gens))**2

    assert f == D*(-2 + x**R.ngens - sum(_**R.ngens for _ in R.gens[1:]))**2
    assert g == D*(+2 + sum(_**R.ngens for _ in R.gens))**2
    assert h == D

    R, x, y, z, t = ring('x y z t', ZZ)

    f, g, h = R.fateman_poly_F_3()
    D = (1 + sum(_**R.ngens for _ in R.gens))**2

    assert f == D*(-2 + x**R.ngens - sum(_**R.ngens for _ in R.gens[1:]))**2
    assert g == D*(+2 + sum(_**R.ngens for _ in R.gens))**2
    assert h == D
