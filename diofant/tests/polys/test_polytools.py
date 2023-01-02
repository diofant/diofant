"""Tests for user-friendly public interface to polynomial functions."""

import functools
import math

import pytest

from diofant import (CC, EX, FF, LC, LM, LT, QQ, RR, ZZ, CoercionFailedError,
                     ComputationFailedError, Derivative, DomainError, E, Eq,
                     ExactQuotientFailedError, Expr, FlagError, Float,
                     GeneratorsError, GeneratorsNeededError, GroebnerBasis, I,
                     Integer, Integral, MatrixSymbol, Mul,
                     MultivariatePolynomialError, O, OptionError, Piecewise,
                     PolificationFailedError, Poly, PolynomialError, PurePoly,
                     Rational, RealField, RootOf, Sum, Symbol, Tuple,
                     UnificationFailedError, cancel, cofactors, compose,
                     content, cos, count_roots, decompose, degree, diff,
                     discriminant, div, exp, expand, exquo, factor,
                     factor_list, false, gcd, gcdex, grevlex, grlex, groebner,
                     half_gcdex, im, invert, lcm, lex, log, monic, nroots, oo,
                     parallel_poly_from_expr, pi, primitive, quo, re,
                     real_roots, reduced, rem, resultant, sin, sqf, sqf_list,
                     sqf_norm, sqf_part, sqrt, subresultants, symbols, sympify,
                     tan, tanh, terms_gcd, true, trunc)
from diofant.abc import a, b, c, d, p, q, t, w, x, y, z
from diofant.core.mul import _keep_coeff
from diofant.polys.polytools import to_rational_coeffs


__all__ = ()


def _epsilon_eq(a, b):
    for x, y in zip(a, b):
        if abs(x - y) > 1e-10:
            return False
    return True


def test_Poly_from_dict():
    K = FF(3)

    assert Poly.from_dict({0: 1, 1: 2}, gens=x,
                          domain=K).rep.all_coeffs() == [K(1), K(2)]
    assert Poly.from_dict({0: 1, 1: 5}, gens=x,
                          domain=K).rep.all_coeffs() == [K(1), K(2)]

    assert Poly.from_dict({(0,): 1, (1,): 2}, gens=x,
                          domain=K).rep.all_coeffs() == [K(1), K(2)]
    assert Poly.from_dict({(0,): 1, (1,): 5}, gens=x,
                          domain=K).rep.all_coeffs() == [K(1), K(2)]

    assert dict(Poly.from_dict({(0, 0): 1, (1, 1): 2}, gens=(x, y),
                               domain=K).rep) == {(0, 0): K(1), (1, 1): K(2)}

    assert Poly.from_dict({0: 1, 1: 2}, gens=x).rep.all_coeffs() == [ZZ(1), ZZ(2)]
    assert Poly.from_dict({0: 1, 1: 2}, gens=x,
                          field=True).rep.all_coeffs() == [QQ(1), QQ(2)]

    assert Poly.from_dict({0: 1, 1: 2}, gens=x,
                          domain=ZZ).rep.all_coeffs() == [ZZ(1), ZZ(2)]
    assert Poly.from_dict({0: 1, 1: 2}, gens=x,
                          domain=QQ).rep.all_coeffs() == [QQ(1), QQ(2)]

    assert Poly.from_dict({(0,): 1, (1,): 2},
                          gens=x).rep.all_coeffs() == [ZZ(1), ZZ(2)]
    assert Poly.from_dict({(0,): 1, (1,): 2}, gens=x,
                          field=True).rep.all_coeffs() == [QQ(1), QQ(2)]

    assert Poly.from_dict({(0,): 1, (1,): 2}, gens=x,
                          domain=ZZ).rep.all_coeffs() == [ZZ(1), ZZ(2)]
    assert Poly.from_dict({(0,): 1, (1,): 2}, gens=x,
                          domain=QQ).rep.all_coeffs() == [QQ(1), QQ(2)]

    assert Poly.from_dict({(1,): sin(y)}, gens=x, composite=False) == \
        (sin(y)*x).as_poly(x, domain=EX)
    assert Poly.from_dict({(1,): y}, gens=x, composite=False) == \
        (y*x).as_poly(x, domain=EX)
    assert Poly.from_dict({(1, 1): 1}, gens=(x, y), composite=False) == \
        (x*y).as_poly(x, y, domain=ZZ)
    assert Poly.from_dict({(1, 0): y}, gens=(x, z), composite=False) == \
        (y*x).as_poly(x, z, domain=EX)

    pytest.raises(GeneratorsError,
                  lambda: Poly.from_dict({(1,): x, (0,): 1}, gens=(x,)))


def test_Poly_from_list():
    K = FF(3)

    assert Poly.from_list([2, 1], gens=x, domain=K).rep.all_coeffs() == [K(2), K(1)]
    assert Poly.from_list([5, 1], gens=x, domain=K).rep.all_coeffs() == [K(2), K(1)]

    assert Poly.from_list([2, 1], gens=x).rep.all_coeffs() == [ZZ(2), ZZ(1)]
    assert Poly.from_list([2, 1], gens=x, field=True).rep.all_coeffs() == [QQ(2), QQ(1)]

    assert Poly.from_list([2, 1], gens=x, domain=ZZ).rep.all_coeffs() == [ZZ(2), ZZ(1)]
    assert Poly.from_list([2, 1], gens=x, domain=QQ).rep.all_coeffs() == [QQ(2), QQ(1)]

    assert Poly.from_list([0, 1.0], gens=x).rep.all_coeffs() == [RR(0), RR(1.0)]
    assert Poly.from_list([1.0, 0], gens=x).rep.all_coeffs() == [RR(1.0)]

    pytest.raises(MultivariatePolynomialError, lambda: Poly.from_list([[]], gens=(x, y)))
    pytest.raises(GeneratorsError, lambda: Poly.from_list([x, 1], gens=(x,)))


def test_Poly_from_poly():
    f = (x + 7).as_poly(x, domain=ZZ)
    g = (x + 2).as_poly(x, modulus=3)
    h = (x + y).as_poly(x, y, domain=ZZ)

    K = FF(3)

    assert Poly.from_poly(f) == f
    assert Poly.from_poly(f, domain=K).rep.all_coeffs() == [K(1), K(1)]
    assert Poly.from_poly(f, domain=ZZ).rep.all_coeffs() == [ZZ(7), ZZ(1)]
    assert Poly.from_poly(f, domain=QQ).rep.all_coeffs() == [QQ(7), QQ(1)]

    assert Poly.from_poly(f, gens=x) == f
    assert Poly.from_poly(f, gens=x, domain=K).rep.all_coeffs() == [K(1), K(1)]
    assert Poly.from_poly(f, gens=x, domain=ZZ).rep.all_coeffs() == [ZZ(7), ZZ(1)]
    assert Poly.from_poly(f, gens=x, domain=QQ).rep.all_coeffs() == [QQ(7), QQ(1)]

    assert Poly.from_poly(f, gens=y) == (x + 7).as_poly(y, domain=ZZ.inject(x))
    pytest.raises(CoercionFailedError, lambda: Poly.from_poly(f, gens=y, domain=K))
    pytest.raises(CoercionFailedError, lambda: Poly.from_poly(f, gens=y, domain=ZZ))
    pytest.raises(CoercionFailedError, lambda: Poly.from_poly(f, gens=y, domain=QQ))

    assert Poly.from_poly(f, gens=(x, y)) == (x + 7).as_poly(x, y, domain=ZZ)
    assert Poly.from_poly(
        f, gens=(x, y), domain=ZZ) == (x + 7).as_poly(x, y, domain=ZZ)
    assert Poly.from_poly(
        f, gens=(x, y), domain=QQ) == (x + 7).as_poly(x, y, domain=QQ)
    assert Poly.from_poly(
        f, gens=(x, y), modulus=3) == (x + 7).as_poly(x, y, domain=FF(3))

    K = FF(2)

    assert Poly.from_poly(g) == g
    assert Poly.from_poly(g, domain=ZZ).rep.all_coeffs() == [ZZ(2), ZZ(1)]
    pytest.raises(CoercionFailedError, lambda: Poly.from_poly(g, domain=QQ))
    assert Poly.from_poly(g, domain=K).rep.all_coeffs() == [K(0), K(1)]

    assert Poly.from_poly(g, gens=x) == g
    assert Poly.from_poly(g, gens=x, domain=ZZ).rep.all_coeffs() == [ZZ(2), ZZ(1)]
    pytest.raises(CoercionFailedError, lambda: Poly.from_poly(g, gens=x, domain=QQ))
    assert Poly.from_poly(g, gens=x, domain=K).rep.all_coeffs() == [K(0), K(1)]

    K = FF(3)

    assert Poly.from_poly(h) == h
    assert dict(Poly.from_poly(h, domain=ZZ).rep) == {(1, 0): ZZ(1), (0, 1): ZZ(1)}
    assert dict(Poly.from_poly(h, domain=QQ).rep) == {(1, 0): QQ(1), (0, 1): QQ(1)}
    assert dict(Poly.from_poly(h, domain=K).rep) == {(1, 0): K(1), (0, 1): K(1)}

    assert Poly.from_poly(h, gens=x) == (x + y).as_poly(x, domain=ZZ.inject(y))
    pytest.raises(CoercionFailedError, lambda: Poly.from_poly(h, gens=x, domain=ZZ))
    assert Poly.from_poly(
        h, gens=x, domain=ZZ.inject(y)) == (x + y).as_poly(x, domain=ZZ.inject(y))
    pytest.raises(CoercionFailedError, lambda: Poly.from_poly(h, gens=x, domain=QQ))
    assert Poly.from_poly(
        h, gens=x, domain=QQ.inject(y)) == (x + y).as_poly(x, domain=QQ.inject(y))
    pytest.raises(CoercionFailedError, lambda: Poly.from_poly(h, gens=x, modulus=3))

    assert Poly.from_poly(h, gens=y) == (x + y).as_poly(y, domain=ZZ.inject(x))
    pytest.raises(CoercionFailedError, lambda: Poly.from_poly(h, gens=y, domain=ZZ))
    assert Poly.from_poly(
        h, gens=y, domain=ZZ.inject(x)) == (x + y).as_poly(y, domain=ZZ.inject(x))
    pytest.raises(CoercionFailedError, lambda: Poly.from_poly(h, gens=y, domain=QQ))
    assert Poly.from_poly(
        h, gens=y, domain=QQ.inject(x)) == (x + y).as_poly(y, domain=QQ.inject(x))
    pytest.raises(CoercionFailedError, lambda: Poly.from_poly(h, gens=y, modulus=3))

    assert Poly.from_poly(h, gens=(x, y)) == h
    assert dict(Poly.from_poly(h, gens=(x, y),
                               domain=ZZ).rep) == {(1, 0): ZZ(1), (0, 1): ZZ(1)}
    assert dict(Poly.from_poly(h, gens=(x, y),
                               domain=QQ).rep) == {(1, 0): QQ(1), (0, 1): QQ(1)}
    assert dict(Poly.from_poly(h, gens=(x, y),
                               domain=K).rep) == {(1, 0): K(1), (0, 1): K(1)}

    assert dict(Poly.from_poly(h, gens=(y, x)).rep) == {(1, 0): ZZ(1), (0, 1): ZZ(1)}
    assert dict(Poly.from_poly(h, gens=(y, x),
                               domain=ZZ).rep) == {(1, 0): ZZ(1), (0, 1): ZZ(1)}
    assert dict(Poly.from_poly(h, gens=(y, x),
                               domain=QQ).rep) == {(1, 0): QQ(1), (0, 1): QQ(1)}
    assert dict(Poly.from_poly(h, gens=(y, x),
                               domain=K).rep) == {(1, 0): K(1), (0, 1): K(1)}

    assert dict(Poly.from_poly(h, gens=(x, y),
                               field=True).rep) == {(1, 0): QQ(1), (0, 1): QQ(1)}
    assert dict(Poly.from_poly(h, gens=(x, y),
                               field=True).rep) == {(1, 0): QQ(1), (0, 1): QQ(1)}


def test_Poly_from_expr():
    pytest.raises(GeneratorsNeededError, lambda: Poly.from_expr(Integer(0)))
    pytest.raises(GeneratorsNeededError, lambda: Poly.from_expr(Integer(7)))

    F3 = FF(3)

    assert Poly.from_expr(x + 5, domain=F3).rep.all_coeffs() == [F3(2), F3(1)]
    assert Poly.from_expr(y + 5, domain=F3).rep.all_coeffs() == [F3(2), F3(1)]

    assert Poly.from_expr(x + 5, x, domain=F3).rep.all_coeffs() == [F3(2), F3(1)]
    assert Poly.from_expr(y + 5, y, domain=F3).rep.all_coeffs() == [F3(2), F3(1)]

    assert dict(Poly.from_expr(x + y, domain=F3).rep) == {(1, 0): F3(1), (0, 1): F3(1)}
    assert dict(Poly.from_expr(x + y, x, y, domain=F3).rep) == {(1, 0): F3(1), (0, 1): F3(1)}

    assert Poly.from_expr(x + 5).rep.all_coeffs() == [ZZ(5), ZZ(1)]
    assert Poly.from_expr(y + 5).rep.all_coeffs() == [ZZ(5), ZZ(1)]

    assert Poly.from_expr(x + 5, x).rep.all_coeffs() == [ZZ(5), ZZ(1)]
    assert Poly.from_expr(y + 5, y).rep.all_coeffs() == [ZZ(5), ZZ(1)]

    assert Poly.from_expr(x + 5, domain=ZZ).rep.all_coeffs() == [ZZ(5), ZZ(1)]
    assert Poly.from_expr(y + 5, domain=ZZ).rep.all_coeffs() == [ZZ(5), ZZ(1)]

    assert Poly.from_expr(x + 5, x, domain=ZZ).rep.all_coeffs() == [ZZ(5), ZZ(1)]
    assert Poly.from_expr(y + 5, y, domain=ZZ).rep.all_coeffs() == [ZZ(5), ZZ(1)]

    assert dict(Poly.from_expr(x + 5, x, y, domain=ZZ).rep) == {(1, 0): ZZ(1), (0, 0): ZZ(5)}
    assert dict(Poly.from_expr(y + 5, x, y, domain=ZZ).rep) == {(0, 1): ZZ(1), (0, 0): ZZ(5)}


def test_Poly__new__():
    pytest.raises(GeneratorsError, lambda: Poly(x + 1, x, x))

    pytest.raises(GeneratorsError, lambda: Poly(x + y, x, y, domain=ZZ.inject(x)))
    pytest.raises(GeneratorsError, lambda: Poly(x + y, x, y, domain=ZZ.inject(y)))

    pytest.raises(OptionError, lambda: Poly(x + 2, x, modulus=3, domain=QQ))

    pytest.raises(OptionError, lambda: Poly(x + 2, x, domain=ZZ, gaussian=True))
    pytest.raises(OptionError, lambda: Poly(x + 2, x, modulus=3, gaussian=True))

    pytest.raises(OptionError, lambda: Poly(x + 2, x, domain=ZZ, extension=[sqrt(3)]))
    pytest.raises(OptionError, lambda: Poly(x + 2, x, modulus=3, extension=[sqrt(3)]))

    pytest.raises(OptionError, lambda: Poly(x + 2, x, domain=ZZ, extension=True))
    pytest.raises(OptionError, lambda: Poly(x + 2, x, modulus=3, extension=True))

    pytest.raises(OptionError, lambda: Poly(x + 2, x, domain=ZZ, greedy=True))
    pytest.raises(OptionError, lambda: Poly(x + 2, x, domain=QQ, field=True))

    pytest.raises(OptionError, lambda: Poly(x + 2, x, domain=ZZ, greedy=False))
    pytest.raises(OptionError, lambda: Poly(x + 2, x, domain=QQ, field=False))

    f = (x + 1).as_poly(x, modulus=3, order='grlex')

    assert f.get_modulus() == 3
    assert f.rep.ring.order == grlex

    f = (x + 1).as_poly(x, order='grlex')

    assert f.rep.ring.order == grlex

    pytest.raises(GeneratorsNeededError, lambda: Poly({1: 2, 0: 1}))
    pytest.raises(GeneratorsNeededError, lambda: Poly([2, 1]))
    pytest.raises(GeneratorsNeededError, lambda: Poly((2, 1)))

    pytest.raises(GeneratorsNeededError, lambda: Poly(1))

    f = a*x**2 + b*x + c

    assert Poly({2: a, 1: b, 0: c}, x) == f
    assert Poly(iter([c, b, a]), x) == f
    assert Poly([c, b, a], x) == f
    assert Poly((c, b, a), x) == f

    f = Poly({}, x, y, z)

    assert f.gens == (x, y, z)
    assert f.as_expr() == 0

    assert ((a*x + b*y).as_poly(x, y)).as_poly(x) == (a*x + b*y).as_poly(x)

    assert (3*x**2 + 2*x + 1).as_poly(domain=ZZ).all_coeffs() == [1, 2, 3]
    assert (3*x**2 + 2*x + 1).as_poly(domain=QQ).all_coeffs() == [1, 2, 3]
    assert (3*x**2 + 2*x + 1).as_poly(domain=RR).all_coeffs() == [1.0, 2.0, 3.0]

    pytest.raises(CoercionFailedError, lambda: Poly(3*x**2/5 + 2*x/5 + 1, x, domain=ZZ))
    assert (3*x**2/5 + 2*x/5 + 1).as_poly(domain=QQ).all_coeffs() == [1, Rational(2, 5), Rational(3, 5)]
    assert _epsilon_eq(
        (3*x**2/5 + 2*x/5 + 1).as_poly(domain=RR).all_coeffs(), [1.0, 0.4, 0.6])

    assert (3.0*x**2 + 2.0*x + 1).as_poly(domain=ZZ).all_coeffs() == [1, 2, 3]
    assert (3.0*x**2 + 2.0*x + 1).as_poly(domain=QQ).all_coeffs() == [1, 2, 3]
    assert (3.0*x**2 + 2.0*x + 1).as_poly(domain=RR).all_coeffs() == [1.0, 2.0, 3.0]

    pytest.raises(CoercionFailedError, lambda: Poly(3.1*x**2 + 2.1*x + 1, x, domain=ZZ))
    assert (3.1*x**2 + 2.1*x + 1).as_poly(domain=QQ).all_coeffs() == [1, Rational(21, 10), Rational(31, 10)]
    assert (3.1*x**2 + 2.1*x + 1).as_poly(domain=RR).all_coeffs() == [1.0, 2.1, 3.1]

    assert Poly({(2, 1): 1, (1, 2): 2, (1, 1): 3}, x, y) == \
        (x**2*y + 2*x*y**2 + 3*x*y).as_poly(x, y)

    assert (x**2 + 1).as_poly(extension=I).domain == QQ.algebraic_field(I)

    f = 3*x**5 - x**4 + x**3 - x**2 + 65538

    assert f.as_poly(x, modulus=65537) == \
        (3*x**5 + 65536*x**4 + x**3 + 65536*x**2 + 1).as_poly(x, modulus=65537)

    assert isinstance((x**2 + x + 1.0).as_poly().domain, RealField)


def test_Poly_new():
    pytest.raises(PolynomialError, lambda: Poly.new([1], x))
    pytest.raises(PolynomialError, lambda: Poly.new((x + 1).as_poly(x, y).rep, x))


def test_Poly__args():
    assert (x**2 + 1).as_poly().args == (x**2 + 1, x)


def test_Poly_is_number():
    assert Integer(1).as_poly(x).is_number
    assert x.as_poly().is_number is False


def test_Poly__gens():
    assert ((x - p)*(x - q)).as_poly(x).gens == (x,)
    assert ((x - p)*(x - q)).as_poly(p).gens == (p,)
    assert ((x - p)*(x - q)).as_poly(q).gens == (q,)

    assert ((x - p)*(x - q)).as_poly(x, p).gens == (x, p)
    assert ((x - p)*(x - q)).as_poly(x, q).gens == (x, q)

    assert ((x - p)*(x - q)).as_poly(x, p, q).gens == (x, p, q)
    assert ((x - p)*(x - q)).as_poly(p, x, q).gens == (p, x, q)
    assert ((x - p)*(x - q)).as_poly(p, q, x).gens == (p, q, x)

    assert ((x - p)*(x - q)).as_poly().gens == (x, p, q)

    assert ((x - p)*(x - q)).as_poly(sort='x > p > q').gens == (x, p, q)
    assert ((x - p)*(x - q)).as_poly(sort='p > x > q').gens == (p, x, q)
    assert ((x - p)*(x - q)).as_poly(sort='p > q > x').gens == (p, q, x)

    assert ((x - p)*(x - q)).as_poly(x, p, q, sort='p > q > x').gens == (x, p, q)

    assert ((x - p)*(x - q)).as_poly(wrt='x').gens == (x, p, q)
    assert ((x - p)*(x - q)).as_poly(wrt='p').gens == (p, x, q)
    assert ((x - p)*(x - q)).as_poly(wrt='q').gens == (q, x, p)

    assert ((x - p)*(x - q)).as_poly(wrt=x).gens == (x, p, q)
    assert ((x - p)*(x - q)).as_poly(wrt=p).gens == (p, x, q)
    assert ((x - p)*(x - q)).as_poly(wrt=q).gens == (q, x, p)

    assert ((x - p)*(x - q)).as_poly(x, p, q, wrt='p').gens == (x, p, q)

    assert ((x - p)*(x - q)).as_poly(wrt='p', sort='q > x').gens == (p, q, x)
    assert ((x - p)*(x - q)).as_poly(wrt='q', sort='p > x').gens == (q, p, x)


def test_Poly_unify():
    pytest.raises(UnificationFailedError, lambda: x.as_poly().unify(y))
    pytest.raises(UnificationFailedError, lambda: PurePoly(x).unify(y))
    pytest.raises(UnificationFailedError, lambda: PurePoly(x).unify(Poly(x, x, y)))

    assert x.as_poly(modulus=3).unify(y.as_poly(modulus=3)) == \
        (x.as_poly(x, y, modulus=3), y.as_poly(x, y, modulus=3))
    pytest.raises(NotImplementedError,
                  lambda: x.as_poly(modulus=3).unify(y.as_poly(modulus=5)))

    assert y.as_poly(x, y).unify(x.as_poly(modulus=3)) == \
        (y.as_poly(x, y, modulus=3), x.as_poly(x, y, modulus=3))
    assert x.as_poly(modulus=3).unify(y.as_poly(x, y)) == \
        (x.as_poly(x, y, modulus=3), y.as_poly(x, y, modulus=3))

    assert (x + 1).as_poly().unify((x + 2).as_poly()) == \
        ((x + 1).as_poly(domain=ZZ), (x + 2).as_poly(domain=ZZ))
    assert (x + 1).as_poly(domain=QQ).unify((x + 2).as_poly()) == \
        ((x + 1).as_poly(domain=QQ), (x + 2).as_poly(domain=QQ))
    assert (x + 1).as_poly().unify((x + 2).as_poly(domain=QQ)) == \
        ((x + 1).as_poly(domain=QQ), (x + 2).as_poly(domain=QQ))

    assert (x + 1).as_poly().unify((x + 2).as_poly(x, y)) == \
        ((x + 1).as_poly(x, y, domain=ZZ), (x + 2).as_poly(x, y, domain=ZZ))
    assert (x + 1).as_poly(domain=QQ).unify((x + 2).as_poly(x, y)) == \
        ((x + 1).as_poly(x, y, domain=QQ), (x + 2).as_poly(x, y, domain=QQ))
    assert (x + 1).as_poly().unify((x + 2).as_poly(x, y, domain=QQ)) == \
        ((x + 1).as_poly(x, y, domain=QQ), (x + 2).as_poly(x, y, domain=QQ))

    assert (x + 1).as_poly(x, y).unify((x + 2).as_poly()) == \
        ((x + 1).as_poly(x, y, domain=ZZ), (x + 2).as_poly(x, y, domain=ZZ))
    assert (x + 1).as_poly(x, y, domain=QQ).unify((x + 2).as_poly()) == \
        ((x + 1).as_poly(x, y, domain=QQ), (x + 2).as_poly(x, y, domain=QQ))
    assert (x + 1).as_poly(x, y).unify((x + 2).as_poly(domain=QQ)) == \
        ((x + 1).as_poly(x, y, domain=QQ), (x + 2).as_poly(x, y, domain=QQ))

    assert (x + 1).as_poly(x, y).unify((x + 2).as_poly(x, y)) == \
        ((x + 1).as_poly(x, y, domain=ZZ), (x + 2).as_poly(x, y, domain=ZZ))
    assert (x + 1).as_poly(x, y, domain=QQ).unify((x + 2).as_poly(x, y)) == \
        ((x + 1).as_poly(x, y, domain=QQ), (x + 2).as_poly(x, y, domain=QQ))
    assert (x + 1).as_poly(x, y).unify((x + 2).as_poly(x, y, domain=QQ)) == \
        ((x + 1).as_poly(x, y, domain=QQ), (x + 2).as_poly(x, y, domain=QQ))

    assert (x + 1).as_poly().unify((x + 2).as_poly(y, x)) == \
        ((x + 1).as_poly(y, x, domain=ZZ), (x + 2).as_poly(y, x, domain=ZZ))
    assert (x + 1).as_poly(domain=QQ).unify((x + 2).as_poly(y, x)) == \
        ((x + 1).as_poly(y, x, domain=QQ), (x + 2).as_poly(y, x, domain=QQ))
    assert (x + 1).as_poly().unify((x + 2).as_poly(y, x, domain=QQ)) == \
        ((x + 1).as_poly(y, x, domain=QQ), (x + 2).as_poly(y, x, domain=QQ))

    assert (x + 1).as_poly(y, x).unify((x + 2).as_poly()) == \
        ((x + 1).as_poly(y, x, domain=ZZ), (x + 2).as_poly(y, x, domain=ZZ))
    assert (x + 1).as_poly(y, x, domain=QQ).unify((x + 2).as_poly()) == \
        ((x + 1).as_poly(y, x, domain=QQ), (x + 2).as_poly(y, x, domain=QQ))
    assert (x + 1).as_poly(y, x).unify((x + 2).as_poly(domain=QQ)) == \
        ((x + 1).as_poly(y, x, domain=QQ), (x + 2).as_poly(y, x, domain=QQ))

    assert (x + 1).as_poly(x, y).unify((x + 2).as_poly(y, x)) == \
        ((x + 1).as_poly(x, y, domain=ZZ), (x + 2).as_poly(x, y, domain=ZZ))
    assert (x + 1).as_poly(x, y, domain=QQ).unify((x + 2).as_poly(y, x)) == \
        ((x + 1).as_poly(x, y, domain=QQ), (x + 2).as_poly(x, y, domain=QQ))
    assert (x + 1).as_poly(x, y).unify((x + 2).as_poly(y, x, domain=QQ)) == \
        ((x + 1).as_poly(x, y, domain=QQ), (x + 2).as_poly(x, y, domain=QQ))

    assert (x + 1).as_poly(y, x).unify((x + 2).as_poly(x, y)) == \
        ((x + 1).as_poly(y, x, domain=ZZ), (x + 2).as_poly(y, x, domain=ZZ))
    assert (x + 1).as_poly(y, x, domain=QQ).unify((x + 2).as_poly(x, y)) == \
        ((x + 1).as_poly(y, x, domain=QQ), (x + 2).as_poly(y, x, domain=QQ))
    assert (x + 1).as_poly(y, x).unify((x + 2).as_poly(x, y, domain=QQ)) == \
        ((x + 1).as_poly(y, x, domain=QQ), (x + 2).as_poly(y, x, domain=QQ))

    assert (a*x).as_poly(x, domain=ZZ.inject(a)).unify((a*b*x).as_poly(x, domain=ZZ.inject(a, b).field)) == \
        ((a*x).as_poly(x, domain=ZZ.inject(a, b).field), (a*b*x).as_poly(x, domain=ZZ.inject(a, b).field))

    assert (a*x).as_poly(x, domain=ZZ.inject(a).field).unify((a*b*x).as_poly(x, domain=ZZ.inject(a, b).field)) == \
        ((a*x).as_poly(x, domain=ZZ.inject(a, b).field), (a*b*x).as_poly(x, domain=ZZ.inject(a, b).field))

    pytest.raises(CoercionFailedError, lambda: Poly((x**2 + x**2*z).as_poly(y, field=True),
                                                    domain=ZZ.inject(x).field))

    f = (t**2 + t/3 + x).as_poly(t, domain=QQ.inject(x).field)
    g = (t**2 + t/3 + x).as_poly(t, domain=QQ.inject(x))

    assert f.unify(g) == (f, f)


def test_Poly_free_symbols():
    assert (x**2 + 1).as_poly().free_symbols == {x}
    assert (x**2 + y*z).as_poly().free_symbols == {x, y, z}
    assert (x**2 + y*z).as_poly(x).free_symbols == {x, y, z}
    assert (x**2 + sin(y*z)).as_poly().free_symbols == {x, y, z}
    assert (x**2 + sin(y*z)).as_poly(x).free_symbols == {x, y, z}
    assert (x**2 + sin(y*z)).as_poly(x, domain=EX).free_symbols == {x, y, z}


def test_PurePoly_free_symbols():
    assert PurePoly(x**2 + 1).free_symbols == set()
    assert PurePoly(x**2 + y*z).free_symbols == set()
    assert PurePoly(x**2 + y*z, x).free_symbols == {y, z}
    assert PurePoly(x**2 + sin(y*z)).free_symbols == set()
    assert PurePoly(x**2 + sin(y*z), x).free_symbols == {y, z}
    assert PurePoly(x**2 + sin(y*z), x, domain=EX).free_symbols == {y, z}


def test_Poly__eq__():
    assert (x.as_poly() == x.as_poly()) is True
    assert (x.as_poly(domain=QQ) == x.as_poly()) is True
    assert (x.as_poly() == x.as_poly(domain=QQ)) is True

    assert (x.as_poly(domain=ZZ.inject(a)) == x.as_poly()) is True
    assert (x.as_poly() == x.as_poly(domain=ZZ.inject(a))) is True

    assert ((x*y).as_poly(x, y) == x.as_poly()) is False

    assert (x.as_poly(x, y) == x.as_poly()) is False
    assert (x.as_poly() == x.as_poly(x, y)) is False

    assert ((x**2 + 1).as_poly() == (y**2 + 1).as_poly()) is False
    assert ((y**2 + 1).as_poly() == (x**2 + 1).as_poly()) is False

    f = x.as_poly(domain=ZZ)
    g = x.as_poly(domain=QQ)

    assert (f == g) is True
    assert (f != g) is False

    t0 = Symbol('t0')

    f = ((t0/2 + x**2)*t**2 - x**2*t).as_poly(t, domain=QQ.inject(x, t0))
    g = ((t0/2 + x**2)*t**2 - x**2*t).as_poly(t, domain=ZZ.inject(x, t0).field)

    assert (f == g) is True

    assert (x.as_poly() == y.as_poly()) is False

    with pytest.raises(NotImplementedError):
        assert x.as_poly(modulus=2) == x.as_poly(modulus=3)


def test_PurePoly__eq__():
    assert (PurePoly(x, x) == PurePoly(x, x)) is True
    assert (PurePoly(x, x, domain=QQ) == PurePoly(x, x)) is True
    assert (PurePoly(x, x) == PurePoly(x, x, domain=QQ)) is True

    assert (PurePoly(x, x, domain=ZZ.inject(a)) == PurePoly(x, x)) is True
    assert (PurePoly(x, x) == PurePoly(x, x, domain=ZZ.inject(a))) is True

    assert (PurePoly(x*y, x, y) == PurePoly(x, x)) is False

    assert (PurePoly(x, x, y) == PurePoly(x, x)) is False
    assert (PurePoly(x, x) == PurePoly(x, x, y)) is False

    assert (PurePoly(x**2 + 1, x) == PurePoly(y**2 + 1, y)) is True
    assert (PurePoly(y**2 + 1, y) == PurePoly(x**2 + 1, x)) is True

    f = PurePoly(x, x, domain=ZZ)
    g = PurePoly(x, x, domain=QQ)

    assert (f == g) is True
    assert (f != g) is False

    f = PurePoly(x, x, domain=ZZ)
    g = PurePoly(y, y, domain=QQ)

    assert (f == g) is True
    assert (f != g) is False

    assert (f == 1) is False

    assert (f == sin(x)) is False

    with pytest.raises(NotImplementedError):
        assert PurePoly(x, modulus=2) == PurePoly(x, modulus=3)


def test_PurePoly_Poly():
    assert isinstance(PurePoly((x**2 + 1).as_poly()), PurePoly) is True
    assert isinstance(PurePoly(x**2 + 1).as_poly(), Poly) is True


def test_Poly_domain():
    assert (2*x).as_poly().domain == ZZ

    assert (2*x).as_poly(domain=ZZ).domain == ZZ
    assert (2*x).as_poly(domain=QQ).domain == QQ

    assert (x/2).as_poly().domain == QQ

    assert (x/2).as_poly(domain=ZZ) == Poly({(1, 1): 1}, x, Rational(1, 2),
                                            domain=ZZ)
    assert (x/2).as_poly(domain=QQ).domain == QQ

    assert isinstance((0.2*x).as_poly().domain, RealField)


def test_Poly_set_domain():
    assert (2*x + 1).as_poly().set_domain(ZZ) == (2*x + 1).as_poly()
    assert (2*x + 1).as_poly().set_domain(QQ) == (2*x + 1).as_poly(domain=QQ)

    assert (Rational(2, 10)*x + Rational(1, 10)).as_poly().set_domain(RR) == (0.2*x + 0.1).as_poly()
    assert (0.2*x + 0.1).as_poly().set_domain(QQ) == (Rational(2, 10)*x + Rational(1, 10)).as_poly()

    pytest.raises(CoercionFailedError, lambda: (x/2 + 1).as_poly().set_domain(ZZ))
    pytest.raises(CoercionFailedError, lambda: (x + 1).as_poly(modulus=2).set_domain(QQ))

    pytest.raises(GeneratorsError, lambda: (x*y).as_poly(x, y).set_domain(ZZ.inject(y)))


def test_Poly_get_modulus():
    assert (x**2 + 1).as_poly(modulus=2).get_modulus() == 2
    assert (x**2 + 1).as_poly(modulus=8).get_modulus() == 8

    pytest.raises(PolynomialError, lambda: (x**2 + 1).as_poly().get_modulus())


def test_Poly_set_modulus():
    assert (x**2 + 1).as_poly(modulus=2).set_modulus(7) == (x**2 + 1).as_poly(modulus=7)
    assert (x**2 + 5).as_poly(modulus=7).set_modulus(2) == (x**2 + 1).as_poly(modulus=2)

    assert (x**2 + 1).as_poly().set_modulus(2) == (x**2 + 1).as_poly(modulus=2)

    assert (x**2 + 1).as_poly(modulus=2).set_modulus(4) == (x**2 + 1).as_poly(modulus=4)
    assert (x**2 + 7*x + 6).as_poly(modulus=4) == (x**2 + 3*x + 2).as_poly(modulus=4)

    pytest.raises(CoercionFailedError, lambda: (x/2 + 1).as_poly().set_modulus(2))


def test_Poly_quo_ground():
    assert (2*x + 4).as_poly().quo_ground(2) == (x + 2).as_poly()
    assert (2*x + 3).as_poly().quo_ground(2) == (x + 1).as_poly()


def test_Poly_exquo_ground():
    assert (2*x + 4).as_poly().exquo_ground(2) == (x + 2).as_poly()
    pytest.raises(ExactQuotientFailedError, lambda: (2*x + 3).as_poly().exquo_ground(2))


def test_Poly_abs():
    assert abs((-x + 1).as_poly()) == (x + 1).as_poly()


def test_Poly_neg():
    assert -(-x + 1).as_poly() == (x - 1).as_poly()


def test_Poly_add():
    assert Integer(0).as_poly(x) + Integer(0).as_poly(x) == Integer(0).as_poly(x)

    assert Integer(1).as_poly(x) + Integer(0).as_poly(x) == Integer(1).as_poly(x)
    assert Integer(1).as_poly(x, y) + Integer(0).as_poly(x) == Integer(1).as_poly(x, y)
    assert Integer(0).as_poly(x, y) + Integer(1).as_poly(x, y) == Integer(1).as_poly(x, y)

    assert (x + 1).as_poly() + 2 == (x + 3).as_poly()
    assert (x**2 + 1).as_poly() + (x - 2).as_poly() == (x**2 + x - 1).as_poly()

    assert Integer(1).as_poly(x) + x == (x + 1).as_poly()
    assert Integer(1).as_poly(x) + sin(x) == 1 + sin(x)
    assert sin(x) + Integer(1).as_poly(x) == sin(x) + 1

    assert x.as_poly() + 1 == (x + 1).as_poly()
    assert 1 + x.as_poly() == (x + 1).as_poly()


def test_Poly_sub():
    assert Integer(0).as_poly(x) - Integer(0).as_poly(x) == Integer(0).as_poly(x)

    assert Integer(1).as_poly(x) - Integer(0).as_poly(x) == Integer(1).as_poly(x)
    assert Integer(1).as_poly(x) - 1 == Integer(0).as_poly(x)
    assert Integer(1).as_poly(x, y) - Integer(0).as_poly(x) == Integer(1).as_poly(x, y)
    assert Integer(0).as_poly(x) - Integer(1).as_poly(x, y) == Integer(-1).as_poly(x, y)
    assert Integer(0).as_poly(x, y) - Integer(1).as_poly(x, y) == Integer(-1).as_poly(x, y)

    assert Integer(1).as_poly(x) - x == (1 - x).as_poly()
    assert Integer(1).as_poly(x) - sin(x) == 1 - sin(x)
    assert sin(x) - Integer(1).as_poly(x) == sin(x) - 1

    assert x.as_poly() - 1 == (x - 1).as_poly()
    assert 1 - x.as_poly() == (1 - x).as_poly()

    assert (x + 1).as_poly() - 2 == (x - 1).as_poly()
    assert (x**2 + 1).as_poly() - (x - 2).as_poly() == (x**2 - x + 3).as_poly()


def test_Poly_mul():
    assert Integer(0).as_poly(x) * Integer(0).as_poly(x) == Integer(0).as_poly(x)

    assert Integer(2).as_poly(x) * Integer(4).as_poly(x) == Integer(8).as_poly(x)
    assert Integer(2).as_poly(x, y) * Integer(4).as_poly(x) == Integer(8).as_poly(x, y)
    assert Integer(4).as_poly(x) * Integer(2).as_poly(x, y) == Integer(8).as_poly(x, y)
    assert Integer(4).as_poly(x, y) * Integer(2).as_poly(x, y) == Integer(8).as_poly(x, y)

    assert Integer(1).as_poly(x) * x == x.as_poly()
    assert Integer(1).as_poly(x) * sin(x) == sin(x)
    assert sin(x) * Integer(1).as_poly(x) == sin(x)

    assert x.as_poly() * 2 == (2*x).as_poly()
    assert 2 * x.as_poly() == (2*x).as_poly()

    assert (x + 1).as_poly() * 2 == (2*x + 2).as_poly()
    assert (x**2 + 1).as_poly() * (x - 2).as_poly() == (x**3 - 2*x**2 + x - 2).as_poly()


def test_Poly_pow():
    assert x.as_poly()**10 == (x**10).as_poly()
    assert (2*y).as_poly(x, y)**4 == (16*y**4).as_poly(x, y)
    assert (7*x*y).as_poly(x, y)**3 == (343*x**3*y**3).as_poly(x, y)
    assert (x*y + 1).as_poly(x, y)**(-1) == (x*y + 1)**(-1)
    assert (x*y + 1).as_poly(x, y)**x == (x*y + 1)**x
    assert (x - 2).as_poly()**3 == (x**3 - 6*x**2 + 12*x - 8).as_poly()
    assert (x*y).as_poly(x, y)**2 == (x**2*y**2).as_poly(x, y)
    assert (x - 2).as_poly()**2 == (x**2 - 4*x + 4).as_poly()

    f, g = (x**3 + x - 1).as_poly(), (x**3 + 1).as_poly()
    r = pow(f, 3, g)

    assert r == f**3 % g
    assert r == (-6*x**2 + 12*x - 9).as_poly()


def test_Poly_divmod():
    f, g = (x**2).as_poly(), x.as_poly()
    q, r = g, Integer(0).as_poly(x)

    assert divmod(f, g) == (q, r)
    assert f // g == q
    assert f % g == r

    assert divmod(f, x) == (q, r)
    assert f // x == q
    assert f % x == r

    q, r = Integer(0).as_poly(x), Integer(2).as_poly(x)

    assert divmod(2, g) == (q, r)
    assert 2 // g == q
    assert 2 % g == r

    assert x.as_poly()/x.as_poly() == 1
    assert (x**2).as_poly()/x.as_poly() == x
    assert x.as_poly()/(x**2).as_poly() == 1/x


def test_Poly_eq_ne():
    assert ((x + y).as_poly(x, y) == (x + y).as_poly(x, y)) is True
    assert ((x + y).as_poly(x) == (x + y).as_poly(x, y)) is False
    assert ((x + y).as_poly(x, y) == (x + y).as_poly(x)) is False
    assert ((x + y).as_poly(x) == (x + y).as_poly(x)) is True
    assert ((x + y).as_poly(y) == (x + y).as_poly(y)) is True

    assert ((x + y).as_poly(x, y) == x + y) is True
    assert ((x + y).as_poly(x) == x + y) is True
    assert ((x + y).as_poly(x, y) == x + y) is True
    assert ((x + y).as_poly(x) == x + y) is True
    assert ((x + y).as_poly(y) == x + y) is True

    assert ((x + y).as_poly(x, y) != (x + y).as_poly(x, y)) is False
    assert ((x + y).as_poly(x) != (x + y).as_poly(x, y)) is True
    assert ((x + y).as_poly(x, y) != (x + y).as_poly(x)) is True
    assert ((x + y).as_poly(x) != (x + y).as_poly(x)) is False
    assert ((x + y).as_poly(y) != (x + y).as_poly(y)) is False

    assert ((x + y).as_poly(x, y) != x + y) is False
    assert ((x + y).as_poly(x) != x + y) is False
    assert ((x + y).as_poly(x, y) != x + y) is False
    assert ((x + y).as_poly(x) != x + y) is False
    assert ((x + y).as_poly(y) != x + y) is False

    assert (x.as_poly() == sin(x)) is False
    assert (x.as_poly() != sin(x)) is True


def test_Poly_nonzero():
    assert not bool(Integer(0).as_poly(x)) is True
    assert not bool(Integer(1).as_poly(x)) is False


def test_Poly_properties():
    assert Integer(0).as_poly(x).is_zero is True
    assert Integer(1).as_poly(x).is_zero is False

    assert Integer(1).as_poly(x).is_one is True
    assert Integer(2).as_poly(x).is_one is False

    assert (x - 1).as_poly().is_squarefree is True
    assert ((x - 1)**2).as_poly().is_squarefree is False

    assert Integer(1).as_poly(x).is_ground is True
    assert x.as_poly().is_ground is False

    assert (x + y + z + 1).as_poly().is_linear is True
    assert (x*y*z + 1).as_poly().is_linear is False

    assert (x*y + z + 1).as_poly().is_quadratic is True
    assert (x*y*z + 1).as_poly().is_quadratic is False

    assert (x*y).as_poly().is_term is True
    assert (x*y + 1).as_poly().is_term is False

    assert (x**2 + x*y).as_poly().is_homogeneous is True
    assert (x**3 + x*y).as_poly().is_homogeneous is False

    assert x.as_poly().is_univariate is True
    assert (x*y).as_poly().is_univariate is False

    assert (x*y).as_poly().is_multivariate is True
    assert x.as_poly().is_multivariate is False

    assert (x**16 + x**14 - x**10 + x**8 - x**6 +
            x**2 + 1).as_poly().is_cyclotomic is False
    assert (x**16 + x**14 - x**10 - x**8 - x**6 +
            x**2 + 1).as_poly().is_cyclotomic is True


def test_Poly_is_irreducible():
    assert (x**2 + x + 1).as_poly().is_irreducible is True
    assert (x**2 + 2*x + 1).as_poly().is_irreducible is False

    assert (7*x + 3).as_poly(modulus=11).is_irreducible is True
    assert (7*x**2 + 3*x + 1).as_poly(modulus=11).is_irreducible is False


def test_Poly_subs():
    assert (x + 1).as_poly().subs({x: 0}) == 1

    assert (x + 1).as_poly().subs({x: x}) == (x + 1).as_poly()
    assert (x + 1).as_poly().subs({x: y}) == (y + 1).as_poly()

    assert (x*y).as_poly(x).subs({y: x}) == x**2
    assert (x*y).as_poly(x).subs({x: y}) == y**2


def test_Poly_replace():
    assert (x + 1).as_poly().replace(x) == (x + 1).as_poly()
    assert (x + 1).as_poly().replace(y) == (y + 1).as_poly()

    pytest.raises(PolynomialError, lambda: (x + y).as_poly().replace(z))

    assert (x + 1).as_poly().replace(x, x) == (x + 1).as_poly()
    assert (x + 1).as_poly().replace(x, y) == (y + 1).as_poly()

    assert (x + y).as_poly().replace(x, x) == (x + y).as_poly()
    assert (x + y).as_poly().replace(x, z) == (z + y).as_poly(z, y)

    assert (x + y).as_poly().replace(y, y) == (x + y).as_poly()
    assert (x + y).as_poly().replace(y, z) == (x + z).as_poly(x, z)

    pytest.raises(PolynomialError, lambda: (x + y).as_poly().replace(x, y))
    pytest.raises(PolynomialError, lambda: (x + y).as_poly().replace(z, t))

    assert (x + y).as_poly(x).replace(x, z) == (z + y).as_poly(z)
    assert (x + y).as_poly(y).replace(y, z) == (x + z).as_poly(z)

    pytest.raises(PolynomialError, lambda: (x + y).as_poly(x).replace(x, y))
    pytest.raises(PolynomialError, lambda: (x + y).as_poly(y).replace(y, x))


def test_Poly_reorder():
    pytest.raises(PolynomialError, lambda: (x + y).as_poly().reorder(x, z))

    assert (x + y).as_poly().reorder(x, y) == (x + y).as_poly(x, y)
    assert (x + y).as_poly().reorder(y, x) == (x + y).as_poly(y, x)

    assert (x + y).as_poly(y, x).reorder(x, y) == (x + y).as_poly()
    assert (x + y).as_poly(y, x).reorder(y, x) == (x + y).as_poly(y, x)

    assert (x + y).as_poly().reorder(wrt=x) == (x + y).as_poly()
    assert (x + y).as_poly(x, y).reorder(wrt=y) == (x + y).as_poly(y, x)


def test_Poly_has_only_gens():
    assert (x*y + 1).as_poly(x, y, z).has_only_gens(x, y) is True
    assert (x*y + z).as_poly(x, y, z).has_only_gens(x, y) is False

    pytest.raises(GeneratorsError, lambda: (x*y**2 + y**2).as_poly(x, y).has_only_gens(t))


def test_Poly_to_ring():
    assert (2*x + 1).as_poly(domain=ZZ).to_ring() == (2*x + 1).as_poly(domain=ZZ)
    assert (2*x + 1).as_poly(domain=QQ).to_ring() == (2*x + 1).as_poly(domain=ZZ)

    pytest.raises(CoercionFailedError, lambda: (x/2 + 1).as_poly().to_ring())
    pytest.raises(AttributeError, lambda: (2*x + 1).as_poly(modulus=3).to_ring())


def test_Poly_to_field():
    assert (2*x + 1).as_poly(domain=ZZ).to_field() == (2*x + 1).as_poly(domain=QQ)
    assert (2*x + 1).as_poly(domain=QQ).to_field() == (2*x + 1).as_poly(domain=QQ)

    assert (x/2 + 1).as_poly(domain=QQ).to_field() == (x/2 + 1).as_poly(domain=QQ)
    assert (2*x + 1).as_poly(modulus=3).to_field() == (2*x + 1).as_poly(modulus=3)

    assert (2.0*x + 1.0).as_poly().to_field() == (2.0*x + 1.0).as_poly()


def test_Poly_to_exact():
    assert (2*x).as_poly().to_exact() == (2*x).as_poly()
    assert (x/2).as_poly().to_exact() == (x/2).as_poly()

    assert (0.1*x).as_poly().to_exact() == (x/10).as_poly()


def test_Poly_retract():
    f = (x**2 + 1).as_poly(domain=QQ.inject(y))

    assert f.retract() == (x**2 + 1).as_poly(domain=ZZ)
    assert f.retract(field=True) == (x**2 + 1).as_poly(domain=QQ)

    assert Integer(0).as_poly(x, y).retract() == Integer(0).as_poly(x, y)


def test_Poly_slice():
    f = (x**3 + 2*x**2 + 3*x + 4).as_poly()

    assert f.slice(0, 0) == Integer(0).as_poly(x)
    assert f.slice(0, 1) == Integer(4).as_poly(x)
    assert f.slice(0, 2) == (3*x + 4).as_poly()
    assert f.slice(0, 3) == (2*x**2 + 3*x + 4).as_poly()
    assert f.slice(0, 4) == (x**3 + 2*x**2 + 3*x + 4).as_poly()

    assert f.slice(x, 0, 0) == Integer(0).as_poly(x)
    assert f.slice(x, 0, 1) == Integer(4).as_poly(x)
    assert f.slice(x, 0, 2) == (3*x + 4).as_poly()
    assert f.slice(x, 0, 3) == (2*x**2 + 3*x + 4).as_poly()
    assert f.slice(x, 0, 4) == (x**3 + 2*x**2 + 3*x + 4).as_poly()


def test_Poly_coeffs():
    assert Integer(0).as_poly(x).coeffs() == []
    assert Integer(1).as_poly(x).coeffs() == [1]

    assert (2*x + 1).as_poly().coeffs() == [2, 1]

    assert (7*x**2 + 2*x + 1).as_poly().coeffs() == [7, 2, 1]
    assert (7*x**4 + 2*x + 1).as_poly().coeffs() == [7, 2, 1]

    assert (x*y**7 + 2*x**2*y**3).as_poly().coeffs('lex') == [2, 1]
    assert (x*y**7 + 2*x**2*y**3).as_poly().coeffs('grlex') == [1, 2]


def test_Poly_monoms():
    assert Integer(0).as_poly(x).monoms() == []
    assert Integer(1).as_poly(x).monoms() == [(0,)]

    assert (2*x + 1).as_poly().monoms() == [(1,), (0,)]

    assert (7*x**2 + 2*x + 1).as_poly().monoms() == [(2,), (1,), (0,)]
    assert (7*x**4 + 2*x + 1).as_poly().monoms() == [(4,), (1,), (0,)]

    assert (x*y**7 + 2*x**2*y**3).as_poly().monoms('lex') == [(2, 3), (1, 7)]
    assert (x*y**7 + 2*x**2*y**3).as_poly().monoms('grlex') == [(1, 7), (2, 3)]


def test_Poly_terms():
    assert Integer(0).as_poly(x).terms() == []
    assert Integer(1).as_poly(x).terms() == [((0,), 1)]

    assert (2*x + 1).as_poly().terms() == [((1,), 2), ((0,), 1)]

    assert (7*x**2 + 2*x + 1).as_poly().terms() == [((2,), 7), ((1,), 2), ((0,), 1)]

    assert (x*y**7 + 2*x**2*y**3).as_poly().terms('lex') == [((2, 3), 2), ((1, 7), 1)]
    assert (x*y**7 + 2*x**2*y**3).as_poly().terms('grlex') == [((1, 7), 1), ((2, 3), 2)]


def test_Poly_all_coeffs():
    assert Integer(0).as_poly(x).all_coeffs() == [0]
    assert Integer(1).as_poly(x).all_coeffs() == [1]

    assert (2*x + 1).as_poly().all_coeffs() == [1, 2]

    assert (7*x**2 + 2*x + 1).as_poly().all_coeffs() == [1, 2, 7]
    assert (7*x**4 + 2*x + 1).as_poly().all_coeffs() == [1, 2, 0, 0, 7]


def test_Poly_termwise():
    f = (x**2 + 20*x + 400).as_poly()
    g = (x**2 + 2*x + 4).as_poly()

    def func(monom, coeff):
        k, = monom
        return coeff//10**(2 - k)

    assert f.termwise(func) == g

    def func2(monom, coeff):
        k, = monom
        return (k,), coeff//10**(2 - k)

    assert f.termwise(func2) == g

    def func3(monom, coeff):
        k, = monom
        return (k,), coeff // 2

    assert f.termwise(func3) == (10*x + 200).as_poly()

    def func4(monom, coeff):
        k, = monom
        return k % 2, coeff

    pytest.raises(PolynomialError, lambda: f.termwise(func4))


def test_Poly_length():
    assert Integer(0).as_poly(x).length() == 0
    assert Integer(1).as_poly(x).length() == 1
    assert x.as_poly().length() == 1

    assert (x + 1).as_poly().length() == 2
    assert (x**2 + 1).as_poly().length() == 2
    assert (x**2 + x + 1).as_poly().length() == 3


def test_Poly_as_dict():
    assert Integer(0).as_poly(x).as_dict() == {}
    assert Integer(0).as_poly(x, y, z).as_dict() == {}

    assert Integer(1).as_poly(x).as_dict() == {(0,): 1}
    assert Integer(1).as_poly(x, y, z).as_dict() == {(0, 0, 0): 1}

    assert (x**2 + 3).as_poly().as_dict() == {(2,): 1, (0,): 3}
    assert (x**2 + 3).as_poly().as_dict(native=True) == {(2,): ZZ(1), (0,): ZZ(3)}
    assert (x**2 + 3).as_poly(x, y, z).as_dict() == {(2, 0, 0): 1, (0, 0, 0): 3}

    assert (3*x**2*y*z**3 + 4*x*y +
            5*x*z).as_poly().as_dict() == {(2, 1, 3): 3, (1, 1, 0): 4,
                                           (1, 0, 1): 5}


def test_Poly_as_expr():
    assert Integer(0).as_poly(x).as_expr() == 0
    assert Integer(0).as_poly(x, y, z).as_expr() == 0

    assert Integer(1).as_poly(x).as_expr() == 1
    assert Integer(1).as_poly(x, y, z).as_expr() == 1

    assert (x**2 + 3).as_poly().as_expr() == x**2 + 3
    assert (x**2 + 3).as_poly(x, y, z).as_expr() == x**2 + 3

    assert (3*x**2*y*z**3 + 4*x*y +
            5*x*z).as_poly().as_expr() == 3*x**2*y*z**3 + 4*x*y + 5*x*z

    f = (x**2 + 2*x*y**2 - y).as_poly()

    assert f.as_expr() == -y + x**2 + 2*x*y**2

    assert f.as_expr({x: 5}) == 25 - y + 10*y**2
    assert f.as_expr({y: 6}) == -6 + 72*x + x**2

    assert f.as_expr({x: 5, y: 6}) == 379
    assert f.as_expr(5, 6) == 379

    pytest.raises(GeneratorsError, lambda: f.as_expr({z: 7}))


def test_Poly_inject():
    f = (x**2*y + x*y**3 + x*y + 1).as_poly(x)

    assert f.inject() == (x**2*y + x*y**3 + x*y + 1).as_poly()
    assert f.inject(front=True) == (y**3*x + y*x**2 + y*x + 1).as_poly(y, x)

    f = (x**2 + 2*x - 1).as_poly()
    assert f.inject() == f

    f = (x**2 - 2*sqrt(3)*x + 4).as_poly(extension=True)
    assert f.inject().replace(f.domain.ext, y) == (x**2 - 2*x*y + 4).as_poly()


def test_Poly_eject():
    f = (x**2*y + x*y**3 + x*y + 1).as_poly()

    assert f.eject(x) == (x*y**3 + (x**2 + x)*y + 1).as_poly(y, domain=ZZ.inject(x))
    assert f.eject(y) == (y*x**2 + (y**3 + y)*x + 1).as_poly(x, domain=ZZ.inject(y))

    ex = x + y + z + t + w
    g = ex.as_poly()

    assert g.eject(x) == ex.as_poly(y, z, t, w, domain=ZZ.inject(x))
    assert g.eject(x, y) == ex.as_poly(z, t, w, domain=ZZ.inject(x, y))
    assert g.eject(x, y, z) == ex.as_poly(t, w, domain=ZZ.inject(x, y, z))
    assert g.eject(w) == ex.as_poly(x, y, z, t, domain=ZZ.inject(w))
    assert g.eject(t, w) == ex.as_poly(x, y, z, domain=ZZ.inject(w, t))
    assert g.eject(z, t, w) == ex.as_poly(x, y, domain=ZZ.inject(w, t, z))

    pytest.raises(DomainError, lambda: x*y.as_poly(x, y, domain=ZZ.inject(z)).eject(y))

    assert (x*y).as_poly(x, y, z).eject(y) == (x*y).as_poly(x, z, domain=ZZ.inject(y))


def test_Poly_exclude():
    assert x.as_poly(x, y).exclude() == x.as_poly()
    assert (x*y).as_poly(x, y).exclude() == (x*y).as_poly(x, y)
    assert Integer(1).as_poly(x, y).exclude() == Integer(1).as_poly(x, y)
    assert (y**2 + y*z**2).as_poly(x, y, z).exclude() == (y**2 + y*z**2).as_poly(y, z)


def test_Poly__gen_to_level():
    f = Integer(1).as_poly(x, y)

    assert f._gen_to_level(-2) == 0
    assert f._gen_to_level(-1) == 1
    assert f._gen_to_level(+0) == 0
    assert f._gen_to_level(+1) == 1

    pytest.raises(PolynomialError, lambda: f._gen_to_level(-3))
    pytest.raises(PolynomialError, lambda: f._gen_to_level(+2))

    assert f._gen_to_level(x) == 0
    assert f._gen_to_level(y) == 1

    assert f._gen_to_level('x') == 0
    assert f._gen_to_level('y') == 1

    pytest.raises(PolynomialError, lambda: f._gen_to_level(z))
    pytest.raises(PolynomialError, lambda: f._gen_to_level('z'))


def test_Poly_degree():
    assert Integer(0).as_poly(x).degree() == -math.inf
    assert Integer(1).as_poly(x).degree() == 0
    assert x.as_poly().degree() == 1

    assert Integer(0).as_poly(x).degree(gen=0) == -math.inf
    assert Integer(1).as_poly(x).degree(gen=0) == 0
    assert x.as_poly().degree(gen=0) == 1

    assert Integer(0).as_poly(x).degree(gen=x) == -math.inf
    assert Integer(1).as_poly(x).degree(gen=x) == 0
    assert x.as_poly().degree(gen=x) == 1

    assert Integer(0).as_poly(x).degree(gen='x') == -math.inf
    assert Integer(1).as_poly(x).degree(gen='x') == 0
    assert x.as_poly().degree(gen='x') == 1

    f = Integer(1).as_poly(x)

    pytest.raises(PolynomialError, lambda: f.degree(gen=1))
    pytest.raises(PolynomialError, lambda: f.degree(gen=y))
    pytest.raises(PolynomialError, lambda: f.degree(gen='y'))

    assert Integer(1).as_poly(x, y).degree() == 0
    assert (2*y).as_poly(x, y).degree() == 0
    assert (x*y).as_poly(x, y).degree() == 1

    assert Integer(1).as_poly(x, y).degree(gen=x) == 0
    assert (2*y).as_poly(x, y).degree(gen=x) == 0
    assert (x*y).as_poly(x, y).degree(gen=x) == 1

    assert Integer(1).as_poly(x, y).degree(gen=y) == 0
    assert (2*y).as_poly(x, y).degree(gen=y) == 1
    assert (x*y).as_poly(x, y).degree(gen=y) == 1

    assert degree(1, x) == 0
    assert degree(x, x) == 1

    assert degree(x*y**2, gen=x) == 1
    assert degree(x*y**2, gen=y) == 2

    assert degree(x*y**2, x, y) == 1
    assert degree(x*y**2, y, x) == 2

    pytest.raises(ComputationFailedError, lambda: degree(1))

    # issue sympy/sympy#20389
    assert degree(x*(x + 1) - x**2 - x, x) == -oo


@pytest.mark.timeout(30)
def test_sympyissue_6322():
    assert degree((1 + x)**10000) == 10000


def test_Poly_degree_list():
    assert [Integer(0).as_poly(x, y).degree(_) for _ in (x, y)] == [-math.inf]*2
    assert [Integer(0).as_poly(x, y, z).degree(_) for _ in (x, y, z)] == [-math.inf]*3

    assert [Integer(1).as_poly(x, y).degree(_) for _ in (x, y)] == [0, 0]
    assert [Integer(1).as_poly(x, y, z).degree(_) for _ in (x, y, z)] == [0, 0, 0]

    assert [(x**2*y + x**3*z**2 + 1).as_poly().degree(_)
            for _ in (x, y, z)] == [3, 1, 2]

    assert [degree(x*y**2, _) for _ in (x, y)] == [1, 2]
    assert [degree(x**2 + y*x + 1, _) for _ in (x, y)] == [2, 1]


def test_Poly_total_degree():
    assert (x**2*y + x**3*z**2 + 1).as_poly().total_degree() == 5
    assert (x**2 + z**3).as_poly().total_degree() == 3
    assert (x*y*z + z**4).as_poly().total_degree() == 4
    assert (x**3 + x + 1).as_poly().total_degree() == 3


def test_Poly_LC():
    assert Integer(0).as_poly(x).LC() == 0
    assert Integer(1).as_poly(x).LC() == 1
    assert (2*x**2 + x).as_poly().LC() == 2

    assert (x*y**7 + 2*x**2*y**3).as_poly().LC('lex') == 2
    assert (x*y**7 + 2*x**2*y**3).as_poly().LC('grlex') == 1

    assert LC(x*y**7 + 2*x**2*y**3, order='lex') == 2
    assert LC(x*y**7 + 2*x**2*y**3, order='grlex') == 1

    pytest.raises(ComputationFailedError, lambda: LC([1, 2]))


def test_Poly_TC():
    assert Integer(0).as_poly(x).TC() == 0
    assert Integer(1).as_poly(x).TC() == 1
    assert (2*x**2 + x).as_poly().TC() == 0


def test_Poly_EC():
    assert Integer(0).as_poly(x).EC() == 0
    assert Integer(1).as_poly(x).EC() == 1
    assert (2*x**2 + x).as_poly().EC() == 1

    assert (x*y**7 + 2*x**2*y**3).as_poly().EC('lex') == 1
    assert (x*y**7 + 2*x**2*y**3).as_poly().EC('grlex') == 2


def test_Poly_coeff():
    f = Integer(0).as_poly(x)

    assert f.coeff_monomial(1) == 0
    assert f.coeff_monomial(x) == 0
    assert f.coeff_monomial((0,)) == 0
    assert f.coeff_monomial((1,)) == 0

    f = Integer(1).as_poly(x)

    assert f.coeff_monomial(1) == 1
    assert f.coeff_monomial(x) == 0
    assert f.coeff_monomial((0,)) == 1
    assert f.coeff_monomial((1,)) == 0

    f = (x**8).as_poly()

    assert f.coeff_monomial(1) == 0
    assert f.coeff_monomial(x**7) == 0
    assert f.coeff_monomial(x**8) == 1
    assert f.coeff_monomial(x**9) == 0
    assert f.coeff_monomial((0,)) == 0
    assert f.coeff_monomial((7,)) == 0
    assert f.coeff_monomial((8,)) == 1
    assert f.coeff_monomial((9,)) == 0

    f = (3*x*y**2 + 1).as_poly()

    assert f.coeff_monomial(1) == 1
    assert f.coeff_monomial(x*y**2) == 3
    assert f.coeff_monomial((0, 0)) == 1
    assert f.coeff_monomial((1, 2)) == 3

    p = (24*x*y*exp(8) + 23*x).as_poly(x, y)

    assert p.coeff_monomial(x) == 23
    assert p.coeff_monomial(y) == 0
    assert p.coeff_monomial(x*y) == 24*exp(8)

    assert p.as_expr().coeff(x) == 24*y*exp(8) + 23

    pytest.raises(NotImplementedError, lambda: p.coeff(x))

    f = (x + 1).as_poly()

    pytest.raises(ValueError, lambda: f.coeff_monomial(0))
    pytest.raises(ValueError, lambda: f.coeff_monomial(3*x))
    pytest.raises(ValueError, lambda: f.coeff_monomial(3*x*y))

    pytest.raises(ValueError, lambda: (x*y + 1).as_poly().coeff_monomial((1,)))

    assert (x**3 + 2*x**2 + 3*x).as_poly().coeff_monomial((2,)) == 2
    assert (x**3 + 2*x*y**2 + y**2).as_poly().coeff_monomial((1, 2)) == 2
    assert (4*sqrt(x)*y).as_poly().coeff_monomial((1, 1)) == 4


def test_Poly_LM():
    assert Integer(0).as_poly(x).LM() == (0,)
    assert Integer(1).as_poly(x).LM() == (0,)
    assert (2*x**2 + x).as_poly().LM() == (2,)

    assert (x*y**7 + 2*x**2*y**3).as_poly().LM('lex') == (2, 3)
    assert (x*y**7 + 2*x**2*y**3).as_poly().LM('grlex') == (1, 7)

    assert LM(x*y**7 + 2*x**2*y**3, order='lex') == x**2*y**3
    assert LM(x*y**7 + 2*x**2*y**3, order='grlex') == x*y**7

    pytest.raises(ComputationFailedError, lambda: LM([1, 2]))


def test_Poly_LM_custom_order():
    f = (x**2*y**3*z + x**2*y*z**3 + x*y*z + 1).as_poly()

    def rev_lex(monom):
        return tuple(reversed(monom))

    assert f.LM(order='lex') == (2, 3, 1)
    assert f.LM(order=rev_lex) == (2, 1, 3)


def test_Poly_EM():
    assert Integer(0).as_poly(x).EM() == (0,)
    assert Integer(1).as_poly(x).EM() == (0,)
    assert (2*x**2 + x).as_poly().EM() == (1,)

    assert (x*y**7 + 2*x**2*y**3).as_poly().EM('lex') == (1, 7)
    assert (x*y**7 + 2*x**2*y**3).as_poly().EM('grlex') == (2, 3)


def test_Poly_LT():
    assert Integer(0).as_poly(x).LT() == ((0,), 0)
    assert Integer(1).as_poly(x).LT() == ((0,), 1)
    assert (2*x**2 + x).as_poly().LT() == ((2,), 2)

    assert (x*y**7 + 2*x**2*y**3).as_poly().LT('lex') == ((2, 3), 2)
    assert (x*y**7 + 2*x**2*y**3).as_poly().LT('grlex') == ((1, 7), 1)

    assert LT(x*y**7 + 2*x**2*y**3, order='lex') == 2*x**2*y**3
    assert LT(x*y**7 + 2*x**2*y**3, order='grlex') == x*y**7

    pytest.raises(ComputationFailedError, lambda: LT([1, 2]))


def test_Poly_ET():
    assert Integer(0).as_poly(x).ET() == ((0,), 0)
    assert Integer(1).as_poly(x).ET() == ((0,), 1)
    assert (2*x**2 + x).as_poly().ET() == ((1,), 1)

    assert (x*y**7 + 2*x**2*y**3).as_poly().ET('lex') == ((1, 7), 1)
    assert (x*y**7 + 2*x**2*y**3).as_poly().ET('grlex') == ((2, 3), 2)


def test_Poly_clear_denoms():
    coeff, poly = (x + 2).as_poly().clear_denoms()
    assert coeff == 1
    assert poly == (x + 2).as_poly(domain=ZZ)
    assert poly.domain == ZZ

    coeff, poly = (x/2 + 1).as_poly().clear_denoms()
    assert coeff == 2
    assert poly == (x + 2).as_poly(domain=QQ)
    assert poly.domain == QQ

    coeff, poly = (x/2 + 1).as_poly().clear_denoms(convert=True)
    assert coeff == 2
    assert poly == (x + 2).as_poly(domain=ZZ)
    assert poly.domain == ZZ

    coeff, poly = (x/y + 1).as_poly(x).clear_denoms(convert=True)
    assert coeff == y
    assert poly == (x + y).as_poly(x, domain=ZZ.inject(y))
    assert poly.domain == ZZ.inject(y)

    coeff, poly = (x/3 + sqrt(2)).as_poly(x, domain=EX).clear_denoms()
    assert coeff == 3
    assert poly == (x + 3*sqrt(2)).as_poly(x, domain=EX)
    assert poly.domain == EX

    coeff, poly = (x/3 + sqrt(2)).as_poly(x, domain=EX).clear_denoms(convert=True)
    assert coeff == 3
    assert poly == (x + 3*sqrt(2)).as_poly(x, domain=EX)
    assert poly.domain == EX


def test_Poly_rat_clear_denoms():
    f = (x**2/y + 1).as_poly(x)
    g = (x**3 + y).as_poly(x)

    assert f.rat_clear_denoms(g) == \
        ((x**2 + y).as_poly(x), (y*x**3 + y**2).as_poly(x))

    f = f.set_domain(EX)
    g = g.set_domain(EX)

    assert f.rat_clear_denoms(g) == (f, g)


def test_Poly_integrate():
    f = (x + 1).as_poly()

    assert f.integrate() == (x**2/2 + x).as_poly()
    assert f.integrate(x) == (x**2/2 + x).as_poly()
    assert f.integrate((x, 1)) == (x**2/2 + x).as_poly()

    assert (2*x + 1).as_poly().integrate(auto=False) == (x**2 + x).as_poly()

    f = (x*y + 1).as_poly()

    assert f.integrate(x) == (x**2*y/2 + x).as_poly()
    assert f.integrate(y) == (x*y**2/2 + y).as_poly()

    assert f.integrate(x, x) == (x**3*y/6 + x**2/2).as_poly()
    assert f.integrate(y, y) == (x*y**3/6 + y**2/2).as_poly()

    assert f.integrate((x, 2)) == (x**3*y/6 + x**2/2).as_poly()
    assert f.integrate((y, 2)) == (x*y**3/6 + y**2/2).as_poly()

    assert f.integrate(x, y) == (x**2*y**2/4 + x*y).as_poly()
    assert f.integrate(y, x) == (x**2*y**2/4 + x*y).as_poly()


def test_Poly_diff():
    f = (x**2 + x).as_poly()

    assert f.diff() == (2*x + 1).as_poly()
    assert f.diff(x) == (2*x + 1).as_poly()
    assert f.diff((x, 1)) == (2*x + 1).as_poly()

    # issue sympy/sympy#9585
    assert diff(f) == (2*x + 1).as_poly()
    assert diff(f, x, evaluate=False) == Derivative(f, x)
    assert Derivative(f, x).doit() == (2*x + 1).as_poly()
    assert f.diff(x, evaluate=False) == Derivative(f, x)

    f = (x**2 + 2*x + 1).as_poly()

    assert f.diff() == (2*x + 2).as_poly()

    f = (x**2*y**2 + x*y).as_poly()

    assert f.diff(x) == (2*x*y**2 + y).as_poly()
    assert f.diff(y) == (2*x**2*y + x).as_poly()

    assert f.diff(x, x) == (2*y**2).as_poly(x, y)
    assert f.diff(y, y) == (2*x**2).as_poly(x, y)

    assert f.diff((x, 2)) == (2*y**2).as_poly(x, y)
    assert f.diff((y, 2)) == (2*x**2).as_poly(x, y)

    assert f.diff(x, y) == (4*x*y + 1).as_poly()
    assert f.diff(y, x) == (4*x*y + 1).as_poly()

    f = (x*y**2 + x).as_poly()

    assert f.diff((x, 0), (y, 1)) == (2*x*y).as_poly()


def test_Poly_eval():
    f = Integer(0).as_poly(x)

    assert f.eval(7) == 0
    assert f.eval(0, 7) == 0
    assert f.eval(x, 7) == 0
    assert f.eval('x', 7) == 0

    f = Integer(1).as_poly(x)

    assert f.eval(7) == 1
    assert f.eval(0, 7) == 1
    assert f.eval(x, 7) == 1
    assert f.eval('x', 7) == 1

    pytest.raises(PolynomialError, lambda: f.eval(1, 7))
    pytest.raises(PolynomialError, lambda: f.eval(y, 7))
    pytest.raises(PolynomialError, lambda: f.eval('y', 7))

    f = x.as_poly()

    assert f.eval(7) == 7
    assert f.eval(0, 7) == 7
    assert f.eval(x, 7) == 7
    assert f.eval('x', 7) == 7

    f = Integer(123).as_poly(x, y)

    assert f.eval(7) == Integer(123).as_poly(y)
    assert f.eval(x, 7) == Integer(123).as_poly(y)
    assert f.eval(y, 7) == Integer(123).as_poly(x)

    f = (2*y).as_poly(x, y)

    assert f.eval(7) == (2*y).as_poly()
    assert f.eval(x, 7) == (2*y).as_poly()
    assert f.eval(y, 7) == Integer(14).as_poly(x)

    f = (x*y).as_poly()

    assert f.eval(7) == (7*y).as_poly()
    assert f.eval(x, 7) == (7*y).as_poly()
    assert f.eval(y, 7) == (7*x).as_poly()

    f = (x*y + y).as_poly()

    assert f.eval({x: 7}) == (8*y).as_poly()
    assert f.eval({y: 7}) == (7*x + 7).as_poly()

    assert f.eval({x: 6, y: 7}) == 49
    assert f.eval({x: 7, y: 6}) == 48

    assert f.eval((6, 7)) == 49
    assert f.eval([6, 7]) == 49

    pytest.raises(ValueError, lambda: f.eval((6, 7, 8)))

    f = (x + 1).as_poly()

    assert f.eval(Rational(1, 2)) == Rational(3, 2)
    assert f.eval(sqrt(2)) == sqrt(2) + 1

    pytest.raises(DomainError, lambda: f.eval(Rational(1, 2), auto=False))

    # issue sympy/sympy#6344
    alpha = Symbol('alpha')
    result = (2*alpha*z - 2*alpha + z**2 + 3)/(z**2 - 2*z + 1)

    f = (x**2 + (alpha - 1)*x - alpha + 1).as_poly(x, domain=ZZ.inject(alpha))

    assert f.eval((z + 1)/(z - 1)) == result

    f = (x**2 + (alpha - 1)*x - alpha + 1).as_poly(x, y, domain=ZZ.inject(alpha))

    assert f.eval((z + 1)/(z - 1)) == result.as_poly(y, domain=ZZ.inject(alpha, z).field)


def test_Poly___call__():
    f = (2*x*y + 3*x + y + 2*z).as_poly()

    assert f(2) == (5*y + 2*z + 6).as_poly()
    assert f(2, 5) == (2*z + 31).as_poly()
    assert f(2, 5, 7) == 45


def test_parallel_poly_from_expr():
    pytest.raises(PolificationFailedError, lambda: parallel_poly_from_expr([]))
    pytest.raises(PolificationFailedError, lambda: parallel_poly_from_expr([[1, 2]]))

    assert parallel_poly_from_expr(
        [x - 1, x**2 - 1], x)[0] == [(x - 1).as_poly(), (x**2 - 1).as_poly()]
    assert parallel_poly_from_expr(
        [(x - 1).as_poly(), x**2 - 1], x)[0] == [(x - 1).as_poly(), (x**2 - 1).as_poly()]
    assert parallel_poly_from_expr(
        [x - 1, (x**2 - 1).as_poly()], x)[0] == [(x - 1).as_poly(), (x**2 - 1).as_poly()]
    assert parallel_poly_from_expr([(
        x - 1).as_poly(), (x**2 - 1).as_poly()], x)[0] == [(x - 1).as_poly(), (x**2 - 1).as_poly()]

    assert parallel_poly_from_expr(
        [x - 1, x**2 - 1], x, y)[0] == [(x - 1).as_poly(x, y), (x**2 - 1).as_poly(x, y)]
    assert parallel_poly_from_expr([(
        x - 1).as_poly(), x**2 - 1], x, y)[0] == [(x - 1).as_poly(x, y), (x**2 - 1).as_poly(x, y)]
    assert parallel_poly_from_expr([x - 1, (
        x**2 - 1).as_poly()], x, y)[0] == [(x - 1).as_poly(x, y), (x**2 - 1).as_poly(x, y)]
    assert parallel_poly_from_expr([(x - 1).as_poly(), (
        x**2 - 1).as_poly()], x, y)[0] == [(x - 1).as_poly(x, y), (x**2 - 1).as_poly(x, y)]

    assert parallel_poly_from_expr(
        [x - 1, x**2 - 1])[0] == [(x - 1).as_poly(), (x**2 - 1).as_poly()]
    assert parallel_poly_from_expr(
        [(x - 1).as_poly(), x**2 - 1])[0] == [(x - 1).as_poly(), (x**2 - 1).as_poly()]
    assert parallel_poly_from_expr(
        [x - 1, (x**2 - 1).as_poly()])[0] == [(x - 1).as_poly(), (x**2 - 1).as_poly()]
    assert parallel_poly_from_expr(
        [(x - 1).as_poly(), (x**2 - 1).as_poly()])[0] == [(x - 1).as_poly(), (x**2 - 1).as_poly()]

    assert parallel_poly_from_expr(
        [1, x**2 - 1])[0] == [Integer(1).as_poly(x), (x**2 - 1).as_poly()]
    assert parallel_poly_from_expr(
        [1, x**2 - 1])[0] == [Integer(1).as_poly(x), (x**2 - 1).as_poly()]
    assert parallel_poly_from_expr(
        [1, (x**2 - 1).as_poly()])[0] == [Integer(1).as_poly(x), (x**2 - 1).as_poly()]
    assert parallel_poly_from_expr(
        [1, (x**2 - 1).as_poly()])[0] == [Integer(1).as_poly(x), (x**2 - 1).as_poly()]

    assert parallel_poly_from_expr(
        [x**2 - 1, 1])[0] == [(x**2 - 1).as_poly(), Integer(1).as_poly(x)]
    assert parallel_poly_from_expr(
        [x**2 - 1, 1])[0] == [(x**2 - 1).as_poly(), Integer(1).as_poly(x)]
    assert parallel_poly_from_expr(
        [(x**2 - 1).as_poly(), 1])[0] == [(x**2 - 1).as_poly(), Integer(1).as_poly(x)]
    assert parallel_poly_from_expr(
        [(x**2 - 1).as_poly(), 1])[0] == [(x**2 - 1).as_poly(), Integer(1).as_poly(x)]

    assert parallel_poly_from_expr([x.as_poly(x, y), y.as_poly(x, y)], x, y, order='lex')[0] == \
        [x.as_poly(x, y, domain=ZZ), y.as_poly(x, y, domain=ZZ)]

    pytest.raises(PolificationFailedError, lambda: parallel_poly_from_expr([0, 1]))

    assert (parallel_poly_from_expr([(x - 1)**2, 1], expand=False) ==
            ([((x - 1)**2).as_poly(x - 1, expand=False), Integer(1).as_poly(x - 1)],
             {'domain': ZZ, 'expand': False, 'gens': (x - 1,),
              'polys': False}))


def test_div():
    f, g = x**2 - y**2, x - y
    q, r = x + y, Integer(0)

    F, G, Q, R = (h.as_poly(x, y) for h in (f, g, q, r))

    assert F.div(G) == (Q, R)
    assert F.rem(G) == R
    assert F.quo(G) == Q
    assert F.exquo(G) == Q

    assert div(f, g) == (q, r)
    assert rem(f, g) == r
    assert quo(f, g) == q
    assert exquo(f, g) == q

    assert div(f, g, x, y) == (q, r)
    assert rem(f, g, x, y) == r
    assert quo(f, g, x, y) == q
    assert exquo(f, g, x, y) == q

    assert div(F, G) == (Q, R)
    assert rem(F, G) == R
    assert quo(F, G) == Q
    assert exquo(F, G) == Q

    assert div(f, g, polys=True) == (Q, R)
    assert rem(f, g, polys=True) == R
    assert quo(f, g, polys=True) == Q
    assert exquo(f, g, polys=True) == Q

    assert div(F, G, polys=False) == (q, r)
    assert rem(F, G, polys=False) == r
    assert quo(F, G, polys=False) == q
    assert exquo(F, G, polys=False) == q

    pytest.raises(ComputationFailedError, lambda: div(4, 2))
    pytest.raises(ComputationFailedError, lambda: rem(4, 2))
    pytest.raises(ComputationFailedError, lambda: quo(4, 2))
    pytest.raises(ComputationFailedError, lambda: exquo(4, 2))

    f, g = x**2 + 1, 2*x - 4

    qz, rz = 0, x**2 + 1
    qq, rq = x/2 + 1, 5

    assert div(f, g) == (qq, rq)
    assert div(f, g, auto=True) == (qq, rq)
    assert div(f, g, auto=False) == (qz, rz)
    assert div(f, g, domain=ZZ) == (qz, rz)
    assert div(f, g, domain=QQ) == (qq, rq)
    assert div(f, g, domain=ZZ, auto=True) == (qq, rq)
    assert div(f, g, domain=ZZ, auto=False) == (qz, rz)
    assert div(f, g, domain=QQ, auto=True) == (qq, rq)
    assert div(f, g, domain=QQ, auto=False) == (qq, rq)

    assert rem(f, g) == rq
    assert rem(f, g, auto=True) == rq
    assert rem(f, g, auto=False) == rz
    assert rem(f, g, domain=ZZ) == rz
    assert rem(f, g, domain=QQ) == rq
    assert rem(f, g, domain=ZZ, auto=True) == rq
    assert rem(f, g, domain=ZZ, auto=False) == rz
    assert rem(f, g, domain=QQ, auto=True) == rq
    assert rem(f, g, domain=QQ, auto=False) == rq

    assert quo(f, g) == qq
    assert quo(f, g, auto=True) == qq
    assert quo(f, g, auto=False) == qz
    assert quo(f, g, domain=ZZ) == qz
    assert quo(f, g, domain=QQ) == qq
    assert quo(f, g, domain=ZZ, auto=True) == qq
    assert quo(f, g, domain=ZZ, auto=False) == qz
    assert quo(f, g, domain=QQ, auto=True) == qq
    assert quo(f, g, domain=QQ, auto=False) == qq

    f, g, q = x**2, 2*x, x/2

    assert exquo(f, g) == q
    assert exquo(f, g, auto=True) == q
    pytest.raises(ExactQuotientFailedError, lambda: exquo(f, g, auto=False))
    pytest.raises(ExactQuotientFailedError, lambda: exquo(f, g, domain=ZZ))
    assert exquo(f, g, domain=QQ) == q
    assert exquo(f, g, domain=ZZ, auto=True) == q
    pytest.raises(ExactQuotientFailedError, lambda: exquo(f, g, domain=ZZ, auto=False))
    assert exquo(f, g, domain=QQ, auto=True) == q
    assert exquo(f, g, domain=QQ, auto=False) == q

    f, g = (x**2).as_poly(), x.as_poly()

    q, r = f.div(g)
    assert q.domain.is_IntegerRing
    assert r.domain.is_IntegerRing
    r = f.rem(g)
    assert r.domain.is_IntegerRing
    q = f.quo(g)
    assert q.domain.is_IntegerRing
    q = f.exquo(g)
    assert q.domain.is_IntegerRing

    f, g, q = x**2 + 1, 2*x - 9, QQ(85, 4)
    assert rem(f, g) == q

    f, g = a*x**2 + b*x + c, 3*x + 2
    assert div(f, g) == (a*x/3 - 2*a/9 + b/3, 4*a/9 - 2*b/3 + c)

    f, g = (x + y).as_poly(x), (2*x + y).as_poly(x)
    q, r = f.div(g)
    assert q.domain.is_FractionField
    assert r.domain.is_FractionField


def test_gcdex():
    f, g = 2*x, x**2 - 16
    s, t, h = x/32, -Rational(1, 16), Integer(1)

    F, G, S, T, H = (u.as_poly(x, domain=QQ) for u in (f, g, s, t, h))

    assert F.half_gcdex(G) == (S, H)
    assert F.gcdex(G) == (S, T, H)
    assert F.invert(G) == S

    assert half_gcdex(f, g) == (s, h)
    assert gcdex(f, g) == (s, t, h)
    assert invert(f, g) == s

    assert half_gcdex(f, g, x) == (s, h)
    assert gcdex(f, g, x) == (s, t, h)
    assert invert(f, g, x) == s

    assert half_gcdex(F, G) == (S, H)
    assert gcdex(F, G) == (S, T, H)
    assert invert(F, G) == S

    assert half_gcdex(f, g, polys=True) == (S, H)
    assert gcdex(f, g, polys=True) == (S, T, H)
    assert invert(f, g, polys=True) == S

    assert half_gcdex(F, G, polys=False) == (s, h)
    assert gcdex(F, G, polys=False) == (s, t, h)
    assert invert(F, G, polys=False) == s

    assert half_gcdex(100, 2004) == (-20, 4)
    assert gcdex(100, 2004) == (-20, 1, 4)
    assert invert(3, 7) == 5

    pytest.raises(DomainError, lambda: half_gcdex(x + 1, 2*x + 1, auto=False))
    pytest.raises(DomainError, lambda: gcdex(x + 1, 2*x + 1, auto=False))
    pytest.raises(DomainError, lambda: invert(x + 1, 2*x + 1, auto=False))


def test_subresultants():
    f, g, h = x**2 - 2*x + 1, x**2 - 1, 2*x - 2
    F, G, H = f.as_poly(), g.as_poly(), h.as_poly()

    assert F.subresultants(G) == [F, G, H]
    assert subresultants(f, g) == [f, g, h]
    assert subresultants(f, g, x) == [f, g, h]
    assert subresultants(F, G) == [F, G, H]
    assert subresultants(f, g, polys=True) == [F, G, H]
    assert subresultants(F, G, polys=False) == [f, g, h]

    pytest.raises(ComputationFailedError, lambda: subresultants(4, 2))


def test_resultant():
    f, g, h = x**2 - 2*x + 1, x**2 - 1, Integer(0)
    F, G = f.as_poly(), g.as_poly()

    assert F.resultant(G) == h
    assert PurePoly(f).resultant(PurePoly(g)) == h
    assert resultant(f, g) == h
    assert resultant(f, g, x) == h
    assert resultant(F, G) == h
    assert resultant(f, g, polys=True) == h
    assert resultant(F, G, polys=False) == h
    assert resultant(f, g, includePRS=True) == (h, [f, g, 2*x - 2])
    assert resultant(f, g, polys=True,
                     includePRS=True) == (0, [f, g, (2*x - 2).as_poly()])

    f, g, h = x - a, x - b, a - b
    F, G, H = f.as_poly(), g.as_poly(), h.as_poly()

    assert F.resultant(G) == H
    assert resultant(f, g) == h
    assert resultant(f, g, x) == h
    assert resultant(F, G) == H
    assert resultant(f, g, polys=True) == H
    assert resultant(F, G, polys=False) == h

    pytest.raises(ComputationFailedError, lambda: resultant(4, 2))

    assert resultant(PurePoly(x**2 + y),
                     PurePoly(x*y - 1)) == PurePoly(y**3 + 1)

    f = (x**4 - x**2 + 1).as_poly()

    assert f.resultant((x**2 - y).as_poly()) == ((y**2 - y + 1)**2).as_poly()
    assert f.resultant((x**3 - y).as_poly()) == ((y**2 + 1)**2).as_poly()
    assert f.resultant((x**4 - y).as_poly()) == ((y**2 + y + 1)**2).as_poly()
    assert f.resultant((x**12 - y).as_poly()) == ((y - 1)**4).as_poly()


def test_discriminant():
    f, g = x**3 + 3*x**2 + 9*x - 13, -11664
    F = f.as_poly()

    assert F.discriminant() == g
    assert discriminant(f) == g
    assert discriminant(f, x) == g
    assert discriminant(F) == g
    assert discriminant(f, polys=True) == g
    assert discriminant(F, polys=False) == g

    f, g = a*x**2 + b*x + c, b**2 - 4*a*c
    F, G = f.as_poly(), g.as_poly()

    assert F.discriminant() == G
    assert discriminant(f) == g
    assert discriminant(f, x, a, b, c) == g
    assert discriminant(F) == G
    assert discriminant(f, polys=True) == G
    assert discriminant(F, polys=False) == g

    pytest.raises(ComputationFailedError, lambda: discriminant(4))


def test_dispersion():
    pytest.raises(AttributeError, lambda: (x*y).as_poly().dispersionset(x.as_poly()))
    pytest.raises(ValueError, lambda: x.as_poly().dispersionset(y.as_poly()))

    fp = Integer(0).as_poly(x)

    assert fp.dispersionset() == {0}

    fp = Integer(2).as_poly(x)

    assert fp.dispersionset() == {0}

    fp = (x + 1).as_poly()

    assert fp.dispersionset() == {0}

    fp = ((x + 1)*(x + 2)).as_poly()

    assert fp.dispersionset() == {0, 1}

    fp = (x**4 - 3*x**2 + 1).as_poly()
    gp = fp.shift(-3)

    assert fp.dispersionset(gp) == {2, 3, 4}
    assert gp.dispersionset(fp) == set()

    fp = (x*(x + 3)).as_poly()

    assert fp.dispersionset() == {0, 3}

    fp = ((x - 3)*(x + 3)).as_poly()

    assert fp.dispersionset() == {0, 6}

    fp = (x**2 + 2*x - 1).as_poly()
    gp = (x**2 + 2*x + 3).as_poly()

    assert fp.dispersionset(gp) == set()

    fp = (x*(3*x**2 + a)*(x - 2536)*(x**3 + a)).as_poly(x)
    gp = fp.as_expr().subs({x: x - 345}).as_poly(x)

    assert fp.dispersionset(gp) == {345, 2881}
    assert gp.dispersionset(fp) == {2191}

    fp = ((x - 2)**2*(x - 3)**3*(x - 5)**3).as_poly()
    gp = (fp + 4)**2

    assert fp.dispersionset() == {0, 1, 2, 3}
    assert fp.dispersionset(gp) == {1, 2}

    fp = (x*(x + 2)*(x - 1)).as_poly()

    assert fp.dispersionset() == {0, 1, 2, 3}

    fp = (x**2 + sqrt(5)*x - 1).as_poly(x)
    gp = (x**2 + (2 + sqrt(5))*x + sqrt(5)).as_poly(x)

    assert fp.dispersionset(gp) == {2}
    assert gp.dispersionset(fp) == {1, 4}

    # There are some difficulties if we compute over Z[a]
    # and alpha happenes to lie in Z[a] instead of simply Z.
    # Hence we can not decide if alpha is indeed integral
    # in general.

    fp = (4*x**4 + (4*a + 8)*x**3 + (a**2 + 6*a + 4)*x**2 + (a**2 + 2*a)*x).as_poly(x)

    assert fp.dispersionset() == {0, 1}

    # For any specific value of a, the dispersion is 3*a
    # but the algorithm can not find this in general.
    # This is the point where the resultant based Ansatz
    # is superior to the current one.
    fp = (a**2*x**3 + (a**3 + a**2 + a + 1)*x).as_poly(x)
    gp = fp.as_expr().subs({x: x - 3*a}).as_poly(x)

    assert fp.dispersionset(gp) == set()

    fpa = fp.as_expr().subs({a: 2}).as_poly(x)
    gpa = gp.as_expr().subs({a: 2}).as_poly(x)

    assert fpa.dispersionset(gpa) == {6}


def test_lcm():
    F = [x**3 - 1, x**2 - 1, x**2 - 3*x + 2]

    assert functools.reduce(lcm, F) == x**5 - x**4 - 2*x**3 - x**2 + x + 2
    assert functools.reduce(lambda x, y: lcm(x, y, polys=True),
                            F) == (x**5 - x**4 - 2*x**3 - x**2 + x + 2).as_poly()

    assert lcm(1, 2) == 2
    assert functools.reduce(lcm, [4, 6, 8]) == 24


def test_gcd():
    f, g = x**3 - 1, x**2 - 1
    s, t = x**2 + x + 1, x + 1
    h, r = x - 1, x**4 + x**3 - x - 1

    F, G, S, T, H, R = (u.as_poly() for u in (f, g, s, t, h, r))

    assert F.cofactors(G) == (H, S, T)
    assert F.gcd(G) == H
    assert F.lcm(G) == R

    assert cofactors(f, g) == (h, s, t)
    assert gcd(f, g) == h
    assert lcm(f, g) == r

    assert cofactors(f, g, x) == (h, s, t)
    assert gcd(f, g, x) == h
    assert lcm(f, g, x) == r

    assert cofactors(F, G) == (H, S, T)
    assert gcd(F, G) == H
    assert lcm(F, G) == R

    assert cofactors(f, g, polys=True) == (H, S, T)
    assert gcd(f, g, polys=True) == H
    assert lcm(f, g, polys=True) == R

    assert cofactors(F, G, polys=False) == (h, s, t)
    assert gcd(F, G, polys=False) == h
    assert lcm(F, G, polys=False) == r

    f, g = 1.0*x**2 - 1.0, 1.0*x - 1.0
    h, s, t = g, 1.0*x + 1.0, 1.0

    assert cofactors(f, g) == (h, s, t)
    assert gcd(f, g) == h
    assert lcm(f, g) == f

    f, g = 1.0*x**2 - 1.0, 1.0*x - 1.0
    h, s, t = g, 1.0*x + 1.0, 1.0

    assert cofactors(f, g) == (h, s, t)
    assert gcd(f, g) == h
    assert lcm(f, g) == f

    assert cofactors(8, 6) == (2, 4, 3)
    assert gcd(8, 6) == 2
    assert lcm(8, 6) == 24

    f, g = x**2 + 8*x + 7, x**3 + 7*x**2 + x + 7
    l = x**4 + 8*x**3 + 8*x**2 + 8*x + 7
    h, s, t = x + 7, x + 1, x**2 + 1

    assert cofactors(f, g, modulus=11) == (h, s, t)
    assert gcd(f, g, modulus=11) == h
    assert lcm(f, g, modulus=11) == l

    pytest.raises(TypeError, lambda: gcd(x))
    pytest.raises(TypeError, lambda: lcm(x))

    F = [x**3 - 1, x**2 - 1, x**2 - 3*x + 2]

    assert functools.reduce(gcd, F) == x - 1
    assert functools.reduce(lambda x, y: gcd(x, y, polys=True), F) == (x - 1).as_poly()

    F = [x**3 - 1, x**2 - 2, x**2 - 3*x + 2]

    assert functools.reduce(gcd, F) == 1


def test_gcd_numbers_vs_polys():
    assert isinstance(gcd(3, 9), Integer)
    assert isinstance(gcd(3*x, 9), Integer)

    assert gcd(3, 9) == 3
    assert gcd(3*x, 9) == 3

    assert isinstance(gcd(Rational(3, 2), Rational(9, 4)), Rational)
    assert isinstance(gcd(Rational(3, 2)*x, Rational(9, 4)), Rational)

    assert gcd(Rational(3, 2), Rational(9, 4)) == Rational(3, 4)
    assert gcd(Rational(3, 2)*x, Rational(9, 4)) == 1

    assert isinstance(gcd(3.0, 9.0), Float)
    assert isinstance(gcd(3.0*x, 9.0), Float)

    assert gcd(3.0, 9.0) == 1.0
    assert gcd(3.0*x, 9.0) == 1.0


def test_terms_gcd():
    assert terms_gcd(1) == 1
    assert terms_gcd(1, x) == 1

    assert terms_gcd(x - 1) == x - 1
    assert terms_gcd(-x - 1) == -x - 1

    assert terms_gcd(2*x + 3) == 2*x + 3
    assert terms_gcd(6*x + 4) == Mul(2, 3*x + 2, evaluate=False)

    assert terms_gcd(x**3*y + x*y**3) == x*y*(x**2 + y**2)
    assert terms_gcd(2*x**3*y + 2*x*y**3) == 2*x*y*(x**2 + y**2)
    assert terms_gcd(x**3*y/2 + x*y**3/2) == x*y/2*(x**2 + y**2)

    assert terms_gcd(x**3*y + 2*x*y**3) == x*y*(x**2 + 2*y**2)
    assert terms_gcd(2*x**3*y + 4*x*y**3) == 2*x*y*(x**2 + 2*y**2)
    assert terms_gcd(2*x**3*y/3 + 4*x*y**3/5) == 2*x*y/15*(5*x**2 + 6*y**2)

    assert terms_gcd(2.0*x**3*y + 4.1*x*y**3) == x*y*(2.0*x**2 + 4.1*y**2)
    assert terms_gcd(2.0*x + 3) == 2.0*x + 3

    assert terms_gcd((3 + 3*x)*(x + x*y), expand=False) == \
        (3*x + 3)*(x*y + x)
    assert terms_gcd((3 + 3*x)*(x + x*sin(3 + 3*y)), expand=False, deep=True) == \
        3*x*(x + 1)*(sin(Mul(3, y + 1, evaluate=False)) + 1)
    assert terms_gcd(sin(x + x*y), deep=True) == \
        sin(x*(y + 1))

    eq = Eq(2*x, 2*y + 2*z*y)
    assert terms_gcd(eq) == eq
    assert terms_gcd(eq, deep=True) == Eq(2*x, 2*y*(z + 1))


def test_trunc():
    f, g = x**5 + 2*x**4 + 3*x**3 + 4*x**2 + 5*x + 6, x**5 - x**4 + x**2 - x
    F, G = f.as_poly(), g.as_poly()

    assert F.trunc(3) == G
    assert trunc(f, 3) == g
    assert trunc(f, 3, x) == g
    assert trunc(F, 3) == G
    assert trunc(f, 3, polys=True) == G
    assert trunc(F, 3, polys=False) == g

    f, g = 6*x**5 + 5*x**4 + 4*x**3 + 3*x**2 + 2*x + 1, -x**4 + x**3 - x + 1
    F, G = f.as_poly(), g.as_poly()

    assert F.trunc(3) == G
    assert trunc(f, 3) == g
    assert trunc(f, 3, x) == g
    assert trunc(F, 3) == G
    assert trunc(f, 3, polys=True) == G
    assert trunc(F, 3, polys=False) == g

    f = (x**2 + 2*x + 3).as_poly(modulus=5)

    assert f.trunc(2) == (x**2 + 1).as_poly(modulus=5)

    pytest.raises(ComputationFailedError, lambda: trunc([1, 2], 2))


def test_monic():
    f, g = 2*x - 1, x - Rational(1, 2)
    F, G = f.as_poly(domain=QQ), g.as_poly()

    assert F.monic() == G
    assert monic(f) == g
    assert monic(f, x) == g
    assert monic(F) == G
    assert monic(f, polys=True) == G
    assert monic(F, polys=False) == g

    pytest.raises(ComputationFailedError, lambda: monic(4))

    assert monic(2*x**2 + 6*x + 4, auto=False) == x**2 + 3*x + 2
    pytest.raises(ExactQuotientFailedError, lambda: monic(2*x + 6*x + 1, auto=False))

    assert monic(2.0*x**2 + 6.0*x + 4.0) == 1.0*x**2 + 3.0*x + 2.0
    assert monic(2*x**2 + 3*x + 4, modulus=5) == x**2 + 4*x + 2

    assert monic(x + 2) == x + 2
    assert monic(2*x + 2) == x + 1

    assert monic(x - 1) == x - 1
    assert monic(2*x - 1) == x - Rational(1, 2)


def test_content():
    f, F = 4*x + 2, (4*x + 2).as_poly()

    assert F.content() == 2
    assert content(f) == 2

    pytest.raises(ComputationFailedError, lambda: content(4))

    f = (2*x).as_poly(modulus=3)

    assert f.content() == 1


def test_primitive():
    f, g = 4*x + 2, 2*x + 1
    F, G = f.as_poly(), g.as_poly()

    assert F.primitive() == (2, G)
    assert primitive(f) == (2, g)
    assert primitive(f, x) == (2, g)
    assert primitive(F) == (2, G)
    assert primitive(f, polys=True) == (2, G)
    assert primitive(F, polys=False) == (2, g)

    pytest.raises(ComputationFailedError, lambda: primitive(4))

    f = (2*x).as_poly(modulus=3)
    g = (2.0*x).as_poly(domain=RR)

    assert f.primitive() == (1, f)
    assert g.primitive() == (1.0, g)

    assert primitive(-3*x/4 + y + Rational(11, 8)) == \
        (Rational(-1, 8), 6*x - 8*y - 11)

    assert primitive(3*x + 2) == (1, 3*x + 2)
    assert primitive(4*x + 2) == (2, 2*x + 1)
    assert primitive(2*x**2 + 6*x + 12) == (2, x**2 + 3*x + 6)
    assert primitive(x**2 + 3*x + 6) == (1, x**2 + 3*x + 6)


def test_compose():
    f = x**12 + 20*x**10 + 150*x**8 + 500*x**6 + 625*x**4 - 2*x**3 - 10*x + 9
    g = x**4 - 2*x + 9
    h = x**3 + 5*x

    F, G, H = map(Poly, (f, g, h))

    assert G.compose(H) == F
    assert compose(g, h) == f
    assert compose(g, h, x) == f
    assert compose(G, H) == F
    assert compose(g, h, polys=True) == F
    assert compose(G, H, polys=False) == f

    assert F.decompose() == [G, H]
    assert decompose(f) == [g, h]
    assert decompose(f, x) == [g, h]
    assert decompose(F) == [G, H]
    assert decompose(f, polys=True) == [G, H]
    assert decompose(F, polys=False) == [g, h]

    pytest.raises(ComputationFailedError, lambda: compose(4, 2))
    pytest.raises(ComputationFailedError, lambda: decompose(4))

    assert compose(x**2 - y**2, x - y, x, y) == x**2 - 2*x*y
    assert compose(x**2 - y**2, x - y, y, x) == -y**2 + 2*x*y


def test_shift():
    assert (x**2 - 2*x + 1).as_poly().shift(2) == (x**2 + 2*x + 1).as_poly()


def test_sqf_norm():
    assert sqf_norm(x**2 - 2, extension=sqrt(3)) == \
        (1, x**2 - 2*sqrt(3)*x + 1, x**4 - 10*x**2 + 1)
    assert sqf_norm(x**2 - 3, extension=sqrt(2)) == \
        (1, x**2 - 2*sqrt(2)*x - 1, x**4 - 10*x**2 + 1)
    assert sqf_norm(x**2 - 3, extension=sqrt(2), polys=True) == \
        (1, (x**2 - 2*sqrt(2)*x - 1).as_poly(extension=sqrt(2)),
         (x**4 - 10*x**2 + 1).as_poly())

    pytest.raises(ComputationFailedError, lambda: sqf_norm([1, 2]))

    assert (x**2 - 2).as_poly(extension=sqrt(3)).sqf_norm() == \
        (1, (x**2 - 2*sqrt(3)*x + 1).as_poly(x, extension=sqrt(3)),
            (x**4 - 10*x**2 + 1).as_poly(x, domain=QQ))

    assert (x**2 - 3).as_poly(extension=sqrt(2)).sqf_norm() == \
        (1, (x**2 - 2*sqrt(2)*x - 1).as_poly(x, extension=sqrt(2)),
            (x**4 - 10*x**2 + 1).as_poly(x, domain=QQ))


def test_sqf():
    f = x**5 - x**3 - x**2 + 1
    g = x**3 + 2*x**2 + 2*x + 1
    h = x - 1

    p = x**4 + x**3 - x - 1

    F, G, H, P = map(Poly, (f, g, h, p))

    assert F.sqf_part() == P
    assert sqf_part(f) == p
    assert sqf_part(f, x) == p
    assert sqf_part(F) == P
    assert sqf_part(f, polys=True) == P
    assert sqf_part(F, polys=False) == p

    assert F.sqf_list() == (1, [(G, 1), (H, 2)])
    assert sqf_list(f) == (1, [(g, 1), (h, 2)])
    assert sqf_list(f, x) == (1, [(g, 1), (h, 2)])
    assert sqf_list(F) == (1, [(G, 1), (H, 2)])
    assert sqf_list(f, polys=True) == (1, [(G, 1), (H, 2)])
    assert sqf_list(F, polys=False) == (1, [(g, 1), (h, 2)])

    pytest.raises(PolynomialError, lambda: sqf_list([1, 2]))
    pytest.raises(ComputationFailedError, lambda: sqf_part(4))

    assert sqf(1) == 1
    assert sqf_list(1) == (1, [])

    assert sqf((2*x**2 + 2)**7) == 128*(x**2 + 1)**7

    assert sqf(f) == g*h**2
    assert sqf(f, x) == g*h**2

    d = x**2 + y**2

    assert sqf(f/d) == (g*h**2)/d
    assert sqf(f/d, x) == (g*h**2)/d

    assert sqf(x - 1) == x - 1
    assert sqf(-x - 1) == -x - 1

    assert sqf(x - 1) == x - 1
    assert sqf(6*x - 10) == Mul(2, 3*x - 5, evaluate=False)

    assert sqf((6*x - 10)/(3*x - 6)) == Rational(2, 3)*((3*x - 5)/(x - 2))
    assert sqf((x**2 - 2*x + 1).as_poly()) == (x - 1)**2

    f = 3 + x - x*(1 + x) + x**2

    assert sqf(f) == 3

    f = (x**2 + 2*x + 1)**20000000000

    assert sqf(f) == (x + 1)**40000000000
    assert sqf_list(f) == (1, [(x + 1, 40000000000)])

    pytest.raises(PolynomialError, lambda: sqf_list(x/(x**2 - 1), frac=False))
    assert sqf_list(x/(x**2 - 1), frac=True) == (1, [(x, 1)], [(x**2 - 1, 1)])


def test_factor():
    f = x**5 - x**3 - x**2 + 1

    u = x + 1
    v = x - 1
    w = x**2 + x + 1

    F, U, V, W = map(Poly, (f, u, v, w))

    assert F.factor_list() == (1, [(U, 1), (V, 2), (W, 1)])
    assert factor_list(f) == (1, [(u, 1), (v, 2), (w, 1)])
    assert factor_list(f, x) == (1, [(u, 1), (v, 2), (w, 1)])
    assert factor_list(F) == (1, [(U, 1), (V, 2), (W, 1)])
    assert factor_list(f, polys=True) == (1, [(U, 1), (V, 2), (W, 1)])
    assert factor_list(F, polys=False) == (1, [(u, 1), (v, 2), (w, 1)])

    assert factor_list(1) == (1, [])
    assert factor_list(6) == (6, [])
    assert factor_list(sqrt(3), x) == (1, [(3, Rational(1, 2))])
    assert factor_list((-1)**x, x) == (1, [(-1, x)])
    assert factor_list((2*x)**y, x) == (1, [(2, y), (x, y)])
    assert factor_list(sqrt(x*y), x) == (1, [(x*y, Rational(1, 2))])

    # issue sympy/sympy#11198
    assert factor_list(sqrt(2)*x) == (1, [(2, Rational(1, 2)), (x, 1)])
    assert factor_list(sqrt(2)*sin(x),
                       sin(x)) == (1, [(2, Rational(1, 2)), (sin(x), 1)])

    assert factor(6) == 6
    assert factor(6).is_Integer

    assert factor_list(3*x) == (3, [(x, 1)])
    assert factor_list(3*x**2) == (3, [(x, 2)])

    assert factor(3*x) == 3*x
    assert factor(3*x**2) == 3*x**2

    assert factor((2*x**2 + 2)**7) == 128*(x**2 + 1)**7

    assert factor(f) == u*v**2*w
    assert factor(f, x) == u*v**2*w

    g, p, q, r = x**2 - y**2, x - y, x + y, x**2 + 1

    assert factor(f/g) == (u*v**2*w)/(p*q)
    assert factor(f/g, x) == (u*v**2*w)/(p*q)

    p = Symbol('p', positive=True)
    i = Symbol('i', integer=True)
    r = Symbol('r', extended_real=True)

    assert factor(sqrt(x*y)).is_Pow is True

    assert factor(sqrt(3*x**2 - 3)) == sqrt(3)*sqrt((x - 1)*(x + 1))
    assert factor(sqrt(3*x**2 + 3)) == sqrt(3)*sqrt(x**2 + 1)

    assert factor((y*x**2 - y)**i) == y**i*(x - 1)**i*(x + 1)**i
    assert factor((y*x**2 + y)**i) == y**i*(x**2 + 1)**i

    assert factor((y*x**2 - y)**t) == (y*(x - 1)*(x + 1))**t
    assert factor((y*x**2 + y)**t) == (y*(x**2 + 1))**t

    f = sqrt(expand((r**2 + 1)*(p + 1)*(p - 1)*(p - 2)**3))
    g = sqrt((p - 2)**3*(p - 1))*sqrt(p + 1)*sqrt(r**2 + 1)

    assert factor(f) == g
    assert factor(g) == g

    g = (x - 1)**5*(r**2 + 1)
    f = sqrt(expand(g))

    assert factor(f) == sqrt(g)

    f = (sin(1)*x + 1).as_poly(x, domain=EX)

    assert f.factor_list() == (1, [(f, 1)])

    f = x**4 + 1

    assert factor(f) == f
    assert factor(f, extension=I) == (x**2 - I)*(x**2 + I)
    assert factor(f, gaussian=True) == (x**2 - I)*(x**2 + I)
    assert factor(f, extension=sqrt(2)) == (x**2 + sqrt(2)*x +
                                            1)*(x**2 - sqrt(2)*x + 1)

    f = x**2 + 2*sqrt(2)*x + 2

    assert factor(f, extension=sqrt(2)) == (x + sqrt(2))**2
    assert factor(f, extension=True) == (x + sqrt(2))**2
    assert factor(f**3, extension=sqrt(2)) == (x + sqrt(2))**6
    assert factor(f**3, extension=True) == (x + sqrt(2))**6

    # issue sympy/sympy#24346
    f = (x**2 - 2)/(x + sqrt(2))
    assert factor(f, extension=sqrt(2)) == x - sqrt(2)
    assert factor(f, extension=True) == x - sqrt(2)

    assert factor(x**2 - 2*y**2,
                  extension=sqrt(2)) == (x + sqrt(2)*y)*(x - sqrt(2)*y)
    assert factor(2*x**2 - 4*y**2,
                  extension=sqrt(2)) == 2*((x + sqrt(2)*y)*(x - sqrt(2)*y))

    assert factor(x - 1) == x - 1
    assert factor(-x - 1) == -x - 1

    assert factor(x - 1) == x - 1

    assert factor(6*x - 10) == Mul(2, 3*x - 5, evaluate=False)

    assert factor(x**11 + x + 1,
                  modulus=65537) == (x**2 + x + 1)*(x**9 + 65536*x**8 +
                                                    x**6 + 65536*x**5 +
                                                    x**3 + 65536*x**2 + 1)

    assert (factor(x**3 + 3*x + 2, modulus=4) ==
            factor((x**3 + 3*x + 2).as_poly(modulus=4)) ==
            (x + 1)*(x**2 + x + 2))
    assert (factor_list((x**3 + 3*x + 2).as_poly(modulus=4)) ==
            (1, [((x + 1).as_poly(modulus=4), 1),
                 ((x**2 + x + 2).as_poly(modulus=4), 1)]))

    f = x/pi + x*sin(x)/pi
    g = y/(pi**2 + 2*pi + 1) + y*sin(x)/(pi**2 + 2*pi + 1)

    assert factor(f) == x*(sin(x) + 1)/pi
    assert factor(g) == y*(sin(x) + 1)/(pi + 1)**2

    assert factor(Eq(x**2 + 2*x + 1, x**3 + 1)) == Eq((x + 1)**2,
                                                      (x + 1)*(x**2 - x + 1))

    f = (x**2 - 1)/(x**2 + 4*x + 4)

    assert factor(f) == (x + 1)*(x - 1)/(x + 2)**2
    assert factor(f, x) == (x + 1)*(x - 1)/(x + 2)**2

    f = 3 + x - x*(1 + x) + x**2

    assert factor(f) == 3
    assert factor(f, x) == 3

    assert factor(1/(x**2 + 2*x + 1/x) - 1) == -((1 - x + 2*x**2 +
                                                  x**3)/(1 + 2*x**2 + x**3))

    assert factor(f, expand=False) == f
    pytest.raises(PolynomialError, lambda: factor(f, x, expand=False))

    pytest.raises(FlagError, lambda: factor(x**2 - 1, polys=True))

    assert (factor([x, Eq(x**2 - y**2,
                          Tuple(x**2 - z**2, 1/x + 1/y))]) ==
            [x, Eq((x - y)*(x + y), Tuple((x - z)*(x + z), (x + y)/x/y))])

    assert not isinstance((x**3 + x + 1).as_poly().factor_list()[1][0][0],
                          PurePoly) is True
    assert isinstance(PurePoly(x**3 + x + 1).factor_list()[1][0][0],
                      PurePoly) is True

    assert factor(sqrt(-x)) == sqrt(-x)

    # issue sympy/sympy#5917
    e = (-2*x*(-x + 1)*(x - 1)*(-x*(-x + 1)*(x - 1) - x*(x - 1)**2) *
         (x**2*(x - 1) - x*(x - 1) - x) -
         (-2*x**2*(x - 1)**2 - x*(-x + 1)*(-x*(-x + 1) + x*(x - 1))) *
         (x**2*(x - 1)**4 - x*(-x*(-x + 1)*(x - 1) - x*(x - 1)**2)))
    assert factor(e) == 0

    # deep option
    assert factor(sin(x**2 + x) + x, deep=True) == sin(x*(x + 1)) + x

    assert factor(sqrt(x**2)) == sqrt(x**2)

    # issue sympy/sympy#7902
    assert (2*Sum(3*x, (x, 1, 9))).factor() == 6*Sum(x, (x, 1, 9))
    assert (2*Sum(x**2, (x, 1, 9))).factor() == 2*Sum(x**2, (x, 1, 9))

    # issue sympy/sympy#24346
    i = Integral(2*x, x, y)
    assert factor(i) == 2*Integral(x, x)*Integral(1, y)
    assert factor(1 + i) == 1 + i
    i2 = Integral(2*x*y, x, y)
    assert factor(i2) == 2*Integral(x, x)*Integral(y, y)

    A, B = symbols('A B', commutative=False)
    f = (x - A)*(y - B)
    assert factor(f.expand()) == f

    assert factor(Sum(4*x, (x, 1, y))) == 4*Sum(x, (x, 1, y))

    # issue sympy/sympy#13149
    assert (factor(expand((0.5*x + 1)*(0.5*y + 1))) ==
            Mul(4.0, 0.25*x + 0.5, 0.25*y + 0.5))
    assert factor(expand((0.5*x + 1)**2)) == 4.0*(0.25*x + 0.5)**2

    assert factor(x**4/2 + 5*x**3/12 - x**2/3) == x**2*(2*x - 1)*(3*x + 4)/12

    assert factor(x**6 - 4*x**4 + 4*x**3 -
                  x**2) == x**2*(x - 1)**2*(x**2 + 2*x - 1)

    # issue sympy/sympy#9607
    assert factor(1e-20*x - 7.292115e-5) == 1e-20*x - 7.292115e-5


def test_factor_large():
    f = (x**2 + 4*x + 4)**10000000*(x**2 + 1)*(x**2 + 2*x + 1)**1234567
    g = ((x**2 + 2*x + 1)**3000*y**2 + (x**2 + 2*x + 1)**3000*2*y +
         (x**2 + 2*x + 1)**3000)

    assert factor(f) == (x + 2)**20000000*(x**2 + 1)*(x + 1)**2469134
    assert factor(g) == (x + 1)**6000*(y + 1)**2

    assert factor_list(f) == (1, [(x**2 + 1, 1), (x + 1, 2469134),
                                  (x + 2, 20000000)])
    assert factor_list(g) == (1, [(y + 1, 2), (x + 1, 6000)])

    f = (x**2 - y**2)**200000*(x**7 + 1)
    g = (x**2 + y**2)**200000*(x**7 + 1)

    assert factor(f) == \
        (x + 1)*(x - y)**200000*(x + y)**200000*(x**6 - x**5 +
                                                 x**4 - x**3 + x**2 - x + 1)
    assert factor(g, gaussian=True) == \
        (x + 1)*(x - I*y)**200000*(x + I*y)**200000*(x**6 - x**5 +
                                                     x**4 - x**3 + x**2 - x + 1)

    assert factor_list(f) == \
        (1, [(x + 1, 1), (x - y, 200000), (x + y, 200000), (x**6 -
                                                            x**5 + x**4 - x**3 + x**2 - x + 1, 1)])
    assert factor_list(g, gaussian=True) == \
        (1, [(x + 1, 1), (x - I*y, 200000), (x + I*y, 200000), (
            x**6 - x**5 + x**4 - x**3 + x**2 - x + 1, 1)])


def test_factor_noeval():
    assert factor(6*x - 10) == Mul(2, 3*x - 5, evaluate=False)
    assert factor((6*x - 10)/(3*x - 6)) == Rational(2, 3)*((3*x - 5)/(x - 2))


def test_count_roots():
    assert count_roots(x**2 - 2) == 2

    assert count_roots(x**2 - 2, inf=-oo) == 2
    assert count_roots(x**2 - 2, sup=+oo) == 2
    assert count_roots(x**2 - 2, inf=-oo, sup=+oo) == 2

    assert count_roots(x**2 - 2, inf=-2) == 2
    assert count_roots(x**2 - 2, inf=-1) == 1

    assert count_roots(x**2 - 2, sup=1) == 1
    assert count_roots(x**2 - 2, sup=2) == 2

    assert count_roots(x**2 - 2, inf=-1, sup=1) == 0
    assert count_roots(x**2 - 2, inf=-2, sup=2) == 2

    assert count_roots(x**2 - 2, inf=-1, sup=1) == 0
    assert count_roots(x**2 - 2, inf=-2, sup=2) == 2

    assert count_roots(x**2 + 2) == 0
    assert count_roots(x**2 + 2, inf=-2*I) == 2
    assert count_roots(x**2 + 2, sup=+2*I) == 2
    assert count_roots(x**2 + 2, inf=-2*I, sup=+2*I) == 2

    assert count_roots(x**2 + 2, inf=0) == 0
    assert count_roots(x**2 + 2, sup=0) == 0

    assert count_roots(x**2 + 2, inf=-I) == 1
    assert count_roots(x**2 + 2, sup=+I) == 1

    assert count_roots(x**2 + 2, inf=+I/2, sup=+I) == 0
    assert count_roots(x**2 + 2, inf=-I, sup=-I/2) == 0

    assert count_roots(x**2 + 1, inf=-I, sup=1) == 1

    assert count_roots(x**4 - 4, inf=0, sup=1 + 3*I) == 1

    pytest.raises(PolynomialError, lambda: count_roots(1))


def test_sympyissue_12602():
    expr = 11355363812949950368319856364342755712460471081301053527133568171268803160551855579764764406412332136789657466300880824616465555590045220022768132246969281371700283178427904690172215428157788636395727*t**14/500000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 + 1182614101238502509994197939548011046110362591360244720959032955996959698293886871005468894084128139099293801809189908060595593758885614886473934547400040763077455747185622083724725964710198605960741*t**13/6250000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 + 12922822504467751142249299122933324184092020356108036157731097049497758652003692943810675925067800052327142015387959211427374009396705154181837176763552511140169473943304565171121276347837419884681487*t**12/1000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 + 2247204646780185719430022273864876084706708953097666720876560045907791931848809022047384483255204086759310635258105261945382035979316693256857004274231432751741774866992749719943120236265693542959121*t**11/10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 + 48361843519832766388699493325560944345496391872452676523436328806727211606456243561884964568166128629309073817110912281835520854136140406763166099011063597874739148993632932049857510934073377073756943*t**10/10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 + 80164052804759260531194459240730834126540153423610579212661871973340813173703351959915044156338949310408408075892534630817446213642840221172696829016781427774802300251456296296549939454943896121381103*t**9/1000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 + 11820298688854366697091577677719451524251323356734720588058286064886307168520540386661816641085576247246659191024148765432674755172844550760156289246049795015707597236348994207857339854655564489007679*t**8/10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 + 14908242237501135776503697900202844008307102016637953004464689811953233096672981556176039254271036473296837438230462049062156575391853236190707958824171141585643979204278060814079164562347710212783143*t**7/1000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 + 4004804534801340363315980630036975745885352263060800540982728634021250371208042415913703863593076428445232409745085269553713246947392078483528154594010983406996758112107177357774230732091992176355051*t**6/25000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 + 14405532019175131023474809758591882760843355517617477976264800133366074549575009123545658907344444270748700666056555232135755778765629022007752913521423634118196613981546114590366304386999628027814879*t**5/10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 + 10574051936471402714931747549059296527795248964344182731742838740900131997560377847217760142735497446739729085272580576056569967115897852649538335744262984346103108139783500273358254849041434565692381*t**4/1000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 + 7621024600438344015600198602526834713218087596509021419939455089811715884880919464619193508300267251654333701818067555976498078593210932512841643006106774171602557804624587725019604901643333157908251*t**3/125000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 + 810596792271328855063337004546095119237768466927680286629788036582515398849767422695474356109018097281604030779219064519000249092915587584571381512608327847454533913096943271752401133736953700148513*t**2/3125000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 + 7250440010579498280072969338541700246700434564503594127678121996192953098480655821398616200867454166215020508242889954318334847876038061179796990134727332078272146610064625333734530143125317393151907*t/10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 + 1
    assert count_roots(expr) == 0


def test_Poly_root():
    f = (2*x**3 - 7*x**2 + 4*x + 4).as_poly()

    assert f.root(0) == -Rational(1, 2)
    assert f.root(1) == 2
    assert f.root(2) == 2
    pytest.raises(IndexError, lambda: f.root(3))

    assert (x**5 + x + 1).as_poly().root(0) == RootOf(x**3 - x**2 + 1, 0)


def test_real_roots():
    assert real_roots(x) == [0]
    assert real_roots(x, multiple=False) == [(0, 1)]

    assert real_roots(x**3) == [0, 0, 0]
    assert real_roots(x**3, multiple=False) == [(0, 3)]

    assert real_roots(x*(x**3 + x + 3)) == [RootOf(x**3 + x + 3, 0), 0]
    assert real_roots(x*(x**3 + x + 3), multiple=False) == [(RootOf(
        x**3 + x + 3, 0), 1), (0, 1)]

    assert real_roots(
        x**3*(x**3 + x + 3)) == [RootOf(x**3 + x + 3, 0), 0, 0, 0]
    assert real_roots(x**3*(x**3 + x + 3), multiple=False) == [(RootOf(
        x**3 + x + 3, 0), 1), (0, 3)]

    f = 2*x**3 - 7*x**2 + 4*x + 4
    g = x**3 + x + 1

    assert f.as_poly().real_roots() == [-Rational(1, 2), 2, 2]
    assert g.as_poly().real_roots() == [RootOf(g, 0)]

    pytest.raises(PolynomialError, lambda: real_roots(1))


def test_all_roots():
    f = 2*x**3 - 7*x**2 + 4*x + 4
    g = x**3 + x + 1

    assert f.as_poly().all_roots() == [-Rational(1, 2), 2, 2]
    assert f.as_poly().all_roots(multiple=False) == [(-Rational(1, 2), 1), (2, 2)]
    assert g.as_poly().all_roots() == [RootOf(g, 0), RootOf(g, 1), RootOf(g, 2)]

    f = (x**7 - x).as_poly(modulus=7)

    # issue sympy/sympy#22673
    assert f.all_roots() == [RootOf(f, i, evaluate=False) for i in range(7)]


def test_nroots():
    assert not Integer(0).as_poly(x).nroots()
    assert not Integer(1).as_poly(x).nroots()

    assert (x**2 - 1).as_poly().nroots() == [-1.0, 1.0]
    assert (x**2 + 1).as_poly().nroots() == [-1.0*I, 1.0*I]

    roots = (x**2 - 1).as_poly().nroots()
    assert roots == [-1.0, 1.0]

    roots = (x**2 + 1).as_poly().nroots()
    assert roots == [-1.0*I, 1.0*I]

    roots = (x**2/3 - Rational(1, 3)).as_poly().nroots()
    assert roots == [-1.0, 1.0]

    roots = (x**2/3 + Rational(1, 3)).as_poly().nroots()
    assert roots == [-1.0*I, 1.0*I]

    assert (x**2 + 2*I).as_poly(x).nroots() == [-1.0 + 1.0*I, 1.0 - 1.0*I]
    assert (
        x**2 + 2*I).as_poly(x, extension=I).nroots() == [-1.0 + 1.0*I, 1.0 - 1.0*I]

    assert (0.2*x + 0.1).as_poly().nroots() == [-0.5]

    roots = nroots(x**5 + x + 1, n=5)
    eps = Float('1e-5')

    assert re(roots[0]).epsilon_eq(-0.75487, eps) is true
    assert im(roots[0]) == 0.0
    assert re(roots[1]) == -0.5
    assert im(roots[1]).epsilon_eq(-0.86602, eps) is true
    assert re(roots[2]) == -0.5
    assert im(roots[2]).epsilon_eq(+0.86602, eps) is true
    assert re(roots[3]).epsilon_eq(+0.87743, eps) is true
    assert im(roots[3]).epsilon_eq(-0.74486, eps) is true
    assert re(roots[4]).epsilon_eq(+0.87743, eps) is true
    assert im(roots[4]).epsilon_eq(+0.74486, eps) is true

    eps = Float('1e-6')

    assert re(roots[0]).epsilon_eq(-0.75487, eps) is false
    assert im(roots[0]) == 0.0
    assert re(roots[1]) == -0.5
    assert im(roots[1]).epsilon_eq(-0.86602, eps) is false
    assert re(roots[2]) == -0.5
    assert im(roots[2]).epsilon_eq(+0.86602, eps) is false
    assert re(roots[3]).epsilon_eq(+0.87743, eps) is false
    assert im(roots[3]).epsilon_eq(-0.74486, eps) is false
    assert re(roots[4]).epsilon_eq(+0.87743, eps) is false
    assert im(roots[4]).epsilon_eq(+0.74486, eps) is false

    pytest.raises(DomainError, lambda: (x + y).as_poly(x).nroots())
    pytest.raises(MultivariatePolynomialError, lambda: (x + y).as_poly().nroots())

    assert nroots(x**2 - 1) == [-1.0, 1.0]

    roots = nroots(x**2 - 1)
    assert roots == [-1.0, 1.0]

    assert nroots(x + I) == [-1.0*I]
    assert nroots(x + 2*I) == [-2.0*I]

    pytest.raises(PolynomialError, lambda: nroots(0))

    # issue sympy/sympy#8296
    f = (x**4 - 1).as_poly()
    assert f.nroots(2) == [w.evalf(2) for w in f.all_roots()]


def test_cancel():
    assert cancel(0) == 0
    assert cancel(7) == 7
    assert cancel(x) == x

    assert cancel(oo) == oo

    assert cancel((2, 3)) == (1, 2, 3)
    assert cancel((0, 1, 2)) == (0, 1, 2)

    assert cancel((1, 0), x) == (1, 1, 0)
    assert cancel((0, 1), x) == (1, 0, 1)

    f, g, p, q = 4*x**2 - 4, 2*x - 2, 2*x + 2, Integer(1)
    F, G, P, Q = (u.as_poly(x) for u in (f, g, p, q))

    assert F.cancel(G) == (1, P, Q)
    assert cancel((f, g)) == (1, p, q)
    assert cancel((f, g), x) == (1, p, q)
    assert cancel((F, G)) == (1, P, Q)
    assert cancel((f, g), polys=True) == (1, P, Q)
    assert cancel((F, G), polys=False) == (1, p, q)

    f = (x**2 - 2)/(x + sqrt(2))

    assert cancel(f) == f
    assert cancel(f, greedy=False) == x - sqrt(2)

    f = (x**2 - 2)/(x - sqrt(2))

    assert cancel(f) == f
    assert cancel(f, greedy=False) == x + sqrt(2)

    assert cancel((x**2/4 - 1, x/2 - 1)) == (Rational(1, 2), x + 2, 1)

    assert cancel((x**2 - y)/(x - y)) == 1/(x - y)*(x**2 - y)

    assert cancel((x**2 - y**2)/(x - y), x) == x + y
    assert cancel((x**2 - y**2)/(x - y), y) == x + y
    assert cancel((x**2 - y**2)/(x - y)) == x + y

    assert cancel((x**3 - 1)/(x**2 - 1)) == (x**2 + x + 1)/(x + 1)
    assert cancel((x**3/2 - Rational(1, 2))/(x**2 - 1)) == (x**2 + x + 1)/(2*x + 2)

    assert cancel((exp(2*x) + 2*exp(x) + 1)/(exp(x) + 1)) == exp(x) + 1

    f = (x**2 - a**2).as_poly(x)
    g = (x - a).as_poly(x)

    F = (x + a).as_poly(x)
    G = Integer(1).as_poly(x)

    assert cancel((f, g)) == (1, F, G)

    f = x**3 + (sqrt(2) - 2)*x**2 - (2*sqrt(2) + 3)*x - 3*sqrt(2)
    g = x**2 - 2

    assert cancel((f, g), extension=True) == (1, x**2 - 2*x - 3, x - sqrt(2))

    f = (-2*x + 3).as_poly()
    g = (-x**9 + x**8 + x**6 - x**5 + 2*x**2 - 3*x + 1).as_poly()

    assert cancel((f, g)) == (1, -f, -g)

    Zx = ZZ.inject(x)
    Zxf = Zx.field

    f = y.as_poly(y, domain=Zxf)
    g = Integer(1).as_poly(y, domain=Zx)

    assert f.cancel(
        g) == (1, y.as_poly(y, domain=Zxf), Integer(1).as_poly(y, domain=Zxf))
    assert f.cancel(g, include=True) == (
        y.as_poly(y, domain=Zxf), Integer(1).as_poly(y, domain=Zxf))

    f = (5*x*y + x).as_poly(y, domain=Zxf)
    g = (2*x**2*y).as_poly(y, domain=Zxf)

    assert f.cancel(g, include=True) == (
        (5*y + 1).as_poly(y, domain=Zxf), (2*x*y).as_poly(y, domain=Zxf))

    f = -(-2*x - 4*y + 0.005*(z - y)**2)/((z - y)*(-z + y + 2))
    assert cancel(f).is_Mul

    P = tanh(x - 3.0)
    Q = tanh(x + 3.0)
    f = ((-2*P**2 + 2)*(-P**2 + 1)*Q**2/2 + (-2*P**2 + 2)*(-2*Q**2 + 2)*P*Q - (-2*P**2 + 2)*P**2*Q**2 + (-2*Q**2 + 2)*(-Q**2 + 1)*P**2/2 - (-2*Q**2 + 2)*P**2*Q**2)/(2*sqrt(P**2*Q**2 + 0.0001)) \
        + (-(-2*P**2 + 2)*P*Q**2/2 - (-2*Q**2 + 2)*P**2*Q/2)*((-2*P**2 + 2)*P*Q**2/2 + (-2*Q**2 + 2)*P**2*Q/2)/(2*(P**2*Q**2 + 0.0001)**Rational(3, 2))
    assert cancel(f).is_Mul

    # issue sympy/sympy#7022
    A = Symbol('A', commutative=False)
    p1 = Piecewise((A*(x**2 - 1)/(x + 1), x > 1), ((x + 2)/(x**2 + 2*x), True))
    p2 = Piecewise((A*(x - 1), x > 1), (1/x, True))
    assert cancel(p1) == p2
    assert cancel(2*p1) == 2*p2
    assert cancel(1 + p1) == 1 + p2
    assert cancel((x**2 - 1)/(x + 1)*p1) == (x - 1)*p2
    assert cancel((x**2 - 1)/(x + 1) + p1) == (x - 1) + p2
    p3 = Piecewise(((x**2 - 1)/(x + 1), x > 1), ((x + 2)/(x**2 + 2*x), True))
    p4 = Piecewise(((x - 1), x > 1), (1/x, True))
    assert cancel(p3) == p4
    assert cancel(2*p3) == 2*p4
    assert cancel(1 + p3) == 1 + p4
    assert cancel((x**2 - 1)/(x + 1)*p3) == (x - 1)*p4
    assert cancel((x**2 - 1)/(x + 1) + p3) == (x - 1) + p4

    # issue sympy/sympy#9363
    M = MatrixSymbol('M', 5, 5)
    assert cancel(M[0, 0] + 7) == M[0, 0] + 7
    expr = (z*sin(M[1, 4] + M[2, 1] * 5 * M[4, 0]) - 5 * M[1, 2])/z
    assert cancel(expr) == expr

    assert cancel(((x - 1)**2/(x - 1), (x + 2*x**2)/x,
                   (x - x**3)/x)) == (x - 1, 2*x + 1, -x**2 + 1)

    # issue sympy/sympy#12531
    e = (x**4/24 - x*(x**3/24 + Rational(7, 8)) +
         13*x/12)/((x**3/24 + Rational(7, 8))*(-x**4/6 - x/3) +
                   (x**3/6 - Rational(1, 2))*(x**4/24 + 13*x/12))
    assert cancel(e) == Rational(-1, 4)


def test_reduced():
    f = 2*x**4 + y**2 - x**2 + y**3
    G = [x**3 - x, y**3 - y]

    Q = [2*x, 1]
    r = x**2 + y**2 + y

    assert reduced(f, G) == (Q, r)
    assert reduced(f, G, x, y) == (Q, r)

    H = groebner(G)

    assert H.reduce(f) == (Q, r)

    Q = [(2*x).as_poly(x, y), Integer(1).as_poly(x, y)]
    r = (x**2 + y**2 + y).as_poly()

    assert reduced(f, G, polys=True) == (Q, r)
    assert reduced(f, G, x, y, polys=True) == (Q, r)

    H = groebner(G, polys=True)

    assert H.reduce(f) == (Q, r)

    f = 2*x**3 + y**3 + 3*y
    G = groebner([x**2 + y**2 - 1, x*y - 2])

    Q = [x**2 - x*y**3/2 + x*y/2 + y**6/4 - y**4/2 + y**2/4, -y**5/4 + y**3/2 + 3*y/4]
    r = 0

    assert reduced(f, G) == (Q, r)
    assert G.reduce(f) == (Q, r)

    assert reduced(f, G, auto=False)[1] != 0
    assert G.reduce(f, auto=False)[1] != 0

    assert G.contains(f) is True
    assert G.contains(f + 1) is False

    assert reduced(1, [1], x) == ([1], 0)
    pytest.raises(ComputationFailedError, lambda: reduced(1, [1]))


def test_groebner():
    assert not groebner([], x, y, z)

    assert groebner([x**2 + 1, y**4*x + x**3], x, y, order='lex') == [1 + x**2, -1 + y**4]
    assert groebner([x**2 + 1, y**4*x + x**3, x*y*z**3], x, y, z, order='grevlex') == [-1 + y**4, z**3, 1 + x**2]

    assert groebner([x**2 + 1, y**4*x + x**3], x, y, order='lex', polys=True) == \
        [(1 + x**2).as_poly(x, y), (-1 + y**4).as_poly(x, y)]
    assert groebner([x**2 + 1, y**4*x + x**3, x*y*z**3], x, y, z, order='grevlex', polys=True) == \
        [(-1 + y**4).as_poly(x, y, z, order='grevlex'), (z**3).as_poly(x, y, z, order='grevlex'),
         (1 + x**2).as_poly(x, y, z, order='grevlex')]

    assert groebner([x**3 - 1, x**2 - 1]) == [x - 1]

    F = [3*x**2 + y*z - 5*x - 1, 2*x + 3*x*y + y**2, x - 3*y + x*z - 2*z**2]
    f = z**9 - x**2*y**3 - 3*x*y**2*z + 11*y*z**2 + x**2*z**2 - 5

    G = groebner(F, x, y, z, modulus=7)

    assert G == [1 + x + y + 3*z + 2*z**2 + 2*z**3 + 6*z**4 + z**5,
                 1 + 3*y + y**2 + 6*z**2 + 3*z**3 + 3*z**4 + 3*z**5 + 4*z**6,
                 1 + 4*y + 4*z + y*z + 4*z**3 + z**4 + z**6,
                 6 + 6*z + z**2 + 4*z**3 + 3*z**4 + 6*z**5 + 3*z**6 + z**7]

    Q, r = reduced(f, G, x, y, z, modulus=7, polys=True)

    assert sum((q*g for q, g in zip(Q, G.polys)), r) == f.as_poly(modulus=7)

    F = [x*y - 2*y, 2*y**2 - x**2]

    assert groebner(F, x, y) == \
        [x**2 - 2*y**2, x*y - 2*y, y**3 - 2*y]
    assert groebner(F, x, y, order='grlex') == \
        [y**3 - 2*y, x**2 - 2*y**2, x*y - 2*y]
    assert groebner(F, y, x, order='grevlex') == \
        [x**3 - 2*x**2, -x**2 + 2*y**2, x*y - 2*y]
    assert groebner(F, order='grevlex', field=True) == \
        [y**3 - 2*y, x**2 - 2*y**2, x*y - 2*y]

    assert groebner([1], x) == [1]

    pytest.raises(ComputationFailedError, lambda: groebner([1]))

    assert groebner([x**2 - 1, x**3 + 1], method='buchberger') == [x + 1]
    assert groebner([x**2 - 1, x**3 + 1], method='f5b') == [x + 1]

    pytest.raises(ValueError, lambda: groebner([x, y], method='unknown'))

    F = [x**2 - x - 1, (2*x - 1) * y - (x**10 - (1 - x)**10)]
    assert groebner(F, x, y, method='buchberger') == [x**2 - x - 1, y - 55]
    assert groebner(F, x, y, method='f5b') == [x**2 - x - 1, y - 55]

    # issue sympy/sympy#11623
    pytest.raises(ValueError,
                  lambda: groebner([0.144*x*y + 0.018*x**2 + 0.05*x - 1.577,
                                    0.072*y**2 + 0.036*x*y + 0.05*y - 1.423],
                                   x, y))


def test_set_order():
    F = [a + b + c + d, a*b + a*d + b*c + b*d, a*b*c + a*b*d + a*c*d + b*c*d, a*b*c*d - 1]
    G = groebner(F, a, b, c, d, order=grlex)

    B = [4*a + 3*d**9 - 4*d**5 - 3*d,
         4*b + 4*c - 3*d**9 + 4*d**5 + 7*d,
         4*c**2 + 3*d**10 - 4*d**6 - 3*d**2,
         4*c*d**4 + 4*c - d**9 + 4*d**5 + 5*d,
         d**12 - d**8 - d**4 + 1]

    assert groebner(F, a, b, c, d, order=lex) == B
    assert G.set_order(lex) == B
    assert G.set_order(grlex) == G

    F = [9*x**8 + 36*x**7 - 32*x**6 - 252*x**5 - 78*x**4 + 468*x**3 + 288*x**2 - 108*x + 9,
         -72*t*x**7 - 252*t*x**6 + 192*t*x**5 + 1260*t*x**4 + 312*t*x**3 - 404*t*x**2 - 576*t*x +
         108*t - 72*x**7 - 256*x**6 + 192*x**5 + 1280*x**4 + 312*x**3 - 576*x + 96]
    G = groebner(F, t, x, order=grlex)

    B = [203577793572507451707*t + 627982239411707112*x**7 - 666924143779443762*x**6 -
         10874593056632447619*x**5 + 5119998792707079562*x**4 + 72917161949456066376*x**3 +
         20362663855832380362*x**2 - 142079311455258371571*x + 183756699868981873194,
         9*x**8 + 36*x**7 - 32*x**6 - 252*x**5 - 78*x**4 + 468*x**3 + 288*x**2 - 108*x + 9]

    assert groebner(F, t, x, order=lex) == B
    assert G.set_order(lex) == B

    F = [x**2 - x - 3*y + 1, -2*x + y**2 + y - 1]
    G = groebner(F, x, y, order=lex)
    B = [x**2 - x - 3*y + 1, y**2 - 2*x + y - 1]

    assert groebner(F, x, y, order=grlex) == B
    assert G.set_order(grlex) == B
    assert G.set_order(grlex).set_order(lex) == G
    assert G == [2*x - y**2 - y + 1, y**4 + 2*y**3 - 3*y**2 - 16*y + 7]

    F = [x**2 - 2*y + 1, x + y/2]
    G = groebner(F, x, y, order=grlex)
    B = [y**2 - 8*y + 4, x + y/2]

    assert G == B
    assert G.set_order(lex) == reversed(B)

    G = groebner([x**3 - y**3], x, y, order='grlex')
    pytest.raises(NotImplementedError, lambda: G.set_order('lex'))


def test_dimension_and_independent_sets():
    assert groebner((x, y)).dimension == 0
    assert groebner((x**3 + y**2,)).dimension == 1
    assert groebner((x, y, z)).dimension == 0
    assert groebner((x, y, z), x, y, z, t).dimension == 1
    assert groebner((x*y - z, y*z - x, x*y - y)).dimension == 0
    assert groebner((x**2 - 2*x*z + 5, x*y**2 + y*z**3, 3*y**2 - 8*z**2)).dimension == 0

    assert groebner((x + y, x - y)).independent_sets == [[]]
    assert groebner((x + y, 2*x + 2*y)).independent_sets == [[y]]
    assert groebner((x**2 + y**2,)).independent_sets == [[y]]
    assert groebner((x**3*y**2 - 1,)).independent_sets == [[y], [x]]
    assert groebner((x**3 - y**3,)).independent_sets == [[y]]
    assert groebner((y - x, y - x - 1)).independent_sets is None
    assert groebner((x*y - z**2 - z, x**2 + x - y*z, x*z - y**2 - y)).independent_sets == [[z]]
    assert groebner((x*y*z,)).independent_sets == [[y, z], [x, z], [x, y]]
    assert groebner((x**2 - 1, (x - 1)*y, (x + 1)*z)).independent_sets == [[z], [y]]
    assert groebner((x**2 + y**2 + z**2, x + y - z, y + z**2)).independent_sets == [[]]
    assert groebner((x*z - 2*y + 1, y*z - 1 + z, y*z + x*y*z + z)).independent_sets == [[]]
    assert groebner((x**3*y*z - x*z**2, x*y**2*z - x*y*z, x**2*y**2 - z)).independent_sets == [[z], [y], [x]]
    assert groebner((x*y**2 - z - z**2, x**2*y - y, y**2 - z**2)).independent_sets == [[x]]
    assert groebner((x*y + z - 1, x - y - z**2, x**2 - 2*y**2 + 1)).independent_sets == [[]]
    assert groebner((z*x - y - x + x*y, y*z - z + x**2 + y*x**2, x - x**2 + y, z)).independent_sets == [[]]
    assert groebner((x*y - x*z + y**2, y*z - x**2 + x**2*y, x - x*y + y)).independent_sets == [[z]]
    assert groebner((y*z + x**2 + z, x*y*z + x*z - y**3, x*z + y**2)).independent_sets == [[]]
    assert groebner((x**2 + z**2*y + y*z, y**2 - z*x + x, x*y + z**2 - 1)).independent_sets == [[]]
    assert groebner((x + y**2*z - 2*y**2 + 4*y - 2*z - 1, -x + y**2*z - 1)).independent_sets == [[z], [y]]
    assert groebner((x, y - 1, z)).independent_sets == [[]]

    # H. Kredel and V. Weispfennig. Computing dimension and independent sets for
    # polynomial ideals. J. Symbolic Computation, 6(1):231–247, November 1988.

    # Ex. 4.1.
    V = A31, A32, A21, B1, B2, B3, C3, C2 = symbols('A31 A32 A21 B1 B2 B3 C3 C2')
    S = (C2 - A21, C3 - A31 - A32, B1 + B2 + B3 - 1,
         B2*C2 + B3*C3 - QQ(1, 2), B2*C2**2 + B3*C3**2 - QQ(1, 3),
         B3*A32*C2 - QQ(1, 6))
    G = groebner(S, *V, domain=QQ)
    assert G.independent_sets == [[C3, C2], [B3, C2], [B2, C3], [B2, B3], [A32, C3], [A32, B2]]
    assert G.dimension == 2

    # Ex. 4.3
    V = B1, A32, B2, B3, A, C3, C2, B = symbols('B1 A32 B2 B3 A C3 C2 B')
    S = (B1 + B2 + B3 - A - B,
         B2*C2 + B3*C3 - QQ(1, 2) - B/2 - B**2 + A*B,
         B2*C2**2 + B3*C3**2 - A/3 - A*B**2 + 4*B/3 + B**2 + B**3,
         B3*A32*C2 - A/6 - A*B/2 - A*B**2 + 2*B/3 + B**2 + B**3,
         B2*C2**3 + B3*C3**3 - QQ(1, 4) - B/4 - 5*B**2/2 - 3*B**3/2 - B**4 + A*B + A*B**3,
         B3*C3*A32*C2 - QQ(1, 8) - 3*B/8 - 7*B**2/4 - 3*B**3/2 - B**4 + A*B/2 + A*B**2/2 + A*B**3,
         B3*A32*C2**2 - QQ(1, 12) - B/12 - 7*B**2/6 - 3*B**3/2 - B**4 + 2*A*B/3 + A*B**2 + A*B**3,
         QQ(1, 24) + 7*B/24 + 13*B**2/12 + 3*B**3/2 + B**4 - A*B/3 - A*B**2 - A*B**3)
    G = groebner(S, *V, domain=QQ)
    assert G.independent_sets == [[B3, C2], [A32, C3, C2], [A32, B2, C3], [A32, B2, B3]]
    assert G.dimension == 3

    # Ex. 4.4
    V = L7, L6, L4, L1, L5, L3, L2 = symbols('L7 L6 L4 L1 L5 L3 L2')
    S = (L1*(L4 - L5/2 + L6),
         (2*L1**2/7 - L4)*(-10*L1 + 5*L2 - L3),
         (2*L1**2/7 - L4)*(3*L4 - L5 + L6),
         (-2*L1**2 + L1*L2 + 2*L1*L3 - L2**2 - 7*L5 + 21*L6)*(-3*L1 + 2*L2) + 21*(7*L7 - 2*L1*L4 + 3*L1**3/7),
         (-2*L1**2 + L1*L2 + 2*L1*L3 - L2**2 - 7*L5 + 21*L6)*(2*L4 - 2*L5) + (7*L7 - 2*L1*L4 + 3*L1**3/7)*(-45*L1 + 15*L2 - 3*L3),
         2*(-2*L1**2 + L1*L2 + 2*L1*L3 - L2**2 - 7*L5 + 21*L6)*L7 + (7*L7 - 2*L1*L4 + 3*L1**3/7)*(12*L4 - 3*L5 + 2*L6),
         (L1*(5*L1 - 3*L2 + L3))*(2*L2 - L1) + 7*(L1*(2*L6 - 4*L4)),
         (L1*(5*L1 - 3*L2 + L3))*L3+7*(L1*(2*L6 - 4*L4)),
         (L1*(5*L1 - 3*L2 + L3))*(-2*L4 - 2*L5) + (L1*(2*L6 - 4*L4))*(2*L2 - 8*L1) + 42*L1*L7,
         (L1*(5*L1 - 3*L2 + L3))*(8*L5/3 + 6*L6) + (L1*(2*L6 - 4*L4))*(11*L1 - 17*L2/3 + 5*L3/3) - 84*L1*L7,
         15*L7*(L1*(5*L1 - 3*L2 + L3)) + (L1*(2*L6 - 4*L4))*(5*L4 - 2*L5) + L1*L7*(-120*L1 + 30*L2 - 6*L3)/2,
         -3*(L1*(5*L1 - 3*L2 + L3))*L7 + (L1*(2*L6 - 4*L4))*(-L4/2 + L5/4 - L6/2) + L1*L7/2*(24*L1 - 6*L2),
         3*(L1*(2*L6 - 4*L4))*L7 + L1*L7*(40*L4 - 8*L5 + 4*L6)/2)
    G = groebner(S, *V, domain=QQ)
    assert G.independent_sets == [[L5, L3, L2], [L6, L3]]
    assert G.dimension == 3

    # Algebraic Solution of Nonlinear Equation Systems in REDUCE, p.7.
    V = ax, bx, cx, gx, jx, lx, mx, nx, q = symbols('ax bx cx gx jx lx mx nx q')
    S = (ax*q - lx*q - mx, ax - gx*q - lx, bx*q**2 + cx*q - jx*q - nx,
         q*(-ax*q + lx*q + mx), q*(-ax + gx*q + lx))
    G = groebner(S, *V, domain=QQ)
    assert G.independent_sets == [[cx, jx, lx, mx, nx, q], [cx, gx, jx, lx, mx, nx], [bx, cx, gx, jx, lx, nx]]
    assert G.dimension == 6


def test_GroebnerBasis():
    F = [x*y - 2*y, 2*y**2 - x**2]

    G = groebner(F, x, y, order='grevlex')
    assert groebner(F + [0], x, y, order='grevlex') == G

    assert G.args == ((y**3 - 2*y, x**2 - 2*y**2, x*y - 2*y), x, y)

    H = [y**3 - 2*y, x**2 - 2*y**2, x*y - 2*y]
    P = [h.as_poly(x, y, order='grevlex') for h in H]

    assert isinstance(G, GroebnerBasis) is True

    assert len(G) == 3

    assert G[0] == H[0]
    assert not G[0].is_Poly
    assert G[1] == H[1]
    assert not G[1].is_Poly
    assert G[2] == H[2]
    assert not G[2].is_Poly

    assert G[1:] == H[1:]
    assert not any(g.is_Poly for g in G[1:])
    assert G[:2] == H[:2]
    assert not any(g.is_Poly for g in G[1:])

    assert G.exprs == H
    assert G.polys == P
    assert G.gens == (x, y)
    assert G.domain == ZZ
    assert G.order == grevlex

    assert G == H
    assert G == tuple(H)
    assert G == P
    assert G == tuple(P)

    assert G

    G = groebner(F, x, y, order='grevlex', polys=True)

    assert G[0] == P[0]
    assert G[0].is_Poly
    assert G[1] == P[1]
    assert G[1].is_Poly
    assert G[2] == P[2]
    assert G[2].is_Poly

    assert G[1:] == P[1:]
    assert all(g.is_Poly for g in G[1:])
    assert G[:2] == P[:2]
    assert all(g.is_Poly for g in G[1:])

    assert tuple(G) == ((y**3 - 2*y).as_poly(x, y, order='grevlex'),
                        (x**2 - 2*y**2).as_poly(x, y, order='grevlex'),
                        (x*y - 2*y).as_poly(x, y, order='grevlex'))

    G = groebner(F, x, y, order='grevlex', polys=True)
    assert hash(G) == hash(groebner([_.as_poly() for _ in F], order='grevlex'))

    assert (G == 1) is False


def test_Poly_from_expr_recursive():
    assert (x*(x**2 + x - 1)**2).as_poly() == (x**5 + 2*x**4 - x**3 -
                                               2*x**2 + x).as_poly()

    assert (x + y).as_poly(wrt=y) == (x + y).as_poly(y, x)
    assert (x + sin(x)).as_poly(wrt=sin(x)) == (x + sin(x)).as_poly(sin(x), x)

    assert (2*(y + z)**2 - 1).as_poly() == (2*y**2 + 4*y*z +
                                            2*z**2 - 1).as_poly()
    assert (x*(y + z)**2 - 1).as_poly() == (x*y**2 + 2*x*y*z +
                                            x*z**2 - 1).as_poly()
    assert (2*x*(y + z)**2 - 1).as_poly() == (2*x*y**2 + 4*x*y*z +
                                              2*x*z**2 - 1).as_poly()

    assert (2*(y + z)**2 - x - 1).as_poly() == (2*y**2 + 4*y*z + 2*z**2 -
                                                x - 1).as_poly()
    assert (x*(y + z)**2 - x - 1).as_poly() == (x*y**2 + 2*x*y*z +
                                                x*z**2 - x - 1).as_poly()
    assert (2*x*(y + z)**2 - x - 1).as_poly() == (2*x*y**2 + 4*x*y*z + 2 *
                                                  x*z**2 - x - 1).as_poly()

    assert (x*y + (x + y)**2 + (x + z)**2).as_poly() == (2*x*z + 3*x*y + y**2 +
                                                         z**2 + 2*x**2).as_poly()
    assert (x*y*(x + y)*(x + z)**2).as_poly() == (x**3*y**2 + x*y**2*z**2 +
                                                  y*x**2*z**2 + 2*z*x**2*y**2 +
                                                  2*y*z*x**3 + y*x**4).as_poly()

    assert ((x + y)**2).as_poly(x) == (x**2 + 2*x*y + y**2).as_poly(x)
    assert ((x + y)**2).as_poly(x, expand=True) == (x**2 + 2*x*y +
                                                    y**2).as_poly(x)
    assert ((x + y)**2).as_poly(y) == (x**2 + 2*x*y + y**2).as_poly(y)
    assert ((x + y)**2 - y**2 - 2*x*y).as_poly() == (x**2).as_poly(x, y)

    e = x**2 + (1 + sqrt(2))*x + 1

    assert (e.as_poly(x, greedy=False) ==
            e.as_poly(x, domain=QQ.algebraic_field(sqrt(2))))

    # issue sympy/sympy#12400
    assert ((1/(1 + sqrt(2))).as_poly(x) ==
            (1/(1 + sqrt(2))).as_poly(x, domain=QQ.algebraic_field(1/(1 + sqrt(2)))))

    # issue sympy/sympy#19755
    assert ((x + (2*x + 3)**2/5 + Rational(6, 5)).as_poly() ==
            (4*x**2/5 + 17*x/5 + 3).as_poly(domain=QQ))
    assert (((x + 1)**2)/2).as_poly() == (x**2/2 + x +
                                          Rational(1, 2)).as_poly(domain=QQ)


def test_keep_coeff():
    u = Mul(2, x + 1, evaluate=False)
    assert _keep_coeff(Integer(1), x) == x
    assert _keep_coeff(Integer(-1), x) == -x
    assert _keep_coeff(Float(1.0), x) == 1.0*x
    assert _keep_coeff(Float(-1.0), x) == -1.0*x
    assert _keep_coeff(Integer(1), 2*x) == 2*x
    assert _keep_coeff(Integer(2), x/2) == x
    assert _keep_coeff(Integer(2), sin(x)) == 2*sin(x)
    assert _keep_coeff(Integer(2), x + 1) == u
    assert _keep_coeff(x, 1/x) == 1
    assert _keep_coeff(x + 1, Integer(2)) == u


def test_poly_matching_consistency():
    # Test for sympy/sympy#5514
    assert I * x.as_poly() == (I*x).as_poly(x)
    assert x.as_poly() * I == (I*x).as_poly(x)


def test_sympyissue_5786():
    f, g = z - I*t, x - I*y
    assert factor(expand(f*g), extension=[I]) == f*g
    assert factor(expand(f**2*g), extension=[I]) == f**2*g
    assert factor(expand(f*g**3), extension=[I]) == f*g**3

    # issue sympy/sympy#18895
    e = (x - 1)*(y - 1)
    assert factor(expand(e)) == e
    assert factor(expand(e), extension=[I]) == e


def test_noncommutative():
    class Foo(Expr):
        is_commutative = False
    e = x/(x + x*y)
    c = 1/(1 + y)
    fe, fc = map(Foo, [e, c])
    assert cancel(fe) == fc
    assert cancel(e + fe) == c + fc
    assert cancel(e*fc) == c*fc


def test_to_rational_coeffs():
    assert to_rational_coeffs((x**3 + y*x**2 + sqrt(y)).as_poly(x, domain=EX)) is None
    assert to_rational_coeffs((((x**2 - 1)*(x - 2)*y).subs({x: x*(1 + sqrt(2))})).as_poly(x, y, domain=EX)) is None
    assert to_rational_coeffs((x**5 + sqrt(2)*x**2 + 1).as_poly(x, domain=EX)) is None


def test_sympyissue_8754():
    z = 0.0001*(x*(x + (4.0*y))) + 0.0001*(y*(x + (4.0*y)))
    w = expand(z)
    v = factor(w)
    assert v == 10000.0*((0.0001*x + 0.0001*y)*(0.0001*x +
                                                0.00040000000000000002*y))
    assert expand(v) == w


def test_factor_terms():
    # issue sympy/sympy#7067
    assert factor_list(x*(x + y)) == (1, [(x, 1), (x + y, 1)])
    assert sqf_list(x*(x + y)) == (1, [(x**2 + x*y, 1)])


def test_sympyissue_8210():
    p = Integer(0).as_poly(x)
    p2 = p.copy()
    assert id(p) != id(p2)
    assert p == p2


def test_sympyissue_11775():
    e = y**4 + x*y**3 + y**2 + x*y
    assert factor_list(e, y) == (1, [(y, 1), (y**2 + 1, 1), (x + y, 1)])


def test_sympyissue_5602():
    (Integral(x, (x, 0, 1))*x + x**2).as_poly(x)


def test_sympyissue_15798():
    o1 = (x + y).as_poly(x, y, z)
    o2 = o1.copy()
    assert o1 == o2
    p = (x + sqrt(2)).as_poly(x)
    assert p == p.copy()


@pytest.mark.timeout(20)
def test_sympyissue_19670():
    (E**100000000).as_poly()


def test_sympyissue_8810():
    e = y**3 + y**2*sqrt(x) + y + x
    p = e.as_poly(y)
    c = e.as_poly(y, composite=True)

    assert c == e.as_poly(y, domain=ZZ.inject(x, sqrt(x)))
    assert p.as_poly(y, composite=True) == c


def test_sympyissue_8695():
    e = (x**2 + 1) * (x - 1)**2 * (x - 2)**3 * (x - 3)**3
    r = (1, [(x**2 + 1, 1), (x - 1, 2), (x**2 - 5*x + 6, 3)])
    assert sqf_list(e) == r
    assert e.as_poly().sqf_list() == r

    # regression test from the issue thread, not related to the issue
    e = (x + 2)**2 * (y + 4)**5
    assert sqf(e) == sqf(e.expand()) == e


def test_sympyissue_19070():
    e = (5*x).as_poly(modulus=19)
    r = e*2
    assert r == (10*x).as_poly(modulus=19)
    assert r.get_modulus() == 19


def test_sympyissue_19161():
    assert sympify('x**2').as_poly().simplify() == (x**2).as_poly()


def test_sympyissue_20484():
    assert (x*y*z).as_poly().eval(x, y*z) == (y**2*z**2).as_poly()


def test_sympyissue_20640():
    p = (x**2 + y).as_poly(field=True)
    p0 = y.as_poly(x, y, field=True)

    assert div(p, p0) == (Integer(1).as_poly(x, y, field=True),
                          (x**2).as_poly(x, y, field=True))
    assert div(p.as_expr(), p0.as_expr(), field=True) == (1, x**2)


def test_sympyissue_20973():
    e = exp(1 + O(x))
    assert cancel(e) == e


def test_sympyissue_20985():
    assert degree(1.0 + I*x/y, domain=CC.frac_field(y)) == 1


def test_sympyissue_21180():
    f = (x**4 + 6*x**3 + 4*x**2 - 30*x - 45).as_poly()
    assert factor(f) == (x + 3)**2*(x**2 - 5)


def test_sympyissue_20444():
    e = 33*log(x) + log(8) + 58

    assert LT(e) == 3*log(2)


def test_sympyissue_13029():
    assert sqf_part(a*(x - 1)**2*(y - 3)**3, x, y) == x*y - 3*x - y + 3


@pytest.mark.timeout(5)
def test_sympyissue_21760():
    _, r = (x**2016 - x**2015 + x**1008 + x**1003 +
            1).as_poly().div((x - 1).as_poly())
    assert r == Integer(3).as_poly(x)


def test_sympyissue_21761():
    t = tan(pi/7)
    assert factor(-exp(x)*t + 1, extension=True) == -((exp(x) - 5*t - t**5/7 +
                                                       3*t**3)*t)


def test_sympyissue_22093():
    expr = ((2*y**3*sin(x/y)**2 + x)**2*(y*(-6*y**2*sin(x/y)**2 +
                                            4*y*x*sin(x/y)*cos(x/y)) /
                                         (2*y**3*sin(x/y)**2 + x)**2 +
                                         1/(2*y**3*sin(x/y)**2 + x)) /
            (4*y*(2*y**2*(3*y*sin(x/y) - 2*x*cos(x/y))**2*sin(x/y)**2 /
             (2*y**3*sin(x/y)**2 + x) - 3*y*sin(x/y)**2 +
             4*x*sin(x/y)*cos(x/y) - (3*y*sin(x/y) - 2*x*cos(x/y))*sin(x/y) +
             x**2*sin(x/y)**2/y - x**2*cos(x/y)**2/y)))
    res = -(4*x**2*y**2*sin(x/y)*cos(x/y) + x**2 +
            8*x*y**5*sin(x/y)**3*cos(x/y) - 2*x*y**3*sin(x/y)**2 -
            8*y**6*sin(x/y)**4)/(-4*x**3*sin(x/y)**2 + 4*x**3*cos(x/y)**2 -
                                 8*x**2*y**3*sin(x/y)**4 -
                                 24*x**2*y**3*sin(x/y)**2*cos(x/y)**2 -
                                 24*x**2*y*sin(x/y)*cos(x/y) +
                                 48*x*y**4*sin(x/y)**3*cos(x/y) +
                                 24*x*y**2*sin(x/y)**2 - 24*y**5*sin(x/y)**4)
    assert cancel(expr).equals(res)


def test_sympyissue_22673():
    e = x**7 - x
    p = e.as_poly(modulus=7)
    f = x*(x + 1)*(x + 2)*(x + 3)*(x + 4)*(x + 5)*(x + 6)
    assert factor(e, modulus=7) == factor(p) == f
    assert factor_list(e, modulus=7) == (1, [(x + i, 1) for i in range(7)])
    assert factor_list(p) == (1, [((x + i).as_poly(modulus=7), 1)
                                  for i in range(7)])


@pytest.mark.timeout(10)
def test_sympyissue_23766():
    assert factor(exp((x + 1)**25 + 1),
                  deep=True) == exp((x + 2)*(x**4 + 3*x**3 + 4*x**2 + 2*x +
                                             1)*(x**20 + 20*x**19 + 190*x**18 +
                                                 1140*x**17 + 4845*x**16 +
                                                 15503*x**15 + 38745*x**14 +
                                                 77415*x**13 + 125515*x**12 +
                                                 166595*x**11 + 181754*x**10 +
                                                 162965*x**9 + 119580*x**8 +
                                                 71205*x**7 + 33965*x**6 +
                                                 12752*x**5 + 3685*x**4 +
                                                 795*x**3 + 120*x**2 +
                                                 10*x + 1))


def test_sympyissue_24461():
    a = (x**4/(z**2 + 3*z + 2)*y).as_poly(y)
    b = ((12*z**12*x**8 + 60*z**12*x**7 + 72*z**12*x**6 + 24*z**12*x**5 +
          108*z**11*x**8 + 540*z**11*x**7 + 660*z**11*x**6 + 240*z**11*x**5 +
          12*z**11*x**4 + 414*z**10*x**8 + 2076*z**10*x**7 + 2580*z**10*x**6 +
          1008*z**10*x**5 + 84*z**10*x**4 + 876*z**9*x**8 + 4608*z**9*x**7 +
          5824*z**9*x**6 + 2400*z**9*x**5 + 240*z**9*x**4 + 1194*z**8*x**8 +
          6876*z**8*x**7 + 8820*z**8*x**6 + 3816*z**8*x**5 + 408*z**8*x**4 +
          1203*z**7*x**8 + 7596*z**7*x**7 + 9912*z**7*x**6 + 4560*z**7*x**5 +
          540*z**7*x**4 + 897*z**6*x**8 + 6348*z**6*x**7 + 8474*z**6*x**6 +
          4176*z**6*x**5 + 564*z**6*x**4 + 438*z**5*x**8 + 3816*z**5*x**7 +
          5236*z**5*x**6 + 2784*z**5*x**5 + 408*z**5*x**4 + 96*z**4*x**8 +
          1488*z**4*x**7 + 2170*z**4*x**6 + 1344*z**4*x**5 + 240*z**4*x**4 +
          288*z**3*x**7 + 484*z**3*x**6 + 384*z**3*x**5 +
          96*z**3*x**4)/(72*z**13 + 864*z**12 + 4680*z**11 + 15408*z**10 +
                         35280*z**9 + 61056*z**8 + 83952*z**7 + 93600*z**6 +
                         84744*z**5 + 61920*z**4 + 35784*z**3 + 15408*z**2 +
                         4320*z + 576)*y**2 +
         (4*z**6*x**4 + 12*z**6*x**3 + 4*z**6*x**2 + 20*z**5*x**4 +
          60*z**5*x**3 + 20*z**5*x**2 + 37*z**4*x**4 + 108*z**4*x**3 +
          36*z**4*x**2 + 36*z**3*x**4 + 108*z**3*x**3 + 36*z**3*x**2 +
          32*z**2*x**4 + 96*z**2*x**3 + 32*z**2*x**2 + 16*z*x**4 + 48*z*x**3 +
          16*z*x**2)/(12*z**7 + 84*z**6 + 240*z**5 + 384*z**4 + 420*z**3 +
                      348*z**2 + 192*z + 48)*y + 1).as_poly(y)
    assert a*b == ((12*x**12*z**12 + 108*x**12*z**11 + 414*x**12*z**10 +
                    876*x**12*z**9 + 1194*x**12*z**8 + 1203*x**12*z**7 +
                    897*x**12*z**6 + 438*x**12*z**5 + 96*x**12*z**4 +
                    60*x**11*z**12 + 540*x**11*z**11 + 2076*x**11*z**10 +
                    4608*x**11*z**9 + 6876*x**11*z**8 + 7596*x**11*z**7 +
                    6348*x**11*z**6 + 3816*x**11*z**5 + 1488*x**11*z**4 +
                    288*x**11*z**3 + 72*x**10*z**12 + 660*x**10*z**11 +
                    2580*x**10*z**10 + 5824*x**10*z**9 + 8820*x**10*z**8 +
                    9912*x**10*z**7 + 8474*x**10*z**6 + 5236*x**10*z**5 +
                    2170*x**10*z**4 + 484*x**10*z**3 + 24*x**9*z**12 +
                    240*x**9*z**11 + 1008*x**9*z**10 + 2400*x**9*z**9 +
                    3816*x**9*z**8 + 4560*x**9*z**7 + 4176*x**9*z**6 +
                    2784*x**9*z**5 + 1344*x**9*z**4 + 384*x**9*z**3 +
                    12*x**8*z**11 + 84*x**8*z**10 + 240*x**8*z**9 +
                    408*x**8*z**8 + 540*x**8*z**7 + 564*x**8*z**6 +
                    408*x**8*z**5 + 240*x**8*z**4 +
                    96*x**8*z**3)/(72*z**15 + 1080*z**14 + 7416*z**13 +
                                   31176*z**12 + 90864*z**11 + 197712*z**10 +
                                   337680*z**9 + 467568*z**8 + 533448*z**7 +
                                   503352*z**6 + 391032*z**5 + 246600*z**4 +
                                   122112*z**3 + 44352*z**2 + 10368*z +
                                   1152)*y**3 +
                   (4*x**8*z**6 + 20*x**8*z**5 + 37*x**8*z**4 + 36*x**8*z**3 +
                    32*x**8*z**2 + 16*x**8*z + 12*x**7*z**6 + 60*x**7*z**5 +
                    108*x**7*z**4 + 108*x**7*z**3 + 96*x**7*z**2 + 48*x**7*z +
                    4*x**6*z**6 + 20*x**6*z**5 + 36*x**6*z**4 + 36*x**6*z**3 +
                    32*x**6*z**2 + 16*x**6*z)/(12*z**9 + 120*z**8 + 516*z**7 +
                                               1272*z**6 + 2052*z**5 + 2376*z**4 +
                                               2076*z**3 + 1320*z**2 +
                                               528*z + 96)*y**2 +
                   x**4/(z**2 + 3*z + 2)*y).as_poly(y)
