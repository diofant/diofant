"""Tests for user-friendly public interface to polynomial functions. """

import pytest

from diofant import (Derivative, Eq, Expr, Float, I, Integer, Integral, Mul,
                     Piecewise, Rational, RootOf, Sum, Symbol, Tuple, diff,
                     exp, expand, false, im, oo, pi, re, root, sin, sqrt,
                     symbols, tanh, true)
from diofant.abc import a, b, c, d, p, q, t, w, x, y, z
from diofant.core.mul import _keep_coeff
from diofant.domains import EX, FF, QQ, RR, ZZ
from diofant.domains.realfield import RealField
from diofant.matrices import MatrixSymbol
from diofant.polys.orderings import grevlex, grlex, lex
from diofant.polys.polyerrors import (CoercionFailed, ComputationFailed,
                                      DomainError, ExactQuotientFailed,
                                      FlagError, GeneratorsError,
                                      GeneratorsNeeded,
                                      MultivariatePolynomialError, OptionError,
                                      PolificationFailed, PolynomialError,
                                      RefinementFailed, UnificationFailed)
from diofant.polys.polytools import (LC, LM, LT, GroebnerBasis, Poly, PurePoly,
                                     _torational_factor_list, cancel,
                                     cofactors, compose, content, count_roots,
                                     decompose, degree, degree_list,
                                     discriminant, div, exquo, factor,
                                     factor_list, gcd, gcd_list, gcdex,
                                     groebner, ground_roots, half_gcdex,
                                     intervals, invert, lcm, lcm_list, monic,
                                     nroots, nth_power_roots_poly,
                                     parallel_poly_from_expr, poly, prem,
                                     primitive, quo, real_roots, reduced,
                                     refine_root, rem, resultant, sqf,
                                     sqf_list, sqf_norm, sqf_part, sturm,
                                     subresultants, terms_gcd,
                                     to_rational_coeffs, trunc)
from diofant.polys.rings import ring


__all__ = ()


def _epsilon_eq(a, b):
    for x, y in zip(a, b):
        if abs(x - y) > 1e-10:
            return False
    return True


def test_Poly_from_dict():
    K = FF(3)

    assert Poly.from_dict({0: 1, 1: 2}, gens=x,
                          domain=K).rep.to_dense() == [K(2), K(1)]
    assert Poly.from_dict({0: 1, 1: 5}, gens=x,
                          domain=K).rep.to_dense() == [K(2), K(1)]

    assert Poly.from_dict({(0,): 1, (1,): 2}, gens=x,
                          domain=K).rep.to_dense() == [K(2), K(1)]
    assert Poly.from_dict({(0,): 1, (1,): 5}, gens=x,
                          domain=K).rep.to_dense() == [K(2), K(1)]

    assert Poly.from_dict({(0, 0): 1, (1, 1): 2}, gens=(x, y),
                          domain=K).rep.to_dense() == [[K(2), K(0)], [K(1)]]

    assert Poly.from_dict({0: 1, 1: 2}, gens=x).rep.to_dense() == [ZZ(2), ZZ(1)]
    assert Poly.from_dict({0: 1, 1: 2}, gens=x,
                          field=True).rep.to_dense() == [QQ(2), QQ(1)]

    assert Poly.from_dict({0: 1, 1: 2}, gens=x,
                          domain=ZZ).rep.to_dense() == [ZZ(2), ZZ(1)]
    assert Poly.from_dict({0: 1, 1: 2}, gens=x,
                          domain=QQ).rep.to_dense() == [QQ(2), QQ(1)]

    assert Poly.from_dict({(0,): 1, (1,): 2},
                          gens=x).rep.to_dense() == [ZZ(2), ZZ(1)]
    assert Poly.from_dict({(0,): 1, (1,): 2}, gens=x,
                          field=True).rep.to_dense() == [QQ(2), QQ(1)]

    assert Poly.from_dict({(0,): 1, (1,): 2}, gens=x,
                          domain=ZZ).rep.to_dense() == [ZZ(2), ZZ(1)]
    assert Poly.from_dict({(0,): 1, (1,): 2}, gens=x,
                          domain=QQ).rep.to_dense() == [QQ(2), QQ(1)]

    assert Poly.from_dict({(1,): sin(y)}, gens=x, composite=False) == \
        Poly(sin(y)*x, x, domain='EX')
    assert Poly.from_dict({(1,): y}, gens=x, composite=False) == \
        Poly(y*x, x, domain='EX')
    assert Poly.from_dict({(1, 1): 1}, gens=(x, y), composite=False) == \
        Poly(x*y, x, y, domain='ZZ')
    assert Poly.from_dict({(1, 0): y}, gens=(x, z), composite=False) == \
        Poly(y*x, x, z, domain='EX')

    pytest.raises(GeneratorsError,
                  lambda: Poly.from_dict({(1,): x, (0,): 1}, gens=(x,)))


def test_Poly_from_list():
    K = FF(3)

    assert Poly.from_list([2, 1], gens=x, domain=K).rep.to_dense() == [K(2), K(1)]
    assert Poly.from_list([5, 1], gens=x, domain=K).rep.to_dense() == [K(2), K(1)]

    assert Poly.from_list([2, 1], gens=x).rep.to_dense() == [ZZ(2), ZZ(1)]
    assert Poly.from_list([2, 1], gens=x, field=True).rep.to_dense() == [QQ(2), QQ(1)]

    assert Poly.from_list([2, 1], gens=x, domain=ZZ).rep.to_dense() == [ZZ(2), ZZ(1)]
    assert Poly.from_list([2, 1], gens=x, domain=QQ).rep.to_dense() == [QQ(2), QQ(1)]

    assert Poly.from_list([0, 1.0], gens=x).rep.to_dense() == [RR(1.0)]
    assert Poly.from_list([1.0, 0], gens=x).rep.to_dense() == [RR(1.0), RR(0.0)]

    pytest.raises(MultivariatePolynomialError, lambda: Poly.from_list([[]], gens=(x, y)))
    pytest.raises(GeneratorsError, lambda: Poly.from_list([x, 1], gens=(x,)))


def test_Poly_from_poly():
    f = Poly(x + 7, x, domain=ZZ)
    g = Poly(x + 2, x, modulus=3)
    h = Poly(x + y, x, y, domain=ZZ)

    K = FF(3)

    assert Poly.from_poly(f) == f
    assert Poly.from_poly(f, domain=K).rep.to_dense() == [K(1), K(1)]
    assert Poly.from_poly(f, domain=ZZ).rep.to_dense() == [ZZ(1), ZZ(7)]
    assert Poly.from_poly(f, domain=QQ).rep.to_dense() == [QQ(1), QQ(7)]

    assert Poly.from_poly(f, gens=x) == f
    assert Poly.from_poly(f, gens=x, domain=K).rep.to_dense() == [K(1), K(1)]
    assert Poly.from_poly(f, gens=x, domain=ZZ).rep.to_dense() == [ZZ(1), ZZ(7)]
    assert Poly.from_poly(f, gens=x, domain=QQ).rep.to_dense() == [QQ(1), QQ(7)]

    assert Poly.from_poly(f, gens=y) == Poly(x + 7, y, domain='ZZ[x]')
    pytest.raises(CoercionFailed, lambda: Poly.from_poly(f, gens=y, domain=K))
    pytest.raises(CoercionFailed, lambda: Poly.from_poly(f, gens=y, domain=ZZ))
    pytest.raises(CoercionFailed, lambda: Poly.from_poly(f, gens=y, domain=QQ))

    assert Poly.from_poly(f, gens=(x, y)) == Poly(x + 7, x, y, domain='ZZ')
    assert Poly.from_poly(
        f, gens=(x, y), domain=ZZ) == Poly(x + 7, x, y, domain='ZZ')
    assert Poly.from_poly(
        f, gens=(x, y), domain=QQ) == Poly(x + 7, x, y, domain='QQ')
    assert Poly.from_poly(
        f, gens=(x, y), modulus=3) == Poly(x + 7, x, y, domain='FF(3)')

    K = FF(2)

    assert Poly.from_poly(g) == g
    assert Poly.from_poly(g, domain=ZZ).rep.to_dense() == [ZZ(1), ZZ(2)]
    pytest.raises(CoercionFailed, lambda: Poly.from_poly(g, domain=QQ))
    assert Poly.from_poly(g, domain=K).rep.to_dense() == [K(1), K(0)]

    assert Poly.from_poly(g, gens=x) == g
    assert Poly.from_poly(g, gens=x, domain=ZZ).rep.to_dense() == [ZZ(1), ZZ(2)]
    pytest.raises(CoercionFailed, lambda: Poly.from_poly(g, gens=x, domain=QQ))
    assert Poly.from_poly(g, gens=x, domain=K).rep.to_dense() == [K(1), K(0)]

    K = FF(3)

    assert Poly.from_poly(h) == h
    assert Poly.from_poly(h, domain=ZZ).rep.to_dense() == [[ZZ(1)], [ZZ(1), ZZ(0)]]
    assert Poly.from_poly(h, domain=QQ).rep.to_dense() == [[QQ(1)], [QQ(1), QQ(0)]]
    assert Poly.from_poly(h, domain=K).rep.to_dense() == [[K(1)], [K(1), K(0)]]

    assert Poly.from_poly(h, gens=x) == Poly(x + y, x, domain=ZZ.poly_ring(y))
    pytest.raises(CoercionFailed, lambda: Poly.from_poly(h, gens=x, domain=ZZ))
    assert Poly.from_poly(
        h, gens=x, domain=ZZ.poly_ring(y)) == Poly(x + y, x, domain=ZZ.poly_ring(y))
    pytest.raises(CoercionFailed, lambda: Poly.from_poly(h, gens=x, domain=QQ))
    assert Poly.from_poly(
        h, gens=x, domain=QQ.poly_ring(y)) == Poly(x + y, x, domain=QQ.poly_ring(y))
    pytest.raises(CoercionFailed, lambda: Poly.from_poly(h, gens=x, modulus=3))

    assert Poly.from_poly(h, gens=y) == Poly(x + y, y, domain=ZZ.poly_ring(x))
    pytest.raises(CoercionFailed, lambda: Poly.from_poly(h, gens=y, domain=ZZ))
    assert Poly.from_poly(
        h, gens=y, domain=ZZ.poly_ring(x)) == Poly(x + y, y, domain=ZZ.poly_ring(x))
    pytest.raises(CoercionFailed, lambda: Poly.from_poly(h, gens=y, domain=QQ))
    assert Poly.from_poly(
        h, gens=y, domain=QQ.poly_ring(x)) == Poly(x + y, y, domain=QQ.poly_ring(x))
    pytest.raises(CoercionFailed, lambda: Poly.from_poly(h, gens=y, modulus=3))

    assert Poly.from_poly(h, gens=(x, y)) == h
    assert Poly.from_poly(h, gens=(x, y),
                          domain=ZZ).rep.to_dense() == [[ZZ(1)], [ZZ(1), ZZ(0)]]
    assert Poly.from_poly(h, gens=(x, y),
                          domain=QQ).rep.to_dense() == [[QQ(1)], [QQ(1), QQ(0)]]
    assert Poly.from_poly(h, gens=(x, y),
                          domain=K).rep.to_dense() == [[K(1)], [K(1), K(0)]]

    assert Poly.from_poly(h, gens=(y, x)).rep.to_dense() == [[ZZ(1)], [ZZ(1), ZZ(0)]]
    assert Poly.from_poly(h, gens=(y, x),
                          domain=ZZ).rep.to_dense() == [[ZZ(1)], [ZZ(1), ZZ(0)]]
    assert Poly.from_poly(h, gens=(y, x),
                          domain=QQ).rep.to_dense() == [[QQ(1)], [QQ(1), QQ(0)]]
    assert Poly.from_poly(h, gens=(y, x),
                          domain=K).rep.to_dense() == [[K(1)], [K(1), K(0)]]

    assert Poly.from_poly(h, gens=(x, y),
                          field=True).rep.to_dense() == [[QQ(1)], [QQ(1), QQ(0)]]
    assert Poly.from_poly(h, gens=(x, y),
                          field=True).rep.to_dense() == [[QQ(1)], [QQ(1), QQ(0)]]


def test_Poly_from_expr():
    pytest.raises(GeneratorsNeeded, lambda: Poly.from_expr(Integer(0)))
    pytest.raises(GeneratorsNeeded, lambda: Poly.from_expr(Integer(7)))

    F3 = FF(3)

    assert Poly.from_expr(x + 5, domain=F3).rep.to_dense() == [F3(1), F3(2)]
    assert Poly.from_expr(y + 5, domain=F3).rep.to_dense() == [F3(1), F3(2)]

    assert Poly.from_expr(x + 5, x, domain=F3).rep.to_dense() == [F3(1), F3(2)]
    assert Poly.from_expr(y + 5, y, domain=F3).rep.to_dense() == [F3(1), F3(2)]

    assert Poly.from_expr(x + y, domain=F3).rep.to_dense() == [[F3(1)], [F3(1), F3(0)]]
    assert Poly.from_expr(x + y, x, y, domain=F3).rep.to_dense() == [[F3(1)], [F3(1), F3(0)]]

    assert Poly.from_expr(x + 5).rep.to_dense() == [ZZ(1), ZZ(5)]
    assert Poly.from_expr(y + 5).rep.to_dense() == [ZZ(1), ZZ(5)]

    assert Poly.from_expr(x + 5, x).rep.to_dense() == [ZZ(1), ZZ(5)]
    assert Poly.from_expr(y + 5, y).rep.to_dense() == [ZZ(1), ZZ(5)]

    assert Poly.from_expr(x + 5, domain=ZZ).rep.to_dense() == [ZZ(1), ZZ(5)]
    assert Poly.from_expr(y + 5, domain=ZZ).rep.to_dense() == [ZZ(1), ZZ(5)]

    assert Poly.from_expr(x + 5, x, domain=ZZ).rep.to_dense() == [ZZ(1), ZZ(5)]
    assert Poly.from_expr(y + 5, y, domain=ZZ).rep.to_dense() == [ZZ(1), ZZ(5)]

    assert Poly.from_expr(x + 5, x, y, domain=ZZ).rep.to_dense() == [[ZZ(1)], [ZZ(5)]]
    assert Poly.from_expr(y + 5, x, y, domain=ZZ).rep.to_dense() == [[ZZ(1), ZZ(5)]]


def test_Poly__new__():
    pytest.raises(GeneratorsError, lambda: Poly(x + 1, x, x))

    pytest.raises(GeneratorsError, lambda: Poly(x + y, x, y, domain=ZZ.poly_ring(x)))
    pytest.raises(GeneratorsError, lambda: Poly(x + y, x, y, domain=ZZ.poly_ring(y)))

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

    pytest.raises(NotImplementedError, lambda: Poly(x + 1, x, modulus=3, order='grlex'))
    pytest.raises(NotImplementedError, lambda: Poly(x + 1, x, order='grlex'))

    pytest.raises(GeneratorsNeeded, lambda: Poly({1: 2, 0: 1}))
    pytest.raises(GeneratorsNeeded, lambda: Poly([2, 1]))
    pytest.raises(GeneratorsNeeded, lambda: Poly((2, 1)))

    pytest.raises(GeneratorsNeeded, lambda: Poly(1))

    f = a*x**2 + b*x + c

    assert Poly({2: a, 1: b, 0: c}, x) == f
    assert Poly(iter([a, b, c]), x) == f
    assert Poly([a, b, c], x) == f
    assert Poly((a, b, c), x) == f

    f = Poly({}, x, y, z)

    assert f.gens == (x, y, z) and f.as_expr() == 0

    assert Poly(Poly(a*x + b*y, x, y), x) == Poly(a*x + b*y, x)

    assert Poly(3*x**2 + 2*x + 1, domain='ZZ').all_coeffs() == [3, 2, 1]
    assert Poly(3*x**2 + 2*x + 1, domain='QQ').all_coeffs() == [3, 2, 1]
    assert Poly(3*x**2 + 2*x + 1, domain='RR').all_coeffs() == [3.0, 2.0, 1.0]

    pytest.raises(CoercionFailed, lambda: Poly(3*x**2/5 + 2*x/5 + 1, domain='ZZ'))
    assert Poly(
        3*x**2/5 + 2*x/5 + 1, domain='QQ').all_coeffs() == [Rational(3, 5), Rational(2, 5), 1]
    assert _epsilon_eq(
        Poly(3*x**2/5 + 2*x/5 + 1, domain='RR').all_coeffs(), [0.6, 0.4, 1.0])

    assert Poly(3.0*x**2 + 2.0*x + 1, domain='ZZ').all_coeffs() == [3, 2, 1]
    assert Poly(3.0*x**2 + 2.0*x + 1, domain='QQ').all_coeffs() == [3, 2, 1]
    assert Poly(
        3.0*x**2 + 2.0*x + 1, domain='RR').all_coeffs() == [3.0, 2.0, 1.0]

    pytest.raises(CoercionFailed, lambda: Poly(3.1*x**2 + 2.1*x + 1, domain='ZZ'))
    assert Poly(3.1*x**2 + 2.1*x + 1, domain='QQ').all_coeffs() == [Rational(31, 10), Rational(21, 10), 1]
    assert Poly(3.1*x**2 + 2.1*x + 1, domain='RR').all_coeffs() == [3.1, 2.1, 1.0]

    assert Poly({(2, 1): 1, (1, 2): 2, (1, 1): 3}, x, y) == \
        Poly(x**2*y + 2*x*y**2 + 3*x*y, x, y)

    assert Poly(x**2 + 1, extension=I).domain == QQ.algebraic_field(I)

    f = 3*x**5 - x**4 + x**3 - x**2 + 65538

    assert Poly(f, x, modulus=65537) == \
        Poly(3*x**5 + 65536*x**4 + x**3 + 65536*x**2 + 1, x, modulus=65537)

    assert isinstance(Poly(x**2 + x + 1.0).domain, RealField)


def test_Poly_new():
    pytest.raises(PolynomialError, lambda: Poly.new([1], x))
    pytest.raises(PolynomialError, lambda: Poly.new(Poly(x + 1, x, y).rep, x))


def test_Poly__args():
    assert Poly(x**2 + 1).args == (x**2 + 1, x)


def test_Poly_is_number():
    assert Poly(1, x).is_number
    assert Poly(x, x).is_number is False


def test_Poly__gens():
    assert Poly((x - p)*(x - q), x).gens == (x,)
    assert Poly((x - p)*(x - q), p).gens == (p,)
    assert Poly((x - p)*(x - q), q).gens == (q,)

    assert Poly((x - p)*(x - q), x, p).gens == (x, p)
    assert Poly((x - p)*(x - q), x, q).gens == (x, q)

    assert Poly((x - p)*(x - q), x, p, q).gens == (x, p, q)
    assert Poly((x - p)*(x - q), p, x, q).gens == (p, x, q)
    assert Poly((x - p)*(x - q), p, q, x).gens == (p, q, x)

    assert Poly((x - p)*(x - q)).gens == (x, p, q)

    assert Poly((x - p)*(x - q), sort='x > p > q').gens == (x, p, q)
    assert Poly((x - p)*(x - q), sort='p > x > q').gens == (p, x, q)
    assert Poly((x - p)*(x - q), sort='p > q > x').gens == (p, q, x)

    assert Poly((x - p)*(x - q), x, p, q, sort='p > q > x').gens == (x, p, q)

    assert Poly((x - p)*(x - q), wrt='x').gens == (x, p, q)
    assert Poly((x - p)*(x - q), wrt='p').gens == (p, x, q)
    assert Poly((x - p)*(x - q), wrt='q').gens == (q, x, p)

    assert Poly((x - p)*(x - q), wrt=x).gens == (x, p, q)
    assert Poly((x - p)*(x - q), wrt=p).gens == (p, x, q)
    assert Poly((x - p)*(x - q), wrt=q).gens == (q, x, p)

    assert Poly((x - p)*(x - q), x, p, q, wrt='p').gens == (x, p, q)

    assert Poly((x - p)*(x - q), wrt='p', sort='q > x').gens == (p, q, x)
    assert Poly((x - p)*(x - q), wrt='q', sort='p > x').gens == (q, p, x)


def test_Poly_zero():
    assert Poly(x).zero == Poly(0, x, domain=ZZ)
    assert Poly(x/2).zero == Poly(0, x, domain=QQ)


def test_Poly_one():
    assert Poly(x).one == Poly(1, x, domain=ZZ)
    assert Poly(x/2).one == Poly(1, x, domain=QQ)


def test_Poly_unify():
    pytest.raises(UnificationFailed, lambda: Poly(x).unify(y))
    pytest.raises(UnificationFailed, lambda: PurePoly(x).unify(y))
    pytest.raises(UnificationFailed, lambda: PurePoly(x).unify(Poly(x, x, y)))

    assert Poly(x, x, modulus=3).unify(Poly(y, y, modulus=3)) == \
        (Poly(x, x, y, modulus=3), Poly(y, x, y, modulus=3))
    pytest.raises(NotImplementedError,
                  lambda: Poly(x, x, modulus=3).unify(Poly(y, y, modulus=5)))

    assert Poly(y, x, y).unify(Poly(x, x, modulus=3)) == \
        (Poly(y, x, y, modulus=3), Poly(x, x, y, modulus=3))
    assert Poly(x, x, modulus=3).unify(Poly(y, x, y)) == \
        (Poly(x, x, y, modulus=3), Poly(y, x, y, modulus=3))

    assert Poly(x + 1, x).unify(Poly(x + 2, x)) == \
        (Poly(x + 1, x, domain=ZZ), Poly(x + 2, x, domain=ZZ))
    assert Poly(x + 1, x, domain=QQ).unify(Poly(x + 2, x)) == \
        (Poly(x + 1, x, domain=QQ), Poly(x + 2, x, domain=QQ))
    assert Poly(x + 1, x).unify(Poly(x + 2, x, domain=QQ)) == \
        (Poly(x + 1, x, domain=QQ), Poly(x + 2, x, domain=QQ))

    assert Poly(x + 1, x).unify(Poly(x + 2, x, y)) == \
        (Poly(x + 1, x, y, domain=ZZ), Poly(x + 2, x, y, domain=ZZ))
    assert Poly(x + 1, x, domain=QQ).unify(Poly(x + 2, x, y)) == \
        (Poly(x + 1, x, y, domain=QQ), Poly(x + 2, x, y, domain=QQ))
    assert Poly(x + 1, x).unify(Poly(x + 2, x, y, domain=QQ)) == \
        (Poly(x + 1, x, y, domain=QQ), Poly(x + 2, x, y, domain=QQ))

    assert Poly(x + 1, x, y).unify(Poly(x + 2, x)) == \
        (Poly(x + 1, x, y, domain=ZZ), Poly(x + 2, x, y, domain=ZZ))
    assert Poly(x + 1, x, y, domain=QQ).unify(Poly(x + 2, x)) == \
        (Poly(x + 1, x, y, domain=QQ), Poly(x + 2, x, y, domain=QQ))
    assert Poly(x + 1, x, y).unify(Poly(x + 2, x, domain=QQ)) == \
        (Poly(x + 1, x, y, domain=QQ), Poly(x + 2, x, y, domain=QQ))

    assert Poly(x + 1, x, y).unify(Poly(x + 2, x, y)) == \
        (Poly(x + 1, x, y, domain=ZZ), Poly(x + 2, x, y, domain=ZZ))
    assert Poly(x + 1, x, y, domain=QQ).unify(Poly(x + 2, x, y)) == \
        (Poly(x + 1, x, y, domain=QQ), Poly(x + 2, x, y, domain=QQ))
    assert Poly(x + 1, x, y).unify(Poly(x + 2, x, y, domain=QQ)) == \
        (Poly(x + 1, x, y, domain=QQ), Poly(x + 2, x, y, domain=QQ))

    assert Poly(x + 1, x).unify(Poly(x + 2, y, x)) == \
        (Poly(x + 1, y, x, domain=ZZ), Poly(x + 2, y, x, domain=ZZ))
    assert Poly(x + 1, x, domain=QQ).unify(Poly(x + 2, y, x)) == \
        (Poly(x + 1, y, x, domain=QQ), Poly(x + 2, y, x, domain=QQ))
    assert Poly(x + 1, x).unify(Poly(x + 2, y, x, domain=QQ)) == \
        (Poly(x + 1, y, x, domain=QQ), Poly(x + 2, y, x, domain=QQ))

    assert Poly(x + 1, y, x).unify(Poly(x + 2, x)) == \
        (Poly(x + 1, y, x, domain=ZZ), Poly(x + 2, y, x, domain=ZZ))
    assert Poly(x + 1, y, x, domain=QQ).unify(Poly(x + 2, x)) == \
        (Poly(x + 1, y, x, domain=QQ), Poly(x + 2, y, x, domain=QQ))
    assert Poly(x + 1, y, x).unify(Poly(x + 2, x, domain=QQ)) == \
        (Poly(x + 1, y, x, domain=QQ), Poly(x + 2, y, x, domain=QQ))

    assert Poly(x + 1, x, y).unify(Poly(x + 2, y, x)) == \
        (Poly(x + 1, x, y, domain=ZZ), Poly(x + 2, x, y, domain=ZZ))
    assert Poly(x + 1, x, y, domain=QQ).unify(Poly(x + 2, y, x)) == \
        (Poly(x + 1, x, y, domain=QQ), Poly(x + 2, x, y, domain=QQ))
    assert Poly(x + 1, x, y).unify(Poly(x + 2, y, x, domain=QQ)) == \
        (Poly(x + 1, x, y, domain=QQ), Poly(x + 2, x, y, domain=QQ))

    assert Poly(x + 1, y, x).unify(Poly(x + 2, x, y)) == \
        (Poly(x + 1, y, x, domain=ZZ), Poly(x + 2, y, x, domain=ZZ))
    assert Poly(x + 1, y, x, domain=QQ).unify(Poly(x + 2, x, y)) == \
        (Poly(x + 1, y, x, domain=QQ), Poly(x + 2, y, x, domain=QQ))
    assert Poly(x + 1, y, x).unify(Poly(x + 2, x, y, domain=QQ)) == \
        (Poly(x + 1, y, x, domain=QQ), Poly(x + 2, y, x, domain=QQ))

    assert Poly(a*x, x, domain=ZZ.poly_ring(a)).unify(Poly(a*b*x, x, domain=ZZ.frac_field(a, b))) == \
        (Poly(a*x, x, domain=ZZ.frac_field(a, b)), Poly(a*b*x, x, domain=ZZ.frac_field(a, b)))

    assert Poly(a*x, x, domain=ZZ.frac_field(a)).unify(Poly(a*b*x, x, domain=ZZ.frac_field(a, b))) == \
        (Poly(a*x, x, domain=ZZ.frac_field(a, b)), Poly(a*b*x, x, domain=ZZ.frac_field(a, b)))

    pytest.raises(CoercionFailed, lambda: Poly(Poly(x**2 + x**2*z, y, field=True),
                                               domain=ZZ.frac_field(x)))

    f = Poly(t**2 + t/3 + x, t, domain=QQ.frac_field(x))
    g = Poly(t**2 + t/3 + x, t, domain=QQ.poly_ring(x))

    assert f.unify(g) == (f, f)


def test_Poly_free_symbols():
    assert Poly(x**2 + 1).free_symbols == {x}
    assert Poly(x**2 + y*z).free_symbols == {x, y, z}
    assert Poly(x**2 + y*z, x).free_symbols == {x, y, z}
    assert Poly(x**2 + sin(y*z)).free_symbols == {x, y, z}
    assert Poly(x**2 + sin(y*z), x).free_symbols == {x, y, z}
    assert Poly(x**2 + sin(y*z), x, domain=EX).free_symbols == {x, y, z}


def test_PurePoly_free_symbols():
    assert PurePoly(x**2 + 1).free_symbols == set()
    assert PurePoly(x**2 + y*z).free_symbols == set()
    assert PurePoly(x**2 + y*z, x).free_symbols == {y, z}
    assert PurePoly(x**2 + sin(y*z)).free_symbols == set()
    assert PurePoly(x**2 + sin(y*z), x).free_symbols == {y, z}
    assert PurePoly(x**2 + sin(y*z), x, domain=EX).free_symbols == {y, z}


def test_Poly__eq__():
    assert (Poly(x, x) == Poly(x, x)) is True
    assert (Poly(x, x, domain=QQ) == Poly(x, x)) is True
    assert (Poly(x, x) == Poly(x, x, domain=QQ)) is True

    assert (Poly(x, x, domain=ZZ.poly_ring(a)) == Poly(x, x)) is True
    assert (Poly(x, x) == Poly(x, x, domain=ZZ.poly_ring(a))) is True

    assert (Poly(x*y, x, y) == Poly(x, x)) is False

    assert (Poly(x, x, y) == Poly(x, x)) is False
    assert (Poly(x, x) == Poly(x, x, y)) is False

    assert (Poly(x**2 + 1, x) == Poly(y**2 + 1, y)) is False
    assert (Poly(y**2 + 1, y) == Poly(x**2 + 1, x)) is False

    f = Poly(x, x, domain=ZZ)
    g = Poly(x, x, domain=QQ)

    assert (f == g) is True
    assert (f != g) is False

    t0 = Symbol('t0')

    f = Poly((t0/2 + x**2)*t**2 - x**2*t, t, domain='QQ[x,t0]')
    g = Poly((t0/2 + x**2)*t**2 - x**2*t, t, domain='ZZ(x,t0)')

    assert (f == g) is True


def test_PurePoly__eq__():
    assert (PurePoly(x, x) == PurePoly(x, x)) is True
    assert (PurePoly(x, x, domain=QQ) == PurePoly(x, x)) is True
    assert (PurePoly(x, x) == PurePoly(x, x, domain=QQ)) is True

    assert (PurePoly(x, x, domain=ZZ.poly_ring(a)) == PurePoly(x, x)) is True
    assert (PurePoly(x, x) == PurePoly(x, x, domain=ZZ.poly_ring(a))) is True

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


def test_PurePoly_Poly():
    assert isinstance(PurePoly(Poly(x**2 + 1)), PurePoly) is True
    assert isinstance(Poly(PurePoly(x**2 + 1)), Poly) is True


def test_Poly_domain():
    assert Poly(2*x).domain == ZZ

    assert Poly(2*x, domain='ZZ').domain == ZZ
    assert Poly(2*x, domain='QQ').domain == QQ

    assert Poly(x/2).domain == QQ

    pytest.raises(CoercionFailed, lambda: Poly(x/2, domain='ZZ'))
    assert Poly(x/2, domain='QQ').domain == QQ

    assert isinstance(Poly(0.2*x).domain, RealField)


def test_Poly_set_domain():
    assert Poly(2*x + 1).set_domain(ZZ) == Poly(2*x + 1)
    assert Poly(2*x + 1).set_domain('ZZ') == Poly(2*x + 1)

    assert Poly(2*x + 1).set_domain(QQ) == Poly(2*x + 1, domain='QQ')
    assert Poly(2*x + 1).set_domain('QQ') == Poly(2*x + 1, domain='QQ')

    assert Poly(Rational(2, 10)*x + Rational(1, 10)).set_domain('RR') == Poly(0.2*x + 0.1)
    assert Poly(0.2*x + 0.1).set_domain('QQ') == Poly(Rational(2, 10)*x + Rational(1, 10))

    pytest.raises(CoercionFailed, lambda: Poly(x/2 + 1).set_domain(ZZ))
    pytest.raises(CoercionFailed, lambda: Poly(x + 1, modulus=2).set_domain(QQ))

    pytest.raises(GeneratorsError, lambda: Poly(x*y, x, y).set_domain(ZZ.poly_ring(y)))


def test_Poly_get_modulus():
    assert Poly(x**2 + 1, modulus=2).get_modulus() == 2
    pytest.raises(PolynomialError, lambda: Poly(x**2 + 1).get_modulus())


def test_Poly_set_modulus():
    assert Poly(
        x**2 + 1, modulus=2).set_modulus(7) == Poly(x**2 + 1, modulus=7)
    assert Poly(
        x**2 + 5, modulus=7).set_modulus(2) == Poly(x**2 + 1, modulus=2)

    assert Poly(x**2 + 1).set_modulus(2) == Poly(x**2 + 1, modulus=2)

    pytest.raises(CoercionFailed, lambda: Poly(x/2 + 1).set_modulus(2))


def test_Poly_quo_ground():
    assert Poly(2*x + 4).quo_ground(2) == Poly(x + 2)
    assert Poly(2*x + 3).quo_ground(2) == Poly(x + 1)


def test_Poly_exquo_ground():
    assert Poly(2*x + 4).exquo_ground(2) == Poly(x + 2)
    pytest.raises(ExactQuotientFailed, lambda: Poly(2*x + 3).exquo_ground(2))


def test_Poly_abs():
    assert abs(Poly(-x + 1, x)) == Poly(x + 1, x)


def test_Poly_neg():
    assert -Poly(-x + 1, x) == Poly(x - 1, x)


def test_Poly_add():
    assert Poly(0, x) + Poly(0, x) == Poly(0, x)

    assert Poly(1, x) + Poly(0, x) == Poly(1, x)
    assert Poly(1, x, y) + Poly(0, x) == Poly(1, x, y)
    assert Poly(0, x, y) + Poly(1, x, y) == Poly(1, x, y)

    assert Poly(x + 1) + 2 == Poly(x + 3)
    assert Poly(x**2 + 1, x) + Poly(x - 2, x) == Poly(x**2 + x - 1, x)

    assert Poly(1, x) + x == Poly(x + 1, x)
    assert Poly(1, x) + sin(x) == 1 + sin(x)
    assert sin(x) + Poly(1, x) == sin(x) + 1

    assert Poly(x, x) + 1 == Poly(x + 1, x)
    assert 1 + Poly(x, x) == Poly(x + 1, x)


def test_Poly_sub():
    assert Poly(0, x) - Poly(0, x) == Poly(0, x)

    assert Poly(1, x) - Poly(0, x) == Poly(1, x)
    assert Poly(1, x) - 1 == Poly(0, x)
    assert Poly(1, x, y) - Poly(0, x) == Poly(1, x, y)
    assert Poly(0, x) - Poly(1, x, y) == Poly(-1, x, y)
    assert Poly(0, x, y) - Poly(1, x, y) == Poly(-1, x, y)

    assert Poly(1, x) - x == Poly(1 - x, x)
    assert Poly(1, x) - sin(x) == 1 - sin(x)
    assert sin(x) - Poly(1, x) == sin(x) - 1

    assert Poly(x, x) - 1 == Poly(x - 1, x)
    assert 1 - Poly(x, x) == Poly(1 - x, x)

    assert Poly(x + 1) - 2 == Poly(x - 1)
    assert Poly(x**2 + 1, x) - Poly(x - 2, x) == Poly(x**2 - x + 3)


def test_Poly_mul():
    assert Poly(0, x) * Poly(0, x) == Poly(0, x)

    assert Poly(2, x) * Poly(4, x) == Poly(8, x)
    assert Poly(2, x, y) * Poly(4, x) == Poly(8, x, y)
    assert Poly(4, x) * Poly(2, x, y) == Poly(8, x, y)
    assert Poly(4, x, y) * Poly(2, x, y) == Poly(8, x, y)

    assert Poly(1, x) * x == Poly(x, x)
    assert Poly(1, x) * sin(x) == sin(x)
    assert sin(x) * Poly(1, x) == sin(x)

    assert Poly(x, x) * 2 == Poly(2*x, x)
    assert 2 * Poly(x, x) == Poly(2*x, x)

    assert Poly(x + 1) * 2 == Poly(2*x + 2)
    assert Poly(x**2 + 1, x) * Poly(x - 2, x) == Poly(x**3 - 2*x**2 + x - 2)


def test_Poly_pow():
    assert Poly(x, x)**10 == Poly(x**10, x)
    assert Poly(2*y, x, y)**4 == Poly(16*y**4, x, y)
    assert Poly(7*x*y, x, y)**3 == Poly(343*x**3*y**3, x, y)
    assert Poly(x*y + 1, x, y)**(-1) == (x*y + 1)**(-1)
    assert Poly(x*y + 1, x, y)**x == (x*y + 1)**x
    assert Poly(x - 2, x)**3 == Poly(x**3 - 6*x**2 + 12*x - 8)
    assert Poly(x*y, x, y)**2 == Poly(x**2*y**2, x, y)
    assert Poly(x - 2, x)**2 == Poly(x**2 - 4*x + 4, x)


def test_Poly_divmod():
    f, g = Poly(x**2), Poly(x)
    q, r = g, Poly(0, x)

    assert divmod(f, g) == (q, r)
    assert f // g == q
    assert f % g == r

    assert divmod(f, x) == (q, r)
    assert f // x == q
    assert f % x == r

    q, r = Poly(0, x), Poly(2, x)

    assert divmod(2, g) == (q, r)
    assert 2 // g == q
    assert 2 % g == r

    assert Poly(x)/Poly(x) == 1
    assert Poly(x**2)/Poly(x) == x
    assert Poly(x)/Poly(x**2) == 1/x


def test_Poly_eq_ne():
    assert (Poly(x + y, x, y) == Poly(x + y, x, y)) is True
    assert (Poly(x + y, x) == Poly(x + y, x, y)) is False
    assert (Poly(x + y, x, y) == Poly(x + y, x)) is False
    assert (Poly(x + y, x) == Poly(x + y, x)) is True
    assert (Poly(x + y, y) == Poly(x + y, y)) is True

    assert (Poly(x + y, x, y) == x + y) is True
    assert (Poly(x + y, x) == x + y) is True
    assert (Poly(x + y, x, y) == x + y) is True
    assert (Poly(x + y, x) == x + y) is True
    assert (Poly(x + y, y) == x + y) is True

    assert (Poly(x + y, x, y) != Poly(x + y, x, y)) is False
    assert (Poly(x + y, x) != Poly(x + y, x, y)) is True
    assert (Poly(x + y, x, y) != Poly(x + y, x)) is True
    assert (Poly(x + y, x) != Poly(x + y, x)) is False
    assert (Poly(x + y, y) != Poly(x + y, y)) is False

    assert (Poly(x + y, x, y) != x + y) is False
    assert (Poly(x + y, x) != x + y) is False
    assert (Poly(x + y, x, y) != x + y) is False
    assert (Poly(x + y, x) != x + y) is False
    assert (Poly(x + y, y) != x + y) is False

    assert (Poly(x, x) == sin(x)) is False
    assert (Poly(x, x) != sin(x)) is True


def test_Poly_nonzero():
    assert not bool(Poly(0, x)) is True
    assert not bool(Poly(1, x)) is False


def test_Poly_properties():
    assert Poly(0, x).is_zero is True
    assert Poly(1, x).is_zero is False

    assert Poly(1, x).is_one is True
    assert Poly(2, x).is_one is False

    assert Poly(x - 1, x).is_squarefree is True
    assert Poly((x - 1)**2, x).is_squarefree is False

    assert Poly(x - 1, x).is_monic is True
    assert Poly(2*x - 1, x).is_monic is False

    assert Poly(3*x + 2, x).is_primitive is True
    assert Poly(4*x + 2, x).is_primitive is False

    assert Poly(1, x).is_ground is True
    assert Poly(x, x).is_ground is False

    assert Poly(x + y + z + 1).is_linear is True
    assert Poly(x*y*z + 1).is_linear is False

    assert Poly(x*y + z + 1).is_quadratic is True
    assert Poly(x*y*z + 1).is_quadratic is False

    assert Poly(x*y).is_term is True
    assert Poly(x*y + 1).is_term is False

    assert Poly(x**2 + x*y).is_homogeneous is True
    assert Poly(x**3 + x*y).is_homogeneous is False

    assert Poly(x).is_univariate is True
    assert Poly(x*y).is_univariate is False

    assert Poly(x*y).is_multivariate is True
    assert Poly(x).is_multivariate is False

    assert Poly(
        x**16 + x**14 - x**10 + x**8 - x**6 + x**2 + 1).is_cyclotomic is False
    assert Poly(
        x**16 + x**14 - x**10 - x**8 - x**6 + x**2 + 1).is_cyclotomic is True


def test_Poly_is_irreducible():
    assert Poly(x**2 + x + 1).is_irreducible is True
    assert Poly(x**2 + 2*x + 1).is_irreducible is False

    assert Poly(7*x + 3, modulus=11).is_irreducible is True
    assert Poly(7*x**2 + 3*x + 1, modulus=11).is_irreducible is False


def test_Poly_subs():
    assert Poly(x + 1).subs({x: 0}) == 1

    assert Poly(x + 1).subs({x: x}) == Poly(x + 1)
    assert Poly(x + 1).subs({x: y}) == Poly(y + 1)

    assert Poly(x*y, x).subs({y: x}) == x**2
    assert Poly(x*y, x).subs({x: y}) == y**2


def test_Poly_replace():
    assert Poly(x + 1).replace(x) == Poly(x + 1)
    assert Poly(x + 1).replace(y) == Poly(y + 1)

    pytest.raises(PolynomialError, lambda: Poly(x + y).replace(z))

    assert Poly(x + 1).replace(x, x) == Poly(x + 1)
    assert Poly(x + 1).replace(x, y) == Poly(y + 1)

    assert Poly(x + y).replace(x, x) == Poly(x + y)
    assert Poly(x + y).replace(x, z) == Poly(z + y, z, y)

    assert Poly(x + y).replace(y, y) == Poly(x + y)
    assert Poly(x + y).replace(y, z) == Poly(x + z, x, z)

    pytest.raises(PolynomialError, lambda: Poly(x + y).replace(x, y))
    pytest.raises(PolynomialError, lambda: Poly(x + y).replace(z, t))

    assert Poly(x + y, x).replace(x, z) == Poly(z + y, z)
    assert Poly(x + y, y).replace(y, z) == Poly(x + z, z)

    pytest.raises(PolynomialError, lambda: Poly(x + y, x).replace(x, y))
    pytest.raises(PolynomialError, lambda: Poly(x + y, y).replace(y, x))


def test_Poly_reorder():
    pytest.raises(PolynomialError, lambda: Poly(x + y).reorder(x, z))

    assert Poly(x + y, x, y).reorder(x, y) == Poly(x + y, x, y)
    assert Poly(x + y, x, y).reorder(y, x) == Poly(x + y, y, x)

    assert Poly(x + y, y, x).reorder(x, y) == Poly(x + y, x, y)
    assert Poly(x + y, y, x).reorder(y, x) == Poly(x + y, y, x)

    assert Poly(x + y, x, y).reorder(wrt=x) == Poly(x + y, x, y)
    assert Poly(x + y, x, y).reorder(wrt=y) == Poly(x + y, y, x)


def test_Poly_ltrim():
    f = Poly(y**2 + y*z**2, x, y, z).ltrim(y)
    assert f.as_expr() == y**2 + y*z**2 and f.gens == (y, z)

    pytest.raises(PolynomialError, lambda: Poly(x*y**2 + y**2, x, y).ltrim(y))


def test_Poly_has_only_gens():
    assert Poly(x*y + 1, x, y, z).has_only_gens(x, y) is True
    assert Poly(x*y + z, x, y, z).has_only_gens(x, y) is False

    pytest.raises(GeneratorsError, lambda: Poly(x*y**2 + y**2, x, y).has_only_gens(t))


def test_Poly_to_ring():
    assert Poly(2*x + 1, domain='ZZ').to_ring() == Poly(2*x + 1, domain='ZZ')
    assert Poly(2*x + 1, domain='QQ').to_ring() == Poly(2*x + 1, domain='ZZ')

    pytest.raises(CoercionFailed, lambda: Poly(x/2 + 1).to_ring())
    pytest.raises(AttributeError, lambda: Poly(2*x + 1, modulus=3).to_ring())


def test_Poly_to_field():
    assert Poly(2*x + 1, domain='ZZ').to_field() == Poly(2*x + 1, domain='QQ')
    assert Poly(2*x + 1, domain='QQ').to_field() == Poly(2*x + 1, domain='QQ')

    assert Poly(x/2 + 1, domain='QQ').to_field() == Poly(x/2 + 1, domain='QQ')
    assert Poly(2*x + 1, modulus=3).to_field() == Poly(2*x + 1, modulus=3)

    assert Poly(2.0*x + 1.0).to_field() == Poly(2.0*x + 1.0)


def test_Poly_to_exact():
    assert Poly(2*x).to_exact() == Poly(2*x)
    assert Poly(x/2).to_exact() == Poly(x/2)

    assert Poly(0.1*x).to_exact() == Poly(x/10)


def test_Poly_retract():
    f = Poly(x**2 + 1, x, domain=QQ.poly_ring(y))

    assert f.retract() == Poly(x**2 + 1, x, domain='ZZ')
    assert f.retract(field=True) == Poly(x**2 + 1, x, domain='QQ')

    assert Poly(0, x, y).retract() == Poly(0, x, y)


def test_Poly_slice():
    f = Poly(x**3 + 2*x**2 + 3*x + 4)

    assert f.slice(0, 0) == Poly(0, x)
    assert f.slice(0, 1) == Poly(4, x)
    assert f.slice(0, 2) == Poly(3*x + 4, x)
    assert f.slice(0, 3) == Poly(2*x**2 + 3*x + 4, x)
    assert f.slice(0, 4) == Poly(x**3 + 2*x**2 + 3*x + 4, x)

    assert f.slice(x, 0, 0) == Poly(0, x)
    assert f.slice(x, 0, 1) == Poly(4, x)
    assert f.slice(x, 0, 2) == Poly(3*x + 4, x)
    assert f.slice(x, 0, 3) == Poly(2*x**2 + 3*x + 4, x)
    assert f.slice(x, 0, 4) == Poly(x**3 + 2*x**2 + 3*x + 4, x)


def test_Poly_coeffs():
    assert Poly(0, x).coeffs() == []
    assert Poly(1, x).coeffs() == [1]

    assert Poly(2*x + 1, x).coeffs() == [2, 1]

    assert Poly(7*x**2 + 2*x + 1, x).coeffs() == [7, 2, 1]
    assert Poly(7*x**4 + 2*x + 1, x).coeffs() == [7, 2, 1]

    assert Poly(x*y**7 + 2*x**2*y**3).coeffs('lex') == [2, 1]
    assert Poly(x*y**7 + 2*x**2*y**3).coeffs('grlex') == [1, 2]


def test_Poly_monoms():
    assert Poly(0, x).monoms() == []
    assert Poly(1, x).monoms() == [(0,)]

    assert Poly(2*x + 1, x).monoms() == [(1,), (0,)]

    assert Poly(7*x**2 + 2*x + 1, x).monoms() == [(2,), (1,), (0,)]
    assert Poly(7*x**4 + 2*x + 1, x).monoms() == [(4,), (1,), (0,)]

    assert Poly(x*y**7 + 2*x**2*y**3).monoms('lex') == [(2, 3), (1, 7)]
    assert Poly(x*y**7 + 2*x**2*y**3).monoms('grlex') == [(1, 7), (2, 3)]


def test_Poly_terms():
    assert Poly(0, x).terms() == []
    assert Poly(1, x).terms() == [((0,), 1)]

    assert Poly(2*x + 1, x).terms() == [((1,), 2), ((0,), 1)]

    assert Poly(7*x**2 + 2*x + 1, x).terms() == [((2,), 7), ((1,), 2), ((0,), 1)]
    assert Poly(7*x**4 + 2*x + 1, x).terms() == [((4,), 7), ((1,), 2), ((0,), 1)]

    assert Poly(
        x*y**7 + 2*x**2*y**3).terms('lex') == [((2, 3), 2), ((1, 7), 1)]
    assert Poly(
        x*y**7 + 2*x**2*y**3).terms('grlex') == [((1, 7), 1), ((2, 3), 2)]


def test_Poly_all_coeffs():
    assert Poly(0, x).all_coeffs() == [0]
    assert Poly(1, x).all_coeffs() == [1]

    assert Poly(2*x + 1, x).all_coeffs() == [2, 1]

    assert Poly(7*x**2 + 2*x + 1, x).all_coeffs() == [7, 2, 1]
    assert Poly(7*x**4 + 2*x + 1, x).all_coeffs() == [7, 0, 0, 2, 1]


def test_Poly_termwise():
    f = Poly(x**2 + 20*x + 400)
    g = Poly(x**2 + 2*x + 4)

    def func(monom, coeff):
        k, = monom
        return coeff//10**(2 - k)

    assert f.termwise(func) == g

    def func(monom, coeff):
        k, = monom
        return (k,), coeff//10**(2 - k)

    assert f.termwise(func) == g

    def func(monom, coeff):
        k, = monom
        return (k,), coeff // 2

    assert f.termwise(func) == Poly(10*x + 200)

    def func(monom, coeff):
        k, = monom
        return k % 2, coeff

    pytest.raises(PolynomialError, lambda: f.termwise(func))


def test_Poly_length():
    assert Poly(0, x).length() == 0
    assert Poly(1, x).length() == 1
    assert Poly(x, x).length() == 1

    assert Poly(x + 1, x).length() == 2
    assert Poly(x**2 + 1, x).length() == 2
    assert Poly(x**2 + x + 1, x).length() == 3


def test_Poly_as_dict():
    assert Poly(0, x).as_dict() == {}
    assert Poly(0, x, y, z).as_dict() == {}

    assert Poly(1, x).as_dict() == {(0,): 1}
    assert Poly(1, x, y, z).as_dict() == {(0, 0, 0): 1}

    assert Poly(x**2 + 3, x).as_dict() == {(2,): 1, (0,): 3}
    assert Poly(x**2 + 3, x, y, z).as_dict() == {(2, 0, 0): 1, (0, 0, 0): 3}

    assert Poly(3*x**2*y*z**3 + 4*x*y + 5*x*z).as_dict() == {(2, 1, 3): 3,
                                                             (1, 1, 0): 4, (1, 0, 1): 5}


def test_Poly_as_expr():
    assert Poly(0, x).as_expr() == 0
    assert Poly(0, x, y, z).as_expr() == 0

    assert Poly(1, x).as_expr() == 1
    assert Poly(1, x, y, z).as_expr() == 1

    assert Poly(x**2 + 3, x).as_expr() == x**2 + 3
    assert Poly(x**2 + 3, x, y, z).as_expr() == x**2 + 3

    assert Poly(
        3*x**2*y*z**3 + 4*x*y + 5*x*z).as_expr() == 3*x**2*y*z**3 + 4*x*y + 5*x*z

    f = Poly(x**2 + 2*x*y**2 - y, x, y)

    assert f.as_expr() == -y + x**2 + 2*x*y**2

    assert f.as_expr({x: 5}) == 25 - y + 10*y**2
    assert f.as_expr({y: 6}) == -6 + 72*x + x**2

    assert f.as_expr({x: 5, y: 6}) == 379
    assert f.as_expr(5, 6) == 379

    pytest.raises(GeneratorsError, lambda: f.as_expr({z: 7}))


def test_Poly_deflate():
    assert Poly(0, x).deflate() == ((1,), Poly(0, x))
    assert Poly(1, x).deflate() == ((1,), Poly(1, x))
    assert Poly(x, x).deflate() == ((1,), Poly(x, x))

    assert Poly(x**2, x).deflate() == ((2,), Poly(x, x))
    assert Poly(x**17, x).deflate() == ((17,), Poly(x, x))

    assert Poly(
        x**2*y*z**11 + x**4*z**11).deflate() == ((2, 1, 11), Poly(x*y*z + x**2*z))


def test_Poly_inject():
    f = Poly(x**2*y + x*y**3 + x*y + 1, x)

    assert f.inject() == Poly(x**2*y + x*y**3 + x*y + 1, x, y)
    assert f.inject(front=True) == Poly(y**3*x + y*x**2 + y*x + 1, y, x)

    f = Poly(x**2 + 2*x - 1)
    assert f.inject() == f


def test_Poly_eject():
    f = Poly(x**2*y + x*y**3 + x*y + 1, x, y)

    assert f.eject(x) == Poly(x*y**3 + (x**2 + x)*y + 1, y, domain='ZZ[x]')
    assert f.eject(y) == Poly(y*x**2 + (y**3 + y)*x + 1, x, domain='ZZ[y]')

    ex = x + y + z + t + w
    g = Poly(ex, x, y, z, t, w)

    assert g.eject(x) == Poly(ex, y, z, t, w, domain='ZZ[x]')
    assert g.eject(x, y) == Poly(ex, z, t, w, domain='ZZ[x, y]')
    assert g.eject(x, y, z) == Poly(ex, t, w, domain='ZZ[x, y, z]')
    assert g.eject(w) == Poly(ex, x, y, z, t, domain='ZZ[w]')
    assert g.eject(t, w) == Poly(ex, x, y, z, domain='ZZ[w, t]')
    assert g.eject(z, t, w) == Poly(ex, x, y, domain='ZZ[w, t, z]')

    pytest.raises(DomainError, lambda: Poly(x*y, x, y, domain=ZZ.poly_ring(z)).eject(y))

    assert Poly(x*y, x, y, z).eject(y) == Poly(x*y, x, z, domain='ZZ[y]')


def test_Poly_exclude():
    assert Poly(x, x, y).exclude() == Poly(x, x)
    assert Poly(x*y, x, y).exclude() == Poly(x*y, x, y)
    assert Poly(1, x, y).exclude() == Poly(1, x, y)


def test_Poly__gen_to_level():
    assert Poly(1, x, y)._gen_to_level(-2) == 0
    assert Poly(1, x, y)._gen_to_level(-1) == 1
    assert Poly(1, x, y)._gen_to_level(+0) == 0
    assert Poly(1, x, y)._gen_to_level(+1) == 1

    pytest.raises(PolynomialError, lambda: Poly(1, x, y)._gen_to_level(-3))
    pytest.raises(PolynomialError, lambda: Poly(1, x, y)._gen_to_level(+2))

    assert Poly(1, x, y)._gen_to_level(x) == 0
    assert Poly(1, x, y)._gen_to_level(y) == 1

    assert Poly(1, x, y)._gen_to_level('x') == 0
    assert Poly(1, x, y)._gen_to_level('y') == 1

    pytest.raises(PolynomialError, lambda: Poly(1, x, y)._gen_to_level(z))
    pytest.raises(PolynomialError, lambda: Poly(1, x, y)._gen_to_level('z'))


def test_Poly_degree():
    assert Poly(0, x).degree() == -oo
    assert Poly(1, x).degree() == 0
    assert Poly(x, x).degree() == 1

    assert Poly(0, x).degree(gen=0) == -oo
    assert Poly(1, x).degree(gen=0) == 0
    assert Poly(x, x).degree(gen=0) == 1

    assert Poly(0, x).degree(gen=x) == -oo
    assert Poly(1, x).degree(gen=x) == 0
    assert Poly(x, x).degree(gen=x) == 1

    assert Poly(0, x).degree(gen='x') == -oo
    assert Poly(1, x).degree(gen='x') == 0
    assert Poly(x, x).degree(gen='x') == 1

    pytest.raises(PolynomialError, lambda: Poly(1, x).degree(gen=1))
    pytest.raises(PolynomialError, lambda: Poly(1, x).degree(gen=y))
    pytest.raises(PolynomialError, lambda: Poly(1, x).degree(gen='y'))

    assert Poly(1, x, y).degree() == 0
    assert Poly(2*y, x, y).degree() == 0
    assert Poly(x*y, x, y).degree() == 1

    assert Poly(1, x, y).degree(gen=x) == 0
    assert Poly(2*y, x, y).degree(gen=x) == 0
    assert Poly(x*y, x, y).degree(gen=x) == 1

    assert Poly(1, x, y).degree(gen=y) == 0
    assert Poly(2*y, x, y).degree(gen=y) == 1
    assert Poly(x*y, x, y).degree(gen=y) == 1

    assert degree(1, x) == 0
    assert degree(x, x) == 1

    assert degree(x*y**2, gen=x) == 1
    assert degree(x*y**2, gen=y) == 2

    assert degree(x*y**2, x, y) == 1
    assert degree(x*y**2, y, x) == 2

    pytest.raises(ComputationFailed, lambda: degree(1))


def test_Poly_degree_list():
    assert Poly(0, x).degree_list() == (-oo,)
    assert Poly(0, x, y).degree_list() == (-oo, -oo)
    assert Poly(0, x, y, z).degree_list() == (-oo, -oo, -oo)

    assert Poly(1, x).degree_list() == (0,)
    assert Poly(1, x, y).degree_list() == (0, 0)
    assert Poly(1, x, y, z).degree_list() == (0, 0, 0)

    assert Poly(x**2*y + x**3*z**2 + 1).degree_list() == (3, 1, 2)

    assert degree_list(1, x) == (0,)
    assert degree_list(x, x) == (1,)

    assert degree_list(x*y**2) == (1, 2)

    pytest.raises(ComputationFailed, lambda: degree_list(1))


def test_Poly_total_degree():
    assert Poly(x**2*y + x**3*z**2 + 1).total_degree() == 5
    assert Poly(x**2 + z**3).total_degree() == 3
    assert Poly(x*y*z + z**4).total_degree() == 4
    assert Poly(x**3 + x + 1).total_degree() == 3


def test_Poly_LC():
    assert Poly(0, x).LC() == 0
    assert Poly(1, x).LC() == 1
    assert Poly(2*x**2 + x, x).LC() == 2

    assert Poly(x*y**7 + 2*x**2*y**3).LC('lex') == 2
    assert Poly(x*y**7 + 2*x**2*y**3).LC('grlex') == 1

    assert LC(x*y**7 + 2*x**2*y**3, order='lex') == 2
    assert LC(x*y**7 + 2*x**2*y**3, order='grlex') == 1

    pytest.raises(ComputationFailed, lambda: LC([1, 2]))


def test_Poly_TC():
    assert Poly(0, x).TC() == 0
    assert Poly(1, x).TC() == 1
    assert Poly(2*x**2 + x, x).TC() == 0


def test_Poly_EC():
    assert Poly(0, x).EC() == 0
    assert Poly(1, x).EC() == 1
    assert Poly(2*x**2 + x, x).EC() == 1

    assert Poly(x*y**7 + 2*x**2*y**3).EC('lex') == 1
    assert Poly(x*y**7 + 2*x**2*y**3).EC('grlex') == 2


def test_Poly_coeff():
    assert Poly(0, x).coeff_monomial(1) == 0
    assert Poly(0, x).coeff_monomial(x) == 0

    assert Poly(1, x).coeff_monomial(1) == 1
    assert Poly(1, x).coeff_monomial(x) == 0

    assert Poly(x**8, x).coeff_monomial(1) == 0
    assert Poly(x**8, x).coeff_monomial(x**7) == 0
    assert Poly(x**8, x).coeff_monomial(x**8) == 1
    assert Poly(x**8, x).coeff_monomial(x**9) == 0

    assert Poly(3*x*y**2 + 1, x, y).coeff_monomial(1) == 1
    assert Poly(3*x*y**2 + 1, x, y).coeff_monomial(x*y**2) == 3

    p = Poly(24*x*y*exp(8) + 23*x, x, y)

    assert p.coeff_monomial(x) == 23
    assert p.coeff_monomial(y) == 0
    assert p.coeff_monomial(x*y) == 24*exp(8)

    assert p.as_expr().coeff(x) == 24*y*exp(8) + 23
    pytest.raises(NotImplementedError, lambda: p.coeff(x))

    pytest.raises(ValueError, lambda: Poly(x + 1).coeff_monomial(0))
    pytest.raises(ValueError, lambda: Poly(x + 1).coeff_monomial(3*x))
    pytest.raises(ValueError, lambda: Poly(x + 1).coeff_monomial(3*x*y))

    assert Poly(0, x).coeff_monomial((0,)) == 0
    assert Poly(0, x).coeff_monomial((1,)) == 0

    assert Poly(1, x).coeff_monomial((0,)) == 1
    assert Poly(1, x).coeff_monomial((1,)) == 0

    assert Poly(x**8, x).coeff_monomial((0,)) == 0
    assert Poly(x**8, x).coeff_monomial((7,)) == 0
    assert Poly(x**8, x).coeff_monomial((8,)) == 1
    assert Poly(x**8, x).coeff_monomial((9,)) == 0

    assert Poly(3*x*y**2 + 1, x, y).coeff_monomial((0, 0)) == 1
    assert Poly(3*x*y**2 + 1, x, y).coeff_monomial((1, 2)) == 3

    pytest.raises(ValueError, lambda: Poly(x*y + 1, x, y).coeff_monomial((1,)))

    assert Poly(x**3 + 2*x**2 + 3*x, x).coeff_monomial((2,)) == 2
    assert Poly(x**3 + 2*x*y**2 + y**2, x, y).coeff_monomial((1, 2)) == 2
    assert Poly(4*sqrt(x)*y).coeff_monomial((1, 1)) == 4


def test_Poly_LM():
    assert Poly(0, x).LM() == (0,)
    assert Poly(1, x).LM() == (0,)
    assert Poly(2*x**2 + x, x).LM() == (2,)

    assert Poly(x*y**7 + 2*x**2*y**3).LM('lex') == (2, 3)
    assert Poly(x*y**7 + 2*x**2*y**3).LM('grlex') == (1, 7)

    assert LM(x*y**7 + 2*x**2*y**3, order='lex') == x**2*y**3
    assert LM(x*y**7 + 2*x**2*y**3, order='grlex') == x*y**7

    pytest.raises(ComputationFailed, lambda: LM([1, 2]))


def test_Poly_LM_custom_order():
    f = Poly(x**2*y**3*z + x**2*y*z**3 + x*y*z + 1)

    def rev_lex(monom):
        return tuple(reversed(monom))

    assert f.LM(order='lex') == (2, 3, 1)
    assert f.LM(order=rev_lex) == (2, 1, 3)


def test_Poly_EM():
    assert Poly(0, x).EM() == (0,)
    assert Poly(1, x).EM() == (0,)
    assert Poly(2*x**2 + x, x).EM() == (1,)

    assert Poly(x*y**7 + 2*x**2*y**3).EM('lex') == (1, 7)
    assert Poly(x*y**7 + 2*x**2*y**3).EM('grlex') == (2, 3)


def test_Poly_LT():
    assert Poly(0, x).LT() == ((0,), 0)
    assert Poly(1, x).LT() == ((0,), 1)
    assert Poly(2*x**2 + x, x).LT() == ((2,), 2)

    assert Poly(x*y**7 + 2*x**2*y**3).LT('lex') == ((2, 3), 2)
    assert Poly(x*y**7 + 2*x**2*y**3).LT('grlex') == ((1, 7), 1)

    assert LT(x*y**7 + 2*x**2*y**3, order='lex') == 2*x**2*y**3
    assert LT(x*y**7 + 2*x**2*y**3, order='grlex') == x*y**7

    pytest.raises(ComputationFailed, lambda: LT([1, 2]))


def test_Poly_ET():
    assert Poly(0, x).ET() == ((0,), 0)
    assert Poly(1, x).ET() == ((0,), 1)
    assert Poly(2*x**2 + x, x).ET() == ((1,), 1)

    assert Poly(x*y**7 + 2*x**2*y**3).ET('lex') == ((1, 7), 1)
    assert Poly(x*y**7 + 2*x**2*y**3).ET('grlex') == ((2, 3), 2)


def test_Poly_max_norm():
    assert Poly(-1, x).max_norm() == 1
    assert Poly(+0, x).max_norm() == 0
    assert Poly(+1, x).max_norm() == 1


def test_Poly_l1_norm():
    assert Poly(-1, x).l1_norm() == 1
    assert Poly(+0, x).l1_norm() == 0
    assert Poly(+1, x).l1_norm() == 1


def test_Poly_clear_denoms():
    coeff, poly = Poly(x + 2, x).clear_denoms()
    assert coeff == 1 and poly == Poly(
        x + 2, x, domain='ZZ') and poly.domain == ZZ

    coeff, poly = Poly(x/2 + 1, x).clear_denoms()
    assert coeff == 2 and poly == Poly(
        x + 2, x, domain='QQ') and poly.domain == QQ

    coeff, poly = Poly(x/2 + 1, x).clear_denoms(convert=True)
    assert coeff == 2 and poly == Poly(
        x + 2, x, domain='ZZ') and poly.domain == ZZ

    coeff, poly = Poly(x/y + 1, x).clear_denoms(convert=True)
    assert coeff == y and poly == Poly(
        x + y, x, domain='ZZ[y]') and poly.domain == ZZ.poly_ring(y)

    coeff, poly = Poly(x/3 + sqrt(2), x, domain='EX').clear_denoms()
    assert coeff == 3 and poly == Poly(
        x + 3*sqrt(2), x, domain='EX') and poly.domain == EX

    coeff, poly = Poly(
        x/3 + sqrt(2), x, domain='EX').clear_denoms(convert=True)
    assert coeff == 3 and poly == Poly(
        x + 3*sqrt(2), x, domain='EX') and poly.domain == EX


def test_Poly_rat_clear_denoms():
    f = Poly(x**2/y + 1, x)
    g = Poly(x**3 + y, x)

    assert f.rat_clear_denoms(g) == \
        (Poly(x**2 + y, x), Poly(y*x**3 + y**2, x))

    f = f.set_domain(EX)
    g = g.set_domain(EX)

    assert f.rat_clear_denoms(g) == (f, g)


def test_Poly_integrate():
    assert Poly(x + 1).integrate() == Poly(x**2/2 + x)
    assert Poly(x + 1).integrate(x) == Poly(x**2/2 + x)
    assert Poly(x + 1).integrate((x, 1)) == Poly(x**2/2 + x)

    assert Poly(2*x + 1).integrate(auto=False) == Poly(x**2 + x)

    assert Poly(x*y + 1).integrate(x) == Poly(x**2*y/2 + x)
    assert Poly(x*y + 1).integrate(y) == Poly(x*y**2/2 + y)

    assert Poly(x*y + 1).integrate(x, x) == Poly(x**3*y/6 + x**2/2)
    assert Poly(x*y + 1).integrate(y, y) == Poly(x*y**3/6 + y**2/2)

    assert Poly(x*y + 1).integrate((x, 2)) == Poly(x**3*y/6 + x**2/2)
    assert Poly(x*y + 1).integrate((y, 2)) == Poly(x*y**3/6 + y**2/2)

    assert Poly(x*y + 1).integrate(x, y) == Poly(x**2*y**2/4 + x*y)
    assert Poly(x*y + 1).integrate(y, x) == Poly(x**2*y**2/4 + x*y)


def test_Poly_diff():
    assert Poly(x**2 + x).diff() == Poly(2*x + 1)
    assert Poly(x**2 + x).diff(x) == Poly(2*x + 1)
    assert Poly(x**2 + x).diff((x, 1)) == Poly(2*x + 1)

    assert Poly(x**2*y**2 + x*y).diff(x) == Poly(2*x*y**2 + y)
    assert Poly(x**2*y**2 + x*y).diff(y) == Poly(2*x**2*y + x)

    assert Poly(x**2*y**2 + x*y).diff(x, x) == Poly(2*y**2, x, y)
    assert Poly(x**2*y**2 + x*y).diff(y, y) == Poly(2*x**2, x, y)

    assert Poly(x**2*y**2 + x*y).diff((x, 2)) == Poly(2*y**2, x, y)
    assert Poly(x**2*y**2 + x*y).diff((y, 2)) == Poly(2*x**2, x, y)

    assert Poly(x**2*y**2 + x*y).diff(x, y) == Poly(4*x*y + 1)
    assert Poly(x**2*y**2 + x*y).diff(y, x) == Poly(4*x*y + 1)


def test_sympyissue_9585():
    assert diff(Poly(x**2 + x)) == Poly(2*x + 1)
    assert diff(Poly(x**2 + x), x, evaluate=False) == \
        Derivative(Poly(x**2 + x), x)
    assert Derivative(Poly(x**2 + x), x).doit() == Poly(2*x + 1)
    assert Poly(x**2 + x).diff(x, evaluate=False) == \
        Derivative(Poly(x**2 + x), x)


def test_Poly_eval():
    assert Poly(0, x).eval(7) == 0
    assert Poly(1, x).eval(7) == 1
    assert Poly(x, x).eval(7) == 7

    assert Poly(0, x).eval(0, 7) == 0
    assert Poly(1, x).eval(0, 7) == 1
    assert Poly(x, x).eval(0, 7) == 7

    assert Poly(0, x).eval(x, 7) == 0
    assert Poly(1, x).eval(x, 7) == 1
    assert Poly(x, x).eval(x, 7) == 7

    assert Poly(0, x).eval('x', 7) == 0
    assert Poly(1, x).eval('x', 7) == 1
    assert Poly(x, x).eval('x', 7) == 7

    pytest.raises(PolynomialError, lambda: Poly(1, x).eval(1, 7))
    pytest.raises(PolynomialError, lambda: Poly(1, x).eval(y, 7))
    pytest.raises(PolynomialError, lambda: Poly(1, x).eval('y', 7))

    assert Poly(123, x, y).eval(7) == Poly(123, y)
    assert Poly(2*y, x, y).eval(7) == Poly(2*y, y)
    assert Poly(x*y, x, y).eval(7) == Poly(7*y, y)

    assert Poly(123, x, y).eval(x, 7) == Poly(123, y)
    assert Poly(2*y, x, y).eval(x, 7) == Poly(2*y, y)
    assert Poly(x*y, x, y).eval(x, 7) == Poly(7*y, y)

    assert Poly(123, x, y).eval(y, 7) == Poly(123, x)
    assert Poly(2*y, x, y).eval(y, 7) == Poly(14, x)
    assert Poly(x*y, x, y).eval(y, 7) == Poly(7*x, x)

    assert Poly(x*y + y, x, y).eval({x: 7}) == Poly(8*y, y)
    assert Poly(x*y + y, x, y).eval({y: 7}) == Poly(7*x + 7, x)

    assert Poly(x*y + y, x, y).eval({x: 6, y: 7}) == 49
    assert Poly(x*y + y, x, y).eval({x: 7, y: 6}) == 48

    assert Poly(x*y + y, x, y).eval((6, 7)) == 49
    assert Poly(x*y + y, x, y).eval([6, 7]) == 49

    assert Poly(x + 1, domain='ZZ').eval(Rational(1, 2)) == Rational(3, 2)
    assert Poly(x + 1, domain='ZZ').eval(sqrt(2)) == sqrt(2) + 1

    pytest.raises(ValueError, lambda: Poly(x*y + y, x, y).eval((6, 7, 8)))
    pytest.raises(DomainError, lambda: Poly(x + 1, domain='ZZ').eval(Rational(1, 2), auto=False))

    # issue sympy/sympy#6344
    alpha = Symbol('alpha')
    result = (2*alpha*z - 2*alpha + z**2 + 3)/(z**2 - 2*z + 1)

    f = Poly(x**2 + (alpha - 1)*x - alpha + 1, x, domain='ZZ[alpha]')
    assert f.eval((z + 1)/(z - 1)) == result

    g = Poly(x**2 + (alpha - 1)*x - alpha + 1, x, y, domain='ZZ[alpha]')
    assert g.eval((z + 1)/(z - 1)) == Poly(result, y, domain='ZZ(alpha,z)')


def test_Poly___call__():
    f = Poly(2*x*y + 3*x + y + 2*z)

    assert f(2) == Poly(5*y + 2*z + 6)
    assert f(2, 5) == Poly(2*z + 31)
    assert f(2, 5, 7) == 45


def test_parallel_poly_from_expr():
    pytest.raises(PolificationFailed, lambda: parallel_poly_from_expr([[1, 2]]))

    assert parallel_poly_from_expr(
        [x - 1, x**2 - 1], x)[0] == [Poly(x - 1, x), Poly(x**2 - 1, x)]
    assert parallel_poly_from_expr(
        [Poly(x - 1, x), x**2 - 1], x)[0] == [Poly(x - 1, x), Poly(x**2 - 1, x)]
    assert parallel_poly_from_expr(
        [x - 1, Poly(x**2 - 1, x)], x)[0] == [Poly(x - 1, x), Poly(x**2 - 1, x)]
    assert parallel_poly_from_expr([Poly(
        x - 1, x), Poly(x**2 - 1, x)], x)[0] == [Poly(x - 1, x), Poly(x**2 - 1, x)]

    assert parallel_poly_from_expr(
        [x - 1, x**2 - 1], x, y)[0] == [Poly(x - 1, x, y), Poly(x**2 - 1, x, y)]
    assert parallel_poly_from_expr([Poly(
        x - 1, x), x**2 - 1], x, y)[0] == [Poly(x - 1, x, y), Poly(x**2 - 1, x, y)]
    assert parallel_poly_from_expr([x - 1, Poly(
        x**2 - 1, x)], x, y)[0] == [Poly(x - 1, x, y), Poly(x**2 - 1, x, y)]
    assert parallel_poly_from_expr([Poly(x - 1, x), Poly(
        x**2 - 1, x)], x, y)[0] == [Poly(x - 1, x, y), Poly(x**2 - 1, x, y)]

    assert parallel_poly_from_expr(
        [x - 1, x**2 - 1])[0] == [Poly(x - 1, x), Poly(x**2 - 1, x)]
    assert parallel_poly_from_expr(
        [Poly(x - 1, x), x**2 - 1])[0] == [Poly(x - 1, x), Poly(x**2 - 1, x)]
    assert parallel_poly_from_expr(
        [x - 1, Poly(x**2 - 1, x)])[0] == [Poly(x - 1, x), Poly(x**2 - 1, x)]
    assert parallel_poly_from_expr(
        [Poly(x - 1, x), Poly(x**2 - 1, x)])[0] == [Poly(x - 1, x), Poly(x**2 - 1, x)]

    assert parallel_poly_from_expr(
        [1, x**2 - 1])[0] == [Poly(1, x), Poly(x**2 - 1, x)]
    assert parallel_poly_from_expr(
        [1, x**2 - 1])[0] == [Poly(1, x), Poly(x**2 - 1, x)]
    assert parallel_poly_from_expr(
        [1, Poly(x**2 - 1, x)])[0] == [Poly(1, x), Poly(x**2 - 1, x)]
    assert parallel_poly_from_expr(
        [1, Poly(x**2 - 1, x)])[0] == [Poly(1, x), Poly(x**2 - 1, x)]

    assert parallel_poly_from_expr(
        [x**2 - 1, 1])[0] == [Poly(x**2 - 1, x), Poly(1, x)]
    assert parallel_poly_from_expr(
        [x**2 - 1, 1])[0] == [Poly(x**2 - 1, x), Poly(1, x)]
    assert parallel_poly_from_expr(
        [Poly(x**2 - 1, x), 1])[0] == [Poly(x**2 - 1, x), Poly(1, x)]
    assert parallel_poly_from_expr(
        [Poly(x**2 - 1, x), 1])[0] == [Poly(x**2 - 1, x), Poly(1, x)]

    assert parallel_poly_from_expr([Poly(x, x, y), Poly(y, x, y)], x, y, order='lex')[0] == \
        [Poly(x, x, y, domain='ZZ'), Poly(y, x, y, domain='ZZ')]

    pytest.raises(PolificationFailed, lambda: parallel_poly_from_expr([0, 1]))

    assert (parallel_poly_from_expr([(x - 1)**2, 1], expand=False) ==
            ([Poly((x - 1)**2, x - 1, expand=False), Poly(1, x - 1)],
             {'domain': ZZ, 'expand': False, 'gens': (x - 1,),
              'polys': False}))


def test_prem():
    f, g = x**2 - y**2, x - y
    q, r = x + y, 0

    F, G, Q, R = [Poly(h, x, y) for h in (f, g, q, r)]

    assert F.prem(G) == R

    assert prem(f, g) == r
    assert prem(f, g, x, y) == r
    assert prem(f, g, (x, y)) == r
    assert prem(F, G) == R
    assert prem(f, g, polys=True) == R
    assert prem(F, G, polys=False) == r

    pytest.raises(ComputationFailed, lambda: prem(4, 2))


def test_div():
    f, g = x**2 - y**2, x - y
    q, r = x + y, 0

    F, G, Q, R = [Poly(h, x, y) for h in (f, g, q, r)]

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

    assert div(f, g, (x, y)) == (q, r)
    assert rem(f, g, (x, y)) == r
    assert quo(f, g, (x, y)) == q
    assert exquo(f, g, (x, y)) == q

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

    pytest.raises(ComputationFailed, lambda: div(4, 2))
    pytest.raises(ComputationFailed, lambda: rem(4, 2))
    pytest.raises(ComputationFailed, lambda: quo(4, 2))
    pytest.raises(ComputationFailed, lambda: exquo(4, 2))

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
    pytest.raises(ExactQuotientFailed, lambda: exquo(f, g, auto=False))
    pytest.raises(ExactQuotientFailed, lambda: exquo(f, g, domain=ZZ))
    assert exquo(f, g, domain=QQ) == q
    assert exquo(f, g, domain=ZZ, auto=True) == q
    pytest.raises(ExactQuotientFailed, lambda: exquo(f, g, domain=ZZ, auto=False))
    assert exquo(f, g, domain=QQ, auto=True) == q
    assert exquo(f, g, domain=QQ, auto=False) == q

    f, g = Poly(x**2), Poly(x)

    q, r = f.div(g)
    assert q.domain.is_IntegerRing and r.domain.is_IntegerRing
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

    f, g = Poly(x + y, x), Poly(2*x + y, x)
    q, r = f.div(g)
    assert q.domain.is_FractionField and r.domain.is_FractionField


def test_gcdex():
    f, g = 2*x, x**2 - 16
    s, t, h = x/32, -Rational(1, 16), 1

    F, G, S, T, H = [Poly(u, x, domain='QQ') for u in (f, g, s, t, h)]

    assert F.half_gcdex(G) == (S, H)
    assert F.gcdex(G) == (S, T, H)
    assert F.invert(G) == S

    assert half_gcdex(f, g) == (s, h)
    assert gcdex(f, g) == (s, t, h)
    assert invert(f, g) == s

    assert half_gcdex(f, g, x) == (s, h)
    assert gcdex(f, g, x) == (s, t, h)
    assert invert(f, g, x) == s

    assert half_gcdex(f, g, (x,)) == (s, h)
    assert gcdex(f, g, (x,)) == (s, t, h)
    assert invert(f, g, (x,)) == s

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

    R, y, z = ring('y,z', QQ)
    pytest.raises(MultivariatePolynomialError, lambda: (y + z).half_gcdex(y - z))
    pytest.raises(MultivariatePolynomialError, lambda: (y + z).gcdex(y - z))


def test_subresultants():
    f, g, h = x**2 - 2*x + 1, x**2 - 1, 2*x - 2
    F, G, H = Poly(f), Poly(g), Poly(h)

    assert F.subresultants(G) == [F, G, H]
    assert subresultants(f, g) == [f, g, h]
    assert subresultants(f, g, x) == [f, g, h]
    assert subresultants(f, g, (x,)) == [f, g, h]
    assert subresultants(F, G) == [F, G, H]
    assert subresultants(f, g, polys=True) == [F, G, H]
    assert subresultants(F, G, polys=False) == [f, g, h]

    pytest.raises(ComputationFailed, lambda: subresultants(4, 2))


def test_resultant():
    f, g, h = x**2 - 2*x + 1, x**2 - 1, 0
    F, G = Poly(f), Poly(g)

    assert F.resultant(G) == h
    assert PurePoly(f).resultant(PurePoly(g)) == h
    assert resultant(f, g) == h
    assert resultant(f, g, x) == h
    assert resultant(f, g, (x,)) == h
    assert resultant(F, G) == h
    assert resultant(f, g, polys=True) == h
    assert resultant(F, G, polys=False) == h
    assert resultant(f, g, includePRS=True) == (h, [f, g, 2*x - 2])
    assert resultant(f, g, polys=True,
                     includePRS=True) == (0, [f, g, Poly(2*x - 2)])

    f, g, h = x - a, x - b, a - b
    F, G, H = Poly(f), Poly(g), Poly(h)

    assert F.resultant(G) == H
    assert resultant(f, g) == h
    assert resultant(f, g, x) == h
    assert resultant(f, g, (x,)) == h
    assert resultant(F, G) == H
    assert resultant(f, g, polys=True) == H
    assert resultant(F, G, polys=False) == h

    pytest.raises(ComputationFailed, lambda: resultant(4, 2))

    assert resultant(PurePoly(x**2 + y),
                     PurePoly(x*y - 1)) == PurePoly(y**3 + 1)


def test_discriminant():
    f, g = x**3 + 3*x**2 + 9*x - 13, -11664
    F = Poly(f)

    assert F.discriminant() == g
    assert discriminant(f) == g
    assert discriminant(f, x) == g
    assert discriminant(f, (x,)) == g
    assert discriminant(F) == g
    assert discriminant(f, polys=True) == g
    assert discriminant(F, polys=False) == g

    f, g = a*x**2 + b*x + c, b**2 - 4*a*c
    F, G = Poly(f), Poly(g)

    assert F.discriminant() == G
    assert discriminant(f) == g
    assert discriminant(f, x, a, b, c) == g
    assert discriminant(f, (x, a, b, c)) == g
    assert discriminant(F) == G
    assert discriminant(f, polys=True) == G
    assert discriminant(F, polys=False) == g

    pytest.raises(ComputationFailed, lambda: discriminant(4))


def test_dispersion():
    # We test only the API here. For more mathematical
    # tests see the dedicated test file.
    fp = poly((x + 1)*(x + 2), x)
    assert sorted(fp.dispersionset()) == [0, 1]
    assert fp.dispersion() == 1

    fp = poly(x**4 - 3*x**2 + 1, x)
    gp = fp.shift(-3)
    assert sorted(fp.dispersionset(gp)) == [2, 3, 4]
    assert fp.dispersion(gp) == 4


def test_gcd_list():
    F = [x**3 - 1, x**2 - 1, x**2 - 3*x + 2]

    assert gcd_list(F) == x - 1
    assert gcd_list(F, polys=True) == Poly(x - 1)

    assert gcd_list([]) == 0
    assert gcd_list([1, 2]) == 1
    assert gcd_list([4, 6, 8]) == 2

    assert gcd_list([x*(y + 42) - x*y - x*42]) == 0

    gcd = gcd_list([], x)
    assert gcd.is_Number and gcd is Integer(0)

    gcd = gcd_list([], x, polys=True)
    assert gcd.is_Poly and gcd.is_zero

    pytest.raises(ComputationFailed, lambda: gcd_list([], polys=True))

    assert gcd_list([x**3 - 1, x**2 - 2, x**2 - 3*x + 2]) == 1


def test_lcm_list():
    F = [x**3 - 1, x**2 - 1, x**2 - 3*x + 2]

    assert lcm_list(F) == x**5 - x**4 - 2*x**3 - x**2 + x + 2
    assert lcm_list(F, polys=True) == Poly(x**5 - x**4 - 2*x**3 - x**2 + x + 2)

    assert lcm_list([]) == 1
    assert lcm_list([1, 2]) == 2
    assert lcm_list([4, 6, 8]) == 24

    assert lcm_list([x*(y + 42) - x*y - x*42]) == 0

    lcm = lcm_list([], x)
    assert lcm.is_Number and lcm is Integer(1)

    lcm = lcm_list([], x, polys=True)
    assert lcm.is_Poly and lcm.is_one

    pytest.raises(ComputationFailed, lambda: lcm_list([], polys=True))


def test_gcd():
    f, g = x**3 - 1, x**2 - 1
    s, t = x**2 + x + 1, x + 1
    h, r = x - 1, x**4 + x**3 - x - 1

    F, G, S, T, H, R = [Poly(u) for u in (f, g, s, t, h, r)]

    assert F.cofactors(G) == (H, S, T)
    assert F.gcd(G) == H
    assert F.lcm(G) == R

    assert cofactors(f, g) == (h, s, t)
    assert gcd(f, g) == h
    assert lcm(f, g) == r

    assert cofactors(f, g, x) == (h, s, t)
    assert gcd(f, g, x) == h
    assert lcm(f, g, x) == r

    assert cofactors(f, g, (x,)) == (h, s, t)
    assert gcd(f, g, (x,)) == h
    assert lcm(f, g, (x,)) == r

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
    F, G = Poly(f), Poly(g)

    assert F.trunc(3) == G
    assert trunc(f, 3) == g
    assert trunc(f, 3, x) == g
    assert trunc(f, 3, (x,)) == g
    assert trunc(F, 3) == G
    assert trunc(f, 3, polys=True) == G
    assert trunc(F, 3, polys=False) == g

    f, g = 6*x**5 + 5*x**4 + 4*x**3 + 3*x**2 + 2*x + 1, -x**4 + x**3 - x + 1
    F, G = Poly(f), Poly(g)

    assert F.trunc(3) == G
    assert trunc(f, 3) == g
    assert trunc(f, 3, x) == g
    assert trunc(f, 3, (x,)) == g
    assert trunc(F, 3) == G
    assert trunc(f, 3, polys=True) == G
    assert trunc(F, 3, polys=False) == g

    f = Poly(x**2 + 2*x + 3, modulus=5)

    assert f.trunc(2) == Poly(x**2 + 1, modulus=5)

    pytest.raises(ComputationFailed, lambda: trunc([1, 2], 2))


def test_monic():
    f, g = 2*x - 1, x - Rational(1, 2)
    F, G = Poly(f, domain='QQ'), Poly(g)

    assert F.monic() == G
    assert monic(f) == g
    assert monic(f, x) == g
    assert monic(f, (x,)) == g
    assert monic(F) == G
    assert monic(f, polys=True) == G
    assert monic(F, polys=False) == g

    pytest.raises(ComputationFailed, lambda: monic(4))

    assert monic(2*x**2 + 6*x + 4, auto=False) == x**2 + 3*x + 2
    pytest.raises(ExactQuotientFailed, lambda: monic(2*x + 6*x + 1, auto=False))

    assert monic(2.0*x**2 + 6.0*x + 4.0) == 1.0*x**2 + 3.0*x + 2.0
    assert monic(2*x**2 + 3*x + 4, modulus=5) == x**2 + 4*x + 2


def test_content():
    f, F = 4*x + 2, Poly(4*x + 2)

    assert F.content() == 2
    assert content(f) == 2

    pytest.raises(ComputationFailed, lambda: content(4))

    f = Poly(2*x, modulus=3)

    assert f.content() == 1


def test_primitive():
    f, g = 4*x + 2, 2*x + 1
    F, G = Poly(f), Poly(g)

    assert F.primitive() == (2, G)
    assert primitive(f) == (2, g)
    assert primitive(f, x) == (2, g)
    assert primitive(f, (x,)) == (2, g)
    assert primitive(F) == (2, G)
    assert primitive(f, polys=True) == (2, G)
    assert primitive(F, polys=False) == (2, g)

    pytest.raises(ComputationFailed, lambda: primitive(4))

    f = Poly(2*x, modulus=3)
    g = Poly(2.0*x, domain=RR)

    assert f.primitive() == (1, f)
    assert g.primitive() == (1.0, g)

    assert primitive(-3*x/4 + y + Rational(11, 8)) == \
        (Rational(-1, 8), 6*x - 8*y - 11)


def test_compose():
    f = x**12 + 20*x**10 + 150*x**8 + 500*x**6 + 625*x**4 - 2*x**3 - 10*x + 9
    g = x**4 - 2*x + 9
    h = x**3 + 5*x

    F, G, H = map(Poly, (f, g, h))

    assert G.compose(H) == F
    assert compose(g, h) == f
    assert compose(g, h, x) == f
    assert compose(g, h, (x,)) == f
    assert compose(G, H) == F
    assert compose(g, h, polys=True) == F
    assert compose(G, H, polys=False) == f

    assert F.decompose() == [G, H]
    assert decompose(f) == [g, h]
    assert decompose(f, x) == [g, h]
    assert decompose(f, (x,)) == [g, h]
    assert decompose(F) == [G, H]
    assert decompose(f, polys=True) == [G, H]
    assert decompose(F, polys=False) == [g, h]

    pytest.raises(ComputationFailed, lambda: compose(4, 2))
    pytest.raises(ComputationFailed, lambda: decompose(4))

    assert compose(x**2 - y**2, x - y, x, y) == x**2 - 2*x*y
    assert compose(x**2 - y**2, x - y, y, x) == -y**2 + 2*x*y


def test_shift():
    assert Poly(x**2 - 2*x + 1, x).shift(2) == Poly(x**2 + 2*x + 1, x)


def test_sturm():
    f, F = x, Poly(x, domain='QQ')
    g, G = 1, Poly(1, x, domain='QQ')

    assert F.sturm() == [F, G]
    assert sturm(f) == [f, g]
    assert sturm(f, x) == [f, g]
    assert sturm(f, (x,)) == [f, g]
    assert sturm(F) == [F, G]
    assert sturm(f, polys=True) == [F, G]
    assert sturm(F, polys=False) == [f, g]

    pytest.raises(ComputationFailed, lambda: sturm(4))
    pytest.raises(DomainError, lambda: sturm(f, auto=False))

    f = Poly(1024/(15625*pi**8)*x**5
             - 4096/(625*pi**8)*x**4
             + 32/(15625*pi**4)*x**3
             - 128/(625*pi**4)*x**2
             + Rational(1, 62500)*x
             - Rational(1, 625), x, domain='ZZ(pi)')

    assert sturm(f) == \
        [Poly(x**3 - 100*x**2 + pi**4/64*x - 25*pi**4/16, x, domain='ZZ(pi)'),
         Poly(3*x**2 - 200*x + pi**4/64, x, domain='ZZ(pi)'),
         Poly((Rational(20000, 9) - pi**4/96)*x + 25*pi**4/18, x, domain='ZZ(pi)'),
         Poly((-3686400000000*pi**4 - 11520000*pi**8 - 9*pi**12)/(26214400000000 - 245760000*pi**4 + 576*pi**8), x, domain='ZZ(pi)')]


def test_sqf_norm():
    assert sqf_norm(x**2 - 2, extension=sqrt(3)) == \
        (1, x**2 - 2*sqrt(3)*x + 1, x**4 - 10*x**2 + 1)
    assert sqf_norm(x**2 - 3, extension=sqrt(2)) == \
        (1, x**2 - 2*sqrt(2)*x - 1, x**4 - 10*x**2 + 1)
    assert sqf_norm(x**2 - 3, extension=sqrt(2), polys=True) == \
        (1, Poly(x**2 - 2*sqrt(2)*x - 1, extension=sqrt(2)),
         Poly(x**4 - 10*x**2 + 1))

    pytest.raises(ComputationFailed, lambda: sqf_norm([1, 2]))

    assert Poly(x**2 - 2, extension=sqrt(3)).sqf_norm() == \
        (1, Poly(x**2 - 2*sqrt(3)*x + 1, x, extension=sqrt(3)),
            Poly(x**4 - 10*x**2 + 1, x, domain='QQ'))

    assert Poly(x**2 - 3, extension=sqrt(2)).sqf_norm() == \
        (1, Poly(x**2 - 2*sqrt(2)*x - 1, x, extension=sqrt(2)),
            Poly(x**4 - 10*x**2 + 1, x, domain='QQ'))


def test_sqf():
    f = x**5 - x**3 - x**2 + 1
    g = x**3 + 2*x**2 + 2*x + 1
    h = x - 1

    p = x**4 + x**3 - x - 1

    F, G, H, P = map(Poly, (f, g, h, p))

    assert F.sqf_part() == P
    assert sqf_part(f) == p
    assert sqf_part(f, x) == p
    assert sqf_part(f, (x,)) == p
    assert sqf_part(F) == P
    assert sqf_part(f, polys=True) == P
    assert sqf_part(F, polys=False) == p

    assert F.sqf_list() == (1, [(G, 1), (H, 2)])
    assert sqf_list(f) == (1, [(g, 1), (h, 2)])
    assert sqf_list(f, x) == (1, [(g, 1), (h, 2)])
    assert sqf_list(f, (x,)) == (1, [(g, 1), (h, 2)])
    assert sqf_list(F) == (1, [(G, 1), (H, 2)])
    assert sqf_list(f, polys=True) == (1, [(G, 1), (H, 2)])
    assert sqf_list(F, polys=False) == (1, [(g, 1), (h, 2)])

    pytest.raises(PolynomialError, lambda: sqf_list([1, 2]))
    pytest.raises(ComputationFailed, lambda: sqf_part(4))

    assert sqf(1) == 1
    assert sqf_list(1) == (1, [])

    assert sqf((2*x**2 + 2)**7) == 128*(x**2 + 1)**7

    assert sqf(f) == g*h**2
    assert sqf(f, x) == g*h**2
    assert sqf(f, (x,)) == g*h**2

    d = x**2 + y**2

    assert sqf(f/d) == (g*h**2)/d
    assert sqf(f/d, x) == (g*h**2)/d
    assert sqf(f/d, (x,)) == (g*h**2)/d

    assert sqf(x - 1) == x - 1
    assert sqf(-x - 1) == -x - 1

    assert sqf(x - 1) == x - 1
    assert sqf(6*x - 10) == Mul(2, 3*x - 5, evaluate=False)

    assert sqf((6*x - 10)/(3*x - 6)) == Rational(2, 3)*((3*x - 5)/(x - 2))
    assert sqf(Poly(x**2 - 2*x + 1)) == (x - 1)**2

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
    assert factor_list(f, (x,)) == (1, [(u, 1), (v, 2), (w, 1)])
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
    assert factor_list(sqrt(2)*sin(x), sin(x)) == (1, [(2, 1/2), (sin(x), 1)])

    assert factor(6) == 6 and factor(6).is_Integer

    assert factor_list(3*x) == (3, [(x, 1)])
    assert factor_list(3*x**2) == (3, [(x, 2)])

    assert factor(3*x) == 3*x
    assert factor(3*x**2) == 3*x**2

    assert factor((2*x**2 + 2)**7) == 128*(x**2 + 1)**7

    assert factor(f) == u*v**2*w
    assert factor(f, x) == u*v**2*w
    assert factor(f, (x,)) == u*v**2*w

    g, p, q, r = x**2 - y**2, x - y, x + y, x**2 + 1

    assert factor(f/g) == (u*v**2*w)/(p*q)
    assert factor(f/g, x) == (u*v**2*w)/(p*q)
    assert factor(f/g, (x,)) == (u*v**2*w)/(p*q)

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

    f = Poly(sin(1)*x + 1, x, domain=EX)

    assert f.factor_list() == (1, [(f, 1)])

    f = x**4 + 1

    assert factor(f) == f
    assert factor(f, extension=I) == (x**2 - I)*(x**2 + I)
    assert factor(f, gaussian=True) == (x**2 - I)*(x**2 + I)
    assert factor(
        f, extension=sqrt(2)) == (x**2 + sqrt(2)*x + 1)*(x**2 - sqrt(2)*x + 1)

    f = x**2 + 2*sqrt(2)*x + 2

    assert factor(f, extension=sqrt(2)) == (x + sqrt(2))**2
    assert factor(f**3, extension=sqrt(2)) == (x + sqrt(2))**6

    assert factor(x**2 - 2*y**2, extension=sqrt(2)) == \
        (x + sqrt(2)*y)*(x - sqrt(2)*y)
    assert factor(2*x**2 - 4*y**2, extension=sqrt(2)) == \
        2*((x + sqrt(2)*y)*(x - sqrt(2)*y))

    assert factor(x - 1) == x - 1
    assert factor(-x - 1) == -x - 1

    assert factor(x - 1) == x - 1

    assert factor(6*x - 10) == Mul(2, 3*x - 5, evaluate=False)

    assert factor(x**11 + x + 1, modulus=65537) == \
        (x**2 + x + 1)*(x**9 + 65536*x**8 + x**6 + 65536*x**5 +
                        x**3 + 65536*x**2 + 1)

    f = x/pi + x*sin(x)/pi
    g = y/(pi**2 + 2*pi + 1) + y*sin(x)/(pi**2 + 2*pi + 1)

    assert factor(f) == x*(sin(x) + 1)/pi
    assert factor(g) == y*(sin(x) + 1)/(pi + 1)**2

    assert factor(Eq(
        x**2 + 2*x + 1, x**3 + 1)) == Eq((x + 1)**2, (x + 1)*(x**2 - x + 1))

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

    assert factor([x, Eq(x**2 - y**2, Tuple(x**2 - z**2, 1/x + 1/y))]) == \
        [x, Eq((x - y)*(x + y), Tuple((x - z)*(x + z), (x + y)/x/y))]

    assert not isinstance(
        Poly(x**3 + x + 1).factor_list()[1][0][0], PurePoly) is True
    assert isinstance(
        PurePoly(x**3 + x + 1).factor_list()[1][0][0], PurePoly) is True

    assert factor(sqrt(-x)) == sqrt(-x)

    # issue sympy/sympy#5917
    e = (-2*x*(-x + 1)*(x - 1)*(-x*(-x + 1)*(x - 1) - x*(x - 1)**2)*(x**2*(x -
                                                                           1) - x*(x - 1) - x) - (-2*x**2*(x - 1)**2 - x*(-x + 1)*(-x*(-x + 1) +
                                                                                                                                   x*(x - 1)))*(x**2*(x - 1)**4 - x*(-x*(-x + 1)*(x - 1) - x*(x - 1)**2)))
    assert factor(e) == 0

    # deep option
    assert factor(sin(x**2 + x) + x, deep=True) == sin(x*(x + 1)) + x

    assert factor(sqrt(x**2)) == sqrt(x**2)

    # issue sympy/sympy#7902
    assert (2*Sum(3*x, (x, 1, 9))).factor() == 6*Sum(x, (x, 1, 9))
    assert (2*Sum(x**2, (x, 1, 9))).factor() == 2*Sum(x**2, (x, 1, 9))

    A, B = symbols('A B', commutative=False)
    f = (x - A)*(y - B)
    assert factor(f.expand()) == f

    assert factor(Sum(4*x, (x, 1, y))) == 4*Sum(x, (x, 1, y))

    # issue sympy/sympy#13149
    assert (factor(expand((0.5*x + 1)*(0.5*y + 1))) ==
            Mul(4.0, 0.25*x + 0.5, 0.25*y + 0.5))
    assert factor(expand((0.5*x + 1)**2)) == 4.0*(0.25*x + 0.5)**2

    assert factor(x**4/2 + 5*x**3/12 - x**2/3) == x**2*(2*x - 1)*(3*x + 4)/12


def test_factor_large():
    f = (x**2 + 4*x + 4)**10000000*(x**2 + 1)*(x**2 + 2*x + 1)**1234567
    g = ((x**2 + 2*x + 1)**3000*y**2 + (x**2 + 2*x + 1)**3000*2*y +
         (x**2 + 2*x + 1)**3000)

    assert factor(f) == (x + 2)**20000000*(x**2 + 1)*(x + 1)**2469134
    assert factor(g) == (x + 1)**6000*(y + 1)**2

    assert factor_list(f) == (1, [(x + 1, 2469134), (x + 2, 20000000),
                                  (x**2 + 1, 1)])
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


def test_intervals():
    assert intervals(0) == []
    assert intervals(1) == []

    assert intervals(x, sqf=True) == [(0, 0)]
    assert intervals(x) == [((0, 0), 1)]

    assert intervals(x**128) == [((0, 0), 128)]
    assert intervals([x**2, x**4]) == [((0, 0), {0: 2, 1: 4})]
    assert intervals([x**2, x**4], eps=0.1) == [((0, 0), {0: 2, 1: 4})]

    pytest.raises(MultivariatePolynomialError,
                  lambda: intervals([x**2 - y, x*y + 1]))
    pytest.raises(MultivariatePolynomialError, lambda: Poly(x**2 - y).intervals())
    pytest.raises(ValueError, lambda: intervals([x**2, x**4], eps=-0.1))

    f = Poly((2*x/5 - Rational(17, 3))*(4*x + Rational(1, 257)))

    assert f.intervals(sqf=True) == [(-1, 0), (14, 15)]
    assert f.intervals() == [((-1, 0), 1), ((14, 15), 1)]

    assert f.intervals(sqf=True) == [(-1, 0), (14, 15)]
    assert f.intervals() == [((-1, 0), 1), ((14, 15), 1)]

    assert f.intervals(eps=Rational(1, 10)) == f.intervals(eps=0.1) == \
        [((-Rational(1, 513), 0), 1), ((Rational(85, 6), Rational(85, 6)), 1)]
    assert f.intervals(eps=Rational(1, 100)) == f.intervals(eps=0.01) == \
        [((-Rational(1, 513), 0), 1), ((Rational(85, 6), Rational(85, 6)), 1)]
    assert f.intervals(eps=Rational(1, 1000)) == f.intervals(eps=0.001) == \
        [((-Rational(1, 1025), 0), 1), ((Rational(85, 6), Rational(85, 6)), 1)]
    assert f.intervals(eps=Rational(1, 10000)) == f.intervals(eps=0.0001) == \
        [((-Rational(1, 1025), -Rational(65, 66881)), 1), ((Rational(85, 6), Rational(85, 6)), 1)]

    f = (2*x/5 - Rational(17, 3))*(4*x + Rational(1, 257))

    assert intervals(f, sqf=True) == [(-1, 0), (14, 15)]
    assert intervals(f) == [((-1, 0), 1), ((14, 15), 1)]

    assert intervals(f, eps=Rational(1, 10)) == intervals(f, eps=0.1) == \
        [((-Rational(1, 513), 0), 1), ((Rational(85, 6), Rational(85, 6)), 1)]
    assert intervals(f, eps=Rational(1, 100)) == intervals(f, eps=0.01) == \
        [((-Rational(1, 513), 0), 1), ((Rational(85, 6), Rational(85, 6)), 1)]
    assert intervals(f, eps=Rational(1, 1000)) == intervals(f, eps=0.001) == \
        [((-Rational(1, 1025), 0), 1), ((Rational(85, 6), Rational(85, 6)), 1)]
    assert intervals(f, eps=Rational(1, 10000)) == intervals(f, eps=0.0001) == \
        [((-Rational(1, 1025), -Rational(65, 66881)), 1), ((Rational(85, 6), Rational(85, 6)), 1)]

    f = Poly((x**2 - 2)*(x**2 - 3)**7*(x + 1)*(7*x + 3)**3)

    assert f.intervals() == \
        [((-2, -Rational(3, 2)), 7), ((-Rational(3, 2), -1), 1),
         ((-1, -1), 1), ((-1, 0), 3),
         ((1, Rational(3, 2)), 1), ((Rational(3, 2), 2), 7)]

    assert intervals([x**5 - 200, x**5 - 201]) == \
        [((Rational(75, 26), Rational(101, 35)), {0: 1}), ((Rational(309, 107), Rational(26, 9)), {1: 1})]

    assert intervals([x**5 - 200, x**5 - 201]) == \
        [((Rational(75, 26), Rational(101, 35)), {0: 1}), ((Rational(309, 107), Rational(26, 9)), {1: 1})]

    assert intervals([x**2 - 200, x**2 - 201]) == \
        [((-Rational(71, 5), -Rational(85, 6)), {1: 1}), ((-Rational(85, 6), -14), {0: 1}),
         ((14, Rational(85, 6)), {0: 1}), ((Rational(85, 6), Rational(71, 5)), {1: 1})]

    f, g, h = x**2 - 2, x**4 - 4*x**2 + 4, x - 1

    assert intervals(f, inf=Rational(7, 4), sqf=True) == []
    assert intervals(f, inf=Rational(7, 5), sqf=True) == [(Rational(7, 5), Rational(3, 2))]
    assert intervals(f, sup=Rational(7, 4), sqf=True) == [(-2, -1), (1, Rational(3, 2))]
    assert intervals(f, sup=Rational(7, 5), sqf=True) == [(-2, -1)]

    assert intervals(g, inf=Rational(7, 4)) == []
    assert intervals(g, inf=Rational(7, 5)) == [((Rational(7, 5), Rational(3, 2)), 2)]
    assert intervals(g, sup=Rational(7, 4)) == [((-2, -1), 2), ((1, Rational(3, 2)), 2)]
    assert intervals(g, sup=Rational(7, 5)) == [((-2, -1), 2)]

    assert intervals([g, h], inf=Rational(7, 4)) == []
    assert intervals([g, h], inf=Rational(7, 5)) == [((Rational(7, 5), Rational(3, 2)), {0: 2})]
    assert intervals([g, h], sup=Rational(7, 4)) == [((-2, -1), {0: 2}), ((1, 1), {1: 1}), ((1, Rational(3, 2)), {0: 2})]
    assert intervals(
        [g, h], sup=Rational(7, 5)) == [((-2, -1), {0: 2}), ((1, 1), {1: 1})]

    assert intervals([x + 2, x**2 - 2]) == \
        [((-2, -2), {0: 1}), ((-2, -1), {1: 1}), ((1, 2), {1: 1})]
    assert intervals([x + 2, x**2 - 2], strict=True) == \
        [((-2, -2), {0: 1}), ((-Rational(3, 2), -1), {1: 1}), ((1, 2), {1: 1})]

    f = 7*z**4 - 19*z**3 + 20*z**2 + 17*z + 20

    assert intervals(f) == []

    real_part, complex_part = intervals(f, all=True, sqf=True)

    assert real_part == []
    assert all(re(a) < re(r) < re(b) and im(
        a) < im(r) < im(b) for (a, b), r in zip(complex_part, nroots(f)))

    assert complex_part == [(-Rational(40, 7) - 40*I/7, 0), (-Rational(40, 7), 40*I/7),
                            (-40*I/7, Rational(40, 7)), (0, Rational(40, 7) + 40*I/7)]

    real_part, complex_part = intervals(f, all=True, sqf=True, eps=Rational(1, 10))

    assert real_part == []
    assert all(re(a) < re(r) < re(b) and im(
        a) < im(r) < im(b) for (a, b), r in zip(complex_part, nroots(f)))

    pytest.raises(ValueError, lambda: intervals(x**2 - 2, eps=10**-100000))
    pytest.raises(ValueError, lambda: Poly(x**2 - 2).intervals(eps=10**-100000))
    pytest.raises(
        ValueError, lambda: intervals([x**2 - 2, x**2 - 3], eps=10**-100000))


@pytest.mark.xfail
def test_intervals_xfail():
    assert intervals([x + 1, x + 2, x - 1, x + 1, 1, x - 1, x - 1, (x - 2)**2]) == \
        [((-2, -2), {1: 1}), ((-1, -1), {0: 1, 3: 1}),
         ((1, 1), {2: 1, 5: 1, 6: 1}), ((2, 2), {7: 2})]


def test_refine_root():
    f = Poly(x**2 - 2)

    assert f.refine_root(1, 2, steps=0) == (1, 2)
    assert f.refine_root(-2, -1, steps=0) == (-2, -1)

    assert f.refine_root(1, 2, steps=None) == (1, Rational(3, 2))
    assert f.refine_root(-2, -1, steps=None) == (-Rational(3, 2), -1)

    assert f.refine_root(1, 2, steps=1) == (1, Rational(3, 2))
    assert f.refine_root(-2, -1, steps=1) == (-Rational(3, 2), -1)

    assert f.refine_root(1, 2, steps=1) == (1, Rational(3, 2))
    assert f.refine_root(-2, -1, steps=1) == (-Rational(3, 2), -1)

    assert f.refine_root(1, 2, eps=Rational(1, 100)) == (Rational(24, 17), Rational(17, 12))
    assert f.refine_root(1, 2, eps=1e-2) == (Rational(24, 17), Rational(17, 12))

    pytest.raises(PolynomialError, lambda: (f**2).refine_root(1, 2, check_sqf=True))

    pytest.raises(RefinementFailed, lambda: (f**2).refine_root(1, 2))
    pytest.raises(RefinementFailed, lambda: (f**2).refine_root(2, 3))

    f = x**2 - 2

    assert refine_root(f, 1, 2, steps=1) == (1, Rational(3, 2))
    assert refine_root(f, -2, -1, steps=1) == (-Rational(3, 2), -1)

    assert refine_root(f, 1, 2, steps=1) == (1, Rational(3, 2))
    assert refine_root(f, -2, -1, steps=1) == (-Rational(3, 2), -1)

    assert refine_root(f, 1, 2, eps=Rational(1, 100)) == (Rational(24, 17), Rational(17, 12))
    assert refine_root(f, 1, 2, eps=1e-2) == (Rational(24, 17), Rational(17, 12))

    pytest.raises(PolynomialError, lambda: refine_root(1, 7, 8, eps=Rational(1, 100)))

    pytest.raises(ValueError, lambda: Poly(f).refine_root(1, 2, eps=10**-100000))
    pytest.raises(ValueError, lambda: refine_root(f, 1, 2, eps=10**-100000))


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
    f = Poly(2*x**3 - 7*x**2 + 4*x + 4)

    assert f.root(0) == -Rational(1, 2)
    assert f.root(1) == 2
    assert f.root(2) == 2
    pytest.raises(IndexError, lambda: f.root(3))

    assert Poly(x**5 + x + 1).root(0) == RootOf(x**3 - x**2 + 1, 0)


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

    assert Poly(f).real_roots() == [-Rational(1, 2), 2, 2]
    assert Poly(g).real_roots() == [RootOf(g, 0)]

    pytest.raises(PolynomialError, lambda: real_roots(1))


def test_all_roots():
    f = 2*x**3 - 7*x**2 + 4*x + 4
    g = x**3 + x + 1

    assert Poly(f).all_roots() == [-Rational(1, 2), 2, 2]
    assert Poly(f).all_roots(multiple=False) == [(-Rational(1, 2), 1), (2, 2)]
    assert Poly(g).all_roots() == [RootOf(g, 0), RootOf(g, 1), RootOf(g, 2)]


def test_nroots():
    assert Poly(0, x).nroots() == []
    assert Poly(1, x).nroots() == []

    assert Poly(x**2 - 1, x).nroots() == [-1.0, 1.0]
    assert Poly(x**2 + 1, x).nroots() == [-1.0*I, 1.0*I]

    roots = Poly(x**2 - 1, x).nroots()
    assert roots == [-1.0, 1.0]

    roots = Poly(x**2 + 1, x).nroots()
    assert roots == [-1.0*I, 1.0*I]

    roots = Poly(x**2/3 - Rational(1, 3), x).nroots()
    assert roots == [-1.0, 1.0]

    roots = Poly(x**2/3 + Rational(1, 3), x).nroots()
    assert roots == [-1.0*I, 1.0*I]

    assert Poly(x**2 + 2*I, x).nroots() == [-1.0 + 1.0*I, 1.0 - 1.0*I]
    assert Poly(
        x**2 + 2*I, x, extension=I).nroots() == [-1.0 + 1.0*I, 1.0 - 1.0*I]

    assert Poly(0.2*x + 0.1).nroots() == [-0.5]

    roots = nroots(x**5 + x + 1, n=5)
    eps = Float("1e-5")

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

    eps = Float("1e-6")

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

    pytest.raises(DomainError, lambda: Poly(x + y, x).nroots())
    pytest.raises(MultivariatePolynomialError, lambda: Poly(x + y).nroots())

    assert nroots(x**2 - 1) == [-1.0, 1.0]

    roots = nroots(x**2 - 1)
    assert roots == [-1.0, 1.0]

    assert nroots(x + I) == [-1.0*I]
    assert nroots(x + 2*I) == [-2.0*I]

    pytest.raises(PolynomialError, lambda: nroots(0))

    # issue sympy/sympy#8296
    f = Poly(x**4 - 1)
    assert f.nroots(2) == [w.evalf(2) for w in f.all_roots()]


def test_ground_roots():
    f = x**6 - 4*x**4 + 4*x**3 - x**2

    assert Poly(f).ground_roots() == {1: 2, 0: 2}
    assert ground_roots(f) == {1: 2, 0: 2}

    pytest.raises(MultivariatePolynomialError,
                  lambda: Poly(x + y).ground_roots())
    pytest.raises(ComputationFailed, lambda: ground_roots(1))


def test_nth_power_roots_poly():
    f = x**4 - x**2 + 1

    f_2 = (x**2 - x + 1)**2
    f_3 = (x**2 + 1)**2
    f_4 = (x**2 + x + 1)**2
    f_12 = (x - 1)**4

    assert nth_power_roots_poly(f, 1) == f

    pytest.raises(ValueError, lambda: nth_power_roots_poly(f, 0))
    pytest.raises(ValueError, lambda: nth_power_roots_poly(f, x))

    assert factor(nth_power_roots_poly(f, 2)) == f_2
    assert factor(nth_power_roots_poly(f, 3)) == f_3
    assert factor(nth_power_roots_poly(f, 4)) == f_4
    assert factor(nth_power_roots_poly(f, 12)) == f_12
    assert nth_power_roots_poly(f, 12, polys=True) == Poly(f_12)

    pytest.raises(MultivariatePolynomialError, lambda: nth_power_roots_poly(
        x + y, 2, x, y))
    pytest.raises(ComputationFailed, lambda: nth_power_roots_poly(1, 2))


def test_torational_factor_list():
    p = expand(((x**2 - 1)*(x - 2)).subs({x: x*(1 + sqrt(2))}))
    assert _torational_factor_list(p, x) == (-2, [
        (-x*(1 + sqrt(2))/2 + 1, 1),
        (-x*(1 + sqrt(2)) - 1, 1),
        (-x*(1 + sqrt(2)) + 1, 1)])

    p = expand(((x**2 - 1)*(x - 2)).subs({x: x*(1 + root(2, 4))}))
    assert _torational_factor_list(p, x) is None

    p = expand(((x**2 - 1)*(x - 2)).subs({x: x + sqrt(2)}))
    assert _torational_factor_list(p, x) == (1, [(x - 2 + sqrt(2), 1),
                                                 (x - 1 + sqrt(2), 1),
                                                 (x + 1 + sqrt(2), 1)])


def test_cancel():
    assert cancel(0) == 0
    assert cancel(7) == 7
    assert cancel(x) == x

    assert cancel(oo) == oo

    assert cancel((2, 3)) == (1, 2, 3)
    assert cancel((0, 1, 2)) == (0, 1, 2)

    assert cancel((1, 0), x) == (1, 1, 0)
    assert cancel((0, 1), x) == (1, 0, 1)

    f, g, p, q = 4*x**2 - 4, 2*x - 2, 2*x + 2, 1
    F, G, P, Q = [Poly(u, x) for u in (f, g, p, q)]

    assert F.cancel(G) == (1, P, Q)
    assert cancel((f, g)) == (1, p, q)
    assert cancel((f, g), x) == (1, p, q)
    assert cancel((f, g), (x,)) == (1, p, q)
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

    f = Poly(x**2 - a**2, x)
    g = Poly(x - a, x)

    F = Poly(x + a, x)
    G = Poly(1, x)

    assert cancel((f, g)) == (1, F, G)

    f = x**3 + (sqrt(2) - 2)*x**2 - (2*sqrt(2) + 3)*x - 3*sqrt(2)
    g = x**2 - 2

    assert cancel((f, g), extension=True) == (1, x**2 - 2*x - 3, x - sqrt(2))

    f = Poly(-2*x + 3, x)
    g = Poly(-x**9 + x**8 + x**6 - x**5 + 2*x**2 - 3*x + 1, x)

    assert cancel((f, g)) == (1, -f, -g)

    f = Poly(y, y, domain='ZZ(x)')
    g = Poly(1, y, domain='ZZ[x]')

    assert f.cancel(
        g) == (1, Poly(y, y, domain='ZZ(x)'), Poly(1, y, domain='ZZ(x)'))
    assert f.cancel(g, include=True) == (
        Poly(y, y, domain='ZZ(x)'), Poly(1, y, domain='ZZ(x)'))

    f = Poly(5*x*y + x, y, domain='ZZ(x)')
    g = Poly(2*x**2*y, y, domain='ZZ(x)')

    assert f.cancel(g, include=True) == (
        Poly(5*y + 1, y, domain='ZZ(x)'), Poly(2*x*y, y, domain='ZZ(x)'))

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


def test_reduced():
    f = 2*x**4 + y**2 - x**2 + y**3
    G = [x**3 - x, y**3 - y]

    Q = [2*x, 1]
    r = x**2 + y**2 + y

    assert reduced(f, G) == (Q, r)
    assert reduced(f, G, x, y) == (Q, r)

    H = groebner(G)

    assert H.reduce(f) == (Q, r)

    Q = [Poly(2*x, x, y), Poly(1, x, y)]
    r = Poly(x**2 + y**2 + y, x, y)

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
    pytest.raises(ComputationFailed, lambda: reduced(1, [1]))


def test_groebner():
    assert groebner([], x, y, z) == []

    assert groebner([x**2 + 1, y**4*x + x**3], x, y, order='lex') == [1 + x**2, -1 + y**4]
    assert groebner([x**2 + 1, y**4*x + x**3, x*y*z**3], x, y, z, order='grevlex') == [-1 + y**4, z**3, 1 + x**2]

    assert groebner([x**2 + 1, y**4*x + x**3], x, y, order='lex', polys=True) == \
        [Poly(1 + x**2, x, y), Poly(-1 + y**4, x, y)]
    assert groebner([x**2 + 1, y**4*x + x**3, x*y*z**3], x, y, z, order='grevlex', polys=True) == \
        [Poly(-1 + y**4, x, y, z), Poly(z**3, x, y, z), Poly(1 + x**2, x, y, z)]

    assert groebner([x**3 - 1, x**2 - 1]) == [x - 1]
    assert groebner([Eq(x**3, 1), Eq(x**2, 1)]) == [x - 1]

    F = [3*x**2 + y*z - 5*x - 1, 2*x + 3*x*y + y**2, x - 3*y + x*z - 2*z**2]
    f = z**9 - x**2*y**3 - 3*x*y**2*z + 11*y*z**2 + x**2*z**2 - 5

    G = groebner(F, x, y, z, modulus=7)

    assert G == [1 + x + y + 3*z + 2*z**2 + 2*z**3 + 6*z**4 + z**5,
                 1 + 3*y + y**2 + 6*z**2 + 3*z**3 + 3*z**4 + 3*z**5 + 4*z**6,
                 1 + 4*y + 4*z + y*z + 4*z**3 + z**4 + z**6,
                 6 + 6*z + z**2 + 4*z**3 + 3*z**4 + 6*z**5 + 3*z**6 + z**7]

    Q, r = reduced(f, G, x, y, z, modulus=7, polys=True)

    assert sum((q*g for q, g in zip(Q, G.polys)), r) == Poly(f, modulus=7)

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

    pytest.raises(ComputationFailed, lambda: groebner([1]))

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
                                   [x, y]))


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
    G = groebner(S, V, domain=QQ)
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
    G = groebner(S, V, domain=QQ)
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
    G.independent_sets == [[L5, L3, L2], [L6, L3]]
    assert G.dimension == 3

    # Algebraic Solution of Nonlinear Equation Systems in REDUCE, p.7.
    V = ax, bx, cx, gx, jx, lx, mx, nx, q = symbols('ax bx cx gx jx lx mx nx q')
    S = (ax*q - lx*q - mx, ax - gx*q - lx, bx*q**2 + cx*q - jx*q - nx,
         q*(-ax*q + lx*q + mx), q*(-ax + gx*q + lx))
    G = groebner(S, V, domain=QQ)
    assert G.independent_sets == [[cx, jx, lx, mx, nx, q], [cx, gx, jx, lx, mx, nx], [bx, cx, gx, jx, lx, nx]]
    assert G.dimension == 6


def test_GroebnerBasis():
    F = [x*y - 2*y, 2*y**2 - x**2]

    G = groebner(F, x, y, order='grevlex')
    assert groebner(F + [0], x, y, order='grevlex') == G

    assert G.args == ((y**3 - 2*y, x**2 - 2*y**2, x*y - 2*y), (x, y))

    H = [y**3 - 2*y, x**2 - 2*y**2, x*y - 2*y]
    P = [Poly(h, x, y) for h in H]

    assert isinstance(G, GroebnerBasis) is True

    assert len(G) == 3

    assert G[0] == H[0] and not G[0].is_Poly
    assert G[1] == H[1] and not G[1].is_Poly
    assert G[2] == H[2] and not G[2].is_Poly

    assert G[1:] == H[1:] and not any(g.is_Poly for g in G[1:])
    assert G[:2] == H[:2] and not any(g.is_Poly for g in G[1:])

    assert G.exprs == H
    assert G.polys == P
    assert G.gens == (x, y)
    assert G.domain == ZZ
    assert G.order == grevlex

    assert G == H
    assert G == tuple(H)
    assert G == P
    assert G == tuple(P)

    assert G != []

    G = groebner(F, x, y, order='grevlex', polys=True)

    assert G[0] == P[0] and G[0].is_Poly
    assert G[1] == P[1] and G[1].is_Poly
    assert G[2] == P[2] and G[2].is_Poly

    assert G[1:] == P[1:] and all(g.is_Poly for g in G[1:])
    assert G[:2] == P[:2] and all(g.is_Poly for g in G[1:])

    assert tuple(G) == (Poly(y**3 - 2*y, x, y), Poly(x**2 - 2*y**2, x, y),
                        Poly(x*y - 2*y, x, y))

    G = groebner(F, x, y, order='grevlex', polys=True)
    assert hash(G) == hash(groebner([Poly(_) for _ in F], order='grevlex'))

    assert (G == 1) is False


def test_poly():
    assert poly(x) == Poly(x, x)
    assert poly(y) == Poly(y, y)

    assert poly(x + y) == Poly(x + y, x, y)
    assert poly(x + sin(x)) == Poly(x + sin(x), x, sin(x))

    assert poly(x + y, wrt=y) == Poly(x + y, y, x)
    assert poly(x + sin(x), wrt=sin(x)) == Poly(x + sin(x), sin(x), x)

    assert poly(x*y + 2*x*z**2 + 17) == Poly(x*y + 2*x*z**2 + 17, x, y, z)

    assert poly(2*(y + z)**2 - 1) == Poly(2*y**2 + 4*y*z + 2*z**2 - 1, y, z)
    assert poly(
        x*(y + z)**2 - 1) == Poly(x*y**2 + 2*x*y*z + x*z**2 - 1, x, y, z)
    assert poly(2*x*(
        y + z)**2 - 1) == Poly(2*x*y**2 + 4*x*y*z + 2*x*z**2 - 1, x, y, z)

    assert poly(2*(
        y + z)**2 - x - 1) == Poly(2*y**2 + 4*y*z + 2*z**2 - x - 1, x, y, z)
    assert poly(x*(
        y + z)**2 - x - 1) == Poly(x*y**2 + 2*x*y*z + x*z**2 - x - 1, x, y, z)
    assert poly(2*x*(y + z)**2 - x - 1) == Poly(2*x*y**2 + 4*x*y*z + 2 *
                                                x*z**2 - x - 1, x, y, z)

    assert poly(x*y + (x + y)**2 + (x + z)**2) == \
        Poly(2*x*z + 3*x*y + y**2 + z**2 + 2*x**2, x, y, z)
    assert poly(x*y*(x + y)*(x + z)**2) == \
        Poly(x**3*y**2 + x*y**2*z**2 + y*x**2*z**2 + 2*z*x**2 *
             y**2 + 2*y*z*x**3 + y*x**4, x, y, z)

    assert poly(Poly(x + y + z, y, x, z)) == Poly(x + y + z, y, x, z)

    assert poly((x + y)**2, x) == Poly(x**2 + 2*x*y + y**2, x, domain=ZZ.poly_ring(y))
    assert poly((x + y)**2, x, expand=True) == Poly(x**2 + 2*x*y + y**2,
                                                    x, domain=ZZ.poly_ring(y))
    assert poly((x + y)**2, y) == Poly(x**2 + 2*x*y + y**2, y, domain=ZZ.poly_ring(x))

    assert poly(1, x) == Poly(1, x)
    pytest.raises(GeneratorsNeeded, lambda: poly(1))

    # issue sympy/sympy#6184
    assert poly(x + y, x, y) == Poly(x + y, x, y)
    assert poly(x + y, y, x) == Poly(x + y, y, x)

    # issue sympy/sympy#12400
    assert (poly(1/(1 + sqrt(2)), x) ==
            Poly(1/(1 + sqrt(2)), x,
                 domain=QQ.algebraic_field(1/(1 + sqrt(2)))))


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
    assert I * Poly(x, x) == Poly(I*x, x)
    assert Poly(x, x) * I == Poly(I*x, x)


def test_sympyissue_5786():
    e = (x - I*y)*(z - I*t)
    assert factor(expand(e), extension=[I]) == e


def test_noncommutative():
    class foo(Expr):
        is_commutative = False
    e = x/(x + x*y)
    c = 1/(1 + y)
    assert cancel(foo(e)) == foo(c)
    assert cancel(e + foo(e)) == c + foo(c)
    assert cancel(e*foo(c)) == c*foo(c)


def test_to_rational_coeffs():
    assert to_rational_coeffs(
        Poly(x**3 + y*x**2 + sqrt(y), x, domain='EX')) is None
    assert to_rational_coeffs(Poly(((x**2 - 1)*(x - 2)*y).subs({x: x*(1 + sqrt(2))}),
                                   x, y, domain='EX')) is None
    assert to_rational_coeffs(Poly(x**5 + sqrt(2)*x**2 + 1,
                              x, domain='EX')) is None


def test_sympyissue_9607():
    assert factor(1e-20*x - 7.292115e-5) == 1e-20*x - 7.292115e-5


def test_sympyissue_8754():
    z = 0.0001*(x*(x + (4.0*y))) + 0.0001*(y*(x + (4.0*y)))
    w = expand(z)
    v = factor(w)
    assert v == Mul(Float('10000.0', 15),
                    Float('0.0001', 15)*x + Float('0.0001', 15)*y,
                    Float('0.0001', 15)*x + Float('0.00040000000000000002', 15)*y,
                    evaluate=False)
    assert expand(v) == w


def test_factor_terms():
    # issue sympy/sympy#7067
    assert factor_list(x*(x + y)) == (1, [(x, 1), (x + y, 1)])
    assert sqf_list(x*(x + y)) == (1, [(x, 1), (x + y, 1)])


def test_sympyissue_8210():
    p = Poly(0, x)
    p2 = p.copy()
    assert id(p) != id(p2)
    assert p == p2


def test_sympyissue_11775():
    e = y**4 + x*y**3 + y**2 + x*y
    assert factor_list(e, y) == (1, [(y, 1), (x + y, 1), (y**2 + 1, 1)])


def test_sympyissue_5602():
    Poly(Integral(x, (x, 0, 1))*x + x**2, x)


def test_sympyissue_15798():
    o1 = Poly(x + y, x, y, z)
    o2 = o1.copy()
    assert o1 == o2
    p = Poly(x + sqrt(2), x)
    assert p == p.copy()


@pytest.mark.timeout(20)
def test_sympyissue_16222():
    Poly(x**100000000)
