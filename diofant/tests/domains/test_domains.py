"""Tests for classes defining properties of ground domains, e.g. ZZ, QQ, ZZ[x] ... """

import pytest

from diofant import Float, I, Integer, Poly, Rational, oo, root, sin, sqrt
from diofant.abc import x, y, z
from diofant.domains import CC, EX, FF, GF, QQ, RR, ZZ, QQ_python, ZZ_python
from diofant.domains.algebraicfield import AlgebraicField
from diofant.domains.complexfield import ComplexField
from diofant.domains.domainelement import DomainElement
from diofant.domains.groundtypes import PythonRational
from diofant.domains.realfield import RealField
from diofant.polys import RootOf, field, ring, roots
from diofant.polys.polyerrors import (CoercionFailed, DomainError,
                                      GeneratorsError, GeneratorsNeeded,
                                      NotInvertible, UnificationFailed)


__all__ = ()

ALG = QQ.algebraic_field(sqrt(2), sqrt(3))


def unify(K0, K1):
    return K0.unify(K1)


def test_Domain_interface():
    pytest.raises(TypeError, lambda: DomainElement().parent)

    assert RR(1).parent is RR
    assert CC(1).parent is CC

    assert RR.has_default_precision
    assert CC.has_default_precision

    RR3 = RealField(prec=53, dps=3)
    assert str(RR3(1.7611107002)) == '1.76'

    assert RealField(tol=3).tolerance == 3.0
    assert RealField(tol=0.1).tolerance == 0.1
    assert RealField(tol="0.1").tolerance == 0.1
    pytest.raises(ValueError, lambda: RealField(tol=object()))

    pytest.raises(AttributeError, lambda: CC.ring)
    pytest.raises(DomainError, lambda: CC.get_exact())

    assert str(EX(1)) == 'EX(1)'

    assert EX(1).as_expr() == Integer(1)
    assert bool(EX(1)) is True
    assert bool(EX(0)) is False


def test_Domain_unify():
    F3 = GF(3)

    assert unify(F3, F3) == F3
    assert unify(F3, ZZ) == F3
    assert unify(F3, QQ) == QQ
    assert unify(F3, ALG) == ALG
    assert unify(F3, RR) == RR
    assert unify(F3, CC) == CC
    assert unify(F3, ZZ.poly_ring(x)) == F3.poly_ring(x)
    assert unify(F3, ZZ.frac_field(x)) == F3.frac_field(x)
    assert unify(F3, EX) == EX

    assert unify(ZZ, F3) == F3
    assert unify(ZZ, ZZ) == ZZ
    assert unify(ZZ, QQ) == QQ
    assert unify(ZZ, ALG) == ALG
    assert unify(ZZ, RR) == RR
    assert unify(ZZ, CC) == CC
    assert unify(ZZ, ZZ.poly_ring(x)) == ZZ.poly_ring(x)
    assert unify(ZZ, ZZ.frac_field(x)) == ZZ.frac_field(x)
    assert unify(ZZ, EX) == EX

    assert unify(QQ, F3) == QQ
    assert unify(QQ, ZZ) == QQ
    assert unify(QQ, QQ) == QQ
    assert unify(QQ, ALG) == ALG
    assert unify(QQ, RR) == RR
    assert unify(QQ, CC) == CC
    assert unify(QQ, ZZ.poly_ring(x)) == QQ.poly_ring(x)
    assert unify(QQ, ZZ.frac_field(x)) == QQ.frac_field(x)
    assert unify(QQ, EX) == EX

    assert unify(RR, F3) == RR
    assert unify(RR, ZZ) == RR
    assert unify(RR, QQ) == RR
    assert unify(RR, ALG) == RR
    assert unify(RR, RR) == RR
    assert unify(RR, CC) == CC
    assert unify(RR, ZZ.poly_ring(x)) == RR.poly_ring(x)
    assert unify(RR, ZZ.frac_field(x)) == RR.frac_field(x)
    assert unify(RR, EX) == EX

    assert unify(CC, F3) == CC
    assert unify(CC, ZZ) == CC
    assert unify(CC, QQ) == CC
    assert unify(CC, ALG) == CC
    assert unify(CC, RR) == CC
    assert unify(CC, CC) == CC
    assert unify(CC, ZZ.poly_ring(x)) == CC.poly_ring(x)
    assert unify(CC, ZZ.frac_field(x)) == CC.frac_field(x)
    assert unify(CC, EX) == EX

    CC2 = ComplexField(prec=20)
    assert unify(CC, CC2) == unify(CC2, CC) == ComplexField(prec=CC.precision,
                                                            tol=CC2.tolerance)
    RR2 = RealField(prec=20)
    assert unify(RR, RR2) == unify(RR2, RR) == RealField(prec=RR.precision,
                                                         tol=RR2.tolerance)

    assert unify(ZZ.poly_ring(x), F3) == F3.poly_ring(x)
    assert unify(ZZ.poly_ring(x), ZZ) == ZZ.poly_ring(x)
    assert unify(ZZ.poly_ring(x), QQ) == QQ.poly_ring(x)
    assert unify(ZZ.poly_ring(x), ALG) == ALG.poly_ring(x)
    assert unify(ZZ.poly_ring(x), RR) == RR.poly_ring(x)
    assert unify(ZZ.poly_ring(x), CC) == CC.poly_ring(x)
    assert unify(ZZ.poly_ring(x), ZZ.poly_ring(x)) == ZZ.poly_ring(x)
    assert unify(ZZ.poly_ring(x), ZZ.frac_field(x)) == ZZ.frac_field(x)
    assert unify(ZZ.poly_ring(x), EX) == EX

    assert unify(ZZ.frac_field(x), F3) == F3.frac_field(x)
    assert unify(ZZ.frac_field(x), ZZ) == ZZ.frac_field(x)
    assert unify(ZZ.frac_field(x), QQ) == QQ.frac_field(x)
    assert unify(ZZ.frac_field(x), ALG) == ALG.frac_field(x)
    assert unify(ZZ.frac_field(x), RR) == RR.frac_field(x)
    assert unify(ZZ.frac_field(x), CC) == CC.frac_field(x)
    assert unify(ZZ.frac_field(x), ZZ.poly_ring(x)) == ZZ.frac_field(x)
    assert unify(ZZ.frac_field(x), ZZ.frac_field(x)) == ZZ.frac_field(x)
    assert unify(ZZ.frac_field(x), EX) == EX

    assert unify(EX, F3) == EX
    assert unify(EX, ZZ) == EX
    assert unify(EX, QQ) == EX
    assert unify(EX, ALG) == EX
    assert unify(EX, RR) == EX
    assert unify(EX, CC) == EX
    assert unify(EX, ZZ.poly_ring(x)) == EX
    assert unify(EX, ZZ.frac_field(x)) == EX
    assert unify(EX, EX) == EX


def test_Domain_unify_composite():
    assert unify(ZZ.poly_ring(x), ZZ) == ZZ.poly_ring(x)
    assert unify(ZZ.poly_ring(x), QQ) == QQ.poly_ring(x)
    assert unify(QQ.poly_ring(x), ZZ) == QQ.poly_ring(x)
    assert unify(QQ.poly_ring(x), QQ) == QQ.poly_ring(x)

    assert unify(ZZ, ZZ.poly_ring(x)) == ZZ.poly_ring(x)
    assert unify(QQ, ZZ.poly_ring(x)) == QQ.poly_ring(x)
    assert unify(ZZ, QQ.poly_ring(x)) == QQ.poly_ring(x)
    assert unify(QQ, QQ.poly_ring(x)) == QQ.poly_ring(x)

    assert unify(ZZ.poly_ring(x, y), ZZ) == ZZ.poly_ring(x, y)
    assert unify(ZZ.poly_ring(x, y), QQ) == QQ.poly_ring(x, y)
    assert unify(QQ.poly_ring(x, y), ZZ) == QQ.poly_ring(x, y)
    assert unify(QQ.poly_ring(x, y), QQ) == QQ.poly_ring(x, y)

    assert unify(ZZ, ZZ.poly_ring(x, y)) == ZZ.poly_ring(x, y)
    assert unify(QQ, ZZ.poly_ring(x, y)) == QQ.poly_ring(x, y)
    assert unify(ZZ, QQ.poly_ring(x, y)) == QQ.poly_ring(x, y)
    assert unify(QQ, QQ.poly_ring(x, y)) == QQ.poly_ring(x, y)

    assert unify(ZZ.frac_field(x), ZZ) == ZZ.frac_field(x)
    assert unify(ZZ.frac_field(x), QQ) == QQ.frac_field(x)
    assert unify(QQ.frac_field(x), ZZ) == QQ.frac_field(x)
    assert unify(QQ.frac_field(x), QQ) == QQ.frac_field(x)

    assert unify(ZZ, ZZ.frac_field(x)) == ZZ.frac_field(x)
    assert unify(QQ, ZZ.frac_field(x)) == QQ.frac_field(x)
    assert unify(ZZ, QQ.frac_field(x)) == QQ.frac_field(x)
    assert unify(QQ, QQ.frac_field(x)) == QQ.frac_field(x)

    assert unify(ZZ.frac_field(x, y), ZZ) == ZZ.frac_field(x, y)
    assert unify(ZZ.frac_field(x, y), QQ) == QQ.frac_field(x, y)
    assert unify(QQ.frac_field(x, y), ZZ) == QQ.frac_field(x, y)
    assert unify(QQ.frac_field(x, y), QQ) == QQ.frac_field(x, y)

    assert unify(ZZ, ZZ.frac_field(x, y)) == ZZ.frac_field(x, y)
    assert unify(QQ, ZZ.frac_field(x, y)) == QQ.frac_field(x, y)
    assert unify(ZZ, QQ.frac_field(x, y)) == QQ.frac_field(x, y)
    assert unify(QQ, QQ.frac_field(x, y)) == QQ.frac_field(x, y)

    assert unify(ZZ.poly_ring(x), ZZ.poly_ring(x)) == ZZ.poly_ring(x)
    assert unify(ZZ.poly_ring(x), QQ.poly_ring(x)) == QQ.poly_ring(x)
    assert unify(QQ.poly_ring(x), ZZ.poly_ring(x)) == QQ.poly_ring(x)
    assert unify(QQ.poly_ring(x), QQ.poly_ring(x)) == QQ.poly_ring(x)

    assert unify(ZZ.poly_ring(x, y), ZZ.poly_ring(x)) == ZZ.poly_ring(x, y)
    assert unify(ZZ.poly_ring(x, y), QQ.poly_ring(x)) == QQ.poly_ring(x, y)
    assert unify(QQ.poly_ring(x, y), ZZ.poly_ring(x)) == QQ.poly_ring(x, y)
    assert unify(QQ.poly_ring(x, y), QQ.poly_ring(x)) == QQ.poly_ring(x, y)

    assert unify(ZZ.poly_ring(x), ZZ.poly_ring(x, y)) == ZZ.poly_ring(x, y)
    assert unify(ZZ.poly_ring(x), QQ.poly_ring(x, y)) == QQ.poly_ring(x, y)
    assert unify(QQ.poly_ring(x), ZZ.poly_ring(x, y)) == QQ.poly_ring(x, y)
    assert unify(QQ.poly_ring(x), QQ.poly_ring(x, y)) == QQ.poly_ring(x, y)

    assert unify(ZZ.poly_ring(x, y), ZZ.poly_ring(x, z)) == ZZ.poly_ring(x, y, z)
    assert unify(ZZ.poly_ring(x, y), QQ.poly_ring(x, z)) == QQ.poly_ring(x, y, z)
    assert unify(QQ.poly_ring(x, y), ZZ.poly_ring(x, z)) == QQ.poly_ring(x, y, z)
    assert unify(QQ.poly_ring(x, y), QQ.poly_ring(x, z)) == QQ.poly_ring(x, y, z)

    assert unify(ZZ.frac_field(x), ZZ.frac_field(x)) == ZZ.frac_field(x)
    assert unify(ZZ.frac_field(x), QQ.frac_field(x)) == QQ.frac_field(x)
    assert unify(QQ.frac_field(x), ZZ.frac_field(x)) == QQ.frac_field(x)
    assert unify(QQ.frac_field(x), QQ.frac_field(x)) == QQ.frac_field(x)

    assert unify(ZZ.frac_field(x, y), ZZ.frac_field(x)) == ZZ.frac_field(x, y)
    assert unify(ZZ.frac_field(x, y), QQ.frac_field(x)) == QQ.frac_field(x, y)
    assert unify(QQ.frac_field(x, y), ZZ.frac_field(x)) == QQ.frac_field(x, y)
    assert unify(QQ.frac_field(x, y), QQ.frac_field(x)) == QQ.frac_field(x, y)

    assert unify(ZZ.frac_field(x), ZZ.frac_field(x, y)) == ZZ.frac_field(x, y)
    assert unify(ZZ.frac_field(x), QQ.frac_field(x, y)) == QQ.frac_field(x, y)
    assert unify(QQ.frac_field(x), ZZ.frac_field(x, y)) == QQ.frac_field(x, y)
    assert unify(QQ.frac_field(x), QQ.frac_field(x, y)) == QQ.frac_field(x, y)

    assert unify(ZZ.frac_field(x, y), ZZ.frac_field(x, z)) == ZZ.frac_field(x, y, z)
    assert unify(ZZ.frac_field(x, y), QQ.frac_field(x, z)) == QQ.frac_field(x, y, z)
    assert unify(QQ.frac_field(x, y), ZZ.frac_field(x, z)) == QQ.frac_field(x, y, z)
    assert unify(QQ.frac_field(x, y), QQ.frac_field(x, z)) == QQ.frac_field(x, y, z)

    assert unify(ZZ.poly_ring(x), ZZ.frac_field(x)) == ZZ.frac_field(x)
    assert unify(ZZ.poly_ring(x), QQ.frac_field(x)) == ZZ.frac_field(x)
    assert unify(QQ.poly_ring(x), ZZ.frac_field(x)) == ZZ.frac_field(x)
    assert unify(QQ.poly_ring(x), QQ.frac_field(x)) == QQ.frac_field(x)

    assert unify(ZZ.poly_ring(x, y), ZZ.frac_field(x)) == ZZ.frac_field(x, y)
    assert unify(ZZ.poly_ring(x, y), QQ.frac_field(x)) == ZZ.frac_field(x, y)
    assert unify(QQ.poly_ring(x, y), ZZ.frac_field(x)) == ZZ.frac_field(x, y)
    assert unify(QQ.poly_ring(x, y), QQ.frac_field(x)) == QQ.frac_field(x, y)

    assert unify(ZZ.poly_ring(x), ZZ.frac_field(x, y)) == ZZ.frac_field(x, y)
    assert unify(ZZ.poly_ring(x), QQ.frac_field(x, y)) == ZZ.frac_field(x, y)
    assert unify(QQ.poly_ring(x), ZZ.frac_field(x, y)) == ZZ.frac_field(x, y)
    assert unify(QQ.poly_ring(x), QQ.frac_field(x, y)) == QQ.frac_field(x, y)

    assert unify(ZZ.poly_ring(x, y), ZZ.frac_field(x, z)) == ZZ.frac_field(x, y, z)
    assert unify(ZZ.poly_ring(x, y), QQ.frac_field(x, z)) == ZZ.frac_field(x, y, z)
    assert unify(QQ.poly_ring(x, y), ZZ.frac_field(x, z)) == ZZ.frac_field(x, y, z)
    assert unify(QQ.poly_ring(x, y), QQ.frac_field(x, z)) == QQ.frac_field(x, y, z)

    assert unify(ZZ.frac_field(x), ZZ.poly_ring(x)) == ZZ.frac_field(x)
    assert unify(ZZ.frac_field(x), QQ.poly_ring(x)) == ZZ.frac_field(x)
    assert unify(QQ.frac_field(x), ZZ.poly_ring(x)) == ZZ.frac_field(x)
    assert unify(QQ.frac_field(x), QQ.poly_ring(x)) == QQ.frac_field(x)

    assert unify(ZZ.frac_field(x, y), ZZ.poly_ring(x)) == ZZ.frac_field(x, y)
    assert unify(ZZ.frac_field(x, y), QQ.poly_ring(x)) == ZZ.frac_field(x, y)
    assert unify(QQ.frac_field(x, y), ZZ.poly_ring(x)) == ZZ.frac_field(x, y)
    assert unify(QQ.frac_field(x, y), QQ.poly_ring(x)) == QQ.frac_field(x, y)

    assert unify(ZZ.frac_field(x), ZZ.poly_ring(x, y)) == ZZ.frac_field(x, y)
    assert unify(ZZ.frac_field(x), QQ.poly_ring(x, y)) == ZZ.frac_field(x, y)
    assert unify(QQ.frac_field(x), ZZ.poly_ring(x, y)) == ZZ.frac_field(x, y)
    assert unify(QQ.frac_field(x), QQ.poly_ring(x, y)) == QQ.frac_field(x, y)

    assert unify(ZZ.frac_field(x, y), ZZ.poly_ring(x, z)) == ZZ.frac_field(x, y, z)
    assert unify(ZZ.frac_field(x, y), QQ.poly_ring(x, z)) == ZZ.frac_field(x, y, z)
    assert unify(QQ.frac_field(x, y), ZZ.poly_ring(x, z)) == ZZ.frac_field(x, y, z)
    assert unify(QQ.frac_field(x, y), QQ.poly_ring(x, z)) == QQ.frac_field(x, y, z)


def test_Domain_unify_algebraic():
    sqrt5 = QQ.algebraic_field(sqrt(5))
    sqrt7 = QQ.algebraic_field(sqrt(7))
    sqrt57 = QQ.algebraic_field(sqrt(5), sqrt(7))

    assert sqrt5.unify(sqrt7) == sqrt57

    assert sqrt5.unify(sqrt5.poly_ring(x, y)) == sqrt5.poly_ring(x, y)
    assert sqrt5.poly_ring(x, y).unify(sqrt5) == sqrt5.poly_ring(x, y)

    assert sqrt5.unify(sqrt5.frac_field(x, y)) == sqrt5.frac_field(x, y)
    assert sqrt5.frac_field(x, y).unify(sqrt5) == sqrt5.frac_field(x, y)

    assert sqrt5.unify(sqrt7.poly_ring(x, y)) == sqrt57.poly_ring(x, y)
    assert sqrt5.poly_ring(x, y).unify(sqrt7) == sqrt57.poly_ring(x, y)

    assert sqrt5.unify(sqrt7.frac_field(x, y)) == sqrt57.frac_field(x, y)
    assert sqrt5.frac_field(x, y).unify(sqrt7) == sqrt57.frac_field(x, y)

    sqrt2 = QQ.algebraic_field(sqrt(2))
    r = RootOf(x**7 - x + 1, 0)
    rootof = QQ.algebraic_field(r)
    ans = QQ.algebraic_field(r + sqrt(2))
    assert sqrt2.unify(rootof) == rootof.unify(sqrt2) == ans

    # here domain created from tuple, not Expr
    p = Poly(x**3 - sqrt(2)*x - 1, x)
    sqrt2 = p.domain
    assert sqrt2.unify(rootof) == rootof.unify(sqrt2) == ans


def test_Domain_unify_with_symbols():
    pytest.raises(UnificationFailed, lambda: ZZ.poly_ring(x, y).unify(ZZ, (y, z)))
    pytest.raises(UnificationFailed, lambda: ZZ.unify(ZZ.poly_ring(x, y), (y, z)))


def test_Domain__contains__():
    assert (0 in EX) is True
    assert (0 in ZZ) is True
    assert (0 in QQ) is True
    assert (0 in RR) is True
    assert (0 in CC) is True
    assert (0 in ALG) is True
    assert (0 in ZZ.poly_ring(x, y)) is True
    assert (0 in QQ.poly_ring(x, y)) is True
    assert (0 in RR.poly_ring(x, y)) is True

    assert (-7 in EX) is True
    assert (-7 in ZZ) is True
    assert (-7 in QQ) is True
    assert (-7 in RR) is True
    assert (-7 in CC) is True
    assert (-7 in ALG) is True
    assert (-7 in ZZ.poly_ring(x, y)) is True
    assert (-7 in QQ.poly_ring(x, y)) is True
    assert (-7 in RR.poly_ring(x, y)) is True

    assert (17 in EX) is True
    assert (17 in ZZ) is True
    assert (17 in QQ) is True
    assert (17 in RR) is True
    assert (17 in CC) is True
    assert (17 in ALG) is True
    assert (17 in ZZ.poly_ring(x, y)) is True
    assert (17 in QQ.poly_ring(x, y)) is True
    assert (17 in RR.poly_ring(x, y)) is True

    assert (-Rational(1, 7) in EX) is True
    assert (-Rational(1, 7) in ZZ) is False
    assert (-Rational(1, 7) in QQ) is True
    assert (-Rational(1, 7) in RR) is True
    assert (-Rational(1, 7) in CC) is True
    assert (-Rational(1, 7) in ALG) is True
    assert (-Rational(1, 7) in ZZ.poly_ring(x, y)) is False
    assert (-Rational(1, 7) in QQ.poly_ring(x, y)) is True
    assert (-Rational(1, 7) in RR.poly_ring(x, y)) is True

    assert (Rational(3, 5) in EX) is True
    assert (Rational(3, 5) in ZZ) is False
    assert (Rational(3, 5) in QQ) is True
    assert (Rational(3, 5) in RR) is True
    assert (Rational(3, 5) in CC) is True
    assert (Rational(3, 5) in ALG) is True
    assert (Rational(3, 5) in ZZ.poly_ring(x, y)) is False
    assert (Rational(3, 5) in QQ.poly_ring(x, y)) is True
    assert (Rational(3, 5) in RR.poly_ring(x, y)) is True

    assert (3.0 in EX) is True
    assert (3.0 in ZZ) is True
    assert (3.0 in QQ) is True
    assert (3.0 in RR) is True
    assert (3.0 in CC) is True
    assert (3.0 in ALG) is True
    assert (3.0 in ZZ.poly_ring(x, y)) is True
    assert (3.0 in QQ.poly_ring(x, y)) is True
    assert (3.0 in RR.poly_ring(x, y)) is True

    assert (3.14 in EX) is True
    assert (3.14 in ZZ) is False
    assert (3.14 in QQ) is True
    assert (3.14 in RR) is True
    assert (3.14 in CC) is True
    assert (3.14 in ALG) is True
    assert (3.14 in ZZ.poly_ring(x, y)) is False
    assert (3.14 in QQ.poly_ring(x, y)) is True
    assert (3.14 in RR.poly_ring(x, y)) is True

    assert (oo in EX) is True
    assert (oo in ZZ) is False
    assert (oo in QQ) is False
    assert (oo in RR) is True
    assert (oo in CC) is True
    assert (oo in ALG) is False
    assert (oo in ZZ.poly_ring(x, y)) is False
    assert (oo in QQ.poly_ring(x, y)) is False
    assert (oo in RR.poly_ring(x, y)) is True

    assert (-oo in EX) is True
    assert (-oo in ZZ) is False
    assert (-oo in QQ) is False
    assert (-oo in RR) is True
    assert (-oo in CC) is True
    assert (-oo in ALG) is False
    assert (-oo in ZZ.poly_ring(x, y)) is False
    assert (-oo in QQ.poly_ring(x, y)) is False
    assert (-oo in RR.poly_ring(x, y)) is True

    assert (sqrt(7) in EX) is True
    assert (sqrt(7) in ZZ) is False
    assert (sqrt(7) in QQ) is False
    assert (sqrt(7) in RR) is True
    assert (sqrt(7) in CC) is True
    assert (sqrt(7) in ALG) is False
    assert (sqrt(7) in ZZ.poly_ring(x, y)) is False
    assert (sqrt(7) in QQ.poly_ring(x, y)) is False
    assert (sqrt(7) in RR.poly_ring(x, y)) is True

    assert (2*sqrt(3) + 1 in EX) is True
    assert (2*sqrt(3) + 1 in ZZ) is False
    assert (2*sqrt(3) + 1 in QQ) is False
    assert (2*sqrt(3) + 1 in RR) is True
    assert (2*sqrt(3) + 1 in CC) is True
    assert (2*sqrt(3) + 1 in ALG) is True
    assert (2*sqrt(3) + 1 in ZZ.poly_ring(x, y)) is False
    assert (2*sqrt(3) + 1 in QQ.poly_ring(x, y)) is False
    assert (2*sqrt(3) + 1 in RR.poly_ring(x, y)) is True

    assert (sin(1) in EX) is True
    assert (sin(1) in ZZ) is False
    assert (sin(1) in QQ) is False
    assert (sin(1) in RR) is True
    assert (sin(1) in CC) is True
    assert (sin(1) in ALG) is False
    assert (sin(1) in ZZ.poly_ring(x, y)) is False
    assert (sin(1) in QQ.poly_ring(x, y)) is False
    assert (sin(1) in RR.poly_ring(x, y)) is True

    assert (x**2 + 1 in EX) is True
    assert (x**2 + 1 in ZZ) is False
    assert (x**2 + 1 in QQ) is False
    assert (x**2 + 1 in RR) is False
    assert (x**2 + 1 in CC) is False
    assert (x**2 + 1 in ALG) is False
    assert (x**2 + 1 in ZZ.poly_ring(x)) is True
    assert (x**2 + 1 in QQ.poly_ring(x)) is True
    assert (x**2 + 1 in RR.poly_ring(x)) is True
    assert (x**2 + 1 in ZZ.poly_ring(x, y)) is True
    assert (x**2 + 1 in QQ.poly_ring(x, y)) is True
    assert (x**2 + 1 in RR.poly_ring(x, y)) is True

    assert (x**2 + y**2 in EX) is True
    assert (x**2 + y**2 in ZZ) is False
    assert (x**2 + y**2 in QQ) is False
    assert (x**2 + y**2 in RR) is False
    assert (x**2 + y**2 in CC) is False
    assert (x**2 + y**2 in ALG) is False
    assert (x**2 + y**2 in ZZ.poly_ring(x)) is False
    assert (x**2 + y**2 in QQ.poly_ring(x)) is False
    assert (x**2 + y**2 in RR.poly_ring(x)) is False
    assert (x**2 + y**2 in ZZ.poly_ring(x, y)) is True
    assert (x**2 + y**2 in QQ.poly_ring(x, y)) is True
    assert (x**2 + y**2 in RR.poly_ring(x, y)) is True

    assert (Rational(3, 2)*x/(y + 1) - z in QQ.poly_ring(x, y, z)) is False


def test_Domain_ring():
    assert ZZ.has_assoc_Ring is True
    assert QQ.has_assoc_Ring is True
    assert ZZ.poly_ring(x).has_assoc_Ring is True
    assert QQ.poly_ring(x).has_assoc_Ring is True
    assert ZZ.poly_ring(x, y).has_assoc_Ring is True
    assert QQ.poly_ring(x, y).has_assoc_Ring is True
    assert ZZ.frac_field(x).has_assoc_Ring is True
    assert QQ.frac_field(x).has_assoc_Ring is True
    assert ZZ.frac_field(x, y).has_assoc_Ring is True
    assert QQ.frac_field(x, y).has_assoc_Ring is True

    assert EX.has_assoc_Ring is False
    assert RR.has_assoc_Ring is False
    assert ALG.has_assoc_Ring is False

    assert ZZ.ring == ZZ
    assert QQ.ring == ZZ
    assert ZZ.poly_ring(x).ring == ZZ.poly_ring(x)
    assert QQ.poly_ring(x).ring == QQ.poly_ring(x)
    assert ZZ.poly_ring(x, y).ring == ZZ.poly_ring(x, y)
    assert QQ.poly_ring(x, y).ring == QQ.poly_ring(x, y)
    assert ZZ.frac_field(x).ring == ZZ.poly_ring(x)
    assert QQ.frac_field(x).ring == QQ.poly_ring(x)
    assert ZZ.frac_field(x, y).ring == ZZ.poly_ring(x, y)
    assert QQ.frac_field(x, y).ring == QQ.poly_ring(x, y)

    assert EX.ring == EX

    pytest.raises(AttributeError, lambda: RR.ring)
    pytest.raises(AttributeError, lambda: ALG.ring)


def test_Domain_field():
    assert EX.has_assoc_Field is True
    assert ZZ.has_assoc_Field is True
    assert QQ.has_assoc_Field is True
    assert RR.has_assoc_Field is True
    assert ALG.has_assoc_Field is True
    assert ZZ.poly_ring(x).has_assoc_Field is True
    assert QQ.poly_ring(x).has_assoc_Field is True
    assert ZZ.poly_ring(x, y).has_assoc_Field is True
    assert QQ.poly_ring(x, y).has_assoc_Field is True

    assert EX.field == EX
    assert ZZ.field == QQ
    assert QQ.field == QQ
    assert RR.field == RR
    assert ALG.field == ALG
    assert ZZ.poly_ring(x).field == ZZ.frac_field(x)
    assert QQ.poly_ring(x).field == QQ.frac_field(x)
    assert ZZ.poly_ring(x, y).field == ZZ.frac_field(x, y)
    assert QQ.poly_ring(x, y).field == QQ.frac_field(x, y)


def test_Domain_get_exact():
    assert EX.get_exact() == EX
    assert ZZ.get_exact() == ZZ
    assert QQ.get_exact() == QQ
    assert RR.get_exact() == QQ
    assert ALG.get_exact() == ALG
    assert ZZ.poly_ring(x).get_exact() == ZZ.poly_ring(x)
    assert QQ.poly_ring(x).get_exact() == QQ.poly_ring(x)
    assert ZZ.poly_ring(x, y).get_exact() == ZZ.poly_ring(x, y)
    assert QQ.poly_ring(x, y).get_exact() == QQ.poly_ring(x, y)
    assert ZZ.frac_field(x).get_exact() == ZZ.frac_field(x)
    assert QQ.frac_field(x).get_exact() == QQ.frac_field(x)
    assert ZZ.frac_field(x, y).get_exact() == ZZ.frac_field(x, y)
    assert QQ.frac_field(x, y).get_exact() == QQ.frac_field(x, y)


def test_Domain_convert():
    assert QQ.convert(10e-52) == QQ(1684996666696915, 1684996666696914987166688442938726917102321526408785780068975640576)

    R, x = ring("x", ZZ)
    assert ZZ.convert(x - x) == 0
    assert ZZ.convert(x - x, R) == 0

    F3 = FF(3)
    assert F3.convert(Float(2.0)) == F3.dtype(2)
    assert F3.convert(PythonRational(2, 1)) == F3.dtype(2)
    pytest.raises(CoercionFailed, lambda: F3.convert(PythonRational(1, 2)))
    assert F3.convert(2.0) == F3.dtype(2)
    pytest.raises(CoercionFailed, lambda: F3.convert(2.1))

    assert RR.convert(CC(1)) == RR(1)
    pytest.raises(CoercionFailed, lambda: RR.convert(CC(1, 2)))

    assert QQ.convert(ALG(1), ALG) == QQ(1)
    pytest.raises(CoercionFailed, lambda: QQ.convert(ALG([1, 1]), ALG))

    assert ZZ.convert(ALG(1), ALG) == ZZ(1)
    pytest.raises(CoercionFailed, lambda: ZZ.convert(ALG([1, 1]), ALG))

    assert EX.convert(ALG([1, 1]), ALG) == sqrt(2) + sqrt(3) + 1

    ALG2 = QQ.algebraic_field(sqrt(2))
    a2 = ALG2.convert(sqrt(2))
    a = ALG.convert(a2, ALG2)
    assert a.rep.to_dense() == [QQ(1, 2), 0, -QQ(9, 2), 0]
    assert RR.convert(a) == RR(1.4142135623730951)
    assert CC.convert(a) == CC(1.4142135623730951)

    assert ZZ_python.convert(3.0) == ZZ_python.dtype(3)
    pytest.raises(CoercionFailed, lambda: ZZ_python.convert(3.2))

    assert CC.convert(QQ_python(1, 2)) == CC(0.5)
    CC01 = ComplexField(tol=0.1)
    assert CC.convert(CC01(0.3)) == CC(0.3)

    assert RR.convert(complex(2 + 0j)) == RR(2)
    pytest.raises(CoercionFailed, lambda: RR.convert(complex(2 + 3j)))

    assert ALG.convert(EX(sqrt(2)), EX) == ALG.from_expr(sqrt(2))
    pytest.raises(CoercionFailed, lambda: ALG.convert(EX(sqrt(5)), EX))

    pytest.raises(CoercionFailed, lambda: ALG2.convert(ALG.unit))


def test_arithmetics():
    assert ZZ.rem(ZZ(2), ZZ(3)) == 2
    assert ZZ.div(ZZ(2), ZZ(3)) == (0, 2)
    assert QQ.rem(QQ(2, 3), QQ(4, 7)) == 0
    assert QQ.div(QQ(2, 3), QQ(4, 7)) == (QQ(7, 6), 0)

    assert QQ_python.factorial(QQ_python(7, 2)) == 6

    assert CC.gcd(CC(1), CC(2)) == 1
    assert CC.lcm(CC(1), CC(2)) == 2

    assert EX(Rational(2, 3)).numerator == 2
    assert EX(Rational(2, 3)).denominator == 3

    assert abs(EX(-2)) == 2

    assert -EX(2) == -2
    assert 2 + EX(3) == EX(3) + 2 == 5
    assert 2 - EX(3) == EX(2) - 3 == -1
    assert 2*EX(3) == EX(3)*2 == 6
    assert 2/EX(3) == EX(2)/3 == EX(Rational(2, 3))
    assert EX(2)**2 == 4

    pytest.raises(TypeError, lambda: EX(2) + object())
    pytest.raises(TypeError, lambda: EX(2) - object())
    pytest.raises(TypeError, lambda: EX(2)*object())
    pytest.raises(TypeError, lambda: EX(2)**object())
    pytest.raises(TypeError, lambda: EX(2)/object())

    F11 = FF(11)
    assert +F11(2) == F11(2)
    assert F11(5) + 7 == 7 + F11(5) == F11(1)
    assert F11(5) - 7 == 5 - F11(7) == F11(9)
    assert F11(5)*7 == 7*F11(5) == F11(2)
    assert F11(5)/9 == 5/F11(9) == F11(3)
    assert F11(4) % 9 == 4 % F11(9) == F11(4)

    pytest.raises(TypeError, lambda: F11(2) + object())
    pytest.raises(TypeError, lambda: F11(2) - object())
    pytest.raises(TypeError, lambda: F11(2)*object())
    pytest.raises(TypeError, lambda: F11(2)**object())
    pytest.raises(TypeError, lambda: F11(2)/object())
    pytest.raises(TypeError, lambda: F11(2) % object())
    pytest.raises(TypeError, lambda: object() % F11(2))


def test_Ring():
    assert ZZ(3).numerator == 3
    assert ZZ(3).denominator == 1


def test_PolynomialRing__init():
    pytest.raises(GeneratorsNeeded, lambda: ZZ.poly_ring())


def test_FractionField__init():
    pytest.raises(GeneratorsNeeded, lambda: ZZ.frac_field())


def test_inject():
    assert ZZ.inject(x, y, z) == ZZ.poly_ring(x, y, z)
    assert ZZ.poly_ring(x).inject(y, z) == ZZ.poly_ring(x, y, z)
    assert ZZ.frac_field(x).inject(y, z) == ZZ.frac_field(x, y, z)
    pytest.raises(GeneratorsError, lambda: ZZ.poly_ring(x).inject(x))


def test_Domain___eq__():
    assert (ZZ.poly_ring(x, y) == ZZ.poly_ring(x, y)) is True
    assert (QQ.poly_ring(x, y) == QQ.poly_ring(x, y)) is True

    assert (ZZ.poly_ring(x, y) == QQ.poly_ring(x, y)) is False
    assert (QQ.poly_ring(x, y) == ZZ.poly_ring(x, y)) is False

    assert (ZZ.frac_field(x, y) == ZZ.frac_field(x, y)) is True
    assert (QQ.frac_field(x, y) == QQ.frac_field(x, y)) is True

    assert (ZZ.frac_field(x, y) == QQ.frac_field(x, y)) is False
    assert (QQ.frac_field(x, y) == ZZ.frac_field(x, y)) is False


def test_Domain__algebraic_field():
    alg = ZZ.algebraic_field(sqrt(3))
    assert alg.minpoly == Poly(x**2 - 3)
    assert alg.domain == QQ
    assert alg.from_expr(sqrt(3)).denominator == 1
    assert alg.from_expr(2*sqrt(3)).denominator == 1
    assert alg.from_expr(sqrt(3)/2).denominator == 2
    assert alg([QQ(7, 38), QQ(3, 2)]).denominator == 38

    alg = QQ.algebraic_field(sqrt(2))
    assert alg.minpoly == Poly(x**2 - 2)
    assert alg.domain == QQ

    alg = QQ.algebraic_field(sqrt(2), sqrt(3))
    assert alg.minpoly == Poly(x**4 - 10*x**2 + 1)
    assert alg.domain == QQ

    assert alg(1).numerator == alg(1)
    assert alg.from_expr(sqrt(3)/2).numerator == alg.from_expr(2*sqrt(3))
    assert alg.from_expr(sqrt(3)/2).denominator == 4

    pytest.raises(DomainError, lambda: AlgebraicField(ZZ, sqrt(2)))

    assert alg.characteristic == 0

    assert alg.is_RealAlgebraicField is True

    assert int(alg(2)) == 2
    assert int(alg.from_expr(Rational(3, 2))) == 1

    alg = QQ.algebraic_field(I)
    assert alg.algebraic_field(I) == alg
    assert alg.is_RealAlgebraicField is False
    pytest.raises(TypeError, lambda: int(alg([1, 1])))

    alg = QQ.algebraic_field(sqrt(2)).algebraic_field(sqrt(3))
    assert alg.minpoly == Poly(x**2 - 3, x, domain=QQ.algebraic_field(sqrt(2)))

    # issue sympy/sympy#14476
    assert QQ.algebraic_field(Rational(1, 7)) is QQ

    alg = QQ.algebraic_field(sqrt(2)).algebraic_field(I)
    assert alg.from_expr(2*sqrt(2) + I/3) == alg([alg.domain([1])/3,
                                                  alg.domain([2, 0])])
    alg2 = QQ.algebraic_field(sqrt(2))
    assert alg2.from_expr(sqrt(2)) == alg2.convert(alg.from_expr(sqrt(2)))

    eq = -x**3 + 2*x**2 + 3*x - 2
    rs = roots(eq, multiple=True)
    alg = QQ.algebraic_field(rs[0])
    assert alg.is_RealAlgebraicField

    alg1 = QQ.algebraic_field(I)
    alg2 = QQ.algebraic_field(sqrt(2)).algebraic_field(I)
    assert alg1 != alg2

    alg3 = QQ.algebraic_field(RootOf(4*x**7 + x - 1, 0))
    assert alg3.is_RealAlgebraicField
    assert int(alg3.unit) == 2
    assert 2.772 > alg3.unit > 2.771
    assert int(alg3([3, 17, 11, -1, 2])) == 622
    assert int(alg3([1, QQ(-11, 4), QQ(125326976730518, 44208605852241),
                     QQ(-16742151878022, 12894796053515),
                     QQ(2331359268715, 10459004949272)])) == 18

    alg4 = QQ.algebraic_field(sqrt(2) + I)
    assert alg4.convert(alg2.unit) == alg4.from_expr(I)


def test_PolynomialRing_from_FractionField():
    F,  x, y = field("x,y", ZZ)
    R,  X, Y = ring("x,y", ZZ)

    f = (x**2 + y**2)/(x + 1)
    g = (x**2 + y**2)/4
    h = x**2 + y**2

    pytest.raises(CoercionFailed, lambda: R.convert(f, F))
    pytest.raises(CoercionFailed, lambda: R.convert(g, F))
    assert R.convert(h, F) == X**2 + Y**2

    F,  x, y = field("x,y", QQ)
    R,  X, Y = ring("x,y", QQ)

    f = (x**2 + y**2)/(x + 1)
    g = (x**2 + y**2)/4
    h = x**2 + y**2

    pytest.raises(CoercionFailed, lambda: R.convert(f, F))
    assert R.convert(g, F) == X**2/4 + Y**2/4
    assert R.convert(h, F) == X**2 + Y**2


def test_FractionField_from_PolynomialRing():
    R,  x, y = ring("x,y", QQ)
    F,  X, Y = field("x,y", ZZ)

    f = 3*x**2 + 5*y**2
    g = x**2/3 + y**2/5

    assert F.convert(f, R) == 3*X**2 + 5*Y**2
    assert F.convert(g, R) == (5*X**2 + 3*Y**2)/15

    RALG,  u, v = ring("u,v", ALG)
    pytest.raises(CoercionFailed,
                  lambda: F.convert(3*u**2 + 5*sqrt(2)*v**2))


def test_FractionField_convert():
    F,  X, Y = field("x,y", QQ)
    F.convert(QQ_python(1, 3)) == F.one/3


def test_FF_of_type():
    assert isinstance(FF(3)(1), FF(3).dtype) is True
    assert isinstance(FF(5)(3), FF(5).dtype) is True
    assert isinstance(FF(7)(3), FF(5).dtype) is False


def test___eq__():
    assert not QQ.poly_ring(x) == ZZ.poly_ring(x)
    assert not QQ.frac_field(x) == ZZ.frac_field(x)

    assert EX(1) != EX(2)

    F11 = FF(11)
    assert F11(2) != F11(3)
    assert F11(2) != object()


def test_RealField_from_expr():
    assert RR.convert(Integer(0)) == RR.dtype(0)
    assert RR.convert(Float(0.0)) == RR.dtype(0.0)
    assert RR.convert(Integer(1)) == RR.dtype(1)
    assert RR.convert(Float(1.0)) == RR.dtype(1.0)
    assert RR.convert(sin(1)) == RR.dtype(sin(1).evalf())
    assert RR.convert(oo) == RR("+inf")
    assert RR.convert(-oo) == RR("-inf")
    pytest.raises(CoercionFailed, lambda: RR.convert(x))


def test_AlgebraicElement():
    A = QQ.algebraic_field(I)

    rep = [QQ(1), QQ(1)]
    mod = [QQ(1), QQ(0), QQ(1)]

    f = A(rep)

    assert f.rep.to_dense() == rep
    assert f.mod.to_dense() == mod
    assert f.domain.domain == QQ

    f = A(1)

    assert f.rep.to_dense() == [QQ(1)]
    assert f.mod.to_dense() == mod
    assert f.domain.domain == QQ

    f = A([QQ(3, 2)])

    assert f.rep.to_dense() == [QQ(3, 2)]
    assert f.mod.to_dense() == mod
    assert f.domain.domain == QQ

    B = QQ.algebraic_field(I*sqrt(2))

    a = A([QQ(1), QQ(1)])
    b = B([QQ(1), QQ(1)])

    assert (a == a) is True
    assert (a != a) is False

    assert (a == b) is False
    assert (a != b) is True

    b = A([QQ(1), QQ(2)])

    assert (a == b) is False
    assert (a != b) is True

    assert A([1, 1]) == A([int(1), int(1)])
    assert hash(A([1, 1])) == hash(A([int(1), int(1)]))

    assert a.to_dict() == {(0,): QQ(1), (1,): QQ(1)}

    assert bool(A([])) is False
    assert bool(A([QQ(1)])) is True

    a = A([QQ(1), -QQ(1), QQ(2)])
    assert a.LC() == -1
    assert a.rep.to_dense() == [-1, 1]

    A = QQ.algebraic_field(root(2, 3))

    assert A.unit > 0
    assert A.unit >= 0
    assert (A.unit < 0) is False
    assert (A.unit <= 0) is False
    pytest.raises(TypeError, lambda: A.unit > x)
    pytest.raises(TypeError, lambda: QQ.algebraic_field(I).unit > 0)

    assert abs(+A.unit) == A.unit
    assert abs(-A.unit) == A.unit

    a = A([QQ(2), QQ(-1), QQ(1)])
    b = A([QQ(1), QQ(2)])

    c = A([QQ(-2), QQ(1), QQ(-1)])

    assert +a == a
    assert -a == c

    c = A([QQ(2), QQ(0), QQ(3)])

    assert a + b == c
    assert b + a == c

    assert c + 1 == A([QQ(2), QQ(0), QQ(4)])
    pytest.raises(TypeError, lambda: c + "x")
    pytest.raises(TypeError, lambda: "x" + c)

    c = A([QQ(2), QQ(-2), QQ(-1)])

    assert a - b == c

    c = A([QQ(-2), QQ(2), QQ(1)])

    assert b - a == c

    assert c - 1 == A([QQ(-2), QQ(2), QQ(0)])
    pytest.raises(TypeError, lambda: c - "x")
    pytest.raises(TypeError, lambda: "x" - c)

    c = A([QQ(3), QQ(-1), QQ(6)])

    assert a*b == c
    assert b*a == c

    assert c*2 == A([QQ(6), QQ(-2), QQ(12)])
    pytest.raises(TypeError, lambda: c*"x")
    pytest.raises(TypeError, lambda: "x"*c)

    c = A([QQ(11, 10), -QQ(1, 5), -QQ(3, 5)])

    assert c/2 == A([QQ(11, 20), -QQ(1, 10), -QQ(3, 10)])
    pytest.raises(TypeError, lambda: c/"x")
    pytest.raises(TypeError, lambda: "x"/c)

    c = A([QQ(-1, 43), QQ(9, 43), QQ(5, 43)])

    assert a**0 == A(1)
    assert a**1 == a
    assert a**-1 == c
    pytest.raises(TypeError, lambda: a**QQ(1, 2))

    assert a*a**(-1) == A(1)
    assert 1/a == a**(-1)

    A = QQ.algebraic_field(I)

    a = A([QQ(1, 2), QQ(1), QQ(2)])
    b = A([ZZ(1), ZZ(1), ZZ(2)])
    c = A([QQ(3, 2), QQ(2), QQ(4)])
    assert a + b == b + a == c


def test_ModularInteger():
    F3 = FF(3)

    a = F3(0)
    assert isinstance(a, F3.dtype) and a == 0
    a = F3(1)
    assert isinstance(a, F3.dtype) and a == 1
    a = F3(2)
    assert isinstance(a, F3.dtype) and a == 2
    a = F3(3)
    assert isinstance(a, F3.dtype) and a == 0
    a = F3(4)
    assert isinstance(a, F3.dtype) and a == 1

    a = F3(F3(0))
    assert isinstance(a, F3.dtype) and a == 0
    a = F3(F3(1))
    assert isinstance(a, F3.dtype) and a == 1
    a = F3(F3(2))
    assert isinstance(a, F3.dtype) and a == 2
    a = F3(F3(3))
    assert isinstance(a, F3.dtype) and a == 0
    a = F3(F3(4))
    assert isinstance(a, F3.dtype) and a == 1

    a = -F3(1)
    assert isinstance(a, F3.dtype) and a == 2
    a = -F3(2)
    assert isinstance(a, F3.dtype) and a == 1

    a = 2 + F3(2)
    assert isinstance(a, F3.dtype) and a == 1
    a = F3(2) + 2
    assert isinstance(a, F3.dtype) and a == 1
    a = F3(2) + F3(2)
    assert isinstance(a, F3.dtype) and a == 1
    a = F3(2) + F3(2)
    assert isinstance(a, F3.dtype) and a == 1

    a = 3 - F3(2)
    assert isinstance(a, F3.dtype) and a == 1
    a = F3(3) - 2
    assert isinstance(a, F3.dtype) and a == 1
    a = F3(3) - F3(2)
    assert isinstance(a, F3.dtype) and a == 1
    a = F3(3) - F3(2)
    assert isinstance(a, F3.dtype) and a == 1

    a = 2*F3(2)
    assert isinstance(a, F3.dtype) and a == 1
    a = F3(2)*2
    assert isinstance(a, F3.dtype) and a == 1
    a = F3(2)*F3(2)
    assert isinstance(a, F3.dtype) and a == 1
    a = F3(2)*F3(2)
    assert isinstance(a, F3.dtype) and a == 1

    a = 2/F3(2)
    assert isinstance(a, F3.dtype) and a == 1
    a = F3(2)/2
    assert isinstance(a, F3.dtype) and a == 1
    a = F3(2)/F3(2)
    assert isinstance(a, F3.dtype) and a == 1
    a = F3(2)/F3(2)
    assert isinstance(a, F3.dtype) and a == 1

    a = 1 % F3(2)
    assert isinstance(a, F3.dtype) and a == 1
    a = F3(1) % 2
    assert isinstance(a, F3.dtype) and a == 1
    a = F3(1) % F3(2)
    assert isinstance(a, F3.dtype) and a == 1
    a = F3(1) % F3(2)
    assert isinstance(a, F3.dtype) and a == 1

    a = F3(2)**0
    assert isinstance(a, F3.dtype) and a == 1
    a = F3(2)**1
    assert isinstance(a, F3.dtype) and a == 2
    a = F3(2)**2
    assert isinstance(a, F3.dtype) and a == 1

    assert bool(F3(3)) is False
    assert bool(F3(4)) is True

    F5 = FF(5)

    a = F5(1)**(-1)
    assert isinstance(a, F5.dtype) and a == 1
    a = F5(2)**(-1)
    assert isinstance(a, F5.dtype) and a == 3
    a = F5(3)**(-1)
    assert isinstance(a, F5.dtype) and a == 2
    a = F5(4)**(-1)
    assert isinstance(a, F5.dtype) and a == 4

    pytest.raises(NotInvertible, lambda: F5(0)**(-1))
    pytest.raises(NotInvertible, lambda: F5(5)**(-1))

    pytest.raises(ValueError, lambda: FF(0))
    pytest.raises(ValueError, lambda: FF(2.1))
    pytest.raises(ValueError, lambda: FF(6))
    pytest.raises(ValueError, lambda: FF(9, [1, 0]))
    pytest.raises(ValueError, lambda: FF(9, [1, 1, 1]))

    assert F5.is_positive(a) is True
    assert F5.is_negative(a) is False
    assert F5.is_positive(F5.zero) is False
    assert F5.is_negative(F5.zero) is False

    assert F5 == FF(5, [1, 0])

    F9 = FF(9, [1, 0, 1])

    assert F9.order == 9
    assert F9.characteristic == 3

    assert F9.zero == F9([0])
    assert F9.one == F9([1])
    assert F9(F9([2, 1])) == F9([2, 1])
    assert F9([2, 1]) + F9([1, 2]) == F9([1, 2]) + F9([2, 1]) == F9.zero
    assert F9([2, 1]) + F9([1, 3]) == F9.one
    assert F9([2, 1]) * F9([1, 3]) == F9([1, 1])
    assert sum(F9.one for _ in range(3)) == F9.zero

    assert int(F9.zero) == 0
    assert int(F9.one) == 1
    assert int(F9([1, 1])) == int(F9(4)) == 4
    assert int(F9([2, 0])) == int(F9(6)) == 6

    F81 = FF(3, [1, 0, 0, 1, 2])

    assert F81([1, 2, 1])*F81([2, 2, 2]) == F81([2, 1, 1])
    assert F81([5, 4, 3, 2, 1]) == F81([1, 0, 0, 0])
    assert F81([0, 1, 0]) == F81([1, 0])
    assert 1 + F81([1, 1]) == F81([1, 2])

    F8 = FF(8)

    assert F8.order == 8
    assert F8.characteristic == 2
    assert F8.dtype.mod.to_dense() == [1, 0, 1, 1]
    assert int(F8([1, 0, 1])) == int(F8(5)) == 5
    assert int(F8(-1)) == int(F8(7)) == 7

    F4 = FF(2, [1, 1, 1])

    assert F4.order == 4
    assert F4.characteristic == 2


def test_QQ_int():
    assert int(QQ(2**2000, 3**1250)) == 455431
    assert int(QQ(2**100, 3)) == 422550200076076467165567735125


def test_RR_double():
    assert RR(3.14) > 1e-50
    assert RR(1e-13) > 1e-50
    assert RR(1e-14) > 1e-50
    assert RR(1e-15) > 1e-50
    assert RR(1e-20) > 1e-50
    assert RR(1e-40) > 1e-50


def test_RR_Float():
    f1 = Float("1.01", 15)
    f2 = Float("1.0000000000000000000001")
    assert f1._prec == 53
    assert f2._prec == 80
    assert RR(f1)-1 > 1e-50
    assert RR(f2)-1 < 1e-50  # RR's precision is lower than f2's

    RR2 = RealField(prec=f2._prec)
    assert RR2(f1)-1 > 1e-50
    assert RR2(f2)-1 > 1e-50  # RR's precision is equal to f2's

    a = RR(2.1)
    assert a.numerator == a and a.denominator == 1


def test_CC_double():
    assert CC(3.14).real > 1e-50
    assert CC(1e-13).real > 1e-50
    assert CC(1e-14).real > 1e-50
    assert CC(1e-15).real > 1e-50
    assert CC(1e-20).real > 1e-50
    assert CC(1e-40).real > 1e-50

    assert CC(3.14j).imag > 1e-50
    assert CC(1e-13j).imag > 1e-50
    assert CC(1e-14j).imag > 1e-50
    assert CC(1e-15j).imag > 1e-50
    assert CC(1e-20j).imag > 1e-50
    assert CC(1e-40j).imag > 1e-50

    a = CC(2.1 + 1j)
    assert a.numerator == a and a.denominator == 1


def test_almosteq():
    assert CC.almosteq(CC(2), 3) is False
    assert CC.almosteq(2, CC(3)) is False
    assert CC.almosteq(2, CC(2.5), 0.1) is False
    assert CC.almosteq(2, CC(2.5), 1.0) is True

    assert RR.almosteq(5, RR(2), 1) is True
    assert RR._context.almosteq(RR(2), 1, None, 1) is True


def test_to_expr():
    assert CC.to_expr(1 - 2j) == 1 - 2*I


def test_EX():
    assert EX.is_positive(EX(2))
    assert EX.is_negative(EX(-1))

    assert (EX(1)/2).numerator == 1
    assert (EX(1)/2).denominator == 2


def test_sympyissue_13545():
    assert Poly(x + 1, x, modulus=2) + 1 == Poly(x, x, modulus=2)
    pytest.raises(NotImplementedError,
                  lambda: Poly(x, modulus=2) + Poly(x, modulus=3))


def test_sympyissue_14294():
    A = QQ.algebraic_field(I)

    a = A.convert(I)
    assert A.convert(a) == a
