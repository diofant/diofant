"""Tests for classes defining properties of ground domains, e.g. ZZ, QQ, ZZ[x] ... """

import pytest

from diofant import Float, I, Integer, Poly, Rational, oo, sin, sqrt
from diofant.abc import x, y, z
from diofant.domains import CC, EX, FF, GF, QQ, RR, ZZ, QQ_python, ZZ_python
from diofant.domains.algebraicfield import AlgebraicField
from diofant.domains.complexfield import ComplexField
from diofant.domains.domainelement import DomainElement
from diofant.domains.groundtypes import PythonRational
from diofant.domains.realfield import RealField
from diofant.polys import RootOf, field, ring
from diofant.polys.polyerrors import (CoercionFailed, DomainError,
                                      GeneratorsError, GeneratorsNeeded,
                                      NotInvertible, UnificationFailed)


__all__ = ()

ALG = QQ.algebraic_field(sqrt(2), sqrt(3))


def unify(K0, K1):
    return K0.unify(K1)


def test_Domain_interface():
    pytest.raises(NotImplementedError, lambda: DomainElement().parent())

    assert RR(1).parent() is RR
    assert CC(1).parent() is CC

    assert RR.has_default_precision
    assert CC.has_default_precision

    RR3 = RealField(prec=53, dps=3)
    assert str(RR3(1.7611107002)) == '1.76'

    assert RealField(tol=3).tolerance == 3.0
    assert RealField(tol=0.1).tolerance == 0.1
    assert RealField(tol="0.1").tolerance == 0.1
    pytest.raises(ValueError, lambda: RealField(tol=object()))

    pytest.raises(DomainError, lambda: CC.get_ring())
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
    assert unify(F3, ZZ[x]) == F3[x]
    assert unify(F3, ZZ.frac_field(x)) == F3.frac_field(x)
    assert unify(F3, EX) == EX

    assert unify(ZZ, F3) == F3
    assert unify(ZZ, ZZ) == ZZ
    assert unify(ZZ, QQ) == QQ
    assert unify(ZZ, ALG) == ALG
    assert unify(ZZ, RR) == RR
    assert unify(ZZ, CC) == CC
    assert unify(ZZ, ZZ[x]) == ZZ[x]
    assert unify(ZZ, ZZ.frac_field(x)) == ZZ.frac_field(x)
    assert unify(ZZ, EX) == EX

    assert unify(QQ, F3) == QQ
    assert unify(QQ, ZZ) == QQ
    assert unify(QQ, QQ) == QQ
    assert unify(QQ, ALG) == ALG
    assert unify(QQ, RR) == RR
    assert unify(QQ, CC) == CC
    assert unify(QQ, ZZ[x]) == QQ[x]
    assert unify(QQ, ZZ.frac_field(x)) == QQ.frac_field(x)
    assert unify(QQ, EX) == EX

    assert unify(RR, F3) == RR
    assert unify(RR, ZZ) == RR
    assert unify(RR, QQ) == RR
    assert unify(RR, ALG) == RR
    assert unify(RR, RR) == RR
    assert unify(RR, CC) == CC
    assert unify(RR, ZZ[x]) == RR[x]
    assert unify(RR, ZZ.frac_field(x)) == RR.frac_field(x)
    assert unify(RR, EX) == EX

    assert unify(CC, F3) == CC
    assert unify(CC, ZZ) == CC
    assert unify(CC, QQ) == CC
    assert unify(CC, ALG) == CC
    assert unify(CC, RR) == CC
    assert unify(CC, CC) == CC
    assert unify(CC, ZZ[x]) == CC[x]
    assert unify(CC, ZZ.frac_field(x)) == CC.frac_field(x)
    assert unify(CC, EX) == EX

    CC2 = ComplexField(prec=20)
    assert unify(CC, CC2) == unify(CC2, CC) == ComplexField(prec=CC.precision,
                                                            tol=CC2.tolerance)
    RR2 = RealField(prec=20)
    assert unify(RR, RR2) == unify(RR2, RR) == RealField(prec=RR.precision,
                                                         tol=RR2.tolerance)

    assert unify(ZZ[x], F3) == F3[x]
    assert unify(ZZ[x], ZZ) == ZZ[x]
    assert unify(ZZ[x], QQ) == QQ[x]
    assert unify(ZZ[x], ALG) == ALG[x]
    assert unify(ZZ[x], RR) == RR[x]
    assert unify(ZZ[x], CC) == CC[x]
    assert unify(ZZ[x], ZZ[x]) == ZZ[x]
    assert unify(ZZ[x], ZZ.frac_field(x)) == ZZ.frac_field(x)
    assert unify(ZZ[x], EX) == EX

    assert unify(ZZ.frac_field(x), F3) == F3.frac_field(x)
    assert unify(ZZ.frac_field(x), ZZ) == ZZ.frac_field(x)
    assert unify(ZZ.frac_field(x), QQ) == QQ.frac_field(x)
    assert unify(ZZ.frac_field(x), ALG) == ALG.frac_field(x)
    assert unify(ZZ.frac_field(x), RR) == RR.frac_field(x)
    assert unify(ZZ.frac_field(x), CC) == CC.frac_field(x)
    assert unify(ZZ.frac_field(x), ZZ[x]) == ZZ.frac_field(x)
    assert unify(ZZ.frac_field(x), ZZ.frac_field(x)) == ZZ.frac_field(x)
    assert unify(ZZ.frac_field(x), EX) == EX

    assert unify(EX, F3) == EX
    assert unify(EX, ZZ) == EX
    assert unify(EX, QQ) == EX
    assert unify(EX, ALG) == EX
    assert unify(EX, RR) == EX
    assert unify(EX, CC) == EX
    assert unify(EX, ZZ[x]) == EX
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

    assert sqrt5.unify(sqrt5[x, y]) == sqrt5[x, y]
    assert sqrt5[x, y].unify(sqrt5) == sqrt5[x, y]

    assert sqrt5.unify(sqrt5.frac_field(x, y)) == sqrt5.frac_field(x, y)
    assert sqrt5.frac_field(x, y).unify(sqrt5) == sqrt5.frac_field(x, y)

    assert sqrt5.unify(sqrt7[x, y]) == sqrt57[x, y]
    assert sqrt5[x, y].unify(sqrt7) == sqrt57[x, y]

    assert sqrt5.unify(sqrt7.frac_field(x, y)) == sqrt57.frac_field(x, y)
    assert sqrt5.frac_field(x, y).unify(sqrt7) == sqrt57.frac_field(x, y)

    sqrt2 = QQ.algebraic_field(sqrt(2))
    r = RootOf(x**7 - x + 1, 0)
    rootof = QQ.algebraic_field(r)
    ans = QQ.algebraic_field(r + sqrt(2))
    assert sqrt2.unify(rootof) == rootof.unify(sqrt2) == ans

    # here domain created from tuple, not Expr
    p = Poly(x**3 - sqrt(2)*x - 1, x, extension=True)
    sqrt2 = p.domain
    assert sqrt2.unify(rootof) == rootof.unify(sqrt2) == ans


def test_Domain_unify_with_symbols():
    pytest.raises(UnificationFailed, lambda: ZZ[x, y].unify_with_symbols(ZZ, (y, z)))
    pytest.raises(UnificationFailed, lambda: ZZ.unify_with_symbols(ZZ[x, y], (y, z)))


def test_Domain__contains__():
    assert (0 in EX) is True
    assert (0 in ZZ) is True
    assert (0 in QQ) is True
    assert (0 in RR) is True
    assert (0 in CC) is True
    assert (0 in ALG) is True
    assert (0 in ZZ[x, y]) is True
    assert (0 in QQ[x, y]) is True
    assert (0 in RR[x, y]) is True

    assert (-7 in EX) is True
    assert (-7 in ZZ) is True
    assert (-7 in QQ) is True
    assert (-7 in RR) is True
    assert (-7 in CC) is True
    assert (-7 in ALG) is True
    assert (-7 in ZZ[x, y]) is True
    assert (-7 in QQ[x, y]) is True
    assert (-7 in RR[x, y]) is True

    assert (17 in EX) is True
    assert (17 in ZZ) is True
    assert (17 in QQ) is True
    assert (17 in RR) is True
    assert (17 in CC) is True
    assert (17 in ALG) is True
    assert (17 in ZZ[x, y]) is True
    assert (17 in QQ[x, y]) is True
    assert (17 in RR[x, y]) is True

    assert (-Rational(1, 7) in EX) is True
    assert (-Rational(1, 7) in ZZ) is False
    assert (-Rational(1, 7) in QQ) is True
    assert (-Rational(1, 7) in RR) is True
    assert (-Rational(1, 7) in CC) is True
    assert (-Rational(1, 7) in ALG) is True
    assert (-Rational(1, 7) in ZZ[x, y]) is False
    assert (-Rational(1, 7) in QQ[x, y]) is True
    assert (-Rational(1, 7) in RR[x, y]) is True

    assert (Rational(3, 5) in EX) is True
    assert (Rational(3, 5) in ZZ) is False
    assert (Rational(3, 5) in QQ) is True
    assert (Rational(3, 5) in RR) is True
    assert (Rational(3, 5) in CC) is True
    assert (Rational(3, 5) in ALG) is True
    assert (Rational(3, 5) in ZZ[x, y]) is False
    assert (Rational(3, 5) in QQ[x, y]) is True
    assert (Rational(3, 5) in RR[x, y]) is True

    assert (3.0 in EX) is True
    assert (3.0 in ZZ) is True
    assert (3.0 in QQ) is True
    assert (3.0 in RR) is True
    assert (3.0 in CC) is True
    assert (3.0 in ALG) is True
    assert (3.0 in ZZ[x, y]) is True
    assert (3.0 in QQ[x, y]) is True
    assert (3.0 in RR[x, y]) is True

    assert (3.14 in EX) is True
    assert (3.14 in ZZ) is False
    assert (3.14 in QQ) is True
    assert (3.14 in RR) is True
    assert (3.14 in CC) is True
    assert (3.14 in ALG) is True
    assert (3.14 in ZZ[x, y]) is False
    assert (3.14 in QQ[x, y]) is True
    assert (3.14 in RR[x, y]) is True

    assert (oo in EX) is True
    assert (oo in ZZ) is False
    assert (oo in QQ) is False
    assert (oo in RR) is True
    assert (oo in CC) is True
    assert (oo in ALG) is False
    assert (oo in ZZ[x, y]) is False
    assert (oo in QQ[x, y]) is False
    assert (oo in RR[x, y]) is True

    assert (-oo in EX) is True
    assert (-oo in ZZ) is False
    assert (-oo in QQ) is False
    assert (-oo in RR) is True
    assert (-oo in CC) is True
    assert (-oo in ALG) is False
    assert (-oo in ZZ[x, y]) is False
    assert (-oo in QQ[x, y]) is False
    assert (-oo in RR[x, y]) is True

    assert (sqrt(7) in EX) is True
    assert (sqrt(7) in ZZ) is False
    assert (sqrt(7) in QQ) is False
    assert (sqrt(7) in RR) is True
    assert (sqrt(7) in CC) is True
    assert (sqrt(7) in ALG) is False
    assert (sqrt(7) in ZZ[x, y]) is False
    assert (sqrt(7) in QQ[x, y]) is False
    assert (sqrt(7) in RR[x, y]) is True

    assert (2*sqrt(3) + 1 in EX) is True
    assert (2*sqrt(3) + 1 in ZZ) is False
    assert (2*sqrt(3) + 1 in QQ) is False
    assert (2*sqrt(3) + 1 in RR) is True
    assert (2*sqrt(3) + 1 in CC) is True
    assert (2*sqrt(3) + 1 in ALG) is True
    assert (2*sqrt(3) + 1 in ZZ[x, y]) is False
    assert (2*sqrt(3) + 1 in QQ[x, y]) is False
    assert (2*sqrt(3) + 1 in RR[x, y]) is True

    assert (sin(1) in EX) is True
    assert (sin(1) in ZZ) is False
    assert (sin(1) in QQ) is False
    assert (sin(1) in RR) is True
    assert (sin(1) in CC) is True
    assert (sin(1) in ALG) is False
    assert (sin(1) in ZZ[x, y]) is False
    assert (sin(1) in QQ[x, y]) is False
    assert (sin(1) in RR[x, y]) is True

    assert (x**2 + 1 in EX) is True
    assert (x**2 + 1 in ZZ) is False
    assert (x**2 + 1 in QQ) is False
    assert (x**2 + 1 in RR) is False
    assert (x**2 + 1 in CC) is False
    assert (x**2 + 1 in ALG) is False
    assert (x**2 + 1 in ZZ[x]) is True
    assert (x**2 + 1 in QQ[x]) is True
    assert (x**2 + 1 in RR[x]) is True
    assert (x**2 + 1 in ZZ[x, y]) is True
    assert (x**2 + 1 in QQ[x, y]) is True
    assert (x**2 + 1 in RR[x, y]) is True

    assert (x**2 + y**2 in EX) is True
    assert (x**2 + y**2 in ZZ) is False
    assert (x**2 + y**2 in QQ) is False
    assert (x**2 + y**2 in RR) is False
    assert (x**2 + y**2 in CC) is False
    assert (x**2 + y**2 in ALG) is False
    assert (x**2 + y**2 in ZZ[x]) is False
    assert (x**2 + y**2 in QQ[x]) is False
    assert (x**2 + y**2 in RR[x]) is False
    assert (x**2 + y**2 in ZZ[x, y]) is True
    assert (x**2 + y**2 in QQ[x, y]) is True
    assert (x**2 + y**2 in RR[x, y]) is True

    assert (Rational(3, 2)*x/(y + 1) - z in QQ[x, y, z]) is False


def test_Domain_get_ring():
    assert ZZ.has_assoc_Ring is True
    assert QQ.has_assoc_Ring is True
    assert ZZ[x].has_assoc_Ring is True
    assert QQ[x].has_assoc_Ring is True
    assert ZZ[x, y].has_assoc_Ring is True
    assert QQ[x, y].has_assoc_Ring is True
    assert ZZ.frac_field(x).has_assoc_Ring is True
    assert QQ.frac_field(x).has_assoc_Ring is True
    assert ZZ.frac_field(x, y).has_assoc_Ring is True
    assert QQ.frac_field(x, y).has_assoc_Ring is True

    assert EX.has_assoc_Ring is False
    assert RR.has_assoc_Ring is False
    assert ALG.has_assoc_Ring is False

    assert ZZ.get_ring() == ZZ
    assert QQ.get_ring() == ZZ
    assert ZZ[x].get_ring() == ZZ[x]
    assert QQ[x].get_ring() == QQ[x]
    assert ZZ[x, y].get_ring() == ZZ[x, y]
    assert QQ[x, y].get_ring() == QQ[x, y]
    assert ZZ.frac_field(x).get_ring() == ZZ[x]
    assert QQ.frac_field(x).get_ring() == QQ[x]
    assert ZZ.frac_field(x, y).get_ring() == ZZ[x, y]
    assert QQ.frac_field(x, y).get_ring() == QQ[x, y]

    assert EX.get_ring() == EX

    pytest.raises(DomainError, lambda: RR.get_ring())
    pytest.raises(DomainError, lambda: ALG.get_ring())


def test_Domain_get_field():
    assert EX.has_assoc_Field is True
    assert ZZ.has_assoc_Field is True
    assert QQ.has_assoc_Field is True
    assert RR.has_assoc_Field is True
    assert ALG.has_assoc_Field is True
    assert ZZ[x].has_assoc_Field is True
    assert QQ[x].has_assoc_Field is True
    assert ZZ[x, y].has_assoc_Field is True
    assert QQ[x, y].has_assoc_Field is True

    assert EX.get_field() == EX
    assert ZZ.get_field() == QQ
    assert QQ.get_field() == QQ
    assert RR.get_field() == RR
    assert ALG.get_field() == ALG
    assert ZZ[x].get_field() == ZZ.frac_field(x)
    assert QQ[x].get_field() == QQ.frac_field(x)
    assert ZZ[x, y].get_field() == ZZ.frac_field(x, y)
    assert QQ[x, y].get_field() == QQ.frac_field(x, y)


def test_Domain_get_exact():
    assert EX.get_exact() == EX
    assert ZZ.get_exact() == ZZ
    assert QQ.get_exact() == QQ
    assert RR.get_exact() == QQ
    assert ALG.get_exact() == ALG
    assert ZZ[x].get_exact() == ZZ[x]
    assert QQ[x].get_exact() == QQ[x]
    assert ZZ[x, y].get_exact() == ZZ[x, y]
    assert QQ[x, y].get_exact() == QQ[x, y]
    assert ZZ.frac_field(x).get_exact() == ZZ.frac_field(x)
    assert QQ.frac_field(x).get_exact() == QQ.frac_field(x)
    assert ZZ.frac_field(x, y).get_exact() == ZZ.frac_field(x, y)
    assert QQ.frac_field(x, y).get_exact() == QQ.frac_field(x, y)


def test_Domain_convert():
    assert QQ.convert(10e-52) == QQ(1684996666696915, 1684996666696914987166688442938726917102321526408785780068975640576)

    R, x = ring("x", ZZ)
    assert ZZ.convert(x - x) == 0
    assert ZZ.convert(x - x, R.to_domain()) == 0

    F3 = FF(3)
    assert F3.convert(Float(2.0)) == F3.dtype(2)
    assert F3.convert(PythonRational(2, 1)) == F3.dtype(2)
    pytest.raises(CoercionFailed, lambda: F3.convert(PythonRational(1, 2)))
    assert F3.convert(2.0) == F3.dtype(2)
    pytest.raises(CoercionFailed, lambda: F3.convert(2.1))

    assert RR.convert(CC(1)) == RR(1)
    pytest.raises(CoercionFailed, lambda: RR.convert(CC(1, 2)))

    assert QQ.convert(ALG.new(1), ALG) == QQ(1)
    pytest.raises(CoercionFailed, lambda: QQ.convert(ALG.new([1, 1]), ALG))

    assert ZZ.convert(ALG.new(1), ALG) == ZZ(1)
    pytest.raises(CoercionFailed, lambda: ZZ.convert(ALG.new([1, 1]), ALG))

    assert EX.convert(ALG.new([1, 1]), ALG) == sqrt(2) + sqrt(3) + 1

    ALG2 = QQ.algebraic_field(sqrt(2))
    a2 = ALG2.from_diofant(sqrt(2))
    a = ALG.convert(a2, ALG2)
    assert a.rep == [QQ(1, 2), 0, -QQ(9, 2), 0]

    assert ZZ_python().convert(3.0) == ZZ_python().dtype(3)
    pytest.raises(CoercionFailed, lambda: ZZ_python().convert(3.2))

    assert CC.convert(QQ_python()(1, 2)) == CC(0.5)
    CC01 = ComplexField(tol=0.1)
    assert CC.convert(CC01(0.3)) == CC(0.3)

    assert RR.convert(complex(2 + 0j)) == RR(2)
    pytest.raises(CoercionFailed, lambda: RR.convert(complex(2 + 3j)))


def test_arithmetics():
    assert ZZ.rem(ZZ(2), ZZ(3)) == 2
    assert ZZ.div(ZZ(2), ZZ(3)) == (0, 2)
    assert QQ.rem(QQ(2, 3), QQ(4, 7)) == 0
    assert QQ.div(QQ(2, 3), QQ(4, 7)) == (QQ(7, 6), 0)

    QQp = QQ_python()
    assert QQp.factorial(QQp(7, 2)) == 6

    assert CC.gcd(CC(1), CC(2)) == 1
    assert CC.lcm(CC(1), CC(2)) == 2

    assert EX(Rational(2, 3)).numer() == 2
    assert EX(Rational(2, 3)).denom() == 3

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
    assert ZZ.numer(ZZ(3)) == 3
    assert ZZ.denom(ZZ(3)) == 1


def test_PolynomialRing__init():
    pytest.raises(GeneratorsNeeded, lambda: ZZ.poly_ring())


def test_FractionField__init():
    pytest.raises(GeneratorsNeeded, lambda: ZZ.frac_field())


def test_inject():
    assert ZZ.inject(x, y, z) == ZZ[x, y, z]
    assert ZZ[x].inject(y, z) == ZZ[x, y, z]
    assert ZZ.frac_field(x).inject(y, z) == ZZ.frac_field(x, y, z)
    pytest.raises(GeneratorsError, lambda: ZZ[x].inject(x))


def test_Domain_map():
    seq = ZZ.map([1, 2, 3, 4])

    assert all(ZZ.of_type(elt) for elt in seq)

    seq = ZZ.map([[1, 2, 3, 4]])

    assert all(ZZ.of_type(elt) for elt in seq[0]) and len(seq) == 1


def test_Domain___eq__():
    assert (ZZ[x, y] == ZZ[x, y]) is True
    assert (QQ[x, y] == QQ[x, y]) is True

    assert (ZZ[x, y] == QQ[x, y]) is False
    assert (QQ[x, y] == ZZ[x, y]) is False

    assert (ZZ.frac_field(x, y) == ZZ.frac_field(x, y)) is True
    assert (QQ.frac_field(x, y) == QQ.frac_field(x, y)) is True

    assert (ZZ.frac_field(x, y) == QQ.frac_field(x, y)) is False
    assert (QQ.frac_field(x, y) == ZZ.frac_field(x, y)) is False


def test_Domain__algebraic_field():
    alg = ZZ.algebraic_field(sqrt(2))
    assert alg.ext.minpoly == Poly(x**2 - 2)
    assert alg.domain == QQ

    alg = QQ.algebraic_field(sqrt(2))
    assert alg.ext.minpoly == Poly(x**2 - 2)
    assert alg.domain == QQ

    alg = alg.algebraic_field(sqrt(3))
    assert alg.ext.minpoly == Poly(x**4 - 10*x**2 + 1)
    assert alg.domain == QQ

    assert alg.is_nonpositive(alg([-1, 1])) is True
    assert alg.is_nonnegative(alg([2, -1])) is True

    assert alg.numer(alg(1)) == alg(1)

    pytest.raises(DomainError, lambda: AlgebraicField(ZZ, sqrt(2)))

    assert alg.characteristic() == 0


def test_PolynomialRing_from_FractionField():
    F,  x, y = field("x,y", ZZ)
    R,  X, Y = ring("x,y", ZZ)

    f = (x**2 + y**2)/(x + 1)
    g = (x**2 + y**2)/4
    h = x**2 + y**2

    assert R.to_domain().from_FractionField(f, F.to_domain()) is None
    assert R.to_domain().from_FractionField(g, F.to_domain()) == X**2/4 + Y**2/4
    assert R.to_domain().from_FractionField(h, F.to_domain()) == X**2 + Y**2

    F,  x, y = field("x,y", QQ)
    R,  X, Y = ring("x,y", QQ)

    f = (x**2 + y**2)/(x + 1)
    g = (x**2 + y**2)/4
    h = x**2 + y**2

    assert R.to_domain().from_FractionField(f, F.to_domain()) is None
    assert R.to_domain().from_FractionField(g, F.to_domain()) == X**2/4 + Y**2/4
    assert R.to_domain().from_FractionField(h, F.to_domain()) == X**2 + Y**2


def test_FractionField_from_PolynomialRing():
    R,  x, y = ring("x,y", QQ)
    F,  X, Y = field("x,y", ZZ)

    f = 3*x**2 + 5*y**2
    g = x**2/3 + y**2/5

    assert F.to_domain().from_PolynomialRing(f, R.to_domain()) == 3*X**2 + 5*Y**2
    assert F.to_domain().from_PolynomialRing(g, R.to_domain()) == (5*X**2 + 3*Y**2)/15

    RALG,  u, v = ring("u,v", ALG)
    pytest.raises(CoercionFailed,
                  lambda: F.to_domain().convert(3*u**2 + 5*sqrt(2)*v**2))


def test_FractionField_convert():
    F,  X, Y = field("x,y", QQ)
    F.to_domain().convert(QQ_python()(1, 3)) == F.one/3


def test_FF_of_type():
    assert FF(3).of_type(FF(3)(1)) is True
    assert FF(5).of_type(FF(5)(3)) is True
    assert FF(5).of_type(FF(7)(3)) is False


def test___eq__():
    assert not QQ[x] == ZZ[x]
    assert not QQ.frac_field(x) == ZZ.frac_field(x)

    assert EX(1) != EX(2)

    F11 = FF(11)
    assert F11(2) != F11(3)
    assert F11(2) != object()


def test_RealField_from_diofant():
    assert RR.convert(Integer(0)) == RR.dtype(0)
    assert RR.convert(Float(0.0)) == RR.dtype(0.0)
    assert RR.convert(Integer(1)) == RR.dtype(1)
    assert RR.convert(Float(1.0)) == RR.dtype(1.0)
    assert RR.convert(sin(1)) == RR.dtype(sin(1).evalf())
    assert RR.convert(oo) == RR("+inf")
    assert RR.convert(-oo) == RR("-inf")
    pytest.raises(CoercionFailed, lambda: RR.convert(x))


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

    assert (F5(1) < F5(2)) is True
    assert (F5(1) <= F5(2)) is True
    assert (F5(1) > F5(2)) is False
    assert (F5(1) >= F5(2)) is False

    assert (F5(3) < F5(2)) is False
    assert (F5(3) <= F5(2)) is False
    assert (F5(3) > F5(2)) is True
    assert (F5(3) >= F5(2)) is True

    assert (F5(1) < F5(7)) is True
    assert (F5(1) <= F5(7)) is True
    assert (F5(1) > F5(7)) is False
    assert (F5(1) >= F5(7)) is False

    assert (F5(3) < F5(7)) is False
    assert (F5(3) <= F5(7)) is False
    assert (F5(3) > F5(7)) is True
    assert (F5(3) >= F5(7)) is True

    assert (F5(1) < 2) is True
    assert (F5(1) <= 2) is True
    assert (F5(1) > 2) is False
    assert (F5(1) >= 2) is False

    assert (F5(3) < 2) is False
    assert (F5(3) <= 2) is False
    assert (F5(3) > 2) is True
    assert (F5(3) >= 2) is True

    assert (F5(1) < 7) is True
    assert (F5(1) <= 7) is True
    assert (F5(1) > 7) is False
    assert (F5(1) >= 7) is False

    assert (F5(3) < 7) is False
    assert (F5(3) <= 7) is False
    assert (F5(3) > 7) is True
    assert (F5(3) >= 7) is True

    pytest.raises(NotInvertible, lambda: F5(0)**(-1))
    pytest.raises(NotInvertible, lambda: F5(5)**(-1))

    pytest.raises(ValueError, lambda: FF(0))
    pytest.raises(ValueError, lambda: FF(2.1))


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
    f1 = Float("1.01")
    f2 = Float("1.0000000000000000000001")
    assert f1._prec == 53
    assert f2._prec == 80
    assert RR(f1)-1 > 1e-50
    assert RR(f2)-1 < 1e-50  # RR's precision is lower than f2's

    RR2 = RealField(prec=f2._prec)
    assert RR2(f1)-1 > 1e-50
    assert RR2(f2)-1 > 1e-50  # RR's precision is equal to f2's


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


def test_almosteq():
    assert CC.almosteq(CC(2), 3) is False
    assert CC.almosteq(2, CC(3)) is False
    assert CC.almosteq(2, CC(2.5), 0.1) is False
    assert CC.almosteq(2, CC(2.5), 1.0) is True

    assert RR.almosteq(5, RR(2), 1) is True
    assert RR._context.almosteq(RR(2), 1, None, 1) is True


def test_to_diofant():
    assert CC.to_diofant(1 - 2j) == 1 - 2*I


def test_EX():
    assert EX.is_positive(EX(2))
    assert EX.is_nonnegative(EX(2))
    assert EX.is_negative(EX(-1))
    assert EX.is_nonpositive(EX(-1))

    assert EX.numer(EX(1)/2) == 1
    assert EX.denom(EX(1)/2) == 2


def test_sympyissue_13545():
    assert Poly(x + 1, x, modulus=2) + 1 == Poly(x, x, modulus=2)
    pytest.raises(NotImplementedError,
                  lambda: Poly(x, modulus=2) + Poly(x, modulus=3))
