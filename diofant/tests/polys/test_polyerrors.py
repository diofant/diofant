from diofant.abc import x
from diofant.domains import EX, RR, ZZ
from diofant.polys import Poly
from diofant.polys.densebasic import dmp_normal
from diofant.polys.polyerrors import (ComputationFailed, ExactQuotientFailed,
                                      OperationNotSupported,
                                      PolificationFailed,
                                      PolynomialDivisionFailed)


__all__ = ()


def test_printing():
    f, g = [dmp_normal([], 0, EX)]*2
    e = PolynomialDivisionFailed(f, g, EX)
    assert str(e)[:57] == ("couldn't reduce degree in a polynomial "
                           "division algorithm")
    assert str(e)[-140:][:57] == ("You may want to use a different "
                                  "simplification algorithm.")

    f, g = [dmp_normal([], 0, RR)]*2
    e = PolynomialDivisionFailed(f, g, RR)
    assert str(e)[-139:][:74] == ("Your working precision or tolerance of "
                                  "computations may be set improperly.")

    f, g = [dmp_normal([], 0, ZZ)]*2
    e = PolynomialDivisionFailed(f, g, ZZ)
    assert str(e)[-168:][:80] == ("Zero detection is guaranteed in this "
                                  "coefficient domain. This may indicate a bug")

    e = OperationNotSupported(Poly(x), 'spam')
    assert str(e).find('spam') >= 0
    assert str(e).find('operation not supported') >= 0

    exc = PolificationFailed(1, x, x**2)
    assert str(exc).find("can't construct a polynomial from x") >= 0
    exc = PolificationFailed(1, [x], [x**2], True)
    assert str(exc).find("can't construct polynomials from x") >= 0

    e = ComputationFailed('LT', 1, exc)
    assert str(e).find('failed without generators') >= 0
    assert str(e).find('x**2') >= 0

    e = ExactQuotientFailed(Poly(x), Poly(x**2))
    assert str(e).find('does not divide') >= 0
    assert str(e).find('x**2') >= 0
    assert str(e).find('in ZZ') < 0
    e = ExactQuotientFailed(Poly(x), Poly(x**2), ZZ)
    assert str(e).find('in ZZ') >= 0
