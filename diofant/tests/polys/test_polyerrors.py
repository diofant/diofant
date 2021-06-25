from diofant import (EX, RR, ZZ, ComputationFailed, ExactQuotientFailed,
                     OperationNotSupported, PolificationFailed,
                     PolynomialDivisionFailed)
from diofant.abc import x


__all__ = ()


def test_printing():
    f, g = [[]]*2
    e = PolynomialDivisionFailed(f, g, EX)
    assert str(e)[:57] == ("couldn't reduce degree in a polynomial "
                           'division algorithm')
    assert str(e)[-140:][:57] == ('You may want to use a different '
                                  'simplification algorithm.')

    f, g = [[]]*2
    e = PolynomialDivisionFailed(f, g, RR)
    assert str(e)[-139:][:74] == ('Your working precision or tolerance of '
                                  'computations may be set improperly.')

    f, g = [[]]*2
    e = PolynomialDivisionFailed(f, g, ZZ)
    assert str(e)[-168:][:80] == ('Zero detection is guaranteed in this '
                                  'coefficient domain. This may indicate a bug')

    e = OperationNotSupported(x.as_poly(), 'spam')
    assert str(e).find('spam') >= 0
    assert str(e).find('operation not supported') >= 0

    exc = PolificationFailed(1, x, x**2)
    assert str(exc).find("can't construct a polynomial from x") >= 0
    exc = PolificationFailed(1, [x], [x**2], True)
    assert str(exc).find("can't construct polynomials from x") >= 0

    e = ComputationFailed('LT', 1, exc)
    assert str(e).find('failed without generators') >= 0
    assert str(e).find('x**2') >= 0

    e = ExactQuotientFailed(x.as_poly(), (x**2).as_poly())
    assert str(e).find('does not divide') >= 0
    assert str(e).find('x**2') >= 0
    assert str(e).find('in ZZ') < 0
    e = ExactQuotientFailed(x.as_poly(), (x**2).as_poly(), ZZ)
    assert str(e).find('in ZZ') >= 0
