"""Tests for tools for constructing domains for expressions. """

from diofant import E, Float, GoldenRatio, Integer, Rational, sin, sqrt
from diofant.abc import x, y
from diofant.domains import EX, QQ, RR, ZZ
from diofant.domains.realfield import RealField
from diofant.polys.constructor import construct_domain


__all__ = ()


def test_construct_domain():
    assert construct_domain([1, 2, 3]) == (ZZ, [ZZ(1), ZZ(2), ZZ(3)])
    assert construct_domain([1, 2, 3], field=True) == (QQ, [QQ(1), QQ(2), QQ(3)])

    assert construct_domain([Integer(1), Integer(2), Integer(3)]) == (ZZ, [ZZ(1), ZZ(2), ZZ(3)])
    assert construct_domain([Integer(1), Integer(2), Integer(3)], field=True) == (QQ, [QQ(1), QQ(2), QQ(3)])

    assert construct_domain([Rational(1, 2), Integer(2)]) == (QQ, [QQ(1, 2), QQ(2)])
    result = construct_domain([3.14, 1, Rational(1, 2)])
    assert isinstance(result[0], RealField)
    assert result[1] == [RR(3.14), RR(1.0), RR(0.5)]

    assert construct_domain([3.14, sqrt(2)], extension=None) == (EX, [EX(3.14), EX(sqrt(2))])
    assert construct_domain([3.14, sqrt(2)], extension=True) == (EX, [EX(3.14), EX(sqrt(2))])
    assert construct_domain([sqrt(2), 3.14], extension=True) == (EX, [EX(sqrt(2)), EX(3.14)])

    assert construct_domain([1, sqrt(2)], extension=None) == (EX, [EX(1), EX(sqrt(2))])

    assert construct_domain([x, sqrt(x)]) == (EX, [EX(x), EX(sqrt(x))])
    assert construct_domain([x, sqrt(x), sqrt(y)]) == (EX, [EX(x), EX(sqrt(x)), EX(sqrt(y))])

    alg = QQ.algebraic_field(sqrt(2))

    assert construct_domain([7, Rational(1, 2), sqrt(2)], extension=True) == \
        (alg, [alg.convert(7), alg.convert(Rational(1, 2)), alg.convert(sqrt(2))])

    alg = QQ.algebraic_field(sqrt(2) + sqrt(3))

    assert construct_domain([7, sqrt(2), sqrt(3)], extension=True) == \
        (alg, [alg.convert(7), alg.convert(sqrt(2)), alg.convert(sqrt(3))])

    dom = ZZ[x]

    assert construct_domain([2*x, 3]) == \
        (dom, [dom.convert(2*x), dom.convert(3)])

    dom = ZZ[x, y]

    assert construct_domain([2*x, 3*y]) == \
        (dom, [dom.convert(2*x), dom.convert(3*y)])

    dom = QQ[x]

    assert construct_domain([x/2, 3]) == \
        (dom, [dom.convert(x/2), dom.convert(3)])

    dom = QQ[x, y]

    assert construct_domain([x/2, 3*y]) == \
        (dom, [dom.convert(x/2), dom.convert(3*y)])

    dom = RR[x]

    assert construct_domain([x/2, 3.5]) == \
        (dom, [dom.convert(x/2), dom.convert(3.5)])

    dom = RR[x, y]

    assert construct_domain([x/2, 3.5*y]) == \
        (dom, [dom.convert(x/2), dom.convert(3.5*y)])

    dom = ZZ.frac_field(x)

    assert construct_domain([2/x, 3]) == \
        (dom, [dom.convert(2/x), dom.convert(3)])

    dom = ZZ.frac_field(x, y)

    assert construct_domain([2/x, 3*y]) == \
        (dom, [dom.convert(2/x), dom.convert(3*y)])

    dom = RR.frac_field(x)

    assert construct_domain([2/x, 3.5]) == \
        (dom, [dom.convert(2/x), dom.convert(3.5)])

    dom = RR.frac_field(x, y)

    assert construct_domain([2/x, 3.5*y]) == \
        (dom, [dom.convert(2/x), dom.convert(3.5*y)])

    assert construct_domain(2) == (ZZ, ZZ(2))
    assert construct_domain(Rational(2, 3)) == (QQ, QQ(2, 3))

    assert construct_domain({}) == (ZZ, {})


def test_composite_option():
    assert construct_domain({(1,): sin(y)}, composite=False) == \
        (EX, {(1,): EX(sin(y))})

    assert construct_domain({(1,): y}, composite=False) == \
        (EX, {(1,): EX(y)})

    assert construct_domain({(1, 1): 1}, composite=False) == \
        (ZZ, {(1, 1): 1})

    assert construct_domain({(1, 0): y}, composite=False) == \
        (EX, {(1, 0): EX(y)})


def test_precision():
    f1 = Float("1.01")
    f2 = Float("1.0000000000000000000001")
    for x in [1, 1e-2, 1e-6, 1e-13, 1e-14, 1e-16, 1e-20, 1e-100, 1e-300,
              f1, f2]:
        result = construct_domain([x])
        y = float(result[1][0])
        assert abs(x - y)/x < 1e-14  # Test relative accuracy
        assert abs(x - y)/x < 1e-14  # Test relative accuracy

    result = construct_domain([f1])
    y = result[1][0]
    assert y-1 > 1e-50

    result = construct_domain([f2])
    y = result[1][0]
    assert y-1 > 1e-50


def test_sympyissue_11538():
    assert construct_domain(E)[0] == ZZ[E]

    assert (construct_domain(x**2 + 2*x + E) ==
            (ZZ[x, E], ZZ[x, E].convert(x**2 + 2*x + E)))

    assert (construct_domain(x + y + GoldenRatio) ==
            (EX, EX.convert(x + y + GoldenRatio)))
