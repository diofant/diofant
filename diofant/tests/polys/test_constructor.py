"""Tests for tools for constructing domains for expressions."""

from diofant import (EX, QQ, RR, ZZ, E, Float, GoldenRatio, I, Integer, Poly,
                     Rational, construct_domain, exp, pi, sin, sqrt)
from diofant.abc import x, y


__all__ = ()


def test_construct_domain():
    assert construct_domain([1, 2, 3]) == (ZZ, [ZZ(1), ZZ(2), ZZ(3)])
    assert construct_domain([1, 2, 3], field=True) == (QQ, [QQ(1), QQ(2), QQ(3)])

    assert construct_domain([Integer(1), Integer(2), Integer(3)]) == (ZZ, [ZZ(1), ZZ(2), ZZ(3)])
    assert construct_domain([Integer(1), Integer(2), Integer(3)], field=True) == (QQ, [QQ(1), QQ(2), QQ(3)])

    assert construct_domain([Rational(1, 2), Integer(2)]) == (QQ, [QQ(1, 2), QQ(2)])
    assert construct_domain([3.14, 1, Rational(1, 2)]) == (RR, [RR(3.14), RR(1.0), RR(0.5)])

    assert construct_domain([3.14, sqrt(2)], extension=False) == (EX, [EX(3.14), EX(sqrt(2))])
    assert construct_domain([3.14, sqrt(2)]) == (EX, [EX(3.14), EX(sqrt(2))])
    assert construct_domain([sqrt(2), 3.14]) == (EX, [EX(sqrt(2)), EX(3.14)])

    assert construct_domain([1, sqrt(2)], extension=False) == (EX, [EX(1), EX(sqrt(2))])

    assert construct_domain([x, sqrt(x)]) == (EX, [EX(x), EX(sqrt(x))])
    assert construct_domain([x, sqrt(x), sqrt(y)]) == (EX, [EX(x), EX(sqrt(x)), EX(sqrt(y))])

    alg = QQ.algebraic_field(sqrt(2))

    assert (construct_domain([7, Rational(1, 2), sqrt(2)]) ==
            (alg, [alg([7]), alg([Rational(1, 2)]), alg([0, 1])]))

    alg = QQ.algebraic_field(sqrt(2) + sqrt(3))

    assert (construct_domain([7, sqrt(2), sqrt(3)]) ==
            (alg, [alg([7]), alg.from_expr(sqrt(2)), alg.from_expr(sqrt(3))]))

    dom = ZZ.inject(x)

    assert construct_domain([2*x, 3]) == (dom, [dom(2*x), dom(3)])

    dom = ZZ.inject(x, y)

    assert construct_domain([2*x, 3*y]) == (dom, [dom(2*x), dom(3*y)])

    dom = QQ.inject(x)

    assert construct_domain([x/2, 3]) == (dom, [dom(x/2), dom(3)])

    dom = QQ.inject(x, y)

    assert construct_domain([x/2, 3*y]) == (dom, [dom(x/2), dom(3*y)])

    dom = RR.inject(x)

    assert construct_domain([x/2, 3.5]) == (dom, [dom(x/2), dom(3.5)])

    dom = RR.inject(x, y)

    assert construct_domain([x/2, 3.5*y]) == (dom, [dom(x/2), dom(3.5*y)])

    dom = ZZ.inject(x).field

    assert construct_domain([2/x, 3]) == (dom, [dom(2/x), dom(3)])

    dom = ZZ.inject(x, y).field

    assert construct_domain([2/x, 3*y]) == (dom, [dom(2/x), dom(3*y)])

    dom = RR.inject(x).field

    assert construct_domain([2/x, 3.5]) == (dom, [dom(2/x), dom(3.5)])

    dom = RR.inject(x, y).field

    assert construct_domain([2/x, 3.5*y]) == (dom, [dom(2/x), dom(3.5*y)])

    assert construct_domain(2) == (ZZ, ZZ(2))
    assert construct_domain(Rational(2, 3)) == (QQ, QQ(2, 3))

    assert construct_domain({}) == (ZZ, {})

    assert construct_domain([-x*y + x*(y + 42) -
                             42*x]) == (EX, [EX(-x*y + x*(y + 42) - 42*x)])


def test_composite_option():
    assert construct_domain({(1,): sin(y)}, composite=False) == (EX, {(1,): EX(sin(y))})
    assert construct_domain({(1,): y}, composite=False) == (EX, {(1,): EX(y)})
    assert construct_domain({(1, 1): 1}, composite=False) == (ZZ, {(1, 1): 1})
    assert construct_domain({(1, 0): y}, composite=False) == (EX, {(1, 0): EX(y)})


def test_precision():
    f1 = Float('1.01')
    f2 = Float('1.0000000000000000000001')
    for x in [1, 1e-2, 1e-6, 1e-13, 1e-14, 1e-16, 1e-20, 1e-100, 1e-300,
              f1, f2]:
        result = construct_domain([x])
        y = float(result[1][0])
        assert abs(x - y)/x < 1e-14  # Test relative accuracy

    result = construct_domain([f1])
    assert result[1][0] - 1 > 1e-50

    result = construct_domain([f2])
    assert result[1][0] - 1 > 1e-50


def test_sympyissue_11538():
    assert construct_domain(E)[0] == ZZ.inject(E)
    assert (construct_domain(x**2 + 2*x + E) == (ZZ.inject(x, E), ZZ.inject(x, E)(x**2 + 2*x + E)))
    assert (construct_domain(x + y + GoldenRatio) == (EX, EX(x + y + GoldenRatio)))


def test_sympyissue_5428_14337():
    assert (x**2 + I).as_poly(x).domain == QQ.algebraic_field(I)
    assert (x**2 + sqrt(2)).as_poly(x).domain == QQ.algebraic_field(sqrt(2))


def test_sympyissue_20617():
    dom = QQ.algebraic_field(2*sqrt(3)).algebraic_field(I)

    w = exp(I*2*pi/3)
    l1 = [w**2, w, 1]

    r1 = Poly(l1, x, extension=True)

    w2 = w.expand()
    l2 = [w2**2, w2, 1]

    r2 = Poly(l2, x, extension=True)

    assert r1 == r2
    assert r1.domain == dom
