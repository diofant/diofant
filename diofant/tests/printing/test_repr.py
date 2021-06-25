"""repr() printing tests."""

from __future__ import annotations

import typing

import pytest

from diofant import (FF, QQ, ZZ, Abs, Catalan, Dummy, E, EulerGamma, Float,
                     Function, GoldenRatio, I, ImmutableMatrix, Integer,
                     Matrix, Rational, Symbol, Wild, WildFunction, false,
                     field, grlex, nan, ones, oo, pi, ring, root, sin, sqrt,
                     srepr, true, zoo)
from diofant.abc import x, y
from diofant.core.exprtools import Factors
from diofant.geometry import Ellipse, Point


__all__ = ()


# eval(repr(expr)) == expr has to succeed in the right environment. The right
# environment is the scope of "from diofant import *" for most cases.
ENV: dict[str, typing.Any] = {}
imports = ['from diofant import *',
           'from diofant.domains.integerring import GMPYIntegerRing, PythonIntegerRing',
           'from diofant.domains.rationalfield import GMPYRationalField, PythonRationalField',
           'from diofant.polys.orderings import GradedLexOrder, LexOrder']
exec('\n'.join(imports), ENV)


def sT(expr, string):
    """
    Tests that repr delivers the expected string and that
    the condition eval(repr(expr))==expr holds.
    """
    assert repr(expr) == string
    assert eval(string, ENV) == expr


def test_printmethod():
    class R(Abs):
        def _diofantrepr(self, printer):
            return f'foo({printer._print(self.args[0])})'
    assert repr(R(x)) == "foo(Symbol('x'))"


def test_Add():
    sT(x + y, "Add(Symbol('x'), Symbol('y'))")
    assert srepr(x**2 + 1, order='lex') == ("Add(Pow(Symbol('x'), "
                                            'Integer(2)), Integer(1))')


def test_Function():
    sT(Function('f')(x), "Function('f')(Symbol('x'))")
    # test unapplied Function
    sT(Function('f'), "Function('f')")

    sT(sin(x), "sin(Symbol('x'))")
    sT(sin, 'sin')


def test_Geometry():
    sT(Point(0, 0), 'Point(Integer(0), Integer(0))')
    sT(Ellipse(Point(0, 0), 5, 1),
       'Ellipse(Point(Integer(0), Integer(0)), Integer(5), Integer(1))')
    # TODO more tests


def test_Singletons():
    sT(Catalan, 'Catalan')
    sT(zoo, 'zoo')
    sT(EulerGamma, 'EulerGamma')
    sT(E, 'E')
    sT(GoldenRatio, 'GoldenRatio')
    sT(Rational(1, 2), 'Rational(1, 2)')
    sT(I, 'I')
    sT(oo, 'oo')
    sT(nan, 'nan')
    sT(-oo, '-oo')
    sT(Integer(-1), 'Integer(-1)')
    sT(Integer(1), 'Integer(1)')
    sT(pi, 'pi')
    sT(Integer(0), 'Integer(0)')


def test_Integer():
    sT(Integer(4), 'Integer(4)')


def test_list():
    sT([x, Integer(4)], "[Symbol('x'), Integer(4)]")


def test_Matrix():
    for cls, name in [(Matrix, 'MutableDenseMatrix'), (ImmutableMatrix, 'ImmutableMatrix')]:
        sT(cls([[x**+1, 1], [y, x + y]]),
           f"{name}([[Symbol('x'), Integer(1)], [Symbol('y'), Add(Symbol('x'), Symbol('y'))]])")

        sT(cls(), f'{name}([])')

        sT(cls([[x**+1, 1], [y, x + y]]), f"{name}([[Symbol('x'), Integer(1)], [Symbol('y'), Add(Symbol('x'), Symbol('y'))]])")


def test_empty_Matrix():
    sT(ones(0, 3), 'MutableDenseMatrix(0, 3, [])')
    sT(ones(4, 0), 'MutableDenseMatrix(4, 0, [])')
    sT(ones(0, 0), 'MutableDenseMatrix([])')


def test_Rational():
    sT(Rational(1, 3), 'Rational(1, 3)')
    sT(Rational(-1, 3), 'Rational(-1, 3)')


def test_Factors():
    assert repr(Factors(x*y**2)) == 'Factors({x: 1, y: 2})'


def test_AlgebraicElement():
    K = QQ.algebraic_field(sqrt(2))
    a = K.unit
    sT(a, f'AlgebraicField({QQ!r}, Pow(Integer(2), Rational(1, 2)))([Integer(0), Integer(1)])')
    K = QQ.algebraic_field(root(-2, 3))
    a = K.unit
    sT(a, f'AlgebraicField({QQ!r}, Pow(Integer(-2), Rational(1, 3)))([Integer(0), Integer(1)])')


def test_Float():
    sT(Float('1.23', dps=3), "Float('1.22998', dps=3)")
    sT(Float('1.23456789', dps=9), "Float('1.23456788994', dps=9)")
    sT(Float('1.234567890123456789', dps=19),
       "Float('1.234567890123456789013', dps=19)")
    sT(Float(
        '0.60038617995049726', 15), "Float('0.60038617995049726', dps=15)")


def test_Symbol():
    sT(x, "Symbol('x')")
    sT(y, "Symbol('y')")
    sT(Symbol('x', negative=True), "Symbol('x', negative=True)")


def test_Symbol_two_assumptions():
    x = Symbol('x', negative=0, integer=1)
    # order could vary
    s1 = "Symbol('x', integer=True, negative=False)"
    s2 = "Symbol('x', negative=False, integer=True)"
    assert repr(x) in (s1, s2)
    assert eval(repr(x), ENV) == x


def test_Symbol_no_special_commutative_treatment():
    sT(Symbol('x'), "Symbol('x')")
    sT(Symbol('x', commutative=False), "Symbol('x', commutative=False)")
    sT(Symbol('x', commutative=0), "Symbol('x', commutative=False)")
    sT(Symbol('x', commutative=True), "Symbol('x', commutative=True)")
    sT(Symbol('x', commutative=1), "Symbol('x', commutative=True)")


def test_Wild():
    sT(Wild('x', even=True), "Wild('x', even=True)")


def test_Dummy():
    # cannot use sT here
    d = Dummy('d', nonzero=True)
    assert repr(d) == "Dummy('d', nonzero=True)"


def test_Dummy_from_Symbol():
    # should not get the full dictionary of assumptions
    n = Symbol('n', integer=True)
    d = n.as_dummy()
    assert repr(d) == "Dummy('n', integer=True)"


def test_tuple():
    sT((x,), "(Symbol('x'),)")
    sT((x, y), "(Symbol('x'), Symbol('y'))")


def test_WildFunction():
    sT(WildFunction('w'), "WildFunction('w')")


def test_settings():
    pytest.raises(TypeError, lambda: srepr(x, method='garbage'))


def test_Mul():
    sT(3*x**3*y, "Mul(Integer(3), Pow(Symbol('x'), Integer(3)), Symbol('y'))")


def test_FiniteField():
    sT(FF(2), 'GF(2)')

    F4 = FF(2, [1, 1, 1])
    repr(F4.one)  # not raises


def test_PolynomialRing():
    sT(ZZ.inject('x'), f"UnivarPolynomialRing({ZZ!r}, (Symbol('x'),), LexOrder())")
    sT(QQ.poly_ring('x', 'y', order=grlex),
       f"PolynomialRing({QQ!r}, (Symbol('x'), Symbol('y')), GradedLexOrder())")
    sT(ZZ.inject('x', 'y', 'z', 't').eject('t'),
       f"PolynomialRing(UnivarPolynomialRing({ZZ!r}, (Symbol('t'),), "
       "LexOrder()), (Symbol('x'), Symbol('y'), Symbol('z')), LexOrder())")


def test_FractionField():
    sT(ZZ.inject('x').field, f"FractionField({ZZ!r}, (Symbol('x'),), LexOrder())")
    sT(QQ.frac_field('x', 'y', order=grlex),
       f"FractionField({QQ!r}, (Symbol('x'), Symbol('y')), GradedLexOrder())")
    sT(ZZ.inject('x', 'y', 'z', 't').eject('t').field,
       f"FractionField(UnivarPolynomialRing({ZZ!r}, (Symbol('t'),), LexOrder()), "
       "(Symbol('x'), Symbol('y'), Symbol('z')), LexOrder())")


def test_PolyElement():
    R, x, y = ring('x y', ZZ)
    g = R.domain.dtype
    assert repr(3*x**2*y + 1) == (f"PolyElement(PolynomialRing({ZZ!r}, (Symbol('x'), "
                                  "Symbol('y')), LexOrder()), [((2, 1), "
                                  f'{g(3)!r}), ((0, 0), {g(1)!r})])')


def test_FracElement():
    F, x, y = field('x y', ZZ)
    g = F.domain.dtype
    assert repr((3*x**2*y + 1)/(x - y**2)) == (f"FracElement(FractionField({ZZ!r}, (Symbol('x'), "
                                               f"Symbol('y')), LexOrder()), [((2, 1), {g(3)!r}), "
                                               f'((0, 0), {g(1)!r})], [((1, 0), {g(1)!r}), '
                                               f'((0, 2), {g(-1)!r})])')


def test_BooleanAtom():
    assert repr(true) == 'true'
    assert repr(false) == 'false'


def test_AlgebraicField():
    sT(QQ.algebraic_field(sqrt(2)),
       f'AlgebraicField({QQ!r}, Pow(Integer(2), Rational(1, 2)))')
