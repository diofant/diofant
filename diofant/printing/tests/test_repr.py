import pytest

from diofant import (Abs, AlgebraicNumber, Catalan, Dummy, E, EulerGamma,
                     Float, Function, GoldenRatio, I, ImmutableMatrix, Integer,
                     Matrix, Rational, Symbol, Wild, WildFunction, false, nan,
                     ones, oo, pi, root, sin, sqrt, symbols, true, zoo)
from diofant.domains import QQ, ZZ
from diofant.geometry import Ellipse, Point
from diofant.polys import field, grlex, lex, ring
from diofant.printing.repr import srepr


__all__ = ()

x, y = symbols('x,y')

# eval(repr(expr)) == expr has to succeed in the right environment. The right
# environment is the scope of "from diofant import *" for most cases.
ENV = {}
imports = """
from diofant import *
from diofant.polys.rings import PolyRing
from diofant.polys.fields import FracField
from diofant.polys.orderings import LexOrder, GradedLexOrder
"""
exec(imports, ENV)


def sT(expr, string):
    """
    sT := reprTest

    Tests that repr delivers the expected string and that
    the condition eval(repr(expr))==expr holds.
    """
    assert repr(expr) == string
    assert eval(string, ENV) == expr


def test_printmethod():
    class R(Abs):
        def _diofantrepr(self, printer):
            return "foo(%s)" % printer._print(self.args[0])
    assert repr(R(x)) == "foo(Symbol('x'))"


def test_Add():
    sT(x + y, "Add(Symbol('x'), Symbol('y'))")
    assert srepr(x**2 + 1, order='lex') == ("Add(Pow(Symbol('x'), "
                                            "Integer(2)), Integer(1))")


def test_Function():
    sT(Function("f")(x), "Function('f')(Symbol('x'))")
    # test unapplied Function
    sT(Function('f'), "Function('f')")

    sT(sin(x), "sin(Symbol('x'))")
    sT(sin, "sin")


def test_Geometry():
    sT(Point(0, 0), "Point2D(Integer(0), Integer(0))")
    sT(Ellipse(Point(0, 0), 5, 1),
       "Ellipse(Point2D(Integer(0), Integer(0)), Integer(5), Integer(1))")
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
    sT(Integer(4), "Integer(4)")


def test_list():
    sT([x, Integer(4)], "[Symbol('x'), Integer(4)]")


def test_Matrix():
    for cls, name in [(Matrix, "MutableDenseMatrix"), (ImmutableMatrix, "ImmutableMatrix")]:
        sT(cls([[x**+1, 1], [y, x + y]]),
           "%s([[Symbol('x'), Integer(1)], [Symbol('y'), Add(Symbol('x'), Symbol('y'))]])" % name)

        sT(cls(), "%s([])" % name)

        sT(cls([[x**+1, 1], [y, x + y]]), "%s([[Symbol('x'), Integer(1)], [Symbol('y'), Add(Symbol('x'), Symbol('y'))]])" % name)


def test_empty_Matrix():
    sT(ones(0, 3), "MutableDenseMatrix(0, 3, [])")
    sT(ones(4, 0), "MutableDenseMatrix(4, 0, [])")
    sT(ones(0, 0), "MutableDenseMatrix([])")


def test_Rational():
    sT(Rational(1, 3), "Rational(1, 3)")
    sT(Rational(-1, 3), "Rational(-1, 3)")


def test_AlgebraicNumber():
    a = AlgebraicNumber(sqrt(2))
    sT(a, "AlgebraicNumber(Pow(Integer(2), Rational(1, 2)), Tuple(Integer(1), Integer(0)))")
    a = AlgebraicNumber(sqrt(2), alias='a')
    sT(a, "AlgebraicNumber(Pow(Integer(2), Rational(1, 2)), Tuple(Integer(1), Integer(0)), Symbol('a'))")
    a = AlgebraicNumber(root(-2, 3))
    sT(a, "AlgebraicNumber(Pow(Integer(-2), Rational(1, 3)), Tuple(Integer(1), Integer(0)))")


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
    pytest.raises(TypeError, lambda: srepr(x, method="garbage"))


def test_Mul():
    sT(3*x**3*y, "Mul(Integer(3), Pow(Symbol('x'), Integer(3)), Symbol('y'))")


def test_PolyRing():
    sT(ring("x", ZZ, lex)[0], "PolyRing((Symbol('x'),), "
                              "%s, LexOrder())" % repr(ZZ))
    sT(ring("x,y", QQ, grlex)[0], "PolyRing((Symbol('x'), Symbol('y')), "
                                  "%s, GradedLexOrder())" % repr(QQ))
    sT(ring("x,y,z", ZZ["t"], lex)[0],
       "PolyRing((Symbol('x'), Symbol('y'), Symbol('z')), "
       "PolynomialRing(PolyRing((Symbol('t'),), "
       "%s, LexOrder())), LexOrder())" % repr(ZZ))


def test_FracField():
    sT(field("x", ZZ, lex)[0], "FracField((Symbol('x'),), "
                               "%s, LexOrder())" % repr(ZZ))
    sT(field("x,y", QQ, grlex)[0], "FracField((Symbol('x'), Symbol('y')), "
                                   "%s, GradedLexOrder())" % repr(QQ))
    sT(field("x,y,z", ZZ["t"], lex)[0],
       "FracField((Symbol('x'), Symbol('y'), Symbol('z')), "
       "PolynomialRing(PolyRing((Symbol('t'),), %s, "
       "LexOrder())), LexOrder())" % repr(ZZ))


def test_PolyElement():
    R, x, y = ring("x,y", ZZ)
    g = R.domain.dtype
    assert repr(3*x**2*y + 1) == ("PolyElement(PolyRing((Symbol('x'), "
                                  "Symbol('y')), %s, LexOrder()), [((2, 1), "
                                  "%s), ((0, 0), %s)])" % (repr(ZZ),
                                                           repr(g(3)),
                                                           repr(g(1))))


def test_FracElement():
    F, x, y = field("x,y", ZZ)
    g = F.domain.dtype
    assert repr((3*x**2*y + 1)/(x - y**2)) == ("FracElement(FracField((Symbol('x'), "
                                               "Symbol('y')), %s, LexOrder()), [((2, 1), %s), "
                                               "((0, 0), %s)], [((1, 0), %s), "
                                               "((0, 2), %s)])" % (repr(ZZ),
                                                                   repr(g(3)),
                                                                   repr(g(1)),
                                                                   repr(g(1)),
                                                                   repr(g(-1))))


def test_BooleanAtom():
    assert repr(true) == "true"
    assert repr(false) == "false"
