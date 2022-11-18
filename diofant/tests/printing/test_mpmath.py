from diofant import QQ, GoldenRatio, Rational, RootOf, Sum, oo, pi, sqrt
from diofant.abc import n, x
from diofant.printing.lambdarepr import MpmathPrinter


__all__ = ()


def test_basic():
    p = MpmathPrinter()

    assert p.doprint(GoldenRatio) == 'phi'
    assert p.doprint(Rational(2)) == '2'
    assert p.doprint(Rational(2, 3)) == '2*power(3, -1)'


def test_Pow():
    p = MpmathPrinter()

    assert p.doprint(sqrt(pi)) == 'root(pi, 2)'
    assert p.doprint(pi**Rational(2, 3)) == 'root(pi, 3)**2'
    assert p.doprint(pi**Rational(-2, 3)) == 'power(root(pi, 3), -2)'
    assert p.doprint(pi**pi) == 'pi**pi'


def test_RootOf():
    p = MpmathPrinter()

    e = RootOf(x**3 + x - 1, x, 0)
    r = f"findroot(lambda x: x**3 + x - 1, ({p.doprint(QQ(0))}, {p.doprint(QQ(1))}), method='bisection')"

    assert p.doprint(e) == r

    e = RootOf(x**3 + x - 1, x, 1)
    r = f"findroot(lambda x: x**3 + x - 1, mpc({p.doprint(QQ(-3, 8))}, {p.doprint(QQ(-9, 8))}), method='secant')"

    assert p.doprint(e) == r


def test_Sum():
    p = MpmathPrinter()

    s = Sum(n**(-2), (n, 1, oo))

    assert p.doprint(s) == 'nsum(lambda n: power(n, -2), (1, inf))'
