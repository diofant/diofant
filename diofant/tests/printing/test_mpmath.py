from diofant import QQ, RootOf, Sum, oo
from diofant.abc import n, x
from diofant.printing.lambdarepr import MpmathPrinter


__all__ = ()


def test_RootOf():
    p = MpmathPrinter()
    e = RootOf(x**3 + x - 1, x, 0)
    r = ("findroot(lambda x: x**3 + x - 1, (%s, %s), "
         "method='bisection')" % (p.doprint(QQ(0)), p.doprint(QQ(1))))
    assert p.doprint(e) == r


def test_Sum():
    p = MpmathPrinter()
    s = Sum(n**(-2), (n, 1, oo))
    assert p.doprint(s) == 'nsum(lambda n: n**(-2), (1, inf))'
