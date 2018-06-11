from diofant import RootOf, Sum, oo
from diofant.abc import n, x
from diofant.printing.lambdarepr import MpmathPrinter


__all__ = ()


def test_RootOf():
    p = MpmathPrinter()
    e = RootOf(x**3 + x - 1, x, 0)
    r = "findroot(lambda x: x**3 + x - 1, (mp.mpq(0, 1), mp.mpq(1, 1)), method='bisection')"
    assert p.doprint(e) == r


def test_Sum():
    p = MpmathPrinter()
    s = Sum(n**(-2), (n, 1, oo))
    assert p.doprint(s) == 'nsum(lambda n: n**(-2), (1, inf))'
