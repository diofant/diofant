from diofant import RootOf
from diofant.abc import x
from diofant.printing.lambdarepr import MpmathPrinter


__all__ = ()


def test_RootOf():
    p = MpmathPrinter()
    e = RootOf(x**3 + x - 1, x, 0)
    r = "findroot(lambda x: x**3 + x - 1, (mp.mpq(0, 1), mp.mpq(1, 1)), method='bisection')"
    assert p.doprint(e) == r
