from diofant import cos, sin
from diofant.abc import x
from diofant.printing.pretty.pretty import PrettyPrinter, pretty


__all__ = ()


def test_render():
    ustr = \
        """\
     3     5      \n\
    x     x     ⎛ \n\
x - ── + ─── + O⎝x\n\
    6    120      \n\
\n\
  \n\
6⎞\n\
 ⎠\n\
  \
"""
    assert pretty(sin(x).series(x), num_columns=20,
                  use_unicode=True) == ustr


def test_stringpict():
    p = PrettyPrinter()

    p1 = p._print(sin(x))
    p2 = p._print(sin(x))
    assert p1 == p2
    assert hash(p1) == hash(p2)

    p3 = p._print(cos(x))
    assert p1 != p3
    assert p1 != 0

    assert len(p1) == 1
