from diofant import cos, pretty, sin
from diofant.abc import x
from diofant.printing.pretty.pretty import PrettyPrinter


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


def test_dumb_term(capsys):
    print(pretty(sin(x)))
    assert capsys.readouterr().out == 'sin(x)\n'


def test_ncol():
    res = pretty(sin(x).series(x, n=20), num_columns=0)
    ans = """\
     3     5     7       9        11          13             15               \n\
    x     x     x       x        x           x              x                x\n\
x - ── + ─── - ──── + ────── - ──────── + ────────── - ───────────── + ───────\n\
    6    120   5040   362880   39916800   6227020800   1307674368000   3556874\n\
\n\
17                 19                 \n\
                  x              ⎛ 20⎞\n\
──────── - ────────────────── + O⎝x  ⎠\n\
28096000   121645100408832000         \
"""
    assert res == ans
