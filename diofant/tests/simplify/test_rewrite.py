from diofant import I, Pow, cos, cot, exp, sin
from diofant.abc import x, y


__all__ = ()


def test_has():
    assert cot(x).has(x)
    assert cot(x).has(cot)
    assert not cot(x).has(sin)
    assert sin(x).has(x)
    assert sin(x).has(sin)
    assert not sin(x).has(cot)


def test_sin_exp_rewrite():
    assert sin(x).rewrite(sin, exp) == -I/2*(exp(I*x) - exp(-I*x))
    assert sin(x).rewrite(sin, exp).rewrite(Pow, sin) == sin(x)
    assert cos(x).rewrite(cos, exp).rewrite(Pow, cos) == cos(x)
    assert (sin(5*y) - sin(
        2*x)).rewrite(sin, exp).rewrite(Pow, sin) == sin(5*y) - sin(2*x)
    assert sin(x + y).rewrite(sin, exp).rewrite(Pow, sin) == sin(x + y)
    assert cos(x + y).rewrite(cos, exp).rewrite(Pow, cos) == cos(x + y)
    # This next test currently passes... not clear whether it should or not?
    assert cos(x).rewrite(cos, exp).rewrite(Pow, sin) == cos(x)
