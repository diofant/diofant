from diofant.core.symbol import Symbol
from diofant.core.cache import cacheit, CACHE, print_cache, clear_cache

from diofant.abc import x


@cacheit
def _emptyfn():
    "test docstring"
    pass


@cacheit
def _identity(x):
    return x


def test_cacheit_doc():
    assert _emptyfn.__doc__ == "test docstring"
    assert _emptyfn.__name__ == "_emptyfn"


def test_cacheit():
    assert _identity(1) == 1
    assert _identity(1) == 1


def test_print_cache(capfd):
    clear_cache()
    wrapped = _identity.__wrapped__
    _identity(x)
    item = str(wrapped)
    head = '='*len(item)
    res = (head + "\n" + item + "\n" + head + "\n" +
           "  ((Symbol('x'), <class 'diofant.core.symbol.Symbol'>), (True,)) : x\n")
    print_cache()
    resout, _ = capfd.readouterr()
    assert resout == res
    assert dict(CACHE)[wrapped] == {((x, Symbol), (True,)): x}
