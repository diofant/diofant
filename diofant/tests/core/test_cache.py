import gc
import sys
import weakref

import pytest

from diofant import cacheit, ordered, sstr, symbols
from diofant.abc import x
from diofant.core.cache import clear_cache, print_cache


__all__ = ()


@cacheit
def _emptyfn():
    """Test docstring."""


@cacheit
def _identity(x):
    return x


def test_cacheit_doc():
    assert _emptyfn.__doc__ == 'Test docstring.'
    assert _emptyfn.__name__ == '_emptyfn'


def test_cacheit():
    assert _identity(1) == 1
    assert _identity(1) == 1


def test_print_cache(capfd):
    clear_cache()
    _identity(x)
    info = _identity.cache_info()
    print_cache()
    resout, _ = capfd.readouterr()
    assert resout.find('_identity ' + str(info)) >= 0


@pytest.fixture(scope='function')
def clear_imports(request):
    # Clear namespace
    orig = sys.modules.copy()
    for m in list(sys.modules):
        if m.startswith('diofant'):
            del sys.modules[m]

    def restore_imports():
        for m in list(sys.modules):
            if m.startswith('diofant'):
                del sys.modules[m]

        for m in orig:
            sys.modules[m] = orig[m]

    request.addfinalizer(restore_imports)


def test_nocache(clear_imports, monkeypatch):
    """Regression tests with DIOFANT_USE_CACHE=False."""
    monkeypatch.setenv('DIOFANT_USE_CACHE', 'False')
    from diofant.core.cache import CACHE
    from diofant.core.symbol import Symbol
    from diofant.functions import exp, sin, sinh, sqrt

    # test that we don't use cache
    assert CACHE == []
    x = Symbol('x')
    assert CACHE == []

    # issue sympy/sympy#8840
    (1 + x)*x  # not raises

    # issue sympy/sympy#9413
    (2*x).is_complex  # not raises

    # see commit c459d18
    sin(x + x)  # not raises

    # see commit 53dd1eb
    mx = -Symbol('x', negative=False)
    assert mx.is_positive is not True

    px = 2*Symbol('x', positive=False)
    assert px.is_positive is not True

    # see commit 2eaaba2
    s = 1/sqrt(x**2)
    y = Symbol('y')
    result = s.subs({sqrt(x**2): y})
    assert result == 1/y

    # problem from https://groups.google.com/forum/#!topic/sympy/LkTMQKC_BOw
    # see commit c459d18
    a = Symbol('a', positive=True)
    f = exp(x*(-a - 1))
    g = sinh(x)
    f*g  # not raises


def test_sympyissue_8825():
    t1, t2 = symbols('t1:3')
    d = weakref.WeakKeyDictionary([(t1, 1), (t2, 2)])
    assert sstr(list(ordered(d.items()))) == '[(t1, 1), (t2, 2)]'
    del t1
    clear_cache()
    gc.collect()
    assert sstr(list(ordered(d.items()))) == '[(t2, 2)]'
