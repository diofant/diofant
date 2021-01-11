import collections

import pytest

from diofant.core.compatibility import HAS_GMPY
from diofant.external import import_module


__all__ = ()

# fixes issue that arose in addressing issue sympy/sympy#6533


def test_no_stdlib_collections():
    """
    Make sure we get the right collections when it is not part of a
    larger list.
    """
    matplotlib = import_module('matplotlib',
                               import__kwargs={'fromlist': ['cm', 'collections']},
                               min_module_version='1.1.0', catch=(RuntimeError,))
    if matplotlib:
        assert collections != matplotlib.collections


def test_no_stdlib_collections2():
    """
    Make sure we get the right collections when it is not part of a
    larger list.
    """
    matplotlib = import_module('matplotlib',
                               import__kwargs={'fromlist': ['collections']},
                               min_module_version='1.1.0', catch=(RuntimeError,))
    if matplotlib:
        assert collections != matplotlib.collections


def test_no_stdlib_collections3():
    """Make sure we get the right collections with no catch."""
    matplotlib = import_module('matplotlib',
                               import__kwargs={'fromlist': ['cm', 'collections']},
                               min_module_version='1.1.0')
    if matplotlib:
        assert collections != matplotlib.collections


def test_interface():
    with pytest.warns(UserWarning) as warn:
        import_module('spam_spam_spam', warn_not_installed=True)
    assert len(warn) == 1
    assert warn[0].message.args[0] == 'spam_spam_spam module is not installed'

    assert import_module('spam_spam_spam') is None

    with pytest.warns(UserWarning) as warn:
        import_module('re', warn_old_version=True, min_module_version='10.1')
    assert len(warn) == 1
    assert warn[0].message.args[0] == ('re version is too old to use '
                                       '(10.1 or newer required)')

    assert import_module('re', warn_old_version=False,
                         min_module_version='10.1') is None
    assert import_module('re', min_module_version='0.1') is not None

    with pytest.warns(UserWarning) as warn:
        import_module('re', warn_old_version=True, min_python_version=(20, 10))
    assert len(warn) == 1
    assert warn[0].message.args[0] == ('Python version is too old to use re '
                                       '(20.10 or newer required)')
    assert import_module('re', warn_old_version=False,
                         min_python_version=(20, 10)) is None
    assert import_module('re', warn_old_version=False,
                         min_python_version=(3, 3)) is not None

    if HAS_GMPY:
        assert import_module('gmpy2', min_module_version='2.0.0',
                             module_version_attr='version',
                             module_version_attr_call_args=(),
                             warn_old_version=False) is not None
