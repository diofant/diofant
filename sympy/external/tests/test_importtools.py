import warnings

import pytest

from sympy.external import import_module

# fixes issue that arose in addressing issue 6533


def test_no_stdlib_collections():
    '''
    make sure we get the right collections when it is not part of a
    larger list
    '''
    import collections
    matplotlib = import_module('matplotlib',
        __import__kwargs={'fromlist': ['cm', 'collections']},
        min_module_version='1.1.0', catch=(RuntimeError,))
    if matplotlib:
        assert collections != matplotlib.collections


def test_no_stdlib_collections2():
    '''
    make sure we get the right collections when it is not part of a
    larger list
    '''
    import collections
    matplotlib = import_module('matplotlib',
        __import__kwargs={'fromlist': ['collections']},
        min_module_version='1.1.0', catch=(RuntimeError,))
    if matplotlib:
        assert collections != matplotlib.collections


def test_no_stdlib_collections3():
    '''make sure we get the right collections with no catch'''
    import collections
    matplotlib = import_module('matplotlib',
        __import__kwargs={'fromlist': ['cm', 'collections']},
        min_module_version='1.1.0')
    if matplotlib:
        assert collections != matplotlib.collections


def test_interface():
    with pytest.warns(UserWarning) as warn:
        import_module('spam_spam_spam', warn_not_installed=True)
    assert len(warn) == 1
    assert warn[0].message.args[0] == "spam_spam_spam module is not installed"

    with pytest.warns(UserWarning) as warn:
        import_module('re', warn_old_version=True, min_module_version="10.1")
    assert len(warn) == 1
    assert warn[0].message.args[0] == ("re version is too old to use "
                                       "(10.1 or newer required)")

    assert import_module('re', min_module_version="10.1") is None

    with pytest.warns(UserWarning) as warn:
        import_module('re', warn_old_version=True, min_python_version=(20, 10))
    assert len(warn) == 1
    assert warn[0].message.args[0] == ("Python version is too old to use re "
                                       "(20.10 or newer required)")
