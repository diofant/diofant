"""Tests of tools for setting up interactive IPython sessions. """

from sympy.interactive.session import (init_ipython_session,
    enable_automatic_symbols, enable_automatic_int_sympification)

from sympy.core import Symbol, Rational, Integer
from sympy.external import import_module

ipython = import_module("IPython", min_module_version="2.3.0")
readline = import_module("readline")

if not ipython:
    # py.test will not execute any tests now
    disabled = True


def test_automatic_symbols():
    # this implicitly requires readline
    if not readline:
        return
    # NOTE: Because of the way the hook works, you have to use run_cell(code,
    # True).  This means that the code must have no Out, or it will be printed
    # during the tests.
    app = init_ipython_session()
    app.run_cell("from sympy import *")

    enable_automatic_symbols(app)

    symbol = "verylongsymbolname"
    assert symbol not in app.user_ns
    app.run_cell("a = %s" % symbol, True)
    assert symbol not in app.user_ns
    app.run_cell("a = type(%s)" % symbol, True)
    assert app.user_ns['a'] == Symbol
    app.run_cell("%s = Symbol('%s')" % (symbol, symbol), True)
    assert symbol in app.user_ns

    # Check that built-in names aren't overridden
    app.run_cell("a = all == __builtin__.all", True)
    assert "all" not in app.user_ns
    assert app.user_ns['a'] is True

    # Check that sympy names aren't overridden
    app.run_cell("import sympy")
    app.run_cell("a = factorial == sympy.factorial", True)
    assert app.user_ns['a'] is True


def test_int_to_Integer():
    # XXX: Warning, don't test with == here.  0.5 == Rational(1, 2) is True!
    app = init_ipython_session()
    app.run_cell("from sympy import Integer")
    app.run_cell("a = 1")
    assert isinstance(app.user_ns['a'], int)

    enable_automatic_int_sympification(app)
    app.run_cell("a = 1/2")
    assert isinstance(app.user_ns['a'], Rational)
    app.run_cell("a = 1")
    assert isinstance(app.user_ns['a'], Integer)
    app.run_cell("a = int(1)")
    assert isinstance(app.user_ns['a'], int)
    app.run_cell("a = (1/\n2)")
    assert app.user_ns['a'] == Rational(1, 2)
    # TODO: How can we test that the output of a SyntaxError is the original
    # input, not the transformed input?
