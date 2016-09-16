"""Tests of tools for setting up interactive IPython sessions. """

import ast
import sys

import pytest

from diofant.interactive.session import (init_ipython_session,
                                         IntegerWrapper)
from diofant.core import Symbol, Rational, Integer
from diofant.external import import_module

ipython = import_module("IPython", min_module_version="2.3.0")
readline = import_module("readline")

if not ipython:
    # py.test will not execute any tests now
    disabled = True


@pytest.mark.skipif(sys.version_info >= (3, 5),
                    reason="XXX python3.5 api changes")
def test_IntegerWrapper():
    tree = ast.parse('1/3')
    dump = ("Module(body=[Expr(value=BinOp(left=Call(func=Name(id='Integer', "
            "ctx=Load()), args=[Num(n=1)], keywords=[], starargs=None, "
            "kwargs=None), op=Div(), right=Call(func=Name(id='Integer', "
            "ctx=Load()), args=[Num(n=3)], keywords=[], starargs=None, "
            "kwargs=None)))])")
    tree = IntegerWrapper().visit(tree)
    assert ast.dump(tree) == dump
    tree2 = ast.parse('Integer(1)/Integer(3)')
    tree_new = IntegerWrapper().visit(tree2)
    assert ast.dump(tree_new) == dump
    dump3 = ("Module(body=[Expr(value=Call(func=Name(id='f', ctx=Load()), "
             "args=[Call(func=Name(id='Integer', ctx=Load()), args=[Num(n=1)], "
             "keywords=[], starargs=None, kwargs=None)], keywords=[], "
             "starargs=None, kwargs=None))])")
    tree3 = ast.parse('f(1)')
    tree_new = IntegerWrapper().visit(tree3)
    assert ast.dump(tree_new) == dump3
    tree_new2 = IntegerWrapper().visit(tree_new)
    assert ast.dump(tree_new2) == dump3


def test_automatic_symbols():
    # this implicitly requires readline
    if not readline:
        return
    # NOTE: Because of the way the hook works, you have to use run_cell(code,
    # True).  This means that the code must have no Out, or it will be printed
    # during the tests.
    app = init_ipython_session(auto_symbols=True)
    app.run_cell("from diofant import *")

    symbol = "verylongsymbolname"
    assert symbol not in app.user_ns
    app.run_cell("a = %s" % symbol, True)
#   assert symbol not in app.user_ns
    app.run_cell("a = type(%s)" % symbol, True)
    assert app.user_ns['a'] == Symbol
    app.run_cell("%s = Symbol('%s')" % (symbol, symbol), True)
    assert symbol in app.user_ns

    # Check that built-in names aren't overridden
    app.run_cell("a = all == __builtin__.all", True)
    assert "all" not in app.user_ns
    assert app.user_ns['a'] is True

    # Check that diofant names aren't overridden
    app.run_cell("import diofant")
    app.run_cell("a = factorial == diofant.factorial", True)
    assert app.user_ns['a'] is True


def test_int_to_Integer():
    # XXX: Warning, don't test with == here.  0.5 == Rational(1, 2) is True!
    app = init_ipython_session(auto_int_to_Integer=True)
    app.run_cell("from diofant import Integer")
    app.run_cell("a = 1")
    assert isinstance(app.user_ns['a'], Integer)

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
