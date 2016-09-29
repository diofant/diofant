"""Tests of tools for setting up interactive IPython sessions. """

import ast
import sys

import pytest

from diofant.interactive.session import IntegerWrapper


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
