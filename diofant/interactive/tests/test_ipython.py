"""Tests of tools for setting up interactive IPython sessions. """

import ast

from diofant.interactive.session import IntegerWrapper

__all__ = ()


def test_IntegerWrapper():
    tree = ast.parse('1/3')
    tree2 = ast.parse('Integer(1)/Integer(3)')
    dump = ast.dump(tree2)
    tree_new = IntegerWrapper().visit(tree)
    assert ast.dump(tree_new) == dump

    tree2_new = IntegerWrapper().visit(tree2)
    assert ast.dump(tree2_new) == dump

    tree3 = ast.parse('f(1)')
    dump3 = ast.dump(ast.parse('f(Integer(1))'))
    tree3_new = IntegerWrapper().visit(tree3)
    assert ast.dump(tree3_new) == dump3

    tree4 = ast.parse('sin(1/5).n()')
    dump4 = ast.dump(ast.parse('sin(Integer(1)/Integer(5)).n()'))
    tree4_new = IntegerWrapper().visit(tree4)
    assert ast.dump(tree4_new) == dump4
