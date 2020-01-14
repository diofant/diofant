"""Tests of tools for setting up interactive IPython sessions."""

import ast

from diofant.interactive.session import IntegerDivisionWrapper


__all__ = ()


def test_IntegerDivisionWrapper():
    tree = ast.parse('1/3')
    tree2 = ast.parse('Rational(1, 3)')
    dump = ast.dump(tree2)
    tree_new = IntegerDivisionWrapper().visit(tree)
    assert ast.dump(tree_new) == dump

    tree = ast.parse('1 + 3')
    tree_new = IntegerDivisionWrapper().visit(tree)
    assert ast.dump(tree_new) == ast.dump(tree)

    tree = ast.parse('-1/3')
    tree2 = ast.parse('Rational(-1, 3)')
    dump = ast.dump(tree2)
    tree_new = IntegerDivisionWrapper().visit(tree)
    assert ast.dump(tree_new) == dump

    tree = ast.parse('2**3/7')
    tree2 = ast.parse('Rational(2**3, 7)')
    dump = ast.dump(tree2)
    tree_new = IntegerDivisionWrapper().visit(tree)
    assert ast.dump(tree_new) == ast.dump(tree)

    tree = ast.parse('(3 + 5)/7')
    tree2 = ast.parse('Rational(3 + 5, 7)')
    dump = ast.dump(tree2)
    tree_new = IntegerDivisionWrapper().visit(tree)
    assert ast.dump(tree_new) == ast.dump(tree)

    tree = ast.parse('2**x/3')
    tree_new = IntegerDivisionWrapper().visit(tree)
    assert ast.dump(tree_new) == ast.dump(tree)
