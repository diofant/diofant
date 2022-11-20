"""Tests of tools for setting up interactive IPython sessions."""

import ast
import uuid

from diofant.interactive.session import (IntegerDivisionWrapper,
                                         unicode_identifiers)


__all__ = ()


def test_IntegerDivisionWrapper():
    tree = ast.parse('1/3')
    tree2 = ast.parse('Fraction(1, 3)')
    dump = ast.dump(tree2)
    tree_new = IntegerDivisionWrapper().visit(tree)
    assert ast.dump(tree_new) == dump

    tree = ast.parse('1 + 3')
    tree_new = IntegerDivisionWrapper().visit(tree)
    assert ast.dump(tree_new) == ast.dump(tree)

    tree = ast.parse('-1/3')
    tree2 = ast.parse('Fraction(-1, 3)')
    dump = ast.dump(tree2)
    tree_new = IntegerDivisionWrapper().visit(tree)
    assert ast.dump(tree_new) == dump

    tree = ast.parse('2**3/7')
    tree_new = IntegerDivisionWrapper().visit(tree)
    assert ast.dump(tree_new) == ast.dump(tree)

    tree = ast.parse('(3 + 5)/7')
    tree_new = IntegerDivisionWrapper().visit(tree)
    assert ast.dump(tree_new) == ast.dump(tree)

    tree = ast.parse('2**x/3')
    tree_new = IntegerDivisionWrapper().visit(tree)
    assert ast.dump(tree_new) == ast.dump(tree)


def test_unicode_identifiers(monkeypatch):
    class dummy_uuid4:
        _counter = 0

        @property
        def hex(self):
            self._counter += 1
            return str(self._counter)

    monkeypatch.setattr(uuid, 'uuid4', dummy_uuid4)

    assert unicode_identifiers(['ℕ = 1\n', 'N = 2\n',
                                'N + ℕ']) == ['_1 =1 \n', 'N =2 \n',
                                              'N +_1 ']
