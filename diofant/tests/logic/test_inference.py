"""For more tests on satisfiability, see test_dimacs"""

import pytest

from diofant import numbered_symbols, pi
from diofant.abc import A, B, C, x, y
from diofant.logic.algorithms.dpll import (dpll, dpll_satisfiable,
                                           find_pure_symbol, find_unit_clause,
                                           unit_propagate)
from diofant.logic.algorithms.dpll2 import SATSolver
from diofant.logic.algorithms.dpll2 import \
    dpll_satisfiable as dpll2_satisfiable
from diofant.logic.boolalg import (And, Boolean, Equivalent, Implies, Or,
                                   false, true)
from diofant.logic.inference import (PropKB, entails, literal_symbol, pl_true,
                                     satisfiable, valid)


__all__ = ()


def test_literal():
    assert literal_symbol(True) is True
    assert literal_symbol(False) is False
    assert literal_symbol(A) is A
    assert literal_symbol(~A) is A

    pytest.raises(ValueError, lambda: literal_symbol(A + B))


def test_find_pure_symbol():
    assert find_pure_symbol([1], [{1}]) == (1, True)
    assert find_pure_symbol([1, 2],
                            [{-1, 2}, {-2, 1}]) == (None, None)
    assert find_pure_symbol([1, 2, 3],
                            [{1, -2}, {-2, -3}, {3, 1}]) == (1, True)
    assert find_pure_symbol([1, 2, 3],
                            [{-1, 2}, {2, -3}, {3, 1}]) == (2, True)
    assert find_pure_symbol([1, 2, 3],
                            [{-1, -2}, {-2, -3}, {3, 1}]) == (2, False)
    assert find_pure_symbol([1, 2, 3],
                            [{-1, 2}, {-2, -3}, {3, 1}]) == (None, None)


def test_unit_clause():
    assert find_unit_clause(map(set, [[1]]), {}) == (1, True)
    assert find_unit_clause(map(set, [[1], [-1]]), {}) == (1, True)
    assert find_unit_clause([{1, 2}], {1: True}) == (2, True)
    assert find_unit_clause([{1, 2}], {2: True}) == (1, True)
    assert find_unit_clause(map(set, [[1, 2, 3], [2, -3], [1, -2]]),
                            {1: True}) == (2, False)
    assert find_unit_clause(map(set, [[1, 2, 3], [3, -3], [1, 2]]),
                            {1: True}) == (2, True)


def test_unit_propagate():
    assert unit_propagate([{1, 2}], 1) == []
    assert unit_propagate(map(set, [[1, 2], [-1, 3],
                                    [-3, 2], [1]]), 1) == [{3}, {-3, 2}]


def test_dpll():
    assert dpll([{1, 2, 4}, {-4, -1, 2, 3}, {-3, -2, 1, 4}],
                [1, 2, 3, 4], {1: False}) == {1: False, 4: True}


def test_dpll_satisfiable():
    assert dpll_satisfiable(false) is False
    assert dpll_satisfiable( A & ~A ) is False
    assert dpll_satisfiable( A & ~B ) == {A: True, B: False}
    assert dpll_satisfiable(
        A | B ) in ({A: True}, {B: True}, {A: True, B: True})
    assert dpll_satisfiable(
        (~A | B) & (~B | A) ) in ({A: True, B: True}, {A: False, B: False})
    assert dpll_satisfiable( (A | B) & (~B | C) ) in ({A: True, B: False},
                                                      {A: True, C: True}, {B: True, C: True})
    assert dpll_satisfiable( A & B & C  ) == {A: True, B: True, C: True}
    assert dpll_satisfiable( (A | B) & (A >> B) ) == {B: True}
    assert dpll_satisfiable( Equivalent(A, B) & A ) == {A: True, B: True}
    assert dpll_satisfiable( Equivalent(A, B) & ~A ) == {A: False, B: False}


def test_dpll2_satisfiable():
    assert dpll2_satisfiable( A & ~A ) is False
    assert dpll2_satisfiable( A & ~B ) == {A: True, B: False}
    assert dpll2_satisfiable(
        A | B ) in ({A: True}, {B: True}, {A: True, B: True})
    assert dpll2_satisfiable(
        (~A | B) & (~B | A) ) in ({A: True, B: True}, {A: False, B: False})
    assert dpll2_satisfiable( (A | B) & (~B | C) ) in ({A: True, B: False, C: True},
                                                       {A: True, B: True, C: True})
    assert dpll2_satisfiable( A & B & C  ) == {A: True, B: True, C: True}
    assert dpll2_satisfiable( (A | B) & (A >> B) ) in ({B: True, A: False},
                                                       {B: True, A: True})
    assert dpll2_satisfiable( Equivalent(A, B) & A ) == {A: True, B: True}
    assert dpll2_satisfiable( Equivalent(A, B) & ~A ) == {A: False, B: False}

    l = SATSolver([], set(), set())
    assert l.lit_heap == []
    assert l._vsids_calculate() == 0


def test_satisfiable():
    assert satisfiable(A & (A >> B) & ~B) is False


def test_valid():
    assert valid(A >> (B >> A)) is True
    assert valid((A >> (B >> C)) >> ((A >> B) >> (A >> C))) is True
    assert valid((~B >> ~A) >> (A >> B)) is True
    assert valid(A | B | C) is False
    assert valid(A >> B) is False


def test_pl_true():
    assert pl_true(True) is True
    assert pl_true( A & B, {A: True, B: True}) is True
    assert pl_true( A | B, {A: True}) is True
    assert pl_true( A | B, {B: True}) is True
    assert pl_true( A | B, {A: None, B: True}) is True
    assert pl_true( A >> B, {A: False}) is True
    assert pl_true( A | B | ~C, {A: False, B: True, C: True}) is True
    assert pl_true(Equivalent(A, B), {A: False, B: False}) is True

    # test for false
    assert pl_true(False) is False
    assert pl_true( A & B, {A: False, B: False}) is False
    assert pl_true( A & B, {A: False}) is False
    assert pl_true( A & B, {B: False}) is False
    assert pl_true( A | B, {A: False, B: False}) is False

    # test for None
    assert pl_true(B, {B: None}) is None
    assert pl_true( A & B, {A: True, B: None}) is None
    assert pl_true( A >> B, {A: True, B: None}) is None
    assert pl_true(Equivalent(A, B), {A: None}) is None
    assert pl_true(Equivalent(A, B), {A: True, B: None}) is None

    # Test for deep
    assert pl_true(A | B, {A: False}, deep=True) is None
    assert pl_true(~A & ~B, {A: False}, deep=True) is None
    assert pl_true(A | B, {A: False, B: False}, deep=True) is False
    assert pl_true(A & B & (~A | ~B), {A: True}, deep=True) is False
    assert pl_true((C >> A) >> (B >> A), {C: True}, deep=True) is True


def test_pl_true_wrong_input():
    pytest.raises(ValueError, lambda: pl_true('John Cleese'))
    pytest.raises(ValueError, lambda: pl_true(42 + pi + pi ** 2))
    pytest.raises(ValueError, lambda: pl_true(42))


def test_entails():
    assert entails(A, [A >> B, ~B]) is False
    assert entails(B, [Equivalent(A, B), A]) is True
    assert entails((A >> B) >> (~A >> ~B)) is False
    assert entails((A >> B) >> (~B >> ~A)) is True


def test_PropKB():
    kb = PropKB()
    assert kb.ask(A >> B) is False
    assert kb.ask(A >> (B >> A)) is True
    kb.tell(A >> B)
    kb.tell(B >> C)
    assert kb.ask(A) is False
    assert kb.ask(B) is False
    assert kb.ask(C) is False
    assert kb.ask(~A) is False
    assert kb.ask(~B) is False
    assert kb.ask(~C) is False
    assert kb.ask(A >> C) is True
    kb.tell(A)
    assert kb.ask(A) is True
    assert kb.ask(B) is True
    assert kb.ask(C) is True
    assert kb.ask(~C) is False
    kb.retract(A)
    assert kb.ask(C) is False
    kb = PropKB((A >> B) & (B >> C))
    assert kb.ask(A >> C) is True


def test_propKB_tolerant():
    """tolerant to bad input"""
    kb = PropKB()
    assert kb.ask(B) is False


def test_satisfiable_non_symbols():
    class zero(Boolean):
        pass

    assumptions = zero(x*y)
    facts = Implies(zero(x*y), zero(x) | zero(y))
    query = ~zero(x) & ~zero(y)
    refutations = [
        {zero(x): True, zero(x*y): True},
        {zero(y): True, zero(x*y): True},
        {zero(x): True, zero(y): True, zero(x*y): True},
        {zero(x): True, zero(y): False, zero(x*y): True},
        {zero(x): False, zero(y): True, zero(x*y): True}]
    assert not satisfiable(And(assumptions, facts, query), algorithm='dpll')
    assert satisfiable(And(assumptions, facts, ~query), algorithm='dpll') in refutations
    assert not satisfiable(And(assumptions, facts, query), algorithm='dpll2')
    assert satisfiable(And(assumptions, facts, ~query), algorithm='dpll2') in refutations


def test_satisfiable_bool():
    assert satisfiable(true) == {true: true}
    assert satisfiable(false) is False


def test_satisfiable_all_models():
    assert next(satisfiable(False, all_models=True)) is False
    assert list(satisfiable((A >> ~A) & A, all_models=True)) == [False]
    assert list(satisfiable(True, all_models=True)) == [{true: true}]

    models = [{A: True, B: False}, {A: False, B: True}]
    result = satisfiable(A ^ B, all_models=True)
    models.remove(next(result))
    models.remove(next(result))
    pytest.raises(StopIteration, lambda: next(result))
    assert not models

    assert list(satisfiable(Equivalent(A, B), all_models=True)) == \
        [{A: False, B: False}, {A: True, B: True}]

    models = [{A: False, B: False}, {A: False, B: True}, {A: True, B: True}]
    for model in satisfiable(A >> B, all_models=True):
        models.remove(model)
    assert not models

    # This is a santiy test to check that only the required number
    # of solutions are generated. The expr below has 2**100 - 1 models
    # which would time out the test if all are generated at once.
    sym = numbered_symbols()
    X = [next(sym) for i in range(100)]
    result = satisfiable(Or(*X), all_models=True)
    for i in range(10):
        assert next(result)
