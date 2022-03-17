"""For more tests on satisfiability, see test_dimacs."""

import copy

import pytest

from diofant import Equivalent, Or, Symbol, false, pi, satisfiable, true
from diofant.abc import a, b, c, x, y
from diofant.logic.algorithms.dpll import (dpll, find_pure_symbol,
                                           find_unit_clause, unit_propagate)
from diofant.logic.algorithms.dpll2 import SATSolver
from diofant.logic.boolalg import Boolean
from diofant.logic.inference import PropKB, entails, pl_true, valid


__all__ = ()


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
    assert find_unit_clause([{1}], {}) == (1, True)
    assert find_unit_clause([{1}, {-1}], {}) == (1, True)
    assert find_unit_clause([{1, 2}], {1: True}) == (2, True)
    assert find_unit_clause([{1, 2}], {2: True}) == (1, True)
    assert find_unit_clause([{1, 2, 3}, {2, -3}, {1, -2}],
                            {1: True}) == (2, False)
    assert find_unit_clause([{1, 2, 3}, {3, -3}, {1, 2}],
                            {1: True}) == (2, True)


def test_unit_propagate():
    assert unit_propagate([{1, 2}], 1) == []
    assert unit_propagate([{1, 2}, {3, -1}, {2, -3}, {1}], 1) == [{3}, {-3, 2}]


def test_dpll():
    assert dpll([{1, 2, 4}, {-4, -1, 2, 3}, {-3, -2, 1, 4}],
                [1, 2, 3, 4], {1: False}) == {1: False, 4: True}


def test_dpll2_satisfiable():
    l = SATSolver([], set(), set())
    assert not l.lit_heap
    assert l._vsids_calculate() == 0

    l = SATSolver([{1}, {-1}], {1}, set())
    with pytest.raises(StopIteration):
        next(l._find_model())
    assert l._clause_sat(0) is False
    assert l._clause_sat(1) is True

    l0 = SATSolver([{2, -3}, {1}, {3, -3}, {2, -2},
                    {3, -2}], {1, 2, 3}, set())

    l = copy.deepcopy(l0)
    assert l.num_learned_clauses == 0
    assert l.lit_scores == {-3: -2, -2: -2, -1: 0, 1: 0, 2: -2, 3: -2}
    l._vsids_clause_added({2, -3})
    assert l.num_learned_clauses == 1
    assert l.lit_scores == {-3: -1, -2: -2, -1: 0, 1: 0, 2: -1, 3: -2}

    l = copy.deepcopy(l0)
    assert l.num_learned_clauses == 0
    assert l.clauses == [[2, -3], [1], [3, -3], [2, -2], [3, -2]]
    assert l.sentinels == {-3: {0, 2}, -2: {3, 4}, 2: {0, 3}, 3: {2, 4}}
    l._simple_add_learned_clause([3])
    assert l.clauses == [[2, -3], [1], [3, -3], [2, -2], [3, -2], [3]]
    assert l.sentinels == {-3: {0, 2}, -2: {3, 4}, 2: {0, 3}, 3: {2, 4, 5}}

    l = copy.deepcopy(l0)
    assert l.lit_scores == {-3: -2, -2: -2, -1: 0, 1: 0, 2: -2, 3: -2}
    l._vsids_decay()
    assert l.lit_scores == {-3: -1, -2: -1, -1: 0, 1: 0, 2: -1, 3: -1}

    l = copy.deepcopy(l0)
    assert next(l._find_model()) == {1: True, 2: False, 3: False}
    assert l._current_level.decision == -3
    assert l._current_level.flipped is False
    assert l._current_level.var_settings == {-3, -2}
    assert l._simple_compute_conflict() == [3]
    assert l._is_sentinel(2, 3) is True
    assert l._is_sentinel(-3, 1) is False
    assert l.var_settings == {-3, -2, 1}
    level = l._current_level
    assert (level.decision, level.var_settings,
            level.flipped) == (-3, {-3, -2}, False)
    l._undo()
    level = l._current_level
    assert (level.decision, level.var_settings, level.flipped) == (0, {1}, False)

    l = copy.deepcopy(l0)
    l._assign_literal(-1)
    with pytest.raises(StopIteration):
        next(l._find_model())
    assert l.var_settings == {-1}

    l = copy.deepcopy(l0)
    assert list(l._find_model()) == [{1: True, 2: False, 3: False},
                                     {1: True, 2: True, 3: True}]

    l = copy.deepcopy(l0)
    assert l.variable_set == [False, False, False, False]
    l._simplify()
    assert l.variable_set == [False, True, False, False]
    assert l.sentinels == {-3: {0, 2}, -2: {3, 4}, -1: set(),
                           2: {0, 3}, 3: {2, 4}}

    l = copy.deepcopy(l0)
    assert l.lit_heap == [(-2, -3), (-2, 2), (-2, -2), (0, 1), (-2, 3), (0, -1)]
    assert l._vsids_calculate() == -3
    assert l.lit_heap == [(-2, -2), (-2, 2), (0, -1), (0, 1), (-2, 3)]

    l = copy.deepcopy(l0)
    l._vsids_lit_unset(2)
    assert l.lit_heap == [(-2, -3), (-2, -2), (-2, -2), (-2, 2), (-2, 3),
                          (0, -1), (-2, 2), (0, 1)]


def test_valid():
    assert valid(a >> (b >> a)) is True
    assert valid((a >> (b >> c)) >> ((a >> b) >> (a >> c))) is True
    assert valid((~b >> ~a) >> (a >> b)) is True
    assert valid(a | b | c) is False
    assert valid(a >> b) is False


def test_pl_true():
    assert pl_true(True) is True
    assert pl_true(a & b, {a: True, b: True}) is True
    assert pl_true(a | b, {a: True}) is True
    assert pl_true(a | b, {b: True}) is True
    assert pl_true(a | b, {a: None, b: True}) is True
    assert pl_true(a >> b, {a: False}) is True
    assert pl_true(a | b | ~c, {a: False, b: True, c: True}) is True
    assert pl_true(Equivalent(a, b), {a: False, b: False}) is True

    assert pl_true(False) is False
    assert pl_true(a & b, {a: False, b: False}) is False
    assert pl_true(a & b, {a: False}) is False
    assert pl_true(a & b, {b: False}) is False
    assert pl_true(a | b, {a: False, b: False}) is False

    assert pl_true(b, {b: None}) is None
    assert pl_true(a & b, {a: True, b: None}) is None
    assert pl_true(a >> b, {a: True, b: None}) is None
    assert pl_true(Equivalent(a, b), {a: None}) is None
    assert pl_true(Equivalent(a, b), {a: True, b: None}) is None

    assert pl_true(a | b, {a: False}, deep=True) is None
    assert pl_true(~a & ~b, {a: False}, deep=True) is None
    assert pl_true(a | b, {a: False, b: False}, deep=True) is False
    assert pl_true(a & b & (~a | ~b), {a: True}, deep=True) is False
    assert pl_true((c >> a) >> (b >> a), {c: True}, deep=True) is True

    pytest.raises(ValueError, lambda: pl_true('John Cleese'))
    pytest.raises(ValueError, lambda: pl_true(42 + pi))
    pytest.raises(ValueError, lambda: pl_true(42))


def test_entails():
    assert entails(a, [a >> b, ~b]) is False
    assert entails(b, [Equivalent(a, b), a]) is True
    assert entails((a >> b) >> (~a >> ~b)) is False
    assert entails((a >> b) >> (~b >> ~a)) is True


def test_PropKB():
    kb = PropKB()

    assert not kb.clauses
    assert kb.ask(b) is False
    assert kb.ask(a >> b) is False
    assert kb.ask(a >> (b >> a)) is True

    kb.tell(a >> b)

    assert kb.clauses == [~a | b]

    kb.tell(b >> c)

    assert kb.clauses == [~a | b, ~b | c]
    assert kb.ask(a) is False
    assert kb.ask(b) is False
    assert kb.ask(c) is False
    assert kb.ask(~a) is False
    assert kb.ask(~b) is False
    assert kb.ask(~c) is False
    assert kb.ask(a >> c) is True

    kb.tell(a)

    assert kb.ask(a) is True
    assert kb.ask(b) is True
    assert kb.ask(c) is True
    assert kb.ask(~c) is False

    kb.retract(a)

    assert kb.ask(c) is False

    kb = PropKB((a >> b) & (b >> c))

    assert kb.ask(a >> c) is True


@pytest.mark.parametrize('algorithm', ['dpll', 'dpll2'])
def test_satisfiable(algorithm):
    assert satisfiable(true, algorithm=algorithm) == {true: true}
    assert satisfiable(false, algorithm=algorithm) is False
    assert satisfiable(a & ~a, algorithm=algorithm) is False
    assert satisfiable(a & ~b, algorithm=algorithm) == {a: True, b: False}
    assert satisfiable(a | b, algorithm=algorithm) in ({a: True}, {b: True},
                                                       {a: True, b: True})
    assert satisfiable((~a | b) & (~b | a),
                       algorithm=algorithm) in ({a: True, b: True},
                                                {a: False, b: False})
    assert satisfiable((a | b) & (~b | c),
                       algorithm=algorithm) in ({a: True, b: False},
                                                {a: True, c: True},
                                                {b: True, c: True},
                                                {a: True, c: True, b: False})
    assert satisfiable(a & (a >> b) & ~b, algorithm=algorithm) is False
    assert satisfiable(a & b & c, algorithm=algorithm) == {a: True, b: True,
                                                           c: True}
    assert satisfiable((a | b) & (a >> b),
                       algorithm=algorithm) in ({b: True}, {b: True, a: False})
    assert satisfiable(Equivalent(a, b) & a, algorithm=algorithm) == {a: True,
                                                                      b: True}
    assert satisfiable(Equivalent(a, b) & ~a, algorithm=algorithm) == {a: False,
                                                                       b: False}

    class Zero(Boolean):
        pass

    assumptions = Zero(x*y)
    facts = Zero(x*y) >> (Zero(x) | Zero(y))
    query = ~Zero(x) & ~Zero(y)
    refutations = [{Zero(x): True, Zero(x*y): True},
                   {Zero(y): True, Zero(x*y): True},
                   {Zero(x): True, Zero(y): True, Zero(x*y): True},
                   {Zero(x): True, Zero(y): False, Zero(x*y): True},
                   {Zero(x): False, Zero(y): True, Zero(x*y): True}]
    assert not satisfiable(assumptions & facts & query, algorithm=algorithm)
    assert satisfiable(assumptions & facts & ~query,
                       algorithm=algorithm) in refutations


def test_satisfiable_all_models():
    assert next(satisfiable(False, all_models=True)) is False
    assert list(satisfiable((a >> ~a) & a, all_models=True)) == [False]
    assert list(satisfiable(True, all_models=True)) == [{true: true}]
    assert next(satisfiable(a & ~a, all_models=True)) is False

    models = [{a: True, b: False}, {a: False, b: True}]
    result = satisfiable(a ^ b, all_models=True)
    models.remove(next(result))
    models.remove(next(result))
    pytest.raises(StopIteration, lambda: next(result))
    assert not models

    assert list(satisfiable(Equivalent(a, b),
                            all_models=True)) == [{a: False, b: False},
                                                  {a: True, b: True}]

    models = [{a: False, b: False}, {a: False, b: True}, {a: True, b: True}]
    for model in satisfiable(a >> b, all_models=True):
        models.remove(model)
    assert not models

    # This is a santiy test to check that only the required number
    # of solutions are generated. The expr below has 2**100 - 1 models
    # which would time out the test if all are generated at once.
    X = [Symbol(f'x{i}') for i in range(100)]
    result = satisfiable(Or(*X), all_models=True)
    for _ in range(10):
        assert next(result)
