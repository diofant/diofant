"""Implementation of DPLL algorithm.

Further improvements: eliminate calls to pl_true, implement branching rules,
efficient unit propagation.

References
==========

* https://en.wikipedia.org/wiki/DPLL_algorithm
* https://www.researchgate.net/publication/242384772

"""

from ...utilities import default_sort_key
from ..boolalg import _find_predicates, conjuncts, to_cnf, to_int_repr


def dpll_satisfiable(expr):
    """Check satisfiability of a propositional sentence.

    It returns a model rather than True when it succeeds

    >>> dpll_satisfiable(a & ~b)
    {a: True, b: False}
    >>> dpll_satisfiable(a & ~a)
    False

    """
    clauses = conjuncts(to_cnf(expr))
    if False in clauses:
        return False
    symbols = sorted(_find_predicates(expr), key=default_sort_key)
    symbols_int_repr = set(range(1, len(symbols) + 1))
    clauses_int_repr = to_int_repr(clauses, symbols)
    result = dpll(clauses_int_repr, symbols_int_repr, {})
    if not result:
        return result
    output = {}
    for key in result:
        output.update({symbols[key - 1]: result[key]})
    return output


def dpll(clauses, symbols, model):
    """Compute satisfiability in a partial model.

    Clauses is an array of conjuncts.  Arguments are expected to be
    in integer representation

    >>> dpll([{1}, {2}, {3}], {1, 2}, {3: False})
    False

    """
    # compute DP kernel
    P, value = find_unit_clause(clauses, model)
    while P:
        model.update({P: value})
        symbols.remove(P)
        if not value:
            P = -P
        clauses = unit_propagate(clauses, P)
        P, value = find_unit_clause(clauses, model)
    P, value = find_pure_symbol(symbols, clauses)
    while P:
        model.update({P: value})
        symbols.remove(P)
        if not value:
            P = -P
        clauses = unit_propagate(clauses, P)
        P, value = find_pure_symbol(symbols, clauses)
    # end DP kernel
    unknown_clauses = []
    for c in clauses:
        val = pl_true_int_repr(c, model)
        if val is False:
            return False
        if val is not True:
            unknown_clauses.append(c)
    if not unknown_clauses:
        return model
    P = symbols.pop()
    model_copy = model.copy()
    model.update({P: True})
    model_copy.update({P: False})
    symbols_copy = symbols.copy()
    return (dpll(unit_propagate(unknown_clauses, P), symbols, model) or
            dpll(unit_propagate(unknown_clauses, -P), symbols_copy, model_copy))


def pl_true_int_repr(clause, model={}):
    """Lightweight version of pl_true.

    Argument clause represents the set of args of an Or clause. This is used
    inside dpll, it is not meant to be used directly.

    >>> pl_true_int_repr({1, 2}, {1: False})
    >>> pl_true_int_repr({1, 2}, {1: False, 2: False})
    False

    """
    result = False
    for lit in clause:
        if lit < 0:
            p = model.get(-lit)
            if p is not None:
                p = not p
        else:
            p = model.get(lit)
        if p is True:
            return True
        elif p is None:
            result = None
    return result


def unit_propagate(clauses, s):
    """Returns an equivalent set of clauses.

    If a set of clauses contains the unit clause l, the other clauses are
    simplified by the application of the two following rules:

      1. every clause containing l is removed
      2. in every clause that contains ~l this literal is deleted

    Arguments are expected to be in integer representation.

    >>> unit_propagate([{1, 2}, {3, -2}, {2}], 2)
    [{3}]

    """
    negated = {-s}
    return [clause - negated for clause in clauses if s not in clause]


def find_pure_symbol(symbols, unknown_clauses):
    """
    Find a symbol and its value if it appears only as a positive literal
    (or only as a negative) in clauses.

    Arguments are expected to be in integer representation.

    >>> find_pure_symbol({1, 2, 3}, [{1, -2}, {-2, -3}, {3, 1}])
    (1, True)

    """
    all_symbols = set().union(*unknown_clauses)
    found_pos = all_symbols.intersection(symbols)
    found_neg = all_symbols.intersection([-s for s in symbols])
    for p in found_pos:
        if -p not in found_neg:
            return p, True
    for p in found_neg:
        if -p not in found_pos:
            return -p, False
    return None, None


def find_unit_clause(clauses, model):
    """Find a unit clause has only 1 variable that is not bound in the model.

    Arguments are expected to be in integer representation.

    >>> find_unit_clause([{1, 2, 3}, {2, -3}, {1, -2}], {1: True})
    (2, False)

    """
    bound = set(model) | {-sym for sym in model}
    for clause in clauses:
        unbound = clause - bound
        if len(unbound) == 1:
            p = unbound.pop()
            if p < 0:
                return -p, False
            else:
                return p, True
    return None, None
