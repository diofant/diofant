"""Inference in propositional logic"""

from diofant.logic.boolalg import And, Not, conjuncts, to_cnf
from diofant.core.compatibility import ordered
from diofant.core.sympify import sympify


def literal_symbol(literal):
    """
    The symbol in this literal (without the negation).

    Examples
    ========

    >>> from diofant.abc import A
    >>> from diofant.logic.inference import literal_symbol
    >>> literal_symbol(A)
    A
    >>> literal_symbol(~A)
    A
    """

    if literal is True or literal is False:
        return literal
    try:
        if literal.is_Symbol:
            return literal
        if literal.is_Not:
            return literal_symbol(literal.args[0])
        else:
            raise ValueError
    except (AttributeError, ValueError):
        raise ValueError("Argument must be a boolean literal.")


def satisfiable(expr, algorithm="dpll2", all_models=False):
    """
    Check satisfiability of a propositional sentence.
    Returns a model when it succeeds.
    Returns {true: true} for trivially true expressions.

    On setting all_models to True, if given expr is satisfiable then
    returns a generator of models. However, if expr is unsatisfiable
    then returns a generator containing the single element False.

    Examples
    ========

    >>> from diofant.abc import A, B
    >>> from diofant.logic.inference import satisfiable
    >>> satisfiable(A & ~B) == {A: True, B: False}
    True
    >>> satisfiable(A & ~A)
    False
    >>> satisfiable(True)
    {true: True}
    >>> next(satisfiable(A & ~A, all_models=True))
    False
    >>> models = satisfiable((A >> B) & B, all_models=True)
    >>> next(models) == {A: False, B: True}
    True
    >>> next(models) == {A: True, B: True}
    True
    >>> def use_models(models):
    ...     for model in models:
    ...         if model:
    ...             # Do something with the model.
    ...             return model
    ...         else:
    ...             # Given expr is unsatisfiable.
    ...             print("UNSAT")
    >>> use_models(satisfiable(A >> ~A, all_models=True))
    {A: False}
    >>> use_models(satisfiable(A ^ A, all_models=True))
    UNSAT
    """
    expr = to_cnf(expr)
    if algorithm == "dpll":
        from diofant.logic.algorithms.dpll import dpll_satisfiable
        return dpll_satisfiable(expr)
    elif algorithm == "dpll2":
        from diofant.logic.algorithms.dpll2 import dpll_satisfiable
        return dpll_satisfiable(expr, all_models)
    else:  # pragma: no cover
        raise NotImplementedError


def valid(expr):
    """
    Check validity of a propositional sentence.
    A valid propositional sentence is True under every assignment.

    Examples
    ========

    >>> from diofant.abc import A, B
    >>> from diofant.logic.inference import valid
    >>> valid(A | ~A)
    True
    >>> valid(A | B)
    False

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Validity

    """
    return not satisfiable(Not(expr))


def pl_true(expr, model={}, deep=False):
    """
    Returns whether the given assignment is a model or not.

    If the assignment does not specify the value for every proposition,
    this may return None to indicate 'not obvious'.

    Parameters
    ==========

    model : dict, optional, default: {}
        Mapping of symbols to boolean values to indicate assignment.
    deep: boolean, optional, default: False
        Gives the value of the expression under partial assignments
        correctly. May still return None to indicate 'not obvious'.


    Examples
    ========

    >>> from diofant.abc import A, B, C
    >>> from diofant.logic.inference import pl_true
    >>> pl_true( A & B, {A: True, B: True})
    True
    >>> pl_true(A & B, {A: False})
    False
    >>> pl_true(A & B, {A: True})
    >>> pl_true(A & B, {A: True}, deep=True)
    >>> pl_true(A >> (B >> A))
    >>> pl_true(A >> (B >> A), deep=True)
    True
    >>> pl_true(A & ~A)
    >>> pl_true(A & ~A, deep=True)
    False
    >>> pl_true(A & B & (~A | ~B), {A: True})
    >>> pl_true(A & B & (~A | ~B), {A: True}, deep=True)
    False
    """

    from diofant.core.symbol import Symbol
    from diofant.logic.boolalg import BooleanFunction
    boolean = (True, False)

    def _validate(expr):
        if isinstance(expr, Symbol) or expr in boolean:
            return True
        if not isinstance(expr, BooleanFunction):
            return False
        return all(_validate(arg) for arg in expr.args)

    if expr in boolean:
        return expr
    expr = sympify(expr)
    if not _validate(expr):
        raise ValueError("%s is not a valid boolean expression" % expr)
    model = {k: v for k, v in model.items() if v in boolean}
    result = expr.subs(model)
    if result in boolean:
        return bool(result)
    if deep:
        model = {k: True for k in result.atoms()}
        if pl_true(result, model):
            if valid(result):
                return True
        else:
            if not satisfiable(result):
                return False
    return


def entails(expr, formula_set={}):
    """
    Check whether the given expr_set entail an expr.
    If formula_set is empty then it returns the validity of expr.

    Examples
    ========

    >>> from diofant.abc import A, B, C
    >>> from diofant.logic.inference import entails
    >>> entails(A, [A >> B, B >> C])
    False
    >>> entails(C, [A >> B, B >> C, A])
    True
    >>> entails(A >> B)
    False
    >>> entails(A >> (B >> A))
    True

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Logical_consequence

    """
    formula_set = list(formula_set)
    formula_set.append(Not(expr))
    return not satisfiable(And(*formula_set))


class KB:
    """Base class for all knowledge bases"""
    def __init__(self, sentence=None):
        self.clauses_ = set()
        if sentence:
            self.tell(sentence)

    def tell(self, sentence):
        raise NotImplementedError  # pragma: no cover

    def ask(self, query):
        raise NotImplementedError  # pragma: no cover

    def retract(self, sentence):
        raise NotImplementedError  # pragma: no cover

    @property
    def clauses(self):
        return list(ordered(self.clauses_))


class PropKB(KB):
    """A KB for Propositional Logic.  Inefficient, with no indexing."""

    def tell(self, sentence):
        """Add the sentence's clauses to the KB

        Examples
        ========

        >>> from diofant.logic.inference import PropKB
        >>> from diofant.abc import x, y
        >>> l = PropKB()
        >>> l.clauses
        []

        >>> l.tell(x | y)
        >>> l.clauses
        [Or(x, y)]

        >>> l.tell(y)
        >>> l.clauses
        [y, Or(x, y)]
        """
        for c in conjuncts(to_cnf(sentence)):
            self.clauses_.add(c)

    def ask(self, query):
        """Checks if the query is true given the set of clauses.

        Examples
        ========

        >>> from diofant.logic.inference import PropKB
        >>> from diofant.abc import x, y
        >>> l = PropKB()
        >>> l.tell(x & ~y)
        >>> l.ask(x)
        True
        >>> l.ask(y)
        False
        """
        return entails(query, self.clauses_)

    def retract(self, sentence):
        """Remove the sentence's clauses from the KB

        Examples
        ========

        >>> from diofant.logic.inference import PropKB
        >>> from diofant.abc import x, y
        >>> l = PropKB()
        >>> l.clauses
        []

        >>> l.tell(x | y)
        >>> l.clauses
        [Or(x, y)]

        >>> l.retract(x | y)
        >>> l.clauses
        []
        """
        for c in conjuncts(to_cnf(sentence)):
            self.clauses_.discard(c)
