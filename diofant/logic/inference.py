"""Inference in propositional logic"""

from ..core.sympify import sympify
from ..utilities import ordered
from .boolalg import And, Not, to_cnf


def satisfiable(expr, algorithm='dpll2', all_models=False):
    """
    Check satisfiability of a propositional sentence.
    Returns a model when it succeeds.
    Returns {true: true} for trivially true expressions.

    On setting all_models to True, if given expr is satisfiable then
    returns a generator of models. However, if expr is unsatisfiable
    then returns a generator containing the single element False.

    Examples
    ========

    >>> satisfiable(a & ~b)
    {a: True, b: False}
    >>> satisfiable(a & ~a)
    False
    >>> satisfiable(True)
    {true: True}
    >>> next(satisfiable(a & ~a, all_models=True))
    False
    >>> models = satisfiable((a >> b) & b, all_models=True)
    >>> next(models)
    {a: False, b: True}
    >>> next(models)
    {a: True, b: True}
    >>> def use_models(models):
    ...     for model in models:
    ...         if model:
    ...             # Do something with the model.
    ...             return model
    ...         else:
    ...             # Given expr is unsatisfiable.
    ...             print('UNSAT')
    >>> use_models(satisfiable(a >> ~a, all_models=True))
    {a: False}
    >>> use_models(satisfiable(a ^ a, all_models=True))
    UNSAT

    """
    expr = to_cnf(expr)
    if algorithm == 'dpll':
        from .algorithms.dpll import dpll_satisfiable
        return dpll_satisfiable(expr)
    if algorithm == 'dpll2':
        from .algorithms.dpll2 import dpll_satisfiable
        return dpll_satisfiable(expr, all_models)
    raise NotImplementedError


def valid(expr):
    """
    Check validity of a propositional sentence.
    A valid propositional sentence is True under every assignment.

    Examples
    ========

    >>> valid(a | ~a)
    True
    >>> valid(a | b)
    False

    References
    ==========

    * https://en.wikipedia.org/wiki/Validity

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

    >>> pl_true(a & b, {a: True, b: True})
    True
    >>> pl_true(a & b, {a: False})
    False
    >>> pl_true(a & b, {a: True})
    >>> pl_true(a & b, {a: True}, deep=True)
    >>> pl_true(a >> (b >> a))
    >>> pl_true(a >> (b >> a), deep=True)
    True
    >>> pl_true(a & ~a)
    >>> pl_true(a & ~a, deep=True)
    False
    >>> pl_true(a & b & (~a | ~b), {a: True})
    >>> pl_true(a & b & (~a | ~b), {a: True}, deep=True)
    False

    """
    from ..core import Symbol
    from .boolalg import BooleanFunction
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
        raise ValueError(f'{expr} is not a valid boolean expression')
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


def entails(expr, formula_set={}):
    """
    Check whether the given expr_set entail an expr.
    If formula_set is empty then it returns the validity of expr.

    Examples
    ========

    >>> entails(a, [a >> b, b >> c])
    False
    >>> entails(c, [a >> b, b >> c, a])
    True
    >>> entails(a >> b)
    False
    >>> entails(a >> (b >> a))
    True

    References
    ==========

    * https://en.wikipedia.org/wiki/Logical_consequence

    """
    formula_set = list(formula_set)
    formula_set.append(Not(expr))
    return not satisfiable(And(*formula_set))


class KB:
    """Base class for all knowledge bases."""

    def __init__(self, sentence=None):
        """Initialize self."""
        self.clauses_ = set()
        if sentence:
            self.tell(sentence)

    def tell(self, sentence):
        raise NotImplementedError

    def ask(self, query):
        raise NotImplementedError

    def retract(self, sentence):
        raise NotImplementedError

    @property
    def clauses(self):
        return list(ordered(self.clauses_))


class PropKB(KB):
    """A KB for Propositional Logic.  Inefficient, with no indexing."""

    def tell(self, sentence):
        """Add the sentence's clauses to the KB

        Examples
        ========

        >>> l = PropKB()
        >>> l.clauses
        []

        >>> l.tell(x | y)
        >>> l.clauses
        [x | y]

        >>> l.tell(y)
        >>> l.clauses
        [y, x | y]

        """
        for c in And.make_args(to_cnf(sentence)):
            self.clauses_.add(c)

    def ask(self, query):
        """Checks if the query is true given the set of clauses.

        Examples
        ========

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

        >>> l = PropKB()
        >>> l.clauses
        []

        >>> l.tell(x | y)
        >>> l.clauses
        [x | y]

        >>> l.retract(x | y)
        >>> l.clauses
        []

        """
        for c in And.make_args(to_cnf(sentence)):
            self.clauses_.discard(c)
