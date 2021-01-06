"""Generic Rules for Diofant.

This file assumes knowledge of Basic and little else.
"""

import functools

from ..utilities.iterables import sift
from .basic import Atom, Basic


__all__ = ('arguments', 'operator', 'term', 'rm_id',
           'glom', 'flatten', 'unpack', 'sort')


@functools.singledispatch
def arguments(o):
    """Extract arguments from an expression."""
    return o.args


@arguments.register(int)
@arguments.register(Atom)
def arguments_atomic(o):
    return ()


@functools.singledispatch
def operator(o):
    """Extract the head of an expression."""
    return o.func


@operator.register(int)
@operator.register(Atom)
def operator_atomic(o):
    return o


@functools.singledispatch
def term(op, args):
    """Build an expression from the head and arguments."""
    return op(*args)


@term.register(int)
@term.register(Atom)
def term_atomic(op, args):
    return op


# Functions that create rules


def rm_id(isid):
    """Create a rule to remove identities

    isid - fn :: x -> Bool  --- whether or not this element is an identity

    >>> remove_zeros = rm_id(lambda x: x == 0)
    >>> remove_zeros(Basic(1, 0, 2))
    Basic(1, 2)
    >>> remove_zeros(Basic(0, 0))  # If only identites then we keep one
    Basic(0)

    See Also
    ========

    unpack

    """
    def ident_remove(expr):
        """Remove identities."""
        ids = list(map(isid, arguments(expr)))
        if sum(ids) == 0:           # No identities. Common case
            return expr
        elif sum(ids) != len(ids):  # there is at least one non-identity
            return term(operator(expr),
                        [arg for arg, x in zip(arguments(expr), ids) if not x])
        else:
            return term(operator(expr), [arguments(expr)[0]])

    return ident_remove


def glom(key, count, combine):
    """Create a rule to conglomerate identical args.

    >>> def key(x):
    ...     return x.as_coeff_Mul()[1]
    >>> def count(x):
    ...     return x.as_coeff_Mul()[0]
    >>> def combine(cnt, arg):
    ...     return cnt * arg
    >>> rl = glom(key, count, combine)

    >>> rl(Add(x, -x, 3*x, 2, 3, evaluate=False))
    3*x + 5

    Wait, how are key, count and combine supposed to work?

    >>> key(2*x)
    x
    >>> count(2*x)
    2
    >>> combine(2, x)
    2*x

    """
    def conglomerate(expr):
        """Conglomerate together identical args x + x -> 2x."""
        groups = sift(arguments(expr), key)
        counts = {k: sum(map(count, args)) for k, args in groups.items()}
        newargs = [combine(cnt, mat) for mat, cnt in counts.items()]
        if set(newargs) != set(arguments(expr)):
            return term(operator(expr), newargs)
        else:
            return expr

    return conglomerate


def sort(key):
    """Create a rule to sort by a key function.

    >>> sort_rl = sort(str)
    >>> sort_rl(Basic(3, 1, 2))
    Basic(1, 2, 3)

    """

    def sort_rl(expr):
        return term(operator(expr), sorted(arguments(expr), key=key))
    return sort_rl


# Functions that are rules


def unpack(expr):
    """Rule to unpack singleton args.

    >>> unpack(Basic(2))
    2

    """
    if len(arguments(expr)) == 1:
        return arguments(expr)[0]
    else:
        return expr


def flatten(expr):
    """Flatten T(a, b, T(c, d), T2(e)) to T(a, b, c, d, T2(e))."""
    cls = operator(expr)
    args = []
    for arg in arguments(expr):
        if operator(arg) == cls:
            args.extend(arguments(arg))
        else:
            args.append(arg)
    return term(cls, args)


def identity(x):
    return x


def switch(key, ruledict):
    """Select a rule based on the result of key called on the function."""
    def switch_rl(expr):
        rl = ruledict.get(key(expr), identity)
        return rl(expr)
    return switch_rl


def typed(ruletypes):
    """Apply rules based on the expression type.

    Examples
    ========

    >>> rm_zeros = rm_id(lambda x: x == 0)
    >>> rm_ones = rm_id(lambda x: x == 1)
    >>> remove_idents = typed({Add: rm_zeros, Mul: rm_ones})
    """
    return switch(type, ruletypes)


def treeapply(tree, join, leaf=identity):
    """Apply functions onto recursive containers (tree).

    join - a dictionary mapping container types to functions
      e.g. ``{list: minimize, tuple: chain}``

    Keys are containers/iterables.  Values are functions [a] -> a.

    Examples
    ========

    >>> tree = [(3, 2), (4, 1)]
    >>> treeapply(tree, {list: max, tuple: min})
    2

    >>> def mul(*args):
    ...     total = 1
    ...     for arg in args:
    ...         total *= arg
    ...     return total
    >>> treeapply(tree, {list: mul, tuple: lambda *args: sum(args)})
    25
    """
    for typ in join:
        if isinstance(tree, typ):
            return join[typ](*map(functools.partial(treeapply, join=join, leaf=leaf),
                                  tree))
    return leaf(tree)


def minimize(*rules, objective=identity):
    """Select result of rules that minimizes objective.

    Examples
    ========

    >>> from diofant.core.strategies import minimize

    >>> rl = minimize(lambda x: x + 1, lambda x: x - 1)
    >>> rl(4)
    3
    """
    def minrule(expr):
        return min((rule(expr) for rule in rules), key=objective)
    return minrule


def chain(*rules):
    """Compose a sequence of rules so that they apply to the expr sequentially."""
    def chain_rl(expr):
        for rule in rules:
            expr = rule(expr)
        return expr
    return chain_rl


def greedy(tree, objective=identity, **kwargs):
    """Execute a strategic tree.  Select alternatives greedily,

    Examples
    ========

    >>> tree = [lambda x: x + 1,
    ...         (lambda x: x - 1, lambda x: 2*x)]  # either inc or dec-then-double
    >>> fn = greedy(tree)
    >>> fn(4)  # lowest value comes from the inc
    5
    >>> fn(1)  # lowest value comes from dec then double
    0

    This function selects between options in a tuple.  The result is chosen that
    minimizes the objective function.

    >>> fn = greedy(tree, objective=lambda x: -x)  # maximize
    >>> fn(4)  # highest value comes from the dec then double
    6
    >>> fn(1)  # highest value comes from the inc
    2
    """
    optimize = functools.partial(minimize, objective=objective)
    return treeapply(tree, {list: optimize, tuple: chain}, **kwargs)


def do_one(rules):
    """Try each of the rules until one works. Then stop."""
    def do_one_rl(expr):
        for rl in rules:
            result = rl(expr)
            if result != expr:
                return result
        return expr
    return do_one_rl


def condition(cond, rule):
    """Only apply rule if condition is true."""
    def conditioned_rl(expr):
        if cond(expr):
            return rule(expr)
        else:
            return expr
    return conditioned_rl


def exhaust(rule):
    """Apply a rule repeatedly until it has no effect."""
    def exhaustive_rl(expr):
        new, old = rule(expr), expr
        while new != old:
            new, old = rule(new), new
        return new
    return exhaustive_rl


basic_fns = {'op': type,
             'new': Basic.__new__,
             'leaf': lambda x: not isinstance(x, Basic) or x.is_Atom,
             'children': lambda x: x.args}


def sall(rule, fns=basic_fns):
    """Strategic all - apply rule to args."""
    op, new, children, leaf = map(fns.get, ('op', 'new', 'children', 'leaf'))

    def all_rl(expr):
        if leaf(expr):
            return expr
        else:
            args = map(rule, children(expr))
            return new(op(expr), *args)

    return all_rl


def bottom_up(rule, fns=basic_fns):
    """Apply a rule down a tree running it on the bottom nodes first."""
    return chain(lambda expr: sall(bottom_up(rule, fns), fns)(expr), rule)


def null_safe(rule):
    """Return original expr if rule returns None."""
    def null_safe_rl(expr):
        result = rule(expr)
        if result is None:
            return expr
        else:
            return result
    return null_safe_rl
