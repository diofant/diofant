""" Generic Rules for Diofant

This file assumes knowledge of Basic and little else.
"""

from strategies.dispatch import dispatch

from ..utilities.iterables import sift
from .basic import Atom, Basic


__all__ = ('arguments', 'operator', 'term', 'rm_id',
           'glom', 'flatten', 'unpack', 'sort')


@dispatch(Basic)
def arguments(o):
    return o.args


@arguments.register((int, Atom))
def arguments_atomic(o):
    return ()


@dispatch(Basic)
def operator(o):
    return o.func


@operator.register((int, Atom))
def operator_atomic(o):
    return o


@dispatch(type, (tuple, list))
def term(op, args):
    return op(*args)


@term.register((int, Atom), (tuple, list))
def term_atomic(op, args):
    return op


# Functions that create rules


def rm_id(isid):
    """ Create a rule to remove identities

    isid - fn :: x -> Bool  --- whether or not this element is an identity

    >>> remove_zeros = rm_id(lambda x: x==0)
    >>> remove_zeros(Basic(1, 0, 2))
    Basic(1, 2)
    >>> remove_zeros(Basic(0, 0)) # If only identites then we keep one
    Basic(0)

    See Also
    ========

    unpack
    """
    def ident_remove(expr):
        """ Remove identities """
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
    """ Create a rule to conglomerate identical args

    >>> key     = lambda x: x.as_coeff_Mul()[1]
    >>> count   = lambda x: x.as_coeff_Mul()[0]
    >>> combine = lambda cnt, arg: cnt * arg
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
        """ Conglomerate together identical args x + x -> 2x """
        groups = sift(arguments(expr), key)
        counts = {k: sum(map(count, args)) for k, args in groups.items()}
        newargs = [combine(cnt, mat) for mat, cnt in counts.items()]
        if set(newargs) != set(arguments(expr)):
            return term(operator(expr), newargs)
        else:
            return expr

    return conglomerate


def sort(key):
    """ Create a rule to sort by a key function

    >>> sort_rl = sort(str)
    >>> sort_rl(Basic(3, 1, 2))
    Basic(1, 2, 3)
    """

    def sort_rl(expr):
        return term(operator(expr), sorted(arguments(expr), key=key))
    return sort_rl


# Functions that are rules


def unpack(expr):
    """ Rule to unpack singleton args

    >>> unpack(Basic(2))
    2
    """
    if len(arguments(expr)) == 1:
        return arguments(expr)[0]
    else:
        return expr


def flatten(expr):
    """ Flatten T(a, b, T(c, d), T2(e)) to T(a, b, c, d, T2(e)) """
    cls = operator(expr)
    args = []
    for arg in arguments(expr):
        if operator(arg) == cls:
            args.extend(arguments(arg))
        else:
            args.append(arg)
    return term(cls, args)
