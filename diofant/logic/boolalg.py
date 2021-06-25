"""
Boolean algebra module for Diofant.
"""

from collections import defaultdict
from itertools import combinations, product

from ..core import Atom, cacheit
from ..core.expr import Expr
from ..core.function import Application
from ..core.numbers import Number
from ..core.operations import LatticeOp
from ..core.singleton import S
from ..core.singleton import SingletonWithManagedProperties as Singleton
from ..core.sympify import converter, sympify
from ..utilities import ordered


class Boolean(Expr):
    """A boolean object is an object for which logic operations make sense."""

    def __and__(self, other):
        """Overloading for & operator."""
        return And(self, other)

    __rand__ = __and__

    def __or__(self, other):
        """Overloading for | operator."""
        return Or(self, other)

    __ror__ = __or__

    def __invert__(self):
        """Overloading for ~ operator."""
        return Not(self)

    def __rshift__(self, other):
        """Overloading for >> operator."""
        return Implies(self, other)

    def __lshift__(self, other):
        """Overloading for << operator."""
        return Implies(other, self)

    __rrshift__ = __lshift__
    __rlshift__ = __rshift__

    def __xor__(self, other):
        return Xor(self, other)

    __rxor__ = __xor__

    def equals(self, other, failing_expression=False):
        """
        Returns True if the given formulas have the same truth table.
        For two formulas to be equal they must have the same literals.

        Examples
        ========

        >>> (a >> b).equals(~b >> ~a)
        True
        >>> Not(And(a, b, c)).equals(And(Not(a), Not(b), Not(c)))
        False
        >>> Not(And(a, Not(a))).equals(Or(b, Not(b)))
        False

        """
        from ..core.relational import Relational
        from .inference import satisfiable

        other = sympify(other)

        if self.has(Relational) or other.has(Relational):
            raise NotImplementedError('handling of relationals')
        return self.atoms() == other.atoms() and \
            not satisfiable(Not(Equivalent(self, other)))


class BooleanAtom(Atom, Boolean):
    """Base class of BooleanTrue and BooleanFalse."""

    is_Boolean = True

    @property
    def canonical(self):
        return self

    def __int__(self):
        return int(bool(self))


class BooleanTrue(BooleanAtom, metaclass=Singleton):
    """Diofant version of True, a singleton that can be accessed via ``true``.

    This is the Diofant version of True, for use in the logic module. The
    primary advantage of using true instead of True is that shorthand boolean
    operations like ~ and >> will work as expected on this class, whereas with
    True they act bitwise on 1. Functions in the logic module will return this
    class when they evaluate to true.

    Notes
    =====

    There is liable to be some confusion as to when ``True`` should
    be used and when ``true`` should be used in various contexts
    throughout Diofant. An important thing to remember is that
    ``sympify(True)`` returns ``true``. This means that for the most
    part, you can just use ``True`` and it will automatically be converted
    to ``true`` when necessary, similar to how you can generally use 1
    instead of ``Integer(1)``.

    The rule of thumb is:

    "If the boolean in question can be replaced by an arbitrary symbolic
    ``Boolean``, like ``Or(x, y)`` or ``x > 1``, use ``true``.
    Otherwise, use ``True``".

    In other words, use ``true`` only on those contexts where the
    boolean is being used as a symbolic representation of truth.
    For example, if the object ends up in the ``.args`` of any expression,
    then it must necessarily be ``true`` instead of ``True``, as
    elements of ``.args`` must be ``Basic``. On the other hand,
    ``==`` is not a symbolic operation in Diofant, since it always returns
    ``True`` or ``False``, and does so in terms of structural equality
    rather than mathematical, so it should return ``True``. The assumptions
    system should use ``True`` and ``False``. Aside from not satisfying
    the above rule of thumb, the
    assumptions system uses a three-valued logic (``True``, ``False``, ``None``),
    whereas ``true`` and ``false`` represent a two-valued logic. When in
    doubt, use ``True``.

    "``true == True is True``."

    While "``true is True``" is ``False``, "``true == True``"
    is ``True``, so if there is any doubt over whether a function or
    expression will return ``true`` or ``True``, just use ``==``
    instead of ``is`` to do the comparison, and it will work in either
    case.  Finally, for boolean flags, it's better to just use ``if x``
    instead of ``if x is True``. To quote PEP 8:

    Don't compare boolean values to ``True`` or ``False``
    using ``==``.

    * Yes:   ``if greeting:``
    * No:    ``if greeting == True:``
    * Worse: ``if greeting is True:``

    Examples
    ========

    >>> sympify(True)
    true
    >>> ~true
    false
    >>> ~True
    -2
    >>> Or(True, False)
    true

    See Also
    ========

    diofant.logic.boolalg.BooleanFalse

    """

    def __bool__(self):
        return True

    def __hash__(self):
        return hash(True)

    def as_set(self):
        """
        Rewrite logic operators and relationals in terms of real sets.

        Examples
        ========

        >>> true.as_set()
        UniversalSet()

        """
        return S.UniversalSet


class BooleanFalse(BooleanAtom, metaclass=Singleton):
    """Diofant version of False, a singleton that can be accessed via ``false``.

    This is the Diofant version of False, for use in the logic module. The
    primary advantage of using false instead of False is that shorthand boolean
    operations like ~ and >> will work as expected on this class, whereas with
    False they act bitwise on 0. Functions in the logic module will return this
    class when they evaluate to false.

    Notes
    =====

    See note in :py:class:`~diofant.logic.boolalg.BooleanTrue`.

    Examples
    ========

    >>> sympify(False)
    false
    >>> false >> false
    true
    >>> False >> False
    0
    >>> Or(True, False)
    true

    See Also
    ========

    diofant.logic.boolalg.BooleanTrue

    """

    def __bool__(self):
        return False

    def __hash__(self):
        return hash(False)

    def as_set(self):
        """
        Rewrite logic operators and relationals in terms of real sets.

        Examples
        ========

        >>> false.as_set()
        EmptySet()

        """
        from ..sets import EmptySet
        return EmptySet()


true = BooleanTrue()
false: BooleanFalse = BooleanFalse()
# We want S.true and S.false to work, rather than S.BooleanTrue and
# S.BooleanFalse, but making the class and instance names the same causes some
# major issues (like the inability to import the class directly from this
# file).
S.true = true
S.false = false

converter[bool] = lambda x: true if x else false


class BooleanFunction(Application, Boolean):
    """Boolean function is a function that lives in a boolean space.

    This is used as base class for And, Or, Not, etc.

    """

    is_Boolean = True

    def _eval_simplify(self, ratio, measure):
        return simplify_logic(self)

    def to_nnf(self, simplify=True):
        return self._to_nnf(*self.args, simplify=simplify)

    @classmethod
    def _to_nnf(cls, *args, **kwargs):
        simplify = kwargs.get('simplify', True)
        argset = set()
        for arg in args:
            if not is_literal(arg):
                arg = arg.to_nnf(simplify)
            if simplify:
                if isinstance(arg, cls):
                    arg = arg.args
                else:
                    arg = arg,
                for a in arg:
                    if Not(a) in argset:
                        return cls.zero
                    argset.add(a)
            else:
                argset.add(arg)
        return cls(*argset)


class And(LatticeOp, BooleanFunction):
    """
    Logical AND function.

    It evaluates its arguments in order, giving False immediately
    if any of them are False, and True if they are all True.

    Examples
    ========

    >>> x & y
    x & y

    Notes
    =====

    The ``&`` operator is provided as a convenience, but note that its use
    here is different from its normal use in Python, which is bitwise
    and. Hence, ``And(a, b)`` and ``a & b`` will return different things if
    ``a`` and ``b`` are integers.

    >>> And(x, y).subs({x: 1})
    y

    """

    zero = false
    identity = true

    nargs = None

    @classmethod
    def _new_args_filter(cls, args):
        newargs = []
        rel = []
        for x in reversed(list(args)):
            if isinstance(x, Number) or x in (0, 1):
                newargs.append(True if x else False)
                continue
            if x.is_Relational:
                c = x.canonical
                if c in rel:
                    continue
                nc = (~c).canonical
                if any(r == nc for r in rel):
                    return [false]
                rel.append(c)
            newargs.append(x)
        return LatticeOp._new_args_filter(newargs, And)

    def as_set(self):
        """
        Rewrite logic operators and relationals in terms of real sets.

        Examples
        ========

        >>> And(x < 2, x > -2).as_set()
        (-2, 2)

        """
        from ..sets import Intersection
        if len(self.free_symbols) == 1:
            return Intersection(*[arg.as_set() for arg in self.args])
        else:
            raise NotImplementedError('Sorry, And.as_set has not yet been'
                                      ' implemented for multivariate'
                                      ' expressions')


class Or(LatticeOp, BooleanFunction):
    """
    Logical OR function

    It evaluates its arguments in order, giving True immediately
    if any of them are True, and False if they are all False.

    Examples
    ========

    >>> x | y
    x | y

    Notes
    =====

    The ``|`` operator is provided as a convenience, but note that its use
    here is different from its normal use in Python, which is bitwise
    or. Hence, ``Or(a, b)`` and ``a | b`` will return different things if
    ``a`` and ``b`` are integers.

    >>> Or(x, y).subs({x: 0})
    y

    """

    zero = true
    identity = false

    @classmethod
    def _new_args_filter(cls, args):
        newargs = []
        rel = []
        for x in args:
            if isinstance(x, Number) or x in (0, 1):
                newargs.append(True if x else False)
                continue
            if x.is_Relational:
                c = x.canonical
                if c in rel:
                    continue
                nc = (~c).canonical
                if any(r == nc for r in rel):
                    return [true]
                rel.append(c)
            newargs.append(x)
        return LatticeOp._new_args_filter(newargs, Or)

    def as_set(self):
        """
        Rewrite logic operators and relationals in terms of real sets.

        Examples
        ========

        >>> Or(x > 2, x < -2).as_set()
        [-oo, -2) U (2, oo]

        """
        from ..sets import Union
        if len(self.free_symbols) == 1:
            return Union(*[arg.as_set() for arg in self.args])
        else:
            raise NotImplementedError('Sorry, Or.as_set has not yet been'
                                      ' implemented for multivariate'
                                      ' expressions')


class Not(BooleanFunction):
    """
    Logical Not function (negation).

    Returns True if the statement is False.
    Returns False if the statement is True.

    Examples
    ========

    >>> Not(True)
    false
    >>> Not(False)
    true
    >>> Not(And(True, False))
    true
    >>> Not(Or(True, False))
    false
    >>> Not(And(And(True, x), Or(x, False)))
    ~x
    >>> ~x
    ~x
    >>> Not(And(Or(x, y), Or(~x, ~y)))
    ~((x | y) & (~x | ~y))

    Notes
    =====

    The ``~`` operator is provided as a convenience, but note that its use
    here is different from its normal use in Python, which is bitwise
    not. In particular, ``~a`` and ``Not(a)`` will be different if ``a`` is
    an integer. Furthermore, since bools in Python subclass from ``int``,
    ``~True`` is the same as ``~1`` which is ``-2``, which has a boolean
    value of True.  To avoid this issue, use the Diofant boolean types
    ``true`` and ``false``.

    >>> ~True
    -2
    >>> ~true
    false

    """

    is_Not = True

    @classmethod
    def eval(cls, arg):
        from ..core import (Equality, GreaterThan, LessThan, StrictGreaterThan,
                            StrictLessThan, Unequality)
        if isinstance(arg, Number) or arg in (True, False):
            return false if arg else true
        if arg.is_Not:
            return arg.args[0]
        # Simplify Relational objects.
        if isinstance(arg, Equality):
            return Unequality(*arg.args)
        if isinstance(arg, Unequality):
            return Equality(*arg.args)
        if isinstance(arg, StrictLessThan):
            return GreaterThan(*arg.args)
        if isinstance(arg, StrictGreaterThan):
            return LessThan(*arg.args)
        if isinstance(arg, LessThan):
            return StrictGreaterThan(*arg.args)
        if isinstance(arg, GreaterThan):
            return StrictLessThan(*arg.args)

    def as_set(self):
        """
        Rewrite logic operators and relationals in terms of real sets.

        Examples
        ========

        >>> Not(x > 0, evaluate=False).as_set()
        (-oo, 0]

        """
        if len(self.free_symbols) == 1:
            return self.args[0].as_set().complement(S.Reals)
        else:
            raise NotImplementedError('Sorry, Not.as_set has not yet been'
                                      ' implemented for mutivariate'
                                      ' expressions')

    def to_nnf(self, simplify=True):
        if is_literal(self):
            return self

        expr = self.args[0]

        func, args = expr.func, expr.args

        if func == And:
            return Or._to_nnf(*[~arg for arg in args], simplify=simplify)

        if func == Or:
            return And._to_nnf(*[~arg for arg in args], simplify=simplify)

        if func == Implies:
            a, b = args
            return And._to_nnf(a, ~b, simplify=simplify)

        if func == Equivalent:
            return And._to_nnf(Or(*args), Or(*[~arg for arg in args]), simplify=simplify)

        if func == Xor:
            result = []
            for i in range(1, len(args)+1, 2):
                for neg in combinations(args, i):
                    clause = [~s if s in neg else s for s in args]
                    result.append(Or(*clause))
            return And._to_nnf(*result, simplify=simplify)

        if func == ITE:
            a, b, c = args
            return And._to_nnf(Or(a, ~c), Or(~a, ~b), simplify=simplify)

        raise ValueError(f'Illegal operator {func} in expression')


class Xor(BooleanFunction):
    """
    Logical XOR (exclusive OR) function.

    Returns True if an odd number of the arguments are True and the rest are
    False.

    Returns False if an even number of the arguments are True and the rest are
    False.

    Examples
    ========

    >>> Xor(True, False)
    true
    >>> Xor(True, True)
    false
    >>> Xor(True, False, True, True, False)
    true
    >>> Xor(True, False, True, False)
    false
    >>> x ^ y
    Xor(x, y)

    Notes
    =====

    The ``^`` operator is provided as a convenience, but note that its use
    here is different from its normal use in Python, which is bitwise xor. In
    particular, ``a ^ b`` and ``Xor(a, b)`` will be different if ``a`` and
    ``b`` are integers.

    >>> Xor(x, y).subs({y: 0})
    x

    """

    def __new__(cls, *args, **kwargs):
        argset = set()
        obj = super().__new__(cls, *args, **kwargs)
        for arg in super(Xor, obj).args:
            if isinstance(arg, Number) or arg in (True, False):
                if not arg:
                    continue
                else:
                    arg = true
            if isinstance(arg, Xor):
                for a in arg.args:
                    argset.remove(a) if a in argset else argset.add(a)
            elif arg in argset:
                argset.remove(arg)
            else:
                argset.add(arg)
        rel = [(r, r.canonical, (~r).canonical) for r in argset if r.is_Relational]
        odd = False  # is number of complimentary pairs odd? start 0 -> False
        remove = []
        for i, (r, c, nc) in enumerate(rel):
            for j in range(i + 1, len(rel)):
                rj, cj = rel[j][:2]
                if cj == nc:
                    odd = ~odd
                    break
                elif cj == c:
                    break
            else:
                continue
            remove.append((r, rj))
        if odd:
            argset.remove(true) if true in argset else argset.add(true)
        for a, b in remove:
            argset.remove(a)
            argset.remove(b)
        if len(argset) == 0:
            return false
        elif len(argset) == 1:
            return argset.pop()
        elif True in argset:
            argset.remove(True)
            return Not(Xor(*argset))
        else:
            obj._args = tuple(ordered(argset))
            obj._argset = frozenset(argset)
            return obj

    @property  # type: ignore[misc]
    @cacheit
    def args(self):
        return tuple(ordered(self._argset))

    def to_nnf(self, simplify=True):
        args = []
        for i in range(0, len(self.args)+1, 2):
            for neg in combinations(self.args, i):
                clause = [~s if s in neg else s for s in self.args]
                args.append(Or(*clause))
        return And._to_nnf(*args, simplify=simplify)


class Nand(BooleanFunction):
    """
    Logical NAND function.

    It evaluates its arguments in order, giving True immediately if any
    of them are False, and False if they are all True.

    Returns True if any of the arguments are False.
    Returns False if all arguments are True.

    Examples
    ========

    >>> Nand(False, True)
    true
    >>> Nand(True, True)
    false
    >>> Nand(x, y)
    ~(x & y)

    """

    @classmethod
    def eval(cls, *args):
        return Not(And(*args))


class Nor(BooleanFunction):
    """
    Logical NOR function.

    It evaluates its arguments in order, giving False immediately if any
    of them are True, and True if they are all False.

    Returns False if any argument is True.
    Returns True if all arguments are False.

    Examples
    ========

    >>> Nor(True, False)
    false
    >>> Nor(True, True)
    false
    >>> Nor(False, True)
    false
    >>> Nor(False, False)
    true
    >>> Nor(x, y)
    ~(x | y)

    """

    @classmethod
    def eval(cls, *args):
        return Not(Or(*args))


class Implies(BooleanFunction):
    """
    Logical implication.

    A implies B is equivalent to !A v B

    Accepts two Boolean arguments; A and B.
    Returns False if A is True and B is False.
    Returns True otherwise.

    Examples
    ========

    >>> Implies(True, False)
    false
    >>> Implies(False, False)
    true
    >>> Implies(True, True)
    true
    >>> Implies(False, True)
    true
    >>> x >> y
    Implies(x, y)
    >>> y << x
    Implies(x, y)

    Notes
    =====

    The ``>>`` and ``<<`` operators are provided as a convenience, but note
    that their use here is different from their normal use in Python, which is
    bit shifts. Hence, ``Implies(a, b)`` and ``a >> b`` will return different
    things if ``a`` and ``b`` are integers.  In particular, since Python
    considers ``True`` and ``False`` to be integers, ``True >> True`` will be
    the same as ``1 >> 1``, i.e., 0, which has a truth value of False.  To
    avoid this issue, use the Diofant objects ``true`` and ``false``.

    >>> True >> False
    1
    >>> true >> false
    false

    """

    @classmethod
    def eval(cls, *args):
        try:
            newargs = []
            for x in args:
                if isinstance(x, Number) or x in (0, 1):
                    newargs.append(True if x else False)
                else:
                    newargs.append(x)
            A, B = newargs
        except ValueError:
            raise ValueError(f'{len(args)} operand(s) used for an Implies '
                             f'(pairs are required): {args!s}')
        if A == true or A == false or B == true or B == false:
            return Or(Not(A), B)
        elif A == B:
            return true
        elif A.is_Relational and B.is_Relational:
            if A.canonical == B.canonical:
                return true
            elif (~A).canonical == B.canonical:
                return B
        else:
            return Expr.__new__(cls, *args)

    def to_nnf(self, simplify=True):
        a, b = self.args
        return Or._to_nnf(~a, b, simplify=simplify)


class Equivalent(BooleanFunction):
    """
    Equivalence relation.

    Equivalent(A, B) is True iff A and B are both True or both False.

    Returns True if all of the arguments are logically equivalent.
    Returns False otherwise.

    Examples
    ========

    >>> Equivalent(False, False, False)
    true
    >>> Equivalent(True, False, False)
    false
    >>> Equivalent(x, And(x, True))
    true

    """

    def __new__(cls, *args, **options):
        from ..core.relational import Relational
        args = [sympify(arg, strict=True) for arg in args]

        argset = set(args)
        for x in args:
            if isinstance(x, Number) or x in [True, False]:  # Includes 0, 1
                argset.discard(x)
                argset.add(True if x else False)
        rel = []
        for r in argset:
            if isinstance(r, Relational):
                rel.append((r, r.canonical, (~r).canonical))
        remove = []
        for i, (r, c, nc) in enumerate(rel):
            for j in range(i + 1, len(rel)):
                rj, cj = rel[j][:2]
                if cj == nc:
                    return false
                elif cj == c:
                    remove.append((r, rj))
                    break
        for a, b in remove:
            argset.remove(a)
            argset.remove(b)
            argset.add(True)
        if len(argset) <= 1:
            return true
        if True in argset:
            argset.discard(True)
            return And(*argset)
        if False in argset:
            argset.discard(False)
            return And(*[~arg for arg in argset])
        _args = frozenset(argset)
        obj = super().__new__(cls, _args)
        obj._argset = _args
        return obj

    @property  # type: ignore[misc]
    @cacheit
    def args(self):
        return tuple(ordered(self._argset))

    def to_nnf(self, simplify=True):
        args = []
        for a, b in zip(self.args, self.args[1:]):
            args.append(Or(~a, b))
        args.append(Or(~self.args[-1], self.args[0]))
        return And._to_nnf(*args, simplify=simplify)


class ITE(BooleanFunction):
    """
    If then else clause.

    ITE(A, B, C) evaluates and returns the result of B if A is true
    else it returns the result of C.

    Examples
    ========

    >>> ITE(True, False, True)
    false
    >>> ITE(Or(True, False), And(True, True), Xor(True, True))
    true
    >>> ITE(x, y, z)
    ITE(x, y, z)
    >>> ITE(True, x, y)
    x
    >>> ITE(False, x, y)
    y
    >>> ITE(x, y, y)
    y

    """

    @classmethod
    def eval(cls, *args):
        try:
            a, b, c = args
        except ValueError:
            raise ValueError('ITE expects exactly 3 arguments')
        if a == true:
            return b
        elif a == false:
            return c
        elif b == c:
            return b
        elif b == true and c == false:
            return a
        elif b == false and c == true:
            return Not(a)

    def to_nnf(self, simplify=True):
        a, b, c = self.args
        return And._to_nnf(Or(~a, b), Or(a, c), simplify=simplify)

    def _eval_derivative(self, x):
        return self.func(self.args[0], *[a.diff(x) for a in self.args[1:]])


# end class definitions. Some useful methods


def conjuncts(expr):
    """Return a list of the conjuncts in the expr s.

    Examples
    ========

    >>> conjuncts(a & b) == frozenset([a, b])
    True
    >>> conjuncts(a | b) == frozenset([Or(a, b)])
    True

    """
    return And.make_args(expr)


def disjuncts(expr):
    """Return a list of the disjuncts in the sentence s.

    Examples
    ========

    >>> disjuncts(a | b) == frozenset([a, b])
    True
    >>> disjuncts(a & b) == frozenset([And(a, b)])
    True

    """
    return Or.make_args(expr)


def distribute_and_over_or(expr):
    """
    Given a sentence s consisting of conjunctions and disjunctions
    of literals, return an equivalent sentence in CNF.

    Examples
    ========

    >>> distribute_and_over_or(Or(a, And(Not(b), Not(c))))
    (a | ~b) & (a | ~c)

    """
    return _distribute((expr, And, Or))


def distribute_or_over_and(expr):
    """
    Given a sentence s consisting of conjunctions and disjunctions
    of literals, return an equivalent sentence in DNF.

    Note that the output is NOT simplified.

    Examples
    ========

    >>> distribute_or_over_and(And(Or(Not(a), b), c))
    (b & c) | (c & ~a)

    """
    return _distribute((expr, Or, And))


def _distribute(info):
    """Distributes info[1] over info[2] with respect to info[0]."""
    if isinstance(info[0], info[2]):
        for arg in info[0].args:
            if isinstance(arg, info[1]):
                conj = arg
                break
        else:
            return info[0]
        rest = info[2](*[a for a in info[0].args if a is not conj])
        return info[1](*list(map(_distribute,
                                 ((info[2](c, rest), info[1], info[2]) for c in conj.args))))
    elif isinstance(info[0], info[1]):
        return info[1](*list(map(_distribute,
                                 ((x, info[1], info[2]) for x in info[0].args))))
    else:
        return info[0]


def to_nnf(expr, simplify=True):
    """
    Converts expr to Negation Normal Form.
    A logical expression is in Negation Normal Form (NNF) if it
    contains only And, Or and Not, and Not is applied only to literals.
    If simplify is True, the result contains no redundant clauses.

    Examples
    ========

    >>> to_nnf(Not((~a & ~b) | (c & d)))
    (a | b) & (~c | ~d)
    >>> to_nnf(Equivalent(a >> b, b >> a))
    (a | ~b | (a & ~b)) & (b | ~a | (b & ~a))

    """
    expr = sympify(expr)
    if is_nnf(expr, simplify):
        return expr
    return expr.to_nnf(simplify)


def to_cnf(expr, simplify=False):
    """
    Convert a propositional logical sentence s to conjunctive normal form.
    That is, of the form ((A | ~B | ...) & (B | C | ...) & ...).
    If simplify is True, the expr is evaluated to its simplest CNF form.

    Examples
    ========

    >>> to_cnf(~(a | b) | c)
    (c | ~a) & (c | ~b)
    >>> to_cnf((a | b) & (a | ~a), True)
    a | b

    """
    expr = sympify(expr)
    if not isinstance(expr, BooleanFunction):
        return expr

    if simplify:
        return simplify_logic(expr, 'cnf', True)

    # Don't convert unless we have to
    if is_cnf(expr):
        return expr

    expr = eliminate_implications(expr)
    return distribute_and_over_or(expr)


def to_dnf(expr, simplify=False):
    """
    Convert a propositional logical sentence s to disjunctive normal form.
    That is, of the form ((A & ~B & ...) | (B & C & ...) | ...).
    If simplify is True, the expr is evaluated to its simplest DNF form.

    Examples
    ========

    >>> to_dnf(b & (a | c))
    (a & b) | (b & c)
    >>> to_dnf((a & b) | (a & ~b) | (b & c) | (~b & c), True)
    a | c

    """
    expr = sympify(expr)
    if not isinstance(expr, BooleanFunction):
        return expr

    if simplify:
        return simplify_logic(expr, 'dnf', True)

    # Don't convert unless we have to
    if is_dnf(expr):
        return expr

    expr = eliminate_implications(expr)
    return distribute_or_over_and(expr)


def is_nnf(expr, simplified=True):
    """
    Checks if expr is in Negation Normal Form.
    A logical expression is in Negation Normal Form (NNF) if it
    contains only And, Or and Not, and Not is applied only to literals.
    If simplified is True, checks if result contains no redundant clauses.

    Examples
    ========

    >>> is_nnf(a & b | ~c)
    True
    >>> is_nnf((a | ~a) & (b | c))
    False
    >>> is_nnf((a | ~a) & (b | c), False)
    True
    >>> is_nnf(Not(a & b) | c)
    False
    >>> is_nnf((a >> b) & (b >> a))
    False

    """
    expr = sympify(expr)
    if is_literal(expr):
        return True

    stack = [expr]

    while stack:
        expr = stack.pop()
        if expr.func in (And, Or):
            if simplified:
                args = expr.args
                for arg in args:
                    if Not(arg) in args:
                        return False
            stack.extend(expr.args)

        elif not is_literal(expr):
            return False

    return True


def is_cnf(expr):
    """
    Test whether or not an expression is in conjunctive normal form.

    Examples
    ========

    >>> is_cnf(a | b | c)
    True
    >>> is_cnf(a & b & c)
    True
    >>> is_cnf((a & b) | c)
    False

    """
    return _is_form(expr, And, Or)


def is_dnf(expr):
    """
    Test whether or not an expression is in disjunctive normal form.

    Examples
    ========

    >>> is_dnf(a | b | c)
    True
    >>> is_dnf(a & b & c)
    True
    >>> is_dnf((a & b) | c)
    True
    >>> is_dnf(a & (b | c))
    False

    """
    return _is_form(expr, Or, And)


def _is_form(expr, function1, function2):
    """Test whether or not an expression is of the required form."""
    expr = sympify(expr)

    # Special case of an Atom
    if expr.is_Atom:
        return True

    # Special case of a single expression of function2
    if isinstance(expr, function2):
        for lit in expr.args:
            if isinstance(lit, Not):
                if not lit.args[0].is_Atom:
                    return False
            else:
                if not lit.is_Atom:
                    return False
        return True

    # Special case of a single negation
    if isinstance(expr, Not):
        if not expr.args[0].is_Atom:
            return False

    if not isinstance(expr, function1):
        return False

    for cls in expr.args:
        if cls.is_Atom:
            continue
        if isinstance(cls, Not):
            if not cls.args[0].is_Atom:
                return False
        elif not isinstance(cls, function2):
            return False
        for lit in cls.args:
            if isinstance(lit, Not):
                if not lit.args[0].is_Atom:
                    return False
            else:
                if not lit.is_Atom:
                    return False

    return True


def eliminate_implications(expr):
    """
    Change >>, <<, and Equivalent into &, |, and ~. That is, return an
    expression that is equivalent to s, but has only &, |, and ~ as logical
    operators.

    Examples
    ========

    >>> eliminate_implications(Implies(a, b))
    b | ~a
    >>> eliminate_implications(Equivalent(a, b))
    (a | ~b) & (b | ~a)
    >>> eliminate_implications(Equivalent(a, b, c))
    (a | ~c) & (b | ~a) & (c | ~b)

    """
    return to_nnf(expr)


def is_literal(expr):
    """
    Returns True if expr is a literal, else False.

    Examples
    ========

    >>> is_literal(a)
    True
    >>> is_literal(~a)
    True
    >>> is_literal(a + b)
    True
    >>> is_literal(Or(a, b))
    False

    """
    if isinstance(expr, Not):
        return not isinstance(expr.args[0], BooleanFunction)
    else:
        return not isinstance(expr, BooleanFunction)


def to_int_repr(clauses, symbols):
    """
    Takes clauses in CNF format and puts them into an integer representation.

    Examples
    ========

    >>> to_int_repr([x | y, y], [x, y])
    [{1, 2}, {2}]

    """
    symbols = dict(zip(symbols, range(1, len(symbols) + 1)))

    def append_symbol(arg, symbols):
        if isinstance(arg, Not):
            return -symbols[arg.args[0]]
        else:
            return symbols[arg]

    return [{append_symbol(arg, symbols) for arg in Or.make_args(c)}
            for c in clauses]


def _check_pair(minterm1, minterm2):
    """
    Checks if a pair of minterms differs by only one bit. If yes, returns
    index, else returns -1.

    """
    index = -1
    for x, (i, j) in enumerate(zip(minterm1, minterm2)):
        if i != j:
            if index == -1:
                index = x
            else:
                return -1
    return index


def _convert_to_varsSOP(minterm, variables):
    """
    Converts a term in the expansion of a function from binary to it's
    variable form (for SOP).

    """
    temp = []
    for i, m in enumerate(minterm):
        if m == 0:
            temp.append(Not(variables[i]))
        elif m == 1:
            temp.append(variables[i])
    return And(*temp)


def _convert_to_varsPOS(maxterm, variables):
    """
    Converts a term in the expansion of a function from binary to it's
    variable form (for POS).

    """
    temp = []
    for i, m in enumerate(maxterm):
        if m == 1:
            temp.append(Not(variables[i]))
        elif m == 0:
            temp.append(variables[i])
    return Or(*temp)


def _simplified_pairs(terms):
    """
    Reduces a set of minterms, if possible, to a simplified set of minterms
    with one less variable in the terms using QM method.

    """
    simplified_terms = []
    todo = list(range(len(terms)))
    for i, ti in enumerate(terms[:-1]):
        for j_i, tj in enumerate(terms[(i + 1):]):
            index = _check_pair(ti, tj)
            if index != -1:
                todo[i] = todo[j_i + i + 1] = None
                newterm = ti[:]
                newterm[index] = 3
                if newterm not in simplified_terms:
                    simplified_terms.append(newterm)
    simplified_terms.extend(
        [terms[i] for i in [_ for _ in todo if _ is not None]])
    return simplified_terms


def _compare_term(minterm, term):
    """
    Return True if a binary term is satisfied by the given term. Used
    for recognizing prime implicants.

    """
    for i, x in enumerate(term):
        if x != 3 and x != minterm[i]:
            return False
    return True


def _rem_redundancy(l1, terms):
    """
    After the truth table has been sufficiently simplified, use the prime
    implicant table method to recognize and eliminate redundant pairs,
    and return the essential arguments.

    """
    essential = []
    for x in terms:
        temporary = []
        for y in l1:
            if _compare_term(x, y):
                temporary.append(y)
        if len(temporary) == 1:
            if temporary[0] not in essential:
                essential.append(temporary[0])
    for x in terms:
        for y in essential:
            if _compare_term(x, y):
                break
        else:
            for z in l1:  # pragma: no branch
                if _compare_term(x, z):
                    assert z not in essential
                    essential.append(z)
                    break

    return essential


def SOPform(variables, minterms, dontcares=None):
    """
    The SOPform function uses simplified_pairs and a redundant group-
    eliminating algorithm to convert the list of all input combos that
    generate '1' (the minterms) into the smallest Sum of Products form.

    The variables must be given as the first argument.

    Return a logical Or function (i.e., the "sum of products" or "SOP"
    form) that gives the desired outcome. If there are inputs that can
    be ignored, pass them as a list, too.

    The result will be one of the (perhaps many) functions that satisfy
    the conditions.

    Examples
    ========

    >>> minterms = [[0, 0, 0, 1], [0, 0, 1, 1],
    ...             [0, 1, 1, 1], [1, 0, 1, 1], [1, 1, 1, 1]]
    >>> dontcares = [[0, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 1]]
    >>> SOPform([t, x, y, z], minterms, dontcares)
    (y & z) | (z & ~t)

    References
    ==========

    * https://en.wikipedia.org/wiki/Quine-McCluskey_algorithm

    """
    variables = [sympify(v) for v in variables]
    if minterms == []:
        return false

    minterms = [list(i) for i in minterms]
    dontcares = [list(i) for i in (dontcares or [])]
    for d in dontcares:
        if d in minterms:
            raise ValueError(f'{d} in minterms is also in dontcares')

    old = None
    new = minterms + dontcares
    while new != old:
        old = new
        new = _simplified_pairs(old)
    essential = _rem_redundancy(new, minterms)
    return Or(*[_convert_to_varsSOP(x, variables) for x in essential])


def POSform(variables, minterms, dontcares=None):
    """
    The POSform function uses simplified_pairs and a redundant-group
    eliminating algorithm to convert the list of all input combinations
    that generate '1' (the minterms) into the smallest Product of Sums form.

    The variables must be given as the first argument.

    Return a logical And function (i.e., the "product of sums" or "POS"
    form) that gives the desired outcome. If there are inputs that can
    be ignored, pass them as a list, too.

    The result will be one of the (perhaps many) functions that satisfy
    the conditions.

    Examples
    ========

    >>> minterms = [[0, 0, 0, 1], [0, 0, 1, 1], [0, 1, 1, 1],
    ...             [1, 0, 1, 1], [1, 1, 1, 1]]
    >>> dontcares = [[0, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 1]]
    >>> POSform([t, x, y, z], minterms, dontcares)
    z & (y | ~t)

    References
    ==========

    * https://en.wikipedia.org/wiki/Quine-McCluskey_algorithm

    """
    variables = [sympify(v) for v in variables]
    if minterms == []:
        return false

    minterms = [list(i) for i in minterms]
    dontcares = [list(i) for i in (dontcares or [])]
    for d in dontcares:
        if d in minterms:
            raise ValueError(f'{d} in minterms is also in dontcares')

    maxterms = []
    for t in product([0, 1], repeat=len(variables)):
        t = list(t)
        if (t not in minterms) and (t not in dontcares):
            maxterms.append(t)
    old = None
    new = maxterms + dontcares
    while new != old:
        old = new
        new = _simplified_pairs(old)
    essential = _rem_redundancy(new, maxterms)
    return And(*[_convert_to_varsPOS(x, variables) for x in essential])


def _find_predicates(expr):
    """Helper to find logical predicates in BooleanFunctions.

    A logical predicate is defined here as anything within a BooleanFunction
    that is not a BooleanFunction itself.

    """
    if not isinstance(expr, BooleanFunction):
        return {expr}
    return set().union(*(_find_predicates(i) for i in expr.args))


def simplify_logic(expr, form=None, deep=True):
    """
    This function simplifies a boolean function to its simplified version
    in SOP or POS form. The return type is an Or or And object in Diofant.

    Parameters
    ==========

    expr : string or boolean expression
    form : string ('cnf' or 'dnf') or None (default).
        If 'cnf' or 'dnf', the simplest expression in the corresponding
        normal form is returned; if None, the answer is returned
        according to the form with fewest args (in CNF by default).
    deep : boolean (default True)
        indicates whether to recursively simplify any
        non-boolean functions contained within the input.

    Examples
    ========

    >>> b = (~x & ~y & ~z) | (~x & ~y & z)
    >>> simplify_logic(b)
    ~x & ~y

    >>> sympify(b)
    (z & ~x & ~y) | (~x & ~y & ~z)
    >>> simplify_logic(_)
    ~x & ~y

    """
    if form == 'cnf' or form == 'dnf' or form is None:
        expr = sympify(expr)
        if not isinstance(expr, BooleanFunction):
            return expr
        variables = _find_predicates(expr)
        truthtable = []
        for t in product([0, 1], repeat=len(variables)):
            t = list(t)
            if expr.xreplace(dict(zip(variables, t))):
                truthtable.append(t)
        if deep:
            from ..simplify import simplify
            variables = [simplify(v) for v in variables]
        if form == 'dnf' or \
           (form is None and len(truthtable) >= (2 ** (len(variables) - 1))):
            return SOPform(variables, truthtable)
        elif form == 'cnf' or form is None:  # pragma: no branch
            return POSform(variables, truthtable)
    else:
        raise ValueError('form can be cnf or dnf only')


def _finger(eq):
    """
    Assign a 5-item fingerprint to each symbol in the equation:
    [
    # of times it appeared as a Symbol,
    # of times it appeared as a Not(symbol),
    # of times it appeared as a Symbol in an And or Or,
    # of times it appeared as a Not(Symbol) in an And or Or,
    sum of the number of arguments with which it appeared,
    counting Symbol as 1 and Not(Symbol) as 2
    ]

    >>> eq = Or(And(Not(y), a), And(Not(y), b), And(x, y))
    >>> dict(_finger(eq))
    {(0, 0, 1, 0, 2): [x],
     (0, 0, 1, 0, 3): [a, b],
     (0, 0, 1, 2, 8): [y]}

    So y and x have unique fingerprints, but a and b do not.

    """
    f = eq.free_symbols
    d = {fi: [0] * 5 for fi in f}
    for a in eq.args:
        if a.is_Symbol:
            d[a][0] += 1
        elif a.is_Not:
            d[a.args[0]][1] += 1
        else:
            o = len(a.args) + sum(isinstance(ai, Not) for ai in a.args)
            for ai in a.args:
                if ai.is_Symbol:
                    d[ai][2] += 1
                    d[ai][-1] += o
                else:
                    d[ai.args[0]][3] += 1
                    d[ai.args[0]][-1] += o
    inv = defaultdict(list)
    for k, v in ordered(d.items()):
        inv[tuple(v)].append(k)
    return inv


def bool_map(bool1, bool2):
    """
    Return the simplified version of bool1, and the mapping of variables
    that makes the two expressions bool1 and bool2 represent the same
    logical behaviour for some correspondence between the variables
    of each.
    If more than one mappings of this sort exist, one of them
    is returned.
    For example, And(x, y) is logically equivalent to And(a, b) for
    the mapping {x: a, y:b} or {x: b, y:a}.
    If no such mapping exists, return False.

    Examples
    ========

    >>> function1 = SOPform([x, z, y], [[1, 0, 1], [0, 0, 1]])
    >>> function2 = SOPform([a, b, c], [[1, 0, 1], [1, 0, 0]])
    >>> bool_map(function1, function2)
    (y & ~z, {y: a, z: b})

    The results are not necessarily unique, but they are canonical. Here,
    ``(t, z)`` could be ``(a, d)`` or ``(d, a)``:

    >>> eq1 = Or(And(Not(y), t), And(Not(y), z), And(x, y))
    >>> eq2 = Or(And(Not(c), a), And(Not(c), d), And(b, c))
    >>> bool_map(eq1, eq2)
    ((x & y) | (t & ~y) | (z & ~y), {t: a, x: b, y: c, z: d})
    >>> eq = And(Xor(a, b), c, And(c, d))
    >>> bool_map(eq, eq.subs({c: x}))
    (c & d & (a | b) & (~a | ~b), {a: a, b: b, c: d, d: x})

    """

    def match(function1, function2):
        """Return the mapping that equates variables between two
        simplified boolean expressions if possible.

        By "simplified" we mean that a function has been denested
        and is either an And (or an Or) whose arguments are either
        symbols (x), negated symbols (Not(x)), or Or (or an And) whose
        arguments are only symbols or negated symbols. For example,
        And(x, Not(y), Or(w, Not(z))).

        Basic.match is not robust enough (see issue sympy/sympy#4835) so this is
        a workaround that is valid for simplified boolean expressions.

        """
        # do some quick checks
        if function1.__class__ != function2.__class__:
            return
        if len(function1.args) != len(function2.args):
            return
        if function1.is_Symbol:
            return {function1: function2}

        # get the fingerprint dictionaries
        f1 = _finger(function1)
        f2 = _finger(function2)

        # more quick checks
        if len(f1) != len(f2):
            return

        # assemble the match dictionary if possible
        matchdict = {}
        for k in f1:
            if k not in f2 or len(f1[k]) != len(f2[k]):
                return
            for i, x in enumerate(f1[k]):
                matchdict[x] = f2[k][i]
        return matchdict if matchdict else None

    a = simplify_logic(bool1)
    b = simplify_logic(bool2)
    m = match(a, b)
    if m:
        return a, m
    return m is not None
