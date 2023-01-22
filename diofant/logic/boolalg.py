"""Boolean algebra module for Diofant."""

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
        >>> (~(a & b & c)).equals(~a & ~b & ~c)
        False
        >>> (~(a & ~a)).equals(b | ~b)
        False

        """
        from ..core.relational import Relational
        from .inference import satisfiable

        other = sympify(other)

        if self.has(Relational) or other.has(Relational):
            raise NotImplementedError('handling of relationals')
        return (self.atoms() == other.atoms() and
                not satisfiable(~Equivalent(self, other)))


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

    BooleanFalse

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

    BooleanTrue

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

    >>> (x & y).subs({x: 1})
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
                newargs.append(bool(x))
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

        >>> ((x < 2) & (x > -2)).as_set()
        (-2, 2)

        """
        from ..sets import Intersection
        if len(self.free_symbols) == 1:
            return Intersection(*[arg.as_set() for arg in self.args])
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

    >>> (x | y).subs({x: 0})
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
                newargs.append(bool(x))
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

        >>> ((x > 2) | (x < -2)).as_set()
        [-oo, -2) U (2, oo]

        """
        from ..sets import Union
        if len(self.free_symbols) == 1:
            return Union(*[arg.as_set() for arg in self.args])
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
    >>> ~And(True, False)
    true
    >>> ~Or(True, False)
    false
    >>> ~(And(True, x) & Or(x, False))
    ~x
    >>> ~x
    ~x
    >>> ~((x | y) & (~x | ~y))
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
        raise NotImplementedError('Sorry, Not.as_set has not yet been'
                                  ' implemented for mutivariate'
                                  ' expressions')


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

    >>> (x ^ y).subs({y: 0})
    x

    """

    def __new__(cls, *args, **kwargs):
        argset = set()
        obj = super().__new__(cls, *args, **kwargs)
        for arg in super(Xor, obj).args:
            if isinstance(arg, Number) or arg in (True, False):
                if not arg:
                    continue
                arg = true
            if isinstance(arg, Xor):
                for a in arg.args:
                    if a in argset:
                        argset.remove(a)
                    else:
                        argset.add(a)
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
                if cj == c:
                    break
            else:
                continue
            remove.append((r, rj))
        if odd:
            assert true not in argset
            argset.add(true)
        for a, b in remove:
            argset.remove(a)
            argset.remove(b)
        if len(argset) == 0:
            return false
        if len(argset) == 1:
            return argset.pop()
        if True in argset:
            argset.remove(True)
            return ~Xor(*argset)
        obj._args = tuple(ordered(argset))
        obj._argset = frozenset(argset)
        return obj

    @property  # type: ignore[misc]
    @cacheit
    def args(self):
        return tuple(ordered(self._argset))


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
        return ~And(*args)


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
        return ~Or(*args)


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
                    newargs.append(bool(x))
                else:
                    newargs.append(x)
            A, B = newargs
        except ValueError as exc:
            raise ValueError(f'{len(args)} operand(s) used for an Implies '
                             f'(pairs are required): {args!s}') from exc
        if A == true or A == false or B == true or B == false:
            return Or(Not(A), B)
        if A == B:
            return true
        if A.is_Relational and B.is_Relational:
            if A.canonical == B.canonical:
                return true
            if (~A).canonical == B.canonical:
                return B
        else:
            return Expr.__new__(cls, *args)


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
                argset.add(bool(x))
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
                if cj == c:
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
        except ValueError as exc:
            raise ValueError('ITE expects exactly 3 arguments') from exc
        if a == true:
            return b
        if a == false:
            return c
        if b == c:
            return b
        if b == true and c == false:
            return a
        if b == false and c == true:
            return ~a

    def _eval_derivative(self, x):
        return self.func(self.args[0], *[a.diff(x) for a in self.args[1:]])


def _distribute(expr, a, b):
    """Distributes a over b with respect to expr."""
    if isinstance(expr, b):
        for arg in expr.args:
            if isinstance(arg, a):
                conj = arg
                break
        else:
            return expr
        rest = b(*(a for a in expr.args if a is not conj))
        return a(*(_distribute(b(c, rest), a, b) for c in conj.args))
    if isinstance(expr, a):
        return a(*(_distribute(c, a, b) for c in expr.args))
    return expr


def to_nnf(expr, simplify=True):
    """
    Converts expr to Negation Normal Form (NNF).

    If simplify is True, the result contains no redundant clauses.

    Examples
    ========

    >>> to_nnf(~((~a & ~b) | (c & d)))
    (a | b) & (~c | ~d)
    >>> to_nnf(Equivalent(a >> b, b >> a))
    (a | ~b | (a & ~b)) & (b | ~a | (b & ~a))

    See Also
    ========

    is_nnf

    """
    expr = sympify(expr)

    if is_nnf(expr, simplify):
        return expr

    if expr.is_Not:
        expr = expr.args[0]

        if isinstance(expr, And):
            expr = Or(*[~arg for arg in expr.args])
        elif isinstance(expr, Or):
            expr = And(*[~arg for arg in expr.args])
        elif isinstance(expr, Implies):
            a, b = expr.args
            expr = a & ~b
        elif isinstance(expr, Equivalent):
            expr = Or(*expr.args) & Or(*[~arg for arg in expr.args])
        elif isinstance(expr, Xor):
            args = []
            for i in range(1, len(expr.args) + 1, 2):
                for neg in combinations(expr.args, i):
                    args.append(Or(*[~s if s in neg else s for s in expr.args]))
            expr = And(*args)
        elif isinstance(expr, ITE):
            a, b, c = expr.args
            expr = (a | ~c) & (~a | ~b)
        else:
            raise ValueError(f'Illegal operator {expr.func} in expression')

    if isinstance(expr, Implies):
        a, b = expr.args
        expr = ~a | b

    if isinstance(expr, Equivalent):
        args = []
        for a, b in zip(expr.args, expr.args[1:]):
            args.append(~a | b)
        args.append(~expr.args[-1] | expr.args[0])
        expr = And(*args)

    if isinstance(expr, Xor):
        args = []
        for i in range(0, len(expr.args) + 1, 2):
            for neg in combinations(expr.args, i):
                args.append(Or(*[~s if s in neg else s for s in expr.args]))
        expr = And(*args)

    if isinstance(expr, ITE):
        a, b, c = expr.args
        expr = (~a | b) & (a | c)

    args = []
    for arg in expr.args:
        if not is_literal(arg):
            arg = to_nnf(arg, simplify)
        if simplify:
            arg = arg.args if isinstance(arg, expr.func) else (arg,)
            for a in arg:
                if ~a in args:
                    return expr.func.zero
                args.append(a)
        else:
            args.append(arg)

    return expr.func(*args)


def to_cnf(expr, simplify=False):
    """
    Convert expr to Conjunctive Normal Form (CNF).

    If simplify is True, the expr is evaluated to its simplest CNF form.

    Examples
    ========

    >>> to_cnf(~(a | b) | c)
    (c | ~a) & (c | ~b)
    >>> to_cnf((a | b) & (a | ~a), True)
    a | b

    See Also
    ========

    is_cnf

    """
    expr = sympify(expr)

    if not isinstance(expr, BooleanFunction):
        return expr

    if simplify:
        return simplify_logic(expr, 'cnf')

    if is_cnf(expr):
        return expr

    return _distribute(to_nnf(expr), And, Or)


def to_dnf(expr, simplify=False):
    """
    Convert expr to Disjunctive Normal Form (DNF).

    If simplify is True, the expr is evaluated to its simplest DNF form.

    Examples
    ========

    >>> to_dnf(b & (a | c))
    (a & b) | (b & c)
    >>> to_dnf((a & b) | (a & ~b) | (b & c) | (~b & c), True)
    a | c

    See Also
    ========

    is_dnf

    """
    expr = sympify(expr)

    if not isinstance(expr, BooleanFunction):
        return expr

    if simplify:
        return simplify_logic(expr, 'dnf')

    if is_dnf(expr):
        return expr

    return _distribute(to_nnf(expr), Or, And)


def is_nnf(expr, simplified=True):
    """
    Checks if expr is in Negation Normal Form (NNF).

    A logical expression is in NNF if the negation operator is only
    applied to literals and the only other allowed boolean functions
    are conjunction and disjunction.

    If simplified is True, checks if result contains no redundant clauses.

    Examples
    ========

    >>> is_nnf(a & b | ~c)
    True
    >>> is_nnf((a | ~a) & (b | c))
    False
    >>> is_nnf((a | ~a) & (b | c), False)
    True
    >>> is_nnf(~(a & b) | c)
    False
    >>> is_nnf((a >> b) & (b >> a))
    False

    See Also
    ========

    to_nnf
    is_cnf
    is_dnf

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
                    if ~arg in args:
                        return False
            stack.extend(expr.args)

        elif not is_literal(expr):
            return False

    return True


def is_cnf(expr):
    """
    Checks if expr is in Conjunctive Normal Form (CNF).

    A logical expression is in CNF if it is a conjunction of one or more
    clauses, where a clause is a disjunction of literals.

    Examples
    ========

    >>> is_cnf(a | b | c)
    True
    >>> is_cnf(a & b & c)
    True
    >>> is_cnf((a & b) | c)
    False

    See Also
    ========

    to_cnf
    is_dnf
    is_nnf

    """
    return _is_form(expr, And, Or)


def is_dnf(expr):
    """
    Checks if expr is in Disjunctive Normal Form (DNF).

    A logical expression is in DNF if it is a disjunction of one or more
    clauses, where a clause is a conjunction of literals.

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

    See Also
    ========

    to_dnf
    is_cnf
    is_nnf

    """
    return _is_form(expr, Or, And)


def _is_form(expr, function1, function2):
    expr = sympify(expr)

    if expr.is_Atom:
        return True

    if isinstance(expr, function2):
        for lit in expr.args:
            if isinstance(lit, Not):
                if not lit.args[0].is_Atom:
                    return False
            else:
                if not lit.is_Atom:
                    return False
        return True

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
    >>> is_literal(a | b)
    False

    """
    if isinstance(expr, Not):
        expr = expr.args[0]
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
            temp.append(~variables[i])
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
            temp.append(~variables[i])
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
        if x not in (3, minterm[i]):
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


def _SOPform(variables, minterms, dontcares=[]):
    """
    The _SOPform function uses simplified_pairs and a redundant group-
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
    >>> _SOPform([t, x, y, z], minterms, dontcares)
    (y & z) | (z & ~t)

    References
    ==========

    * https://en.wikipedia.org/wiki/Quine-McCluskey_algorithm

    """
    variables = [sympify(v) for v in variables]
    if minterms == []:
        return false

    minterms = [list(i) for i in minterms]
    dontcares = [list(i) for i in dontcares]
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


def _POSform(variables, minterms, dontcares=[]):
    """
    The _POSform function uses simplified_pairs and a redundant-group
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
    >>> _POSform([t, x, y, z], minterms, dontcares)
    z & (y | ~t)

    References
    ==========

    * https://en.wikipedia.org/wiki/Quine-McCluskey_algorithm

    """
    variables = [sympify(v) for v in variables]
    if minterms == []:
        return false

    minterms = [list(i) for i in minterms]
    dontcares = [list(i) for i in dontcares]
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


def simplify_logic(expr, form='cnf', deep=True):
    """
    This function simplifies a boolean function to its simplified version
    in SOP or POS form. The return type is an Or or And object in Diofant.

    Parameters
    ==========

    expr : string or boolean expression
    form : string ('cnf' or 'dnf'), default to 'cnf'.
        Selects the normal form in which the result is returned.
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
    if form in ('cnf', 'dnf'):
        expr = sympify(expr)
        variables = _find_predicates(expr)
        truthtable = []
        for t in product([0, 1], repeat=len(variables)):
            t = list(t)
            if expr.xreplace(dict(zip(variables, t))):
                truthtable.append(t)
        if deep:
            variables = [v.simplify() for v in variables]
        if form == 'dnf':
            return _SOPform(variables, truthtable)
        return _POSform(variables, truthtable)
    raise ValueError('form can be cnf or dnf only')
