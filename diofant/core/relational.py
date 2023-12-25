from __future__ import annotations

from ..logic.boolalg import Boolean, BooleanAtom, false, true
from ..utilities import ordered
from .evalf import EvalfMixin
from .evaluate import global_evaluate
from .expr import Expr
from .function import _coeff_isneg
from .symbol import Dummy, Symbol
from .sympify import sympify


# Note, see issue sympy/sympy#4986.  Ideally, we wouldn't want to subclass both Boolean
# and Expr.

class Relational(Boolean, Expr, EvalfMixin):
    """Base class for all relation types.

    Subclasses of Relational should generally be instantiated directly, but
    Relational can be instantiated with a valid `rop` value to dispatch to
    the appropriate subclass.

    Parameters
    ==========
    rop : str or None
        Indicates what subclass to instantiate.  Valid values can be found
        in the keys of Relational.ValidRelationalOperator.

    Examples
    ========

    >>> Rel(y, x+x**2, '==')
    Eq(y, x**2 + x)

    """

    is_Relational = True

    ValidRelationOperator: dict[str | None, type[Relational]]

    # ValidRelationOperator - Defined below, because the necessary classes
    #   have not yet been defined

    def __new__(cls, lhs, rhs=0, rop=None, **assumptions):
        # If called by a subclass, do nothing special and pass on to Expr.
        if cls is not Relational:
            return Expr.__new__(cls, lhs, rhs, **assumptions)
        # If called directly with an operator, look up the subclass
        # corresponding to that operator and delegate to it
        try:
            new_cls = cls.ValidRelationOperator[rop]
            return new_cls(lhs, rhs, **assumptions)
        except KeyError as exc:
            raise ValueError('Invalid relational operator '
                             f'symbol: {rop!r}') from exc

    @property
    def lhs(self):
        """The left-hand side of the relation."""
        return self.args[0]

    @property
    def rhs(self):
        """The right-hand side of the relation."""
        return self.args[1]

    @property
    def reversed(self):
        """Return the relationship with sides (and sign) reversed.

        Examples
        ========

        >>> Eq(x, 1)
        Eq(x, 1)
        >>> _.reversed
        Eq(1, x)
        >>> x < 1
        x < 1
        >>> _.reversed
        1 > x

        """
        ops = {Gt: Lt, Ge: Le, Lt: Gt, Le: Ge}
        a, b = self.args
        return ops.get(self.func, self.func)(b, a, evaluate=False)

    def _eval_evalf(self, prec):
        return self.func(*[s._evalf(prec) for s in self.args])

    @property
    def canonical(self):
        """Return a canonical form of the relational.

        The rules for the canonical form, in order of decreasing priority are:
            1) Number on right if left is not a Number;
            2) Symbol on the left;
            3) Gt/Ge changed to Lt/Le;
            4) Lt/Le are unchanged;
            5) Eq and Ne get ordered args.

        """
        r = self
        if r.func in (Ge, Gt):
            r = r.reversed
        elif r.func in (Lt, Le):
            pass
        elif r.func in (Eq, Ne):
            r = r.func(*ordered(r.args), evaluate=False)
        else:
            raise NotImplementedError
        if r.lhs.is_Number and not r.rhs.is_Number:
            r = r.reversed
        elif r.rhs.is_Symbol and not r.lhs.is_Symbol:
            r = r.reversed
        if _coeff_isneg(r.lhs):
            r = r.reversed.func(-r.lhs, -r.rhs, evaluate=False)
        return r

    def equals(self, other, failing_expression=False):
        """Return True if the sides of the relationship are mathematically
        identical and the type of relationship is the same.
        If failing_expression is True, return the expression whose truth value
        was unknown.

        """
        if isinstance(other, Relational):
            if other in (self, self.reversed):
                return True
            a, b = self, other
            if a.func in (Eq, Ne) or b.func in (Eq, Ne):
                if a.func != b.func:
                    return False
                l = a.lhs.equals(b.lhs, failing_expression=failing_expression)
                r = a.rhs.equals(b.rhs, failing_expression=failing_expression)
                if l is True:
                    return r
                if r is True:
                    return l
                lr = a.lhs.equals(b.rhs, failing_expression=failing_expression)
                rl = a.rhs.equals(b.lhs, failing_expression=failing_expression)
                if lr is True:
                    return rl
                if rl is True:
                    return lr
                e = (l, r, lr, rl)
                if all(i is False for i in e):
                    return False
                if failing_expression:
                    return a.lhs - a.rhs - b.lhs + b.rhs
            else:
                if b.func != a.func:
                    b = b.reversed
                if a.func != b.func:
                    return False
                l = a.lhs.equals(b.lhs, failing_expression=failing_expression)
                if l is False:
                    return False
                r = a.rhs.equals(b.rhs, failing_expression=failing_expression)
                if r is False:
                    return False
                if l is True:
                    return r
                return l

    def _eval_simplify(self, ratio, measure):
        r = self.func(self.lhs.simplify(ratio=ratio, measure=measure),
                      self.rhs.simplify(ratio=ratio, measure=measure))

        if r.is_Relational:
            dif = r.lhs - r.rhs
            # We want a Number to compare with zero and be sure to get a
            # True/False answer.  Check if we can deduce that dif is
            # definitively zero or non-zero.
            if not dif.has(Dummy, Symbol):
                know = dif.equals(0)
                if know is False:
                    if isinstance(r, Eq):
                        return False
                    if isinstance(r, Ne):
                        return True

        r = r.canonical
        if measure(r) < ratio*measure(self):
            return r
        return self

    def __bool__(self):
        raise TypeError('cannot determine truth value of Relational')

    def as_set(self):
        """
        Rewrites univariate inequality in terms of real sets

        Examples
        ========

        >>> (x > 0).as_set()
        (0, oo]
        >>> Eq(x, 0).as_set()
        {0}

        """
        from ..solvers.inequalities import solve_univariate_inequality
        syms = self.free_symbols

        if len(syms) == 1:
            sym = syms.pop()
        else:
            raise NotImplementedError('Sorry, Relational.as_set procedure'
                                      ' is not yet implemented for'
                                      ' multivariate expressions')

        return solve_univariate_inequality(self, sym, relational=False)


Rel = Relational


class Equality(Relational):
    """An equal relation between two objects.

    Represents that two objects are equal.  If they can be easily shown
    to be definitively equal (or unequal), this will reduce to True (or
    False).  Otherwise, the relation is maintained as an unevaluated
    Equality object.  Use the ``simplify`` function on this object for
    more nontrivial evaluation of the equality relation.

    As usual, the keyword argument ``evaluate=False`` can be used to
    prevent any evaluation.

    Examples
    ========

    >>> Eq(y, x + x**2)
    Eq(y, x**2 + x)
    >>> Eq(2, 5)
    false
    >>> Eq(2, 5, evaluate=False)
    Eq(2, 5)
    >>> _.doit()
    false
    >>> Eq(exp(x), exp(x).rewrite(cos))
    Eq(E**x, sinh(x) + cosh(x))
    >>> simplify(_)
    true

    See Also
    ========

    diofant.logic.boolalg.Equivalent : for representing equality between two
        boolean expressions

    Notes
    =====

    This class is not the same as the == operator.  The == operator tests
    for exact structural equality between two expressions; this class
    compares expressions mathematically.

    If either object defines an `_eval_Eq` method, it can be used in place of
    the default algorithm.  If `lhs._eval_Eq(rhs)` or `rhs._eval_Eq(lhs)`
    returns anything other than None, that return value will be substituted for
    the Equality.  If None is returned by `_eval_Eq`, an Equality object will
    be created as usual.

    """

    rel_op = '=='

    is_Equality = True

    def __new__(cls, lhs, rhs=0, **options):  # pylint: disable=signature-differs
        lhs = sympify(lhs, strict=True)
        rhs = sympify(rhs, strict=True)

        evaluate = options.pop('evaluate', global_evaluate[0])

        if evaluate:
            # If one expression has an _eval_Eq, return its results.
            if hasattr(lhs, '_eval_Eq'):
                r = lhs._eval_Eq(rhs)
                if r is not None:
                    return r
            if hasattr(rhs, '_eval_Eq'):
                r = rhs._eval_Eq(lhs)
                if r is not None:
                    return r
            # If expressions have the same structure, they must be equal.
            if lhs == rhs:
                return true
            if all(isinstance(i, BooleanAtom) for i in (rhs, lhs)):
                return false  # equal args already evaluated

            # If appropriate, check if the difference evaluates.  Detect
            # incompatibility such as lhs real and rhs not real.
            if isinstance(lhs, Expr) and isinstance(rhs, Expr):
                r = (lhs - rhs).is_zero
                if r is not None:
                    return sympify(r, strict=True)

        return Relational.__new__(cls, lhs, rhs, **options)


Eq = Equality


class Unequality(Relational):
    """An unequal relation between two objects.

    Represents that two objects are not equal.  If they can be shown to be
    definitively equal, this will reduce to False; if definitively unequal,
    this will reduce to True.  Otherwise, the relation is maintained as an
    Unequality object.

    Examples
    ========

    >>> Ne(y, x+x**2)
    Ne(y, x**2 + x)

    See Also
    ========
    Equality

    Notes
    =====
    This class is not the same as the != operator.  The != operator tests
    for exact structural equality between two expressions; this class
    compares expressions mathematically.

    This class is effectively the inverse of Equality.  As such, it uses the
    same algorithms, including any available `_eval_Eq` methods.

    """

    rel_op = '!='

    def __new__(cls, lhs, rhs=0, **options):  # pylint: disable=signature-differs
        lhs = sympify(lhs, strict=True)
        rhs = sympify(rhs, strict=True)

        evaluate = options.pop('evaluate', global_evaluate[0])

        if evaluate:
            is_equal = Equality(lhs, rhs)
            if is_equal in (true, false):
                return ~is_equal

        return Relational.__new__(cls, lhs, rhs, **options)


Ne = Unequality


class _Inequality(Relational):
    """Internal base class for all *Than types.

    Each subclass must implement _eval_relation to provide the method for
    comparing two real numbers.

    """

    def __new__(cls, lhs, rhs=0, **options):  # pylint: disable=signature-differs
        lhs = sympify(lhs, strict=True)
        rhs = sympify(rhs, strict=True)

        evaluate = options.pop('evaluate', global_evaluate[0])

        if evaluate:
            # First we invoke the appropriate inequality method of `lhs`
            # (e.g., `lhs.__lt__`).  That method will try to reduce to
            # boolean or raise an exception.  It may keep calling
            # superclasses until it reaches `Expr` (e.g., `Expr.__lt__`).
            # In some cases, `Expr` will just invoke us again (if neither it
            # nor a subclass was able to reduce to boolean or raise an
            # exception).  In that case, it must call us with
            # `evaluate=False` to prevent infinite recursion.
            return cls._eval_relation(lhs, rhs)

        # make a "non-evaluated" Expr for the inequality
        return Relational.__new__(cls, lhs, rhs, **options)


class _Greater(_Inequality):
    """Not intended for general use

    _Greater is only used so that GreaterThan and StrictGreaterThan may subclass
    it for the .gts and .lts properties.

    """

    @property
    def gts(self):
        """Greater than side argument."""
        return self.args[0]

    @property
    def lts(self):
        """Less than side argument."""
        return self.args[1]


class _Less(_Inequality):
    """Not intended for general use.

    _Less is only used so that LessThan and StrictLessThan may subclass it for
    the .gts and .lts properties.

    """

    @property
    def gts(self):
        """Greater than side argument."""
        return self.args[1]

    @property
    def lts(self):
        """Less than side argument."""
        return self.args[0]


class GreaterThan(_Greater):
    r"""Class representations of inequalities.

    The ``*Than`` classes represent unequal relationships, where the left-hand
    side is generally bigger or smaller than the right-hand side.  For example,
    the GreaterThan class represents an unequal relationship where the
    left-hand side is at least as big as the right side, if not bigger.  In
    mathematical notation:

    lhs >= rhs

    In total, there are four ``*Than`` classes, to represent the four
    inequalities:

    +-----------------+--------+
    |Class Name       | Symbol |
    +=================+========+
    |GreaterThan      | (>=)   |
    +-----------------+--------+
    |LessThan         | (<=)   |
    +-----------------+--------+
    |StrictGreaterThan| (>)    |
    +-----------------+--------+
    |StrictLessThan   | (<)    |
    +-----------------+--------+

    All classes take two arguments, lhs and rhs.

    +----------------------------+-----------------+
    |Signature Example           | Math equivalent |
    +============================+=================+
    |GreaterThan(lhs, rhs)       |   lhs >= rhs    |
    +----------------------------+-----------------+
    |LessThan(lhs, rhs)          |   lhs <= rhs    |
    +----------------------------+-----------------+
    |StrictGreaterThan(lhs, rhs) |   lhs >  rhs    |
    +----------------------------+-----------------+
    |StrictLessThan(lhs, rhs)    |   lhs <  rhs    |
    +----------------------------+-----------------+

    In addition to the normal .lhs and .rhs of Relations, ``*Than`` inequality
    objects also have the .lts and .gts properties, which represent the "less
    than side" and "greater than side" of the operator.  Use of .lts and .gts
    in an algorithm rather than .lhs and .rhs as an assumption of inequality
    direction will make more explicit the intent of a certain section of code,
    and will make it similarly more robust to client code changes:

    >>> e = GreaterThan(x, 1)
    >>> e
    x >= 1
    >>> f'{e.gts} >= {e.lts} is the same as {e.lts} <= {e.gts}'
    'x >= 1 is the same as 1 <= x'

    Examples
    ========

    One generally does not instantiate these classes directly, but uses various
    convenience methods:

    >>> e1 = Ge(x, 2)  # Ge is a convenience wrapper
    >>> print(e1)
    x >= 2

    >>> print(f'{Ge(x, 2)}\n{Gt(x, 2)}\n{Le(x, 2)}\n{Lt(x, 2)}')
    x >= 2
    x > 2
    x <= 2
    x < 2

    Another option is to use the Python inequality operators (>=, >, <=, <)
    directly.  Their main advantage over the Ge, Gt, Le, and Lt counterparts, is
    that one can write a more "mathematical looking" statement rather than
    littering the math with oddball function calls.  However there are certain
    (minor) caveats of which to be aware (search for 'gotcha', below).

    >>> e2 = x >= 2
    >>> print(e2)
    x >= 2
    >>> print(f'e1: {e1},    e2: {e2}')
    e1: x >= 2,    e2: x >= 2
    >>> e1 == e2
    True

    However, it is also perfectly valid to instantiate a ``*Than`` class less
    succinctly and less conveniently:

    >>> print(f"{Rel(x, 1, '>=')}\n{Relational(x, 1, '>=')}\n{GreaterThan(x, 1)}")
    x >= 1
    x >= 1
    x >= 1

    >>> print(f"{Rel(x, 1, '>')}\n{Relational(x, 1, '>')}\n{StrictGreaterThan(x, 1)}")
    x > 1
    x > 1
    x > 1

    >>> print(f"{Rel(x, 1, '<=')}\n{Relational(x, 1, '<=')}\n{LessThan(x, 1)}")
    x <= 1
    x <= 1
    x <= 1

    >>> print(f"{Rel(x, 1, '<')}\n{Relational(x, 1, '<')}\n{StrictLessThan(x, 1)}")
    x < 1
    x < 1
    x < 1

    Notes
    =====

    There are a couple of "gotchas" when using Python's operators.

    The first enters the mix when comparing against a literal number as the lhs
    argument.  Due to the order that Python decides to parse a statement, it may
    not immediately find two objects comparable.  For example, to evaluate the
    statement (1 < x), Python will first recognize the number 1 as a native
    number, and then that x is *not* a native number.  At this point, because a
    native Python number does not know how to compare itself with a Diofant object
    Python will try the reflective operation, (x > 1).  Unfortunately, there is
    no way available to Diofant to recognize this has happened, so the statement
    (1 < x) will turn silently into (x > 1).

    >>> e1 = x > 1
    >>> e2 = x >= 1
    >>> e3 = x < 1
    >>> e4 = x <= 1
    >>> e5 = 1 > x
    >>> e6 = 1 >= x
    >>> e7 = 1 < x
    >>> e8 = 1 <= x
    >>> print('%s     %s\n'*4 % (e1, e2, e3, e4, e5, e6, e7, e8))
    x > 1     x >= 1
    x < 1     x <= 1
    x < 1     x <= 1
    x > 1     x >= 1

    If the order of the statement is important (for visual output to the
    console, perhaps), one can work around this annoyance in a couple ways: (1)
    "sympify" the literal before comparison, (2) use one of the wrappers, or (3)
    use the less succinct methods described above:

    >>> e1 = Integer(1) > x
    >>> e2 = Integer(1) >= x
    >>> e3 = Integer(1) < x
    >>> e4 = Integer(1) <= x
    >>> e5 = Gt(1, x)
    >>> e6 = Ge(1, x)
    >>> e7 = Lt(1, x)
    >>> e8 = Le(1, x)
    >>> print('%s     %s\n'*4 % (e1, e2, e3, e4, e5, e6, e7, e8))
    1 > x     1 >= x
    1 < x     1 <= x
    1 > x     1 >= x
    1 < x     1 <= x

    The other gotcha is with chained inequalities.  Occasionally, one may be
    tempted to write statements like:

    >>> x < y < z
    Traceback (most recent call last):
    ...
    TypeError: symbolic boolean expression has no truth value.

    Due to an implementation detail or decision of Python, to create a
    chained inequality, the only method currently available is to make use of
    And:

    >>> And(x < y, y < z)
    (x < y) & (y < z)

    """

    rel_op = '>='

    @classmethod
    def _eval_relation(cls, lhs, rhs):
        return sympify(lhs >= rhs, strict=True)


Ge = GreaterThan


class LessThan(_Less):  # noqa: D101
    __doc__ = GreaterThan.__doc__

    rel_op = '<='

    @classmethod
    def _eval_relation(cls, lhs, rhs):
        return sympify(lhs <= rhs, strict=True)


Le = LessThan


class StrictGreaterThan(_Greater):  # noqa: D101
    __doc__ = GreaterThan.__doc__

    rel_op = '>'

    @classmethod
    def _eval_relation(cls, lhs, rhs):
        return sympify(lhs > rhs, strict=True)


Gt = StrictGreaterThan


class StrictLessThan(_Less):  # noqa: D101
    __doc__ = GreaterThan.__doc__

    rel_op = '<'

    @classmethod
    def _eval_relation(cls, lhs, rhs):
        return sympify(lhs < rhs, strict=True)


Lt = StrictLessThan


# A class-specific (not object-specific) data item used for a minor speedup.  It
# is defined here, rather than directly in the class, because the classes that
# it references have not been defined until now (e.g. StrictLessThan).
Relational.ValidRelationOperator = {
    None: Equality,
    '==': Equality,
    'eq': Equality,
    '!=': Unequality,
    '<>': Unequality,
    'ne': Unequality,
    '>=': GreaterThan,
    'ge': GreaterThan,
    '<=': LessThan,
    'le': LessThan,
    '>': StrictGreaterThan,
    'gt': StrictGreaterThan,
    '<': StrictLessThan,
    'lt': StrictLessThan,
}
