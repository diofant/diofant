"""
This module contains the machinery handling assumptions.

All symbolic objects have assumption attributes that can be accessed via
.is_<assumption name> attribute, i.e.
:py:attr:`~diofant.core.expr.Expr.is_integer`.  Full set of defined
assumption names are accessible as
``Expr.default_assumptions.rules.defined_facts`` attribute.

Assumptions determine certain properties of symbolic objects and can
have 3 possible values: True, False, None.  True is returned if the
object has the property and False is returned if it doesn't or can't
(i.e. doesn't make sense):

    >>> I.is_algebraic
    True
    >>> I.is_real
    False
    >>> I.is_prime
    False

When the property cannot be determined (or when a method is not
implemented) None will be returned, e.g. a generic symbol, x, may or
may not be positive so a value of None is returned for x.is_positive.

By default, all symbolic values are in the largest set in the given context
without specifying the property. For example, a symbol that has a property
being integer, is also real, complex, etc.

See Also
========

:py:class:`~diofant.core.numbers.ImaginaryUnit`
:py:attr:`~diofant.core.expr.Expr.is_algebraic`
:py:attr:`~diofant.core.expr.Expr.is_real`
:py:attr:`~diofant.core.expr.Expr.is_prime`

Notes
=====

Assumption values are stored in obj._assumptions dictionary or
are returned by getter methods (with property decorators) or are
attributes of objects/classes.

"""

from random import shuffle

from .facts import FactKB, FactRules
from .sympify import sympify


_assume_rules = FactRules([
    'integer        ->  rational',
    'rational       ->  real',
    'real           ==  extended_real & finite',
    'rational       ->  algebraic',
    'algebraic      ->  complex',
    'real           ->  complex',
    'imaginary      ->  complex',
    'complex        ->  finite & commutative',
    'extended_real  ->  commutative',
    'extended_real  ->  real | infinite',

    'odd            ==  integer & ~even',
    'even           ==  integer & ~odd',

    'extended_real  ==  negative | zero | positive',
    'transcendental ==  complex & ~algebraic',

    'negative       ==  nonpositive & nonzero',
    'positive       ==  nonnegative & nonzero',
    'zero           ==  nonnegative & nonpositive',

    'nonpositive    ==  extended_real & ~positive',
    'nonnegative    ==  extended_real & ~negative',

    'zero           ->  even',

    'prime          ->  integer & positive',
    'composite      ->  integer & positive & ~prime',

    'irrational     ==  real & ~rational',

    'imaginary      ->  ~extended_real | zero',

    'infinite       ->  ~finite',
    'noninteger     ==  real & ~integer',
    'nonzero        ==  ~zero',

    'polar          -> commutative',
])

_assume_defined = frozenset(_assume_rules.defined_facts.copy())
_assume_docs = {
    'commutative':
    """Test if self commutes with any other object wrt multiplication operation.""",
    'polar':
    """
Test if self can have values from the Riemann surface of the logarithm.

See Also
========

diofant.functions.elementary.complexes.polar_lift
diofant.functions.elementary.complexes.principal_branch
diofant.functions.elementary.exponential.exp_polar

""",
    'complex':
    """
Test if self can have only values from the set of complex numbers.

See Also
========

is_real

""",
    'real':
    """
Test if self can have only values from the set of real numbers.

See Also
========

is_complex

References
==========

* https://en.wikipedia.org/wiki/Real_number

""",
    'imaginary':
    """
Test if self is an imaginary number.

I.e. that it can be written as a real number multiplied by
the imaginary unit ``I``.

References
==========

* https://en.wikipedia.org/wiki/Imaginary_number

""",
    'extended_real':
    """
Test if self can have only values on the extended real number line.

See Also
========

is_real

References
==========

* https://en.wikipedia.org/wiki/Extended_real_number_line

""",
    'integer':
    """Test if self can have only values from the set of integers.""",
    'noninteger':
    """
Test if self can have only values from the subset of real numbers,
that aren't integers.

""",
    'odd':
    """
Test if self can have only values from the set of odd integers.

See Also
========

is_even

References
==========

* https://en.wikipedia.org/wiki/Parity_%28mathematics%29

""",
    'even':
    """
Test if self can have only values from the set of even integers.

See Also
========

is_odd

References
==========

* https://en.wikipedia.org/wiki/Parity_%28mathematics%29

""",
    'prime':
    """
Test if self is a natural number greater than ``1`` that has
no positive divisors other than ``1`` and itself.

References
==========

* https://en.wikipedia.org/wiki/Prime_number
""",
    'composite':
    """
Test if self is a positive integer that has at least one positive
divisor other than ``1`` or the number itself.

References
==========

* https://en.wikipedia.org/wiki/Composite_number

""",
    'zero':
    """
Test if self is zero.

See Also
========

is_nonzero

""",
    'nonzero':
    """
Test if self is nonzero.

See Also
========

is_zero

""",
    'rational':
    """Test if self can have only values from the set of rationals.""",
    'algebraic':
    """
Test if self can have only values from the set of algebraic numbers.

References
==========

* https://en.wikipedia.org/wiki/Algebraic_number

""",
    'transcendental':
    """
Test if self can have only values from the set of transcendental numbers.

References
==========

* https://en.wikipedia.org/wiki/Transcendental_number

""",
    'irrational':
    """
Test if self value cannot be represented exactly by Rational.

References
==========

* https://en.wikipedia.org/wiki/Irrational_number

""",
    'finite':
    """
Test if self absolute value is bounded.

References
==========

* https://en.wikipedia.org/wiki/Finite

""",
    'infinite':
    """
Test if self absolute value can be arbitrarily large.

References
==========

* :func:`math.isfinite`
* :obj:`numpy.isfinite`

""",
    'negative':
    """
Test if self can have only negative values.

References
==========

* https://en.wikipedia.org/wiki/Negative_number

""",
    'nonnegative':
    """
Test if self can have only nonnegative values.

See Also
========

is_negative

References
==========

* https://en.wikipedia.org/wiki/Negative_number

""",
    'positive':
    """Test if self can have only positive values.""",
    'nonpositive':
    """Test if self can have only nonpositive values.""",
}


class StdFactKB(FactKB):
    """A FactKB specialised for the built-in rules

    This is the only kind of FactKB that Basic objects should use.

    """

    rules = _assume_rules

    def __init__(self, facts=None):
        # save a copy of the facts dict
        if not facts:
            self._generator = {}
        elif not isinstance(facts, FactKB):
            self._generator = facts.copy()
        else:
            self._generator = facts.generator
        if facts:
            self.deduce_all_facts(facts)

    def copy(self):
        return self.__class__(self)

    @property
    def generator(self):
        return self._generator.copy()


def as_property(fact):
    """Convert a fact name to the name of the corresponding property."""
    return f'is_{fact}'


def make_property(fact):
    """Create the automagic property corresponding to a fact."""

    def getit(self):
        try:
            return self._assumptions[fact]
        except KeyError:
            if self._assumptions is self.default_assumptions:
                self._assumptions = self.default_assumptions.copy()
            return _ask(fact, self)

    getit.func_name = as_property(fact)
    return property(getit, doc=_assume_docs[fact])


def check_assumptions(expr, **assumptions):
    r"""Checks if expression `expr` satisfies all assumptions.

    Parameters
    ==========

    expr : Expr
        Expression to test.
    \*\*assumptions : dict
        Keyword arguments to specify assumptions to test.

    Returns
    =======

    True, False or None (if can't conclude).

    Examples
    ========

    >>> check_assumptions(-5, integer=True)
    True
    >>> x = Symbol('x', positive=True)
    >>> check_assumptions(2*x + 1, negative=True)
    False
    >>> check_assumptions(z, real=True) is None
    True

    """
    expr = sympify(expr)

    result = True
    for key, expected in assumptions.items():
        if expected is None:
            continue
        test = getattr(expr, 'is_' + key, None)
        if test is expected:
            continue
        elif test is not None:
            return False
        else:
            result = None
    return result


def _ask(fact, obj):
    """
    Find the truth value for a property of an object.

    This function is called when a request is made to see what a fact
    value is.

    For this we use several techniques:

    First, the fact-evaluation function is tried, if it exists (for
    example _eval_is_integer). Then we try related facts. For example

        rational   -->   integer

    another example is joined rule:

        integer & ~odd  --> even

    so in the latter case if we are looking at what 'even' value is,
    'integer' and 'odd' facts will be asked.

    In all cases, when we settle on some fact value, its implications are
    deduced, and the result is cached in ._assumptions.

    """
    assumptions = obj._assumptions
    handler_map = obj._prop_handler

    # Store None into the assumptions so that recursive attempts at
    # evaluating the same fact don't trigger infinite recursion.
    assumptions._tell(fact, None)

    # First try the assumption evaluation function if it exists
    try:
        evaluate = handler_map[fact]
    except KeyError:
        pass
    else:
        a = evaluate(obj)
        if a is not None:
            assumptions.deduce_all_facts(((fact, a),))
            return a

    # Try assumption's prerequisites
    prereq = list(_assume_rules.prereq[fact])
    shuffle(prereq)
    for pk in prereq:
        if pk in assumptions:
            continue
        if pk in handler_map:
            _ask(pk, obj)

            # we might have found the value of fact
            ret_val = assumptions.get(fact)
            if ret_val is not None:
                return ret_val


class ManagedProperties(type):
    """Metaclass for classes with old-style assumptions."""

    def __init__(cls, *args, **kws):
        local_defs = {}
        for k in _assume_defined:
            attrname = as_property(k)
            v = cls.__dict__.get(attrname, '')
            if isinstance(v, int):
                v = bool(v)
                local_defs[k] = v

        defs = {}
        for base in reversed(cls.__bases__):
            try:
                defs.update(base._explicit_class_assumptions)
            except AttributeError:
                pass
        defs.update(local_defs)

        cls._explicit_class_assumptions = defs
        cls.default_assumptions = StdFactKB(defs)

        cls._prop_handler = {}
        for k in _assume_defined:
            try:
                cls._prop_handler[k] = getattr(cls, f'_eval_is_{k}')
            except AttributeError:
                pass

        # Put definite results directly into the class dict, for speed
        for k, v in cls.default_assumptions.items():
            setattr(cls, as_property(k), v)

        # protection e.g. for Integer.is_even=F <- (Rational.is_integer=F)
        derived_from_bases = set()
        for base in cls.__bases__:
            try:
                derived_from_bases |= set(base.default_assumptions)
            except AttributeError:
                continue  # not an assumption-aware class
        for fact in derived_from_bases - set(cls.default_assumptions):
            pname = as_property(fact)
            setattr(cls, pname, make_property(fact))

        # Finally, add any missing automagic property (e.g. for Basic)
        for fact in _assume_defined:
            pname = as_property(fact)
            if not hasattr(cls, pname):
                setattr(cls, pname, make_property(fact))
