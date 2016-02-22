"""
This module contains the machinery handling assumptions.

All symbolic objects have assumption attributes that can be accessed via
.is_<assumption name> attribute.

Assumptions determine certain properties of symbolic objects and can
have 3 possible values: True, False, None.  True is returned if the
object has the property and False is returned if it doesn't or can't
(i.e. doesn't make sense):

    >>> from sympy import I
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

Here follows a list of possible assumption names:

.. glossary::

    commutative
        object commutes with any other object with
        respect to multiplication operation.

    complex
        object can have only values from the set
        of complex numbers.

    imaginary
        object value is a number that can be written as a real
        number multiplied by the imaginary unit ``I``.  See [3]_.

    real
        object can have only values from the set
        of real numbers [12]_.

    extended_real
        object can have only values on the extended
        real number line [13]_.

    integer
        object can have only values from the set
        of integers.

    odd
    even
        object can have only values from the set of
        odd (even) integers [2]_.

    prime
        object is a natural number greater than ``1`` that has
        no positive divisors other than ``1`` and itself.  See [6]_.

    composite
        object is a positive integer that has at least one positive
        divisor other than ``1`` or the number itself.  See [4]_.

    zero
    nonzero
        object is zero (not zero).

    rational
        object can have only values from the set
        of rationals.

    algebraic
        object can have only values from the set
        of algebraic numbers [11]_.

    transcendental
        object can have only values from the set
        of transcendental numbers [10]_.

    irrational
        object value cannot be represented exactly by Rational, see [5]_.

    finite
    infinite
        object absolute value is bounded (is value is
        arbitrarily large).  See [7]_, [8]_, [9]_.

    negative
    nonnegative
        object can have only negative (only
        nonnegative) values [1]_.

    positive
    nonpositive
        object can have only positive (only
        nonpositive) values.

    hermitian
    antihermitian
        object belongs to the field of hermitian
        (antihermitian) operators.

Examples
========

    >>> from sympy import Symbol
    >>> x = Symbol('x', real=True); x
    x
    >>> x.is_extended_real
    True
    >>> x.is_complex
    True

See Also
========

.. seealso::

    :py:class:`sympy.core.numbers.ImaginaryUnit`
    :py:class:`sympy.core.numbers.Zero`
    :py:class:`sympy.core.numbers.One`

Notes
=====

Assumption values are stored in obj._assumptions dictionary or
are returned by getter methods (with property decorators) or are
attributes of objects/classes.


References
==========

.. [1] http://en.wikipedia.org/wiki/Negative_number
.. [2] http://en.wikipedia.org/wiki/Parity_%28mathematics%29
.. [3] http://en.wikipedia.org/wiki/Imaginary_number
.. [4] http://en.wikipedia.org/wiki/Composite_number
.. [5] http://en.wikipedia.org/wiki/Irrational_number
.. [6] http://en.wikipedia.org/wiki/Prime_number
.. [7] http://en.wikipedia.org/wiki/Finite
.. [8] https://docs.python.org/3/library/math.html#math.isfinite
.. [9] http://docs.scipy.org/doc/numpy/reference/generated/numpy.isfinite.html
.. [10] http://en.wikipedia.org/wiki/Transcendental_number
.. [11] http://en.wikipedia.org/wiki/Algebraic_number
.. [12] http://en.wikipedia.org/wiki/Real_number
.. [13] http://en.wikipedia.org/wiki/Extended_real_number_line

"""

from random import shuffle

from .facts import FactRules, FactKB


_assume_rules = FactRules([
    'integer        ->  rational',
    'rational       ->  real',
    'real           ==  extended_real & finite',
    'rational       ->  algebraic',
    'algebraic      ->  complex',
    'real           ->  complex & hermitian',
    'imaginary      ->  complex & antihermitian',
    'complex        ->  finite & commutative',
    'extended_real  ->  commutative',

    'odd            ==  integer & ~even',
    'even           ==  integer & ~odd',

    'extended_real  ==  negative | zero | positive',
    'transcendental ==  complex & ~algebraic',

    'negative       ==  nonpositive & nonzero',
    'positive       ==  nonnegative & nonzero',
    'zero           ==  nonnegative & nonpositive',

    'nonpositive    ==  extended_real & ~positive',
    'nonnegative    ==  extended_real & ~negative',

    'zero           ->  even & finite',

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
    """Convert a fact name to the name of the corresponding property"""
    return 'is_%s' % fact


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
    return property(getit)


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

    # Note: the result has already been cached
    return


class ManagedProperties(type):
    """Metaclass for classes with old-style assumptions"""

    def __init__(self, *args, **kws):
        local_defs = {}
        for k in _assume_defined:
            attrname = as_property(k)
            v = self.__dict__.get(attrname, '')
            if isinstance(v, (bool, int, type(None))):
                if v is not None:
                    v = bool(v)
                local_defs[k] = v

        defs = {}
        for base in reversed(self.__bases__):
            try:
                defs.update(base._explicit_class_assumptions)
            except AttributeError:
                pass
        defs.update(local_defs)

        self._explicit_class_assumptions = defs
        self.default_assumptions = StdFactKB(defs)

        self._prop_handler = {}
        for k in _assume_defined:
            try:
                self._prop_handler[k] = getattr(self, '_eval_is_%s' % k)
            except AttributeError:
                pass

        # Put definite results directly into the class dict, for speed
        for k, v in self.default_assumptions.items():
            setattr(self, as_property(k), v)

        # protection e.g. for Integer.is_even=F <- (Rational.is_integer=F)
        derived_from_bases = set()
        for base in self.__bases__:
            try:
                derived_from_bases |= set(base.default_assumptions)
            except AttributeError:
                continue  # not an assumption-aware class
        for fact in derived_from_bases - set(self.default_assumptions):
            pname = as_property(fact)
            if pname not in self.__dict__:
                setattr(self, pname, make_property(fact))

        # Finally, add any missing automagic property (e.g. for Basic)
        for fact in _assume_defined:
            pname = as_property(fact)
            if not hasattr(self, pname):
                setattr(self, pname, make_property(fact))
