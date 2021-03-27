"""Generic set theory interfaces."""

import itertools
import typing

from mpmath import mpf, mpi

from ..core import Basic, Eq, Expr, Mul, S, nan, oo, zoo
from ..core.compatibility import iterable
from ..core.decorators import _sympifyit
from ..core.evalf import EvalfMixin
from ..core.evaluate import global_evaluate
from ..core.singleton import Singleton
from ..core.sympify import sympify
from ..logic import And, Not, Or, false, true
from ..utilities import ordered, subsets
from .contains import Contains


class Set(Basic):
    """
    The base class for any kind of set.

    This is not meant to be used directly as a container of items. It does not
    behave like the builtin ``set``; see :class:`FiniteSet` for that.

    Real intervals are represented by the :class:`Interval` class and unions of
    sets by the :class:`Union` class. The empty set is represented by the
    :class:`EmptySet` class and available as a singleton as ``S.EmptySet``.

    """

    is_number = False
    is_iterable = False
    is_interval = False

    is_FiniteSet = False
    is_Interval = False
    is_ProductSet = False
    is_Union = False
    is_Intersection: typing.Optional[bool] = None
    is_EmptySet: typing.Optional[bool] = None
    is_UniversalSet: typing.Optional[bool] = None
    is_Complement: typing.Optional[bool] = None
    is_SymmetricDifference: typing.Optional[bool] = None

    @staticmethod
    def _infimum_key(expr):
        """Return infimum (if possible) else oo."""
        try:
            infimum = expr.inf
            assert infimum.is_comparable
        except (NotImplementedError,
                AttributeError, AssertionError, ValueError):
            infimum = oo
        return infimum

    def union(self, other):
        """
        Returns the union of 'self' and 'other'.

        Examples
        ========

        As a shortcut it is possible to use the '+' operator:

        >>> Interval(0, 1).union(Interval(2, 3))
        [0, 1] U [2, 3]
        >>> Interval(0, 1) + Interval(2, 3)
        [0, 1] U [2, 3]
        >>> Interval(1, 2, True, True) + FiniteSet(2, 3)
        (1, 2] U {3}

        Similarly it is possible to use the '-' operator for set differences:

        >>> Interval(0, 2) - Interval(0, 1)
        (1, 2]
        >>> Interval(1, 3) - FiniteSet(2)
        [1, 2) U (2, 3]

        """
        return Union(self, other)

    def intersection(self, other):
        """
        Returns the intersection of 'self' and 'other'.


        >>> Interval(1, 3).intersection(Interval(1, 2))
        [1, 2]

        """
        return Intersection(self, other)

    def _intersection(self, other):
        """
        This function should only be used internally

        self._intersection(other) returns a new, intersected set if self knows how
        to intersect itself with other, otherwise it returns ``None``

        When making a new set class you can be assured that other will not
        be a :class:`Union`, :class:`FiniteSet`, or :class:`EmptySet`

        Used within the :class:`Intersection` class

        """
        return

    def is_disjoint(self, other):
        """
        Returns True if 'self' and 'other' are disjoint

        Examples
        ========

        >>> Interval(0, 2).is_disjoint(Interval(1, 2))
        False
        >>> Interval(0, 2).is_disjoint(Interval(3, 4))
        True

        References
        ==========

        * https://en.wikipedia.org/wiki/Disjoint_sets

        """
        return self.intersection(other) == S.EmptySet

    def isdisjoint(self, other):
        """Alias for :meth:`is_disjoint()`."""
        return self.is_disjoint(other)

    def _union(self, other):
        """
        This function should only be used internally

        self._union(other) returns a new, joined set if self knows how
        to join itself with other, otherwise it returns ``None``.
        It may also return a python set of Diofant Sets if they are somehow
        simpler. If it does this it must be idempotent i.e. the sets returned
        must return ``None`` with _union'ed with each other

        Used within the :class:`Union` class

        """
        return

    def complement(self, universe):
        r"""
        The complement of 'self' w.r.t the given the universe.

        Examples
        ========

        >>> Interval(0, 1).complement(S.Reals)
        (-oo, 0) U (1, oo)

        >>> Interval(0, 1).complement(S.UniversalSet)
        UniversalSet() \ [0, 1]

        """
        return Complement(universe, self)

    def _complement(self, other):
        # this behaves as other - self
        if other.is_subset(self):
            return S.EmptySet
        elif isinstance(other, ProductSet):
            # For each set consider it or it's complement
            # We need at least one of the sets to be complemented
            # Consider all 2^n combinations.
            # We can conveniently represent these options easily using a
            # ProductSet

            # XXX: this doesn't work if the dimentions of the sets isn't same.
            # A - B is essentially same as A if B has a different
            # dimensionality than A
            switch_sets = ProductSet(FiniteSet(o, o - s) for s, o in
                                     zip(self.sets, other.sets))
            product_sets = (ProductSet(*set) for set in switch_sets)
            # Union of all combinations but this one
            return Union(p for p in product_sets if p != other)

        elif isinstance(other, Interval):
            if isinstance(self, Interval) or isinstance(self, FiniteSet):
                return Intersection(other, self.complement(S.Reals))

        elif isinstance(other, Union):
            return Union(o - self for o in other.args)

        elif isinstance(other, Complement):
            return Complement(other.args[0], Union(other.args[1], self), evaluate=False)

        elif isinstance(other, FiniteSet):
            unks = FiniteSet(*[el for el in other
                               if self.contains(el) not in [true, false]])
            other = FiniteSet(*[el for el in other
                                if self.contains(el) != true])
            ret = FiniteSet(*[el for el in other if self.contains(el) == false])
            if unks:
                ret |= Complement(FiniteSet(*unks), self, evaluate=False)
            return ret

    def symmetric_difference(self, other):
        """
        Returns symmetric difference of ``self`` and ``other``.

        Examples
        ========

        >>> Interval(1, 3).symmetric_difference(Reals)
        (-oo, 1) U (3, oo)

        References
        ==========

        * https://en.wikipedia.org/wiki/Symmetric_difference

        """
        return SymmetricDifference(self, other)

    def _symmetric_difference(self, other):
        return Union(Complement(self, other), Complement(other, self))

    @property
    def inf(self):
        """
        The infimum of 'self'

        Examples
        ========

        >>> Interval(0, 1).inf
        0
        >>> Union(Interval(0, 1), Interval(2, 3)).inf
        0

        """
        return self._inf

    @property
    def _inf(self):
        raise NotImplementedError(f'({self})._inf')

    @property
    def sup(self):
        """
        The supremum of 'self'

        Examples
        ========

        >>> Interval(0, 1).sup
        1
        >>> Union(Interval(0, 1), Interval(2, 3)).sup
        3

        """
        return self._sup

    @property
    def _sup(self):
        raise NotImplementedError(f'({self})._sup')

    @_sympifyit('other', false)
    def contains(self, other):
        """
        Returns True if 'other' is contained in 'self' as an element.

        As a shortcut it is possible to use the 'in' operator:

        Examples
        ========

        >>> Interval(0, 1).contains(0.5)
        true
        >>> 0.5 in Interval(0, 1)
        True

        """
        ret = self._contains(other)
        if ret is None:
            ret = Contains(other, self, evaluate=False)
        return ret

    def _contains(self, other):
        raise NotImplementedError(f'({self})._contains({other})')

    def is_subset(self, other):
        """
        Returns True if 'self' is a subset of 'other'.

        Examples
        ========

        >>> Interval(0, 0.5).is_subset(Interval(0, 1))
        True
        >>> Interval(0, 1).is_subset(Interval(0, 1, left_open=True))
        False

        """
        if isinstance(other, Set):
            return self.intersection(other) == self
        else:
            raise ValueError(f"Unknown argument '{other}'")

    def issubset(self, other):
        """Alias for :meth:`is_subset()`."""
        return self.is_subset(other)

    def is_proper_subset(self, other):
        """
        Returns True if 'self' is a proper subset of 'other'.

        Examples
        ========

        >>> Interval(0, 0.5).is_proper_subset(Interval(0, 1))
        True
        >>> Interval(0, 1).is_proper_subset(Interval(0, 1))
        False

        """
        if isinstance(other, Set):
            return self != other and self.is_subset(other)
        else:
            raise ValueError(f"Unknown argument '{other}'")

    def is_superset(self, other):
        """
        Returns True if 'self' is a superset of 'other'.

        Examples
        ========

        >>> Interval(0, 0.5).is_superset(Interval(0, 1))
        False
        >>> Interval(0, 1).is_superset(Interval(0, 1, left_open=True))
        True

        """
        if isinstance(other, Set):
            return other.is_subset(self)
        else:
            raise ValueError(f"Unknown argument '{other}'")

    def issuperset(self, other):
        """Alias for :meth:`is_superset()`."""
        return self.is_superset(other)

    def is_proper_superset(self, other):
        """
        Returns True if 'self' is a proper superset of 'other'.

        Examples
        ========

        >>> Interval(0, 1).is_proper_superset(Interval(0, 0.5))
        True
        >>> Interval(0, 1).is_proper_superset(Interval(0, 1))
        False

        """
        if isinstance(other, Set):
            return self != other and self.is_superset(other)
        else:
            raise ValueError(f"Unknown argument '{other}'")

    def _eval_powerset(self):
        raise NotImplementedError(f'Power set not defined for: {self.func}')

    def powerset(self):
        """
        Find the Power set of 'self'.

        Examples
        ========

        >>> A = EmptySet()
        >>> A.powerset()
        {EmptySet()}
        >>> A = FiniteSet(1, 2)
        >>> a, b, c = FiniteSet(1), FiniteSet(2), FiniteSet(1, 2)
        >>> A.powerset() == FiniteSet(a, b, c, EmptySet())
        True

        References
        ==========

        * https://en.wikipedia.org/wiki/Power_set

        """
        return self._eval_powerset()

    @property
    def measure(self):
        """
        The (Lebesgue) measure of 'self'

        Examples
        ========

        >>> Interval(0, 1).measure
        1
        >>> Union(Interval(0, 1), Interval(2, 3)).measure
        2

        """
        return self._measure

    @property
    def boundary(self):
        """
        The boundary or frontier of a set

        A point x is on the boundary of a set S if

        1.  x is in the closure of S.
            I.e. Every neighborhood of x contains a point in S.
        2.  x is not in the interior of S.
            I.e. There does not exist an open set centered on x contained
            entirely within S.

        There are the points on the outer rim of S.  If S is open then these
        points need not actually be contained within S.

        For example, the boundary of an interval is its start and end points.
        This is true regardless of whether or not the interval is open.

        Examples
        ========

        >>> Interval(0, 1).boundary
        {0, 1}
        >>> Interval(0, 1, True, False).boundary
        {0, 1}

        """
        return self._boundary

    @property
    def is_open(self):
        """
        Test if a set is open.

        A set is open if it has an empty intersection with its boundary.

        Examples
        ========

        >>> S.Reals.is_open
        True

        See Also
        ========

        boundary

        """
        if not Intersection(self, self.boundary):
            return True

    @property
    def is_closed(self):
        """
        Test if a set is closed.

        Examples
        ========

        >>> Interval(0, 1).is_closed
        True

        """
        return self.boundary.is_subset(self)

    @property
    def closure(self):
        """
        Return the closure of a set.

        Examples
        ========

        >>> Interval(0, 1, right_open=True).closure
        [0, 1]

        """
        return self + self.boundary

    @property
    def interior(self):
        """
        Return the interior of a set.

        The interior of a set consists all points of a set that do not
        belong to its boundary.

        Examples
        ========

        >>> Interval(0, 1).interior
        (0, 1)
        >>> Interval(0, 1).boundary.interior
        EmptySet()

        """
        return self - self.boundary

    @property
    def _boundary(self):
        raise NotImplementedError()

    def _eval_imageset(self, f):
        from .fancysets import ImageSet
        return ImageSet(f, self)

    @property
    def _measure(self):
        raise NotImplementedError(f'({self})._measure')

    def __add__(self, other):
        return self.union(other)

    def __or__(self, other):
        return self.union(other)

    def __and__(self, other):
        return self.intersection(other)

    def __mul__(self, other):
        return ProductSet(self, other)

    def __xor__(self, other):
        return SymmetricDifference(self, other)

    def __pow__(self, exp):
        if not sympify(exp).is_Integer or exp < 0:
            raise ValueError(f'{exp}: Exponent must be a positive Integer')
        return ProductSet([self]*exp)

    def __sub__(self, other):
        return Complement(self, other)

    def __contains__(self, other):
        symb = self.contains(other)
        if symb not in (true, false):
            raise TypeError(f'contains did not evaluate to a bool: {symb!r}')
        return bool(symb)


class ProductSet(Set):
    """
    Represents a Cartesian Product of Sets.

    Returns a Cartesian product given several sets as either an iterable
    or individual arguments.

    Can use '*' operator on any sets for convenient shorthand.

    Examples
    ========

    >>> I = Interval(0, 5)
    >>> S = FiniteSet(1, 2, 3)
    >>> ProductSet(I, S)
    [0, 5] x {1, 2, 3}

    >>> (2, 2) in ProductSet(I, S)
    True

    >>> Interval(0, 1) * Interval(0, 1)  # The unit square
    [0, 1] x [0, 1]

    >>> H, T = Symbol('H'), Symbol('T')
    >>> coin = FiniteSet(H, T)
    >>> set(coin**2)
    {(H, H), (H, T), (T, H), (T, T)}

    Notes
    =====

    - Passes most operations down to the argument sets
    - Flattens Products of ProductSets

    References
    ==========

    * https://en.wikipedia.org/wiki/Cartesian_product

    """

    is_ProductSet = True

    def __new__(cls, *sets, **assumptions):
        def flatten(arg):
            if isinstance(arg, Set):
                if arg.is_ProductSet:
                    return sum(map(flatten, arg.args), [])
                else:
                    return [arg]
            elif iterable(arg):
                return sum(map(flatten, arg), [])
            raise TypeError('Input must be Sets or iterables of Sets')
        sets = flatten(list(sets))

        if EmptySet() in sets or len(sets) == 0:
            return EmptySet()

        if len(sets) == 1:
            return sets[0]

        return Basic.__new__(cls, *sets, **assumptions)

    def _eval_Eq(self, other):
        if not other.is_ProductSet:
            return

        if len(self.args) != len(other.args):
            return false

        return And(*(Eq(x, y) for x, y in zip(self.args, other.args)))

    def _contains(self, element):
        """
        'in' operator for ProductSets

        Examples
        ========

        >>> (2, 3) in Interval(0, 5) * Interval(0, 5)
        True

        >>> (10, 10) in Interval(0, 5) * Interval(0, 5)
        False

        Passes operation on to constituent sets

        """
        try:
            if len(element) != len(self.args):
                return false
        except TypeError:  # maybe element isn't an iterable
            return false
        return And(*[s.contains(i) for s, i in zip(self.sets, element)])

    def _intersection(self, other):
        """
        This function should only be used internally

        See Set._intersection for docstring

        """
        if not other.is_ProductSet:
            return
        if len(other.args) != len(self.args):
            return S.EmptySet
        return ProductSet(a.intersection(b) for a, b in zip(self.sets, other.sets))

    def _union(self, other):
        if not other.is_ProductSet:
            return
        if len(other.args) != len(self.args):
            return
        if self.args[0] == other.args[0]:
            return self.args[0] * Union(ProductSet(self.args[1:]),
                                        ProductSet(other.args[1:]))
        if self.args[-1] == other.args[-1]:
            return Union(ProductSet(self.args[:-1]),
                         ProductSet(other.args[:-1])) * self.args[-1]

    @property
    def sets(self):
        return self.args

    @property
    def _boundary(self):
        return Union(ProductSet(b + b.boundary if i != j else b.boundary
                                for j, b in enumerate(self.sets))
                     for i, a in enumerate(self.sets))

    @property
    def is_iterable(self):
        return all(set.is_iterable for set in self.sets)

    def __iter__(self):
        from ..utilities.iterables import cantor_product
        if self.is_iterable:
            return cantor_product(*self.sets)
        else:
            raise TypeError('Not all constituent sets are iterable')

    @property
    def _measure(self):
        measure = 1
        for set in self.sets:
            measure *= set.measure
        return measure

    def __len__(self):
        return int(Mul(*[len(s) for s in self.args]))


class Interval(Set, EvalfMixin):
    """
    Represents a real interval as a Set.

    Returns an interval with end points "start" and "end".

    For left_open=True (default left_open is False) the interval
    will be open on the left. Similarly, for right_open=True the interval
    will be open on the right.

    Examples
    ========

    >>> Interval(0, 1)
    [0, 1]
    >>> Interval(0, 1, False, True)
    [0, 1)
    >>> Interval.Ropen(0, 1)
    [0, 1)
    >>> Interval.Lopen(0, 1)
    (0, 1]
    >>> Interval.open(0, 1)
    (0, 1)

    >>> a = Symbol('a', real=True)
    >>> Interval(0, a)
    [0, a]

    Notes
    =====

    - Only real end points are supported
    - Interval(a, b) with a > b will return the empty set
    - Use the evalf() method to turn an Interval into an mpmath
      'mpi' interval instance

    References
    ==========

    * https://en.wikipedia.org/wiki/Interval_%28mathematics%29

    """

    is_Interval = True

    def __new__(cls, start, end, left_open=False, right_open=False):

        start = sympify(start, strict=True)
        end = sympify(end, strict=True)
        left_open = sympify(left_open, strict=True)
        right_open = sympify(right_open, strict=True)

        if not all(isinstance(a, (type(true), type(false)))
                   for a in [left_open, right_open]):
            raise NotImplementedError(
                'left_open and right_open can have only true/false values, '
                f'got {left_open} and {right_open}')

        if not all(i.is_extended_real is not False for i in (start, end)):
            raise ValueError('Non-real intervals are not supported')

        if (end - start).is_negative:
            return S.EmptySet

        is_open = left_open or right_open

        if end == start and is_open:
            return S.EmptySet
        if end == start and not is_open:
            return FiniteSet(end)

        return Basic.__new__(cls, start, end, left_open, right_open)

    @property
    def start(self):
        """
        The left end point of 'self'.

        This property takes the same value as the 'inf' property.

        Examples
        ========

        >>> Interval(0, 1).start
        0

        """
        return self.args[0]

    _inf = left = start

    @classmethod
    def open(cls, a, b):
        """Return an interval including neither boundary."""
        return cls(a, b, True, True)

    @classmethod
    def Lopen(cls, a, b):
        """Return an interval not including the left boundary."""
        return cls(a, b, True, False)

    @classmethod
    def Ropen(cls, a, b):
        """Return an interval not including the right boundary."""
        return cls(a, b, False, True)

    @property
    def end(self):
        """
        The right end point of 'self'.

        This property takes the same value as the 'sup' property.

        Examples
        ========

        >>> Interval(0, 1).end
        1

        """
        return self.args[1]

    _sup = right = end

    @property
    def left_open(self):
        """
        True if 'self' is left-open.

        Examples
        ========

        >>> Interval(0, 1, left_open=True).left_open
        true
        >>> Interval(0, 1, left_open=False).left_open
        false

        """
        return self.args[2]

    @property
    def right_open(self):
        """
        True if 'self' is right-open.

        Examples
        ========

        >>> Interval(0, 1, right_open=True).right_open
        true
        >>> Interval(0, 1, right_open=False).right_open
        false

        """
        return self.args[3]

    def _intersection(self, other):
        """
        This function should only be used internally

        See Set._intersection for docstring

        """
        # We only know how to intersect with other intervals
        if not other.is_Interval:
            return

        # handle unbounded self
        if Eq(self, S.Reals) == true and all((abs(_) < oo) is true or
                                             abs(_) == oo
                                             for _ in other.boundary):
            if other.is_left_unbounded and not other.left_open:
                other = Interval(other.start, other.end, True, other.right_open)
            if other.is_right_unbounded and not other.right_open:
                other = Interval(other.start, other.end, other.left_open, True)
            return other
        elif Eq(self, S.ExtendedReals) == true:
            return other

        # We can't intersect [0,3] with [x,6] -- we don't know if x>0 or x<0
        if not self._is_comparable(other):
            return

        empty = False

        if self.start <= other.end and other.start <= self.end:
            # Get topology right.
            if self.start < other.start:
                start = other.start
                left_open = other.left_open
            elif self.start > other.start:
                start = self.start
                left_open = self.left_open
            else:
                start = self.start
                left_open = self.left_open or other.left_open

            if self.end < other.end:
                end = self.end
                right_open = self.right_open
            elif self.end > other.end:
                end = other.end
                right_open = other.right_open
            else:
                end = self.end
                right_open = self.right_open or other.right_open

            if end - start == 0 and (left_open or right_open):
                empty = True
        else:
            empty = True

        if empty:
            return S.EmptySet

        return Interval(start, end, left_open, right_open)

    def _complement(self, other):
        if other in (S.Reals, S.ExtendedReals):
            a = Interval(-oo, self.start,
                         other.left_open, not self.left_open)
            b = Interval(self.end, oo, not self.right_open, other.right_open)
            return Union(a, b)
        return Set._complement(self, other)

    def _union(self, other):
        """
        This function should only be used internally

        See Set._union for docstring

        """
        if other.is_UniversalSet:
            return S.UniversalSet
        if other.is_Interval and self._is_comparable(other):
            from ..functions import Max, Min

            # Non-overlapping intervals
            end = Min(self.end, other.end)
            start = Max(self.start, other.start)
            if (end < start or
                    (end == start and (end not in self and end not in other))):
                return
            else:
                start = Min(self.start, other.start)
                end = Max(self.end, other.end)

                left_open = ((self.start != start or self.left_open) and
                             (other.start != start or other.left_open))
                right_open = ((self.end != end or self.right_open) and
                              (other.end != end or other.right_open))

                return Interval(start, end, left_open, right_open)

        # If I have open end points and these endpoints are contained in other
        if ((self.left_open and other.contains(self.start) is true) or
                (self.right_open and other.contains(self.end) is true)):
            # Fill in my end points and return
            open_left = self.left_open and self.start not in other
            open_right = self.right_open and self.end not in other
            new_self = Interval(self.start, self.end, open_left, open_right)
            return {new_self, other}

    @property
    def _boundary(self):
        return FiniteSet(self.start, self.end)

    def _contains(self, other):
        if not isinstance(other, Expr) or other in (nan, zoo):
            return false

        if other.is_extended_real is False:
            return false

        if self.left_open:
            expr = other > self.start
        else:
            expr = other >= self.start

            if other == self.start:
                return true

        if self.right_open:
            expr = And(expr, other < self.end)
        else:
            expr = And(expr, other <= self.end)

            if other == self.end:
                return true

        return sympify(expr, strict=True)

    def _eval_imageset(self, f):
        from ..calculus.singularities import singularities
        from ..core import Lambda, diff
        from ..functions import Max, Min
        from ..series import limit
        from ..solvers import solve

        # TODO: handle functions with infinitely many solutions (eg, sin, tan)
        # TODO: handle multivariate functions

        expr = f.expr
        if len(expr.free_symbols) > 1 or len(f.variables) != 1:
            return
        var = f.variables[0]

        if expr.is_Piecewise:
            result = S.EmptySet
            domain_set = self
            for (p_expr, p_cond) in expr.args:
                if p_cond == true:
                    intrvl = domain_set
                else:
                    intrvl = p_cond.as_set()
                    intrvl = Intersection(domain_set, intrvl)

                if p_expr.is_Number:
                    image = FiniteSet(p_expr)
                else:
                    image = imageset(Lambda(var, p_expr), intrvl)
                result = Union(result, image)

                # remove the part which has been `imaged`
                domain_set = Complement(domain_set, intrvl)
                if (p_expr, p_cond) != expr.args[-1] and domain_set.is_EmptySet:
                    break
            return result

        try:
            sing = [x for x in singularities(expr, var)
                    if x.is_extended_real and x in self]
        except NotImplementedError:
            return

        if self.left_open:
            _start = limit(expr, var, self.start, dir='+')
        elif self.start not in sing:
            _start = f(self.start)
        if self.right_open:
            _end = limit(expr, var, self.end, dir='-')
        elif self.end not in sing:
            _end = f(self.end)

        if len(sing) == 0:
            solns = solve(diff(expr, var), var)

            extr = [_start, _end] + [f(x[var]) for x in solns
                                     if (x[var].is_extended_real and
                                         x[var] in self)]
            start, end = Min(*extr), Max(*extr)

            left_open, right_open = False, False
            if _start <= _end:
                # the minimum or maximum value can occur simultaneously
                # on both the edge of the interval and in some interior
                # point
                if start == _start and start not in solns:
                    left_open = self.left_open
                if end == _end and end not in solns:
                    right_open = self.right_open
            else:
                if start == _end and start not in solns:
                    left_open = self.right_open
                if end == _start and end not in solns:
                    right_open = self.left_open

            return Interval(start, end, left_open, right_open)
        else:
            return imageset(f, Interval(self.start, sing[0],
                                        self.left_open, True)) + \
                Union(*[imageset(f, Interval(sing[i], sing[i + 1], True, True))
                        for i in range(len(sing) - 1)]) + \
                imageset(f, Interval(sing[-1], self.end, True, self.right_open))

    @property
    def _measure(self):
        return self.end - self.start

    def to_mpi(self, prec=53):
        return mpi(mpf(self.start._eval_evalf(prec)),
                   mpf(self.end._eval_evalf(prec)))

    def _eval_evalf(self, prec):
        return Interval(self.left._eval_evalf(prec),
                        self.right._eval_evalf(prec),
                        left_open=self.left_open,
                        right_open=self.right_open)

    def _is_comparable(self, other):
        is_comparable = self.start.is_comparable
        is_comparable &= self.end.is_comparable
        is_comparable &= other.start.is_comparable
        is_comparable &= other.end.is_comparable

        return is_comparable

    @property
    def is_left_unbounded(self):
        """Return ``True`` if the left endpoint is negative infinity."""
        return self.left == -oo

    @property
    def is_right_unbounded(self):
        """Return ``True`` if the right endpoint is positive infinity."""
        return self.right == oo

    def as_relational(self, x):
        """Rewrite an interval in terms of inequalities and logic operators."""
        x = sympify(x)
        if self.right_open:
            right = x < self.end
        else:
            right = true if self.is_right_unbounded else x <= self.end
        if self.left_open:
            left = self.start < x
        else:
            left = true if self.is_left_unbounded else self.start <= x
        return And(left, right)

    def _eval_Eq(self, other):
        if not other.is_Interval:
            if (other.is_Union or other.is_Complement or
                    other.is_Intersection or other.is_ProductSet):
                return

            return false

        return And(Eq(self.left, other.left),
                   Eq(self.right, other.right),
                   self.left_open == other.left_open,
                   self.right_open == other.right_open)


class Union(Set, EvalfMixin):
    """
    Represents a union of sets as a :class:`Set`.

    Examples
    ========

    >>> Union(Interval(1, 2), Interval(3, 4))
    [1, 2] U [3, 4]

    The Union constructor will always try to merge overlapping intervals,
    if possible. For example:

    >>> Union(Interval(1, 2), Interval(2, 3))
    [1, 3]

    See Also
    ========

    Intersection

    References
    ==========

    * https://en.wikipedia.org/wiki/Union_%28set_theory%29

    """

    is_Union = True

    def __new__(cls, *args, **kwargs):
        evaluate = kwargs.get('evaluate', global_evaluate[0])

        # flatten inputs to merge intersections and iterables
        args = list(args)

        def flatten(arg):
            if isinstance(arg, Set):
                if arg.is_Union:
                    return sum(map(flatten, arg.args), [])
                else:
                    return [arg]
            if iterable(arg):  # and not isinstance(arg, Set) (implicit)
                return sum(map(flatten, arg), [])
            raise TypeError('Input must be Sets or iterables of Sets')
        args = flatten(args)

        # Union of no sets is EmptySet
        if len(args) == 0:
            return S.EmptySet

        # Reduce sets using known rules
        if evaluate:
            return Union.reduce(args)

        args = list(ordered(args, Set._infimum_key))

        return Basic.__new__(cls, *args)

    @staticmethod
    def reduce(args):
        """
        Simplify a :class:`Union` using known rules

        We first start with global rules like
        'Merge all FiniteSets'

        Then we iterate through all pairs and ask the constituent sets if they
        can simplify themselves with any other constituent

        """
        # ===== Global Rules =====
        # Merge all finite sets
        finite_sets = [x for x in args if x.is_FiniteSet]
        if len(finite_sets) > 1:
            a = (x for set in finite_sets for x in set)
            finite_set = FiniteSet(*a)
            args = [finite_set] + [x for x in args if not x.is_FiniteSet]

        # ===== Pair-wise Rules =====
        # Here we depend on rules built into the constituent sets
        args = set(args)
        new_args = True
        while(new_args):
            for s in args:
                new_args = False
                for t in args - {s}:
                    new_set = s._union(t)
                    # This returns None if s does not know how to intersect
                    # with t. Returns the newly intersected set otherwise
                    if new_set is not None:
                        if not isinstance(new_set, set):
                            new_set = {new_set}
                        new_args = (args - {s, t}).union(new_set)
                        break
                if new_args:
                    args = new_args
                    break

        if len(args) == 1:
            return args.pop()
        else:
            return Union(args, evaluate=False)

    def _complement(self, universe):
        # DeMorgan's Law
        return Intersection(s.complement(universe) for s in self.args)

    @property
    def _inf(self):
        # We use Min so that sup is meaningful in combination with symbolic
        # interval end points.
        from ..functions import Min
        return Min(*[set.inf for set in self.args])

    @property
    def _sup(self):
        # We use Max so that sup is meaningful in combination with symbolic
        # end points.
        from ..functions import Max
        return Max(*[set.sup for set in self.args])

    def _contains(self, other):
        or_args = [the_set.contains(other) for the_set in self.args]
        return Or(*or_args)

    @property
    def _measure(self):
        # Measure of a union is the sum of the measures of the sets minus
        # the sum of their pairwise intersections plus the sum of their
        # triple-wise intersections minus ... etc...

        # Sets is a collection of intersections and a set of elementary
        # sets which made up those intersections (called "sos" for set of sets)
        # An example element might of this list might be:
        #    ( {A,B,C}, A.intersection(B).intersection(C) )

        # Start with just elementary sets (  ({A}, A), ({B}, B), ... )
        # Then get and subtract (  ({A,B}, (A int B), ... ) while non-zero
        sets = [(FiniteSet(s), s) for s in self.args]
        measure = 0
        parity = 1
        while sets:
            # Add up the measure of these sets and add or subtract it to total
            measure += parity * sum(inter.measure for sos, inter in sets)

            # For each intersection in sets, compute the intersection with every
            # other set not already part of the intersection.
            sets = ((sos + FiniteSet(newset), newset.intersection(intersection))
                    for sos, intersection in sets for newset in self.args
                    if newset not in sos)

            # Clear out sets with no measure
            sets = [(sos, inter) for sos, inter in sets if inter.measure != 0]

            # Clear out duplicates
            sos_list = []
            sets_list = []
            for set in sets:
                if set[0] in sos_list:
                    continue
                else:
                    sos_list.append(set[0])
                    sets_list.append(set)
            sets = sets_list

            # Flip Parity - next time subtract/add if we added/subtracted here
            parity *= -1
        return measure

    @property
    def _boundary(self):
        def boundary_of_set(i):
            """The boundary of set i minus interior of all other sets."""
            b = self.args[i].boundary
            for j, a in enumerate(self.args):
                if j != i:
                    b = b - a.interior
            return b
        return Union(map(boundary_of_set, range(len(self.args))))

    def _eval_imageset(self, f):
        return Union(imageset(f, arg) for arg in self.args)

    def as_relational(self, symbol):
        """Rewrite a Union in terms of equalities and logic operators."""
        return Or(*[set.as_relational(symbol) for set in self.args])

    @property
    def is_iterable(self):
        return all(arg.is_iterable for arg in self.args)

    def _eval_evalf(self, prec):
        return Union(set._eval_evalf(prec) for set in self.args)

    def __iter__(self):
        # roundrobin recipe taken from itertools documentation:
        # https://docs.python.org/3/library/itertools.html#itertools-recipes
        def roundrobin(*iterables):
            """roundrobin('ABC', 'D', 'EF') --> A D E B F C."""
            sentinel = object()
            it = itertools.chain.from_iterable(itertools.zip_longest(fillvalue=sentinel, *iterables))
            return (i for i in it if i is not sentinel)

        if all(set.is_iterable for set in self.args):
            return roundrobin(*(iter(arg) for arg in self.args))
        else:
            raise TypeError('Not all constituent sets are iterable')


class Intersection(Set):
    """
    Represents an intersection of sets as a :class:`Set`.

    Examples
    ========

    >>> Intersection(Interval(1, 3), Interval(2, 4))
    [2, 3]

    We often use the .intersect method

    >>> Interval(1, 3).intersection(Interval(2, 4))
    [2, 3]

    See Also
    ========

    Union

    References
    ==========

    * https://en.wikipedia.org/wiki/Intersection_%28set_theory%29

    """

    is_Intersection = True

    def __new__(cls, *args, **kwargs):
        evaluate = kwargs.get('evaluate', global_evaluate[0])

        # flatten inputs to merge intersections and iterables
        args = list(args)

        def flatten(arg):
            if isinstance(arg, Set):
                if arg.is_Intersection:
                    return sum(map(flatten, arg.args), [])
                else:
                    return [arg]
            if iterable(arg):  # and not isinstance(arg, Set) (implicit)
                return sum(map(flatten, arg), [])
            raise TypeError('Input must be Sets or iterables of Sets')
        args = flatten(args)

        if len(args) == 0:
            return S.UniversalSet

        args = list(ordered(args, Set._infimum_key))

        # Reduce sets using known rules
        if evaluate:
            return Intersection.reduce(args)

        return Basic.__new__(cls, *args)

    @property
    def is_iterable(self):
        return any(arg.is_iterable for arg in self.args)

    def _eval_imageset(self, f):
        return Intersection(imageset(f, arg) for arg in self.args)

    def _contains(self, other):
        return And(*[set.contains(other) for set in self.args])

    def __iter__(self):
        for s in self.args:
            if s.is_iterable:
                other_sets = set(self.args) - {s}
                other = Intersection(other_sets, evaluate=False)
                return (x for x in s if x in other)

        raise ValueError('None of the constituent sets are iterable')

    @staticmethod
    def reduce(args):
        """
        Simplify an intersection using known rules

        We first start with global rules like
        'if any empty sets return empty set' and 'distribute any unions'

        Then we iterate through all pairs and ask the constituent sets if they
        can simplify themselves with any other constituent

        """
        # ===== Global Rules =====

        # If any FiniteSets see which elements of that finite set occur within
        # all other sets in the intersection
        for s in args:
            if s.is_FiniteSet:
                args = [a for a in args if a != s]
                res = s.func(*[x for x in s
                               if all(other.contains(x) == true
                                      for other in args)])
                unk = [x for x in s
                       if any(other.contains(x) not in (true, false)
                              for other in args)]
                if unk:
                    res += Intersection(*([s.func(*unk)] + args), evaluate=False)
                return res

        # If any of the sets are unions, return a Union of Intersections
        for s in args:
            if s.is_Union:
                other_sets = set(args) - {s}
                if len(other_sets) > 0:
                    other = Intersection(other_sets)
                    return Union(Intersection(arg, other) for arg in s.args)
                else:
                    return Union(arg for arg in s.args)

        for s in args:
            if s.is_Complement:
                other_sets = args + [s.args[0]]
                other_sets.remove(s)
                return Complement(Intersection(*other_sets), s.args[1])

        # At this stage we are guaranteed not to have any
        # EmptySets, FiniteSets, or Unions in the intersection

        # ===== Pair-wise Rules =====
        # Here we depend on rules built into the constituent sets
        args = set(args)
        new_args = True
        while(new_args):
            for s in args:
                new_args = False
                for t in args - {s}:
                    new_set = s._intersection(t)
                    # This returns None if s does not know how to intersect
                    # with t. Returns the newly intersected set otherwise
                    if new_set is not None:
                        new_args = (args - {s, t}).union({new_set})
                        break
                if new_args:
                    args = new_args
                    break

        if len(args) == 1:
            return args.pop()
        else:
            return Intersection(args, evaluate=False)

    def as_relational(self, symbol):
        """Rewrite an Intersection in terms of equalities and logic operators."""
        return And(*[set.as_relational(symbol) for set in self.args])


class Complement(Set, EvalfMixin):
    r"""
    Represents relative complement of a set with another set.

    `A - B = \{x \in A| x \notin B\}`

    Examples
    ========

    >>> Complement(FiniteSet(0, 1, 2), FiniteSet(1))
    {0, 2}

    See Also
    =========

    Intersection, Union

    References
    ==========

    * https://mathworld.wolfram.com/ComplementSet.html

    """

    is_Complement = True

    def __new__(cls, a, b, evaluate=True):
        if evaluate:
            return Complement.reduce(a, b)

        return Basic.__new__(cls, a, b)

    @staticmethod
    def reduce(A, B):
        """Simplify a :class:`Complement`."""
        result = B._complement(A)
        if result is not None:
            return result
        else:
            return Complement(A, B, evaluate=False)

    def _contains(self, other):
        A = self.args[0]
        B = self.args[1]
        return And(A.contains(other), Not(B.contains(other)))


class EmptySet(Set, metaclass=Singleton):
    """
    Represents the empty set.

    The empty set is available as a singleton as S.EmptySet.

    Examples
    ========

    >>> S.EmptySet
    EmptySet()

    >>> Interval(1, 2).intersection(S.EmptySet)
    EmptySet()

    See Also
    ========

    UniversalSet

    References
    ==========

    * https://en.wikipedia.org/wiki/Empty_set

    """

    is_EmptySet = True
    is_FiniteSet = True

    @property
    def _measure(self):
        return 0

    def _contains(self, other):
        return false

    def as_relational(self, symbol):
        return false

    def __len__(self):
        return 0

    def _union(self, other):
        return other

    def __iter__(self):
        return iter([])

    def _eval_imageset(self, f):
        return self

    def _eval_powerset(self):
        return FiniteSet(self)

    @property
    def _boundary(self):
        return self

    def _complement(self, other):
        return other

    def _symmetric_difference(self, other):
        return other


class UniversalSet(Set, metaclass=Singleton):
    """
    Represents the set of all things.

    The universal set is available as a singleton as S.UniversalSet

    Examples
    ========

    >>> S.UniversalSet
    UniversalSet()

    >>> Interval(1, 2).intersection(S.UniversalSet)
    [1, 2]

    See Also
    ========

    EmptySet

    References
    ==========

    * https://en.wikipedia.org/wiki/Universal_set

    """

    is_UniversalSet = True

    def _intersection(self, other):
        return other

    def _complement(self, other):
        return S.EmptySet

    def _symmetric_difference(self, other):
        return other

    @property
    def _measure(self):
        return oo

    def _contains(self, other):
        return true

    def as_relational(self, symbol):
        return true

    @property
    def _boundary(self):
        return EmptySet()


class FiniteSet(Set, EvalfMixin):
    """
    Represents a finite set of discrete numbers

    Examples
    ========

    >>> FiniteSet(1, 2, 3, 4)
    {1, 2, 3, 4}
    >>> 3 in FiniteSet(1, 2, 3, 4)
    True

    References
    ==========

    * https://en.wikipedia.org/wiki/Finite_set

    """

    is_FiniteSet = True
    is_iterable = True

    def __new__(cls, *args, **kwargs):
        evaluate = kwargs.get('evaluate', global_evaluate[0])
        if evaluate:
            args = list(map(sympify, args))

            if len(args) == 0:
                return EmptySet()
        else:
            args = list(map(sympify, args))

        args = list(ordered(frozenset(tuple(args)), Set._infimum_key))
        obj = Basic.__new__(cls, *args)
        obj._elements = frozenset(args)
        return obj

    def _eval_Eq(self, other):
        if not other.is_FiniteSet:
            if (other.is_Union or other.is_Complement or
                    other.is_Intersection or other.is_ProductSet):
                return

            return false

        if len(self) != len(other):
            return false

        return And(*(Eq(x, y) for x, y in zip(self.args, other.args)))

    def __iter__(self):
        return iter(self.args)

    def _complement(self, other):
        if other.is_Interval:
            nums = sorted(m for m in self.args if m.is_number and m in other)
            syms = [m for m in self.args if m.is_Symbol]
            # Intervals cannot contain elements other than numbers and symbols.

            intervals = S.EmptySet  # Build up a list of intervals between the elements
            if nums:
                intervals |= Interval(other.left, nums[0],
                                      other.left_open, True)
                for a, b in zip(nums[:-1], nums[1:]):
                    intervals |= Interval(a, b, True, True)  # both open
                intervals |= Interval(nums[-1], other.right,
                                      True, other.right_open)
            else:
                intervals |= other

            if syms:
                return Complement(intervals, FiniteSet(*syms), evaluate=False)
            else:
                return intervals

        elif other.is_FiniteSet:
            common = FiniteSet(*[el for el in other
                                 if self.contains(el) == true])
            self2 = FiniteSet(*[el for el in self
                                if common.contains(el) != true])
            if self2.is_EmptySet:
                self2 = common
            other = FiniteSet(*[el for el in other
                                if common.contains(el) != true])
            return Set._complement(FiniteSet(*[el for el in self2
                                               if other.contains(el) != false]),
                                   other)

        return Set._complement(self, other)

    def _union(self, other):
        """
        This function should only be used internally

        See Set._union for docstring

        """
        # If other set contains one of my elements, remove it from myself
        if any(other.contains(x) is true for x in self):
            return {FiniteSet(*[x for x in self
                                if other.contains(x) is not true]),
                    other}

    def _contains(self, other):
        """
        Tests whether an element, other, is in the set.

        Relies on Python's set class. This tests for object equality
        All inputs are sympified

        Examples
        ========

        >>> 1 in FiniteSet(1, 2)
        True
        >>> 5 in FiniteSet(1, 2)
        False

        """
        r = false
        for e in self._elements:
            t = Eq(e, other, evaluate=True)
            if isinstance(t, Eq):
                t = t.simplify()
            if t == true:
                return t
            elif t != false:
                r = None
        return r

    def _eval_imageset(self, f):
        return FiniteSet(*map(f, self))

    @property
    def _boundary(self):
        return self

    @property
    def _inf(self):
        from ..functions import Min
        return Min(*self)

    @property
    def _sup(self):
        from ..functions import Max
        return Max(*self)

    @property
    def measure(self):
        return 0

    def __len__(self):
        return len(self.args)

    def as_relational(self, symbol):
        """Rewrite a FiniteSet in terms of equalities and logic operators."""
        return Or(*[Eq(symbol, elem) for elem in self])

    def _eval_evalf(self, prec):
        return FiniteSet(*[elem._eval_evalf(prec) for elem in self])

    def _hashable_content(self):
        return self._elements,

    @property
    def _sorted_args(self):
        return tuple(ordered(self.args, Set._infimum_key))

    def _eval_powerset(self):
        return self.func(*[self.func(*s) for s in subsets(self.args)])

    def __ge__(self, other):
        if not isinstance(other, Set):
            raise TypeError(f'Invalid comparison of set with {other!r}')
        return other.is_subset(self)

    def __gt__(self, other):
        if not isinstance(other, Set):
            raise TypeError(f'Invalid comparison of set with {other!r}')
        return self.is_proper_superset(other)

    def __le__(self, other):
        if not isinstance(other, Set):
            raise TypeError(f'Invalid comparison of set with {other!r}')
        return self.is_subset(other)

    def __lt__(self, other):
        if not isinstance(other, Set):
            raise TypeError(f'Invalid comparison of set with {other!r}')
        return self.is_proper_subset(other)


class SymmetricDifference(Set):
    """
    Represents the symmetric difference of two sets.

    The set of elements which are in either of the
    sets and not in their intersection.

    Examples
    ========

    >>> SymmetricDifference(FiniteSet(1, 2, 3), FiniteSet(3, 4, 5))
    {1, 2, 4, 5}

    See Also
    ========

    Complement, Union

    References
    ==========

    * https://en.wikipedia.org/wiki/Symmetric_difference

    """

    is_SymmetricDifference = True

    def __new__(cls, a, b, evaluate=True):
        if evaluate:
            return SymmetricDifference.reduce(a, b)

        return Basic.__new__(cls, a, b)

    @staticmethod
    def reduce(A, B):
        result = B._symmetric_difference(A)
        if result is not None:
            return result
        else:
            return SymmetricDifference(A, B, evaluate=False)


def imageset(*args):
    r"""
    Image of set under transformation ``f``.

    If this function can't compute the image, it returns an
    unevaluated ImageSet object.

    .. math::
        { f(x) | x \in self }

    Examples
    ========

    >>> imageset(x, 2*x, Interval(0, 2))
    [0, 4]

    >>> imageset(lambda x: 2*x, Interval(0, 2))
    [0, 4]

    >>> imageset(Lambda(x, sin(x)), Interval(-2, 1))
    ImageSet(Lambda(x, sin(x)), [-2, 1])

    See Also
    ========

    diofant.sets.fancysets.ImageSet

    """
    from ..core import Dummy, Lambda
    from .fancysets import ImageSet
    if len(args) == 3:
        f = Lambda(*args[:2])
    else:
        # var and expr are being defined this way to
        # support Python lambda and not just diofant Lambda
        f = args[0]
        if not isinstance(f, Lambda):
            var = Dummy()
            expr = args[0](var)
            f = Lambda(var, expr)
    set = args[-1]

    r = set._eval_imageset(f)
    if isinstance(r, ImageSet):
        f, set = r.args

    if f.variables[0] == f.expr:
        return set

    if isinstance(set, ImageSet):
        if len(set.lamda.variables) == 1 and len(f.variables) == 1:
            return imageset(Lambda(set.lamda.variables[0],
                                   f.expr.subs({f.variables[0]:
                                                set.lamda.expr})),
                            set.base_set)

    if r is not None:
        return r

    return ImageSet(f, set)
