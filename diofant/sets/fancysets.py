"""Special sets."""

from ..core import Basic, Expr, Integer, Lambda, Rational, S, oo
from ..core.compatibility import as_int
from ..core.singleton import Singleton
from ..core.sympify import converter, sympify
from ..logic import false, true
from ..utilities.iterables import cantor_product
from .sets import EmptySet, FiniteSet, Intersection, Interval, Set


class Naturals(Set, metaclass=Singleton):
    """The set of natural numbers.

    Represents the natural numbers (or counting numbers) which are all
    positive integers starting from 1. This set is also available as
    the Singleton, S.Naturals.

    Examples
    ========

    >>> 5 in S.Naturals
    True
    >>> iterable = iter(S.Naturals)
    >>> next(iterable)
    1
    >>> next(iterable)
    2
    >>> next(iterable)
    3
    >>> S.Naturals.intersection(Interval(0, 10))
    Range(1, 11, 1)

    See Also
    ========

    Naturals0 : non-negative integers
    Integers : also includes negative integers

    """

    is_iterable = True
    inf = Integer(1)
    sup = oo

    def _intersection(self, other):
        if other.is_Interval:
            return Intersection(
                S.Integers, other, Interval(self.inf, oo, False, True))

    def _contains(self, other):
        if not isinstance(other, Expr):
            return false
        if other.is_positive and other.is_integer:
            return true
        if other.is_integer is False or other.is_positive is False:
            return false

    def __iter__(self):
        i = self.inf
        while True:
            yield i
            i = i + 1

    @property
    def boundary(self):
        return self


class Naturals0(Naturals):
    """The set of natural numbers, starting from 0.

    Represents the whole numbers which are all the non-negative
    integers, inclusive of zero.

    See Also
    ========

    Naturals : positive integers
    Integers : also includes the negative integers

    """

    inf = Integer(0)

    def _contains(self, other):
        if not isinstance(other, Expr):
            return false
        if other.is_integer and other.is_nonnegative:
            return true
        if other.is_integer is False or other.is_nonnegative is False:
            return false


class Integers(Set, metaclass=Singleton):
    """The set of all integers.

    Represents all integers: positive, negative and zero. This set
    is also available as the Singleton, S.Integers.

    Examples
    ========

    >>> 5 in S.Naturals
    True
    >>> iterable = iter(S.Integers)
    >>> next(iterable)
    0
    >>> next(iterable)
    1
    >>> next(iterable)
    -1
    >>> next(iterable)
    2

    >>> S.Integers.intersection(Interval(-4, 4))
    Range(-4, 5, 1)

    See Also
    ========

    Naturals0 : non-negative integers
    Integers : positive and negative integers and zero

    """

    is_iterable = True

    def _intersection(self, other):
        from ..functions import ceiling, floor
        if other is Interval(-oo, oo, True, True) or other is S.Reals:
            return self
        if other.is_Interval:
            s = Range(ceiling(other.left), floor(other.right) + 1)
            return s.intersection(other)  # take out endpoints if open interval

    def _contains(self, other):
        if not isinstance(other, Expr):
            return false
        if other.is_integer:
            return true
        if other.is_integer is False:
            return false

    def __iter__(self):
        yield Integer(0)
        i = Integer(1)
        while True:
            yield i
            yield -i
            i = i + 1

    @property
    def inf(self):
        return -oo

    @property
    def sup(self):
        return oo

    @property
    def boundary(self):
        return self

    def _eval_imageset(self, f):
        from ..core import Wild
        expr = f.expr
        if len(f.variables) > 1:
            return
        n = f.variables[0]

        a = Wild('a')
        b = Wild('b')

        match = expr.match(a*n + b)
        if match[a].is_negative:
            expr = -expr

        match = expr.match(a*n + b)
        if match[a] == 1 and match[b].is_integer:
            expr = expr - match[b]

        return ImageSet(Lambda(n, expr), S.Integers)


class Rationals(Set, metaclass=Singleton):
    """The set of all rationals."""

    def _contains(self, other):
        if other.is_rational:
            return true
        if other.is_rational is False:
            return false

    @property
    def inf(self):
        return -oo

    @property
    def sup(self):
        return oo

    @property
    def boundary(self):
        return self

    def __iter__(self):
        seen = []
        for n, d in cantor_product(S.Integers, S.Naturals):  # pragma: no branch
            r = Rational(n, d)
            if r not in seen:
                seen.append(r)
                yield r


class Reals(Interval, metaclass=Singleton):
    """The set of all reals."""

    def __new__(cls, start=-oo, end=oo, left_open=False, right_open=False):
        return super().__new__(cls, left_open=True, right_open=True)

    def __eq__(self, other):
        return other == Interval(-oo, oo, True, True)

    def __hash__(self):
        return hash(Interval(-oo, oo, True, True))


class ExtendedReals(Interval, metaclass=Singleton):
    """The set of all extended reals."""

    def __new__(cls, start=-oo, end=oo, left_open=False, right_open=False):
        return super().__new__(cls)

    def __eq__(self, other):
        return other == Interval(-oo, oo)

    def __hash__(self):
        return hash(Interval(-oo, oo))


class ImageSet(Set):
    """Image of a set under a mathematical function.

    Examples
    ========

    >>> squares = ImageSet(Lambda(x, x**2), S.Naturals)
    >>> 4 in squares
    True
    >>> 5 in squares
    False

    >>> FiniteSet(0, 1, 2, 3, 4, 5, 6, 7, 9, 10).intersection(squares)
    {1, 4, 9}

    >>> square_iterable = iter(squares)
    >>> for i in range(4):
    ...     next(square_iterable)
    1
    4
    9
    16

    If you want to get value for `x` = 2, 1/2 etc. (Please check whether the
    `x` value is in `base_set` or not before passing it as args)

    >>> squares.lamda(2)
    4
    >>> squares.lamda(Rational(1, 2))
    1/4

    """

    def __new__(cls, lamda, base_set):
        return Basic.__new__(cls, lamda, base_set)

    lamda = property(lambda self: self.args[0])
    base_set = property(lambda self: self.args[1])

    def __iter__(self):
        already_seen = set()
        for i in self.base_set:
            val = self.lamda(i)
            if val in already_seen:
                continue
            already_seen.add(val)
            yield val

    def _contains(self, other):
        from ..solvers import solve

        L = self.lamda
        if len(self.lamda.variables) > 1:
            return  # pragma: no cover

        solns = solve(L.expr - other, L.variables[0])

        for soln in solns:
            if soln[L.variables[0]] in self.base_set:
                return true
        return false

    @property
    def is_iterable(self):
        return self.base_set.is_iterable

    def _intersection(self, other):
        from ..core import Dummy, expand_complex
        from ..solvers.diophantine import diophantine
        from .sets import imageset
        if self.base_set is S.Integers:
            if isinstance(other, ImageSet) and other.base_set is S.Integers:
                f, g = self.lamda.expr, other.lamda.expr
                n, m = self.lamda.variables[0], other.lamda.variables[0]

                # Diophantine sorts the solutions according to the alphabetic
                # order of the variable names, since the result should not depend
                # on the variable name, they are replaced by the dummy variables
                # below
                a, b = Dummy('a'), Dummy('b')
                f, g = f.subs({n: a}), g.subs({m: b})
                solns_set = diophantine(f - g)
                if solns_set == set():
                    return EmptySet()
                solns = list(diophantine(f - g))
                if len(solns) == 1:
                    t = list(solns[0][0].free_symbols)[0]
                else:
                    return  # pragma: no cover

                # since 'a' < 'b'
                return imageset(Lambda(t, f.subs({a: solns[0][0]})), S.Integers)

        if other == S.Reals:
            if len(self.lamda.variables) > 1 or self.base_set is not S.Integers:
                return  # pragma: no cover

            f = self.lamda.expr
            n = self.lamda.variables[0]

            n_ = Dummy(n.name, integer=True)
            f_ = f.subs({n: n_})

            re, im = map(expand_complex, f_.as_real_imag())

            sols = list(diophantine(im, n_))
            if all(s[0].has(n_) is False for s in sols):
                s = FiniteSet(*[s[0] for s in sols])
            elif len(sols) == 1 and sols[0][0].has(n_):
                s = imageset(Lambda(n_, sols[0][0]), S.Integers)
            else:
                return  # pragma: no cover

            return imageset(Lambda(n_, re), self.base_set.intersection(s))


class Range(Set):
    """Represents a range of integers.

    Examples
    ========

    >>> list(Range(5))
    [0, 1, 2, 3, 4]
    >>> list(Range(10, 15))
    [10, 11, 12, 13, 14]
    >>> list(Range(10, 20, 2))
    [10, 12, 14, 16, 18]
    >>> list(Range(20, 10, -2))
    [12, 14, 16, 18, 20]

    """

    is_iterable = True

    def __new__(cls, *args):
        from ..functions import ceiling
        if len(args) == 1 and isinstance(args[0], range):
            args = args[0].start, args[0].stop, args[0].step

        # expand range
        slc = slice(*args)
        start, stop, step = slc.start or 0, slc.stop, slc.step or 1
        try:
            start, stop, step = (w if w in [-oo, oo] else Integer(as_int(w))
                                 for w in (start, stop, step))
        except ValueError as exc:
            raise ValueError('Inputs to Range must be Integer Valued\n' +
                             'Use ImageSets of Ranges for other cases') from exc

        if not step.is_finite:
            raise ValueError('Infinite step is not allowed')
        if start == stop:
            return S.EmptySet

        n = ceiling((stop - start)/step)
        if n <= 0:
            return S.EmptySet

        # normalize args: regardless of how they are entered they will show
        # canonically as Range(inf, sup, step) with step > 0
        if n.is_finite:
            start, stop = sorted((start, start + (n - 1)*step))
        else:
            start, stop = sorted((start, stop - step))

        step = abs(step)

        return Basic.__new__(cls, start, stop + step, step)

    start = property(lambda self: self.args[0])
    stop = property(lambda self: self.args[1])
    step = property(lambda self: self.args[2])

    def _intersection(self, other):
        from ..functions import Max, Min, ceiling, floor
        if other.is_Interval:
            osup = other.sup
            oinf = other.inf
            # if other is [0, 10) we can only go up to 9
            if osup.is_integer and other.right_open:
                osup -= 1
            if oinf.is_integer and other.left_open:
                oinf += 1

            # Take the most restrictive of the bounds set by the two sets
            # round inwards
            inf = ceiling(Max(self.inf, oinf))
            sup = floor(Min(self.sup, osup))
            # if we are off the sequence, get back on
            if inf.is_finite and self.inf.is_finite:
                off = (inf - self.inf) % self.step
            else:
                off = Integer(0)
            if off:
                inf += self.step - off

            return Range(inf, sup + 1, self.step)

        if other == S.Naturals:
            return self._intersection(Interval(1, oo, False, True))

        if other == S.Integers:
            return self

    def _contains(self, other):
        if (((self.start - other)/self.step).is_integer or
                ((self.stop - other)/self.step).is_integer):
            return sympify(self.inf <= other <= self.sup, strict=True)
        return false

    def __iter__(self):
        if self.start == -oo:
            i = self.stop - self.step
            step = -self.step
        else:
            i = self.start
            step = self.step

        while self.start <= i < self.stop:
            yield i
            i += step

    def __len__(self):
        return int((self.stop - self.start)//self.step)

    def __bool__(self):
        return True

    def _ith_element(self, i):
        return self.start + i*self.step

    @property
    def _last_element(self):
        if self.stop is oo:
            return oo
        if self.start == -oo:
            return self.stop - self.step
        return self._ith_element(len(self) - 1)

    @property
    def inf(self):
        return self.start

    @property
    def sup(self):
        return self.stop - self.step

    @property
    def boundary(self):
        return self


converter[range] = Range
