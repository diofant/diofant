from ...core import (Add, Dummy, Equality, Expr, Integer, Lambda, Mul, Pow,
                     Rational, Tuple, oo, zoo)
from ...core.compatibility import as_int
from ...core.function import Application, ArgumentIndexError
from ...core.logic import fuzzy_and
from ...core.operations import LatticeOp, ShortCircuitError
from ...core.rules import Transform
from ...core.singleton import SingletonWithManagedProperties as Singleton
from ...core.sympify import sympify
from ...logic import And, Or
from .integers import floor


class IdentityFunction(Lambda, metaclass=Singleton):
    """
    The identity function

    Examples
    ========

    >>> Id(x)
    x

    """

    def __new__(cls):
        from ...sets import FiniteSet
        x = Dummy('dummy_for_IdentityFunction')
        # construct "by hand" to avoid infinite loop
        obj = Expr.__new__(cls, Tuple(x), x)
        obj.nargs = FiniteSet(1)
        return obj


Id = IdentityFunction()

###############################################################################
# ########################### ROOT and SQUARE ROOT FUNCTION ################# #
###############################################################################


def sqrt(arg, **kwargs):
    """The square root function

    sqrt(x) -> Returns the principal square root of x.

    Examples
    ========

    >>> sqrt(x)
    sqrt(x)

    >>> sqrt(x)**2
    x

    Note that sqrt(x**2) does not simplify to x.

    >>> sqrt(x**2)
    sqrt(x**2)

    This is because the two are not equal to each other in general.
    For example, consider x == -1:

    >>> Eq(sqrt(x**2), x).subs({x: -1})
    false

    This is because sqrt computes the principal square root, so the square may
    put the argument in a different branch.  This identity does hold if x is
    positive:

    >>> y = Symbol('y', positive=True)
    >>> sqrt(y**2)
    y

    You can force this simplification by using the powdenest() function with
    the force option set to True:

    >>> sqrt(x**2)
    sqrt(x**2)
    >>> powdenest(sqrt(x**2), force=True)
    x

    To get both branches of the square root you can use the RootOf function:

    >>> [RootOf(x**2 - 3, i) for i in (0, 1)]
    [-sqrt(3), sqrt(3)]

    See Also
    ========

    diofant.polys.rootoftools.RootOf
    diofant.functions.elementary.miscellaneous.root
    diofant.functions.elementary.miscellaneous.real_root

    References
    ==========

    * https://en.wikipedia.org/wiki/Square_root
    * https://en.wikipedia.org/wiki/Principal_value

    """
    return Pow(arg, Rational(1, 2), **kwargs)


def cbrt(arg, **kwargs):
    """This function computes the principial cube root of `arg`, so
    it's just a shortcut for `arg**Rational(1, 3)`.

    Examples
    ========

    >>> cbrt(x)
    x**(1/3)

    >>> cbrt(x)**3
    x

    Note that cbrt(x**3) does not simplify to x.

    >>> cbrt(x**3)
    (x**3)**(1/3)

    This is because the two are not equal to each other in general.
    For example, consider `x == -1`:

    >>> Eq(cbrt(x**3), x).subs({x: -1})
    false

    This is because cbrt computes the principal cube root, this
    identity does hold if `x` is positive:

    >>> y = Symbol('y', positive=True)
    >>> cbrt(y**3)
    y

    See Also
    ========

    diofant.polys.rootoftools.RootOf
    diofant.functions.elementary.miscellaneous.root
    diofant.functions.elementary.miscellaneous.real_root

    References
    ==========

    * https://en.wikipedia.org/wiki/Cube_root
    * https://en.wikipedia.org/wiki/Principal_value

    """
    return Pow(arg, Rational(1, 3), **kwargs)


def root(arg, n, k=0, **kwargs):
    """Returns the k-th n-th root of arg, defaulting to the principle root.

    Examples
    ========

    >>> root(x, 2)
    sqrt(x)

    >>> root(x, 3)
    x**(1/3)

    >>> root(x, n)
    x**(1/n)

    >>> root(x, -Rational(2, 3))
    x**(-3/2)

    To get the k-th n-th root, specify k:

    >>> root(-2, 3, 2)
    -(-1)**(2/3)*2**(1/3)

    To get all n n-th roots you can use the RootOf function.
    The following examples show the roots of unity for n
    equal 2, 3 and 4:

    >>> [RootOf(x**2 - 1, i) for i in range(2)]
    [-1, 1]

    >>> [RootOf(x**3 - 1, i) for i in range(3)]
    [1, -1/2 - sqrt(3)*I/2, -1/2 + sqrt(3)*I/2]

    >>> [RootOf(x**4 - 1, i) for i in range(4)]
    [-1, 1, -I, I]

    Diofant, like other symbolic algebra systems, returns the
    complex root of negative numbers. This is the principal
    root and differs from the text-book result that one might
    be expecting. For example, the cube root of -8 does not
    come back as -2:

    >>> root(-8, 3)
    2*(-1)**(1/3)

    The real_root function can be used to either make the principle
    result real (or simply to return the real root directly):

    >>> real_root(_)
    -2
    >>> real_root(-32, 5)
    -2

    Alternatively, the n//2-th n-th root of a negative number can be
    computed with root:

    >>> root(-32, 5, 5//2)
    -2

    See Also
    ========

    diofant.polys.rootoftools.RootOf
    diofant.core.power.integer_nthroot
    diofant.functions.elementary.miscellaneous.sqrt
    diofant.functions.elementary.miscellaneous.real_root

    References
    ==========

    * https://en.wikipedia.org/wiki/Square_root
    * https://en.wikipedia.org/wiki/Real_root
    * https://en.wikipedia.org/wiki/Root_of_unity
    * https://en.wikipedia.org/wiki/Principal_value
    * https://mathworld.wolfram.com/CubeRoot.html

    """
    n = sympify(n)
    if k:
        return Pow(arg, 1/n, **kwargs)*(-1)**(2*k/n)
    return Pow(arg, 1/n, **kwargs)


def real_root(arg, n=None):
    """Return the real nth-root of arg if possible. If n is omitted then
    all instances of (-n)**(1/odd) will be changed to -n**(1/odd); this
    will only create a real root of a principle root -- the presence of
    other factors may cause the result to not be real.

    Examples
    ========

    >>> real_root(-8, 3)
    -2
    >>> root(-8, 3)
    2*(-1)**(1/3)
    >>> real_root(_)
    -2

    If one creates a non-principle root and applies real_root, the
    result will not be real (so use with caution):

    >>> root(-8, 3, 2)
    -2*(-1)**(2/3)
    >>> real_root(_)
    -2*(-1)**(2/3)


    See Also
    ========

    diofant.polys.rootoftools.RootOf
    diofant.core.power.integer_nthroot
    diofant.functions.elementary.miscellaneous.root
    diofant.functions.elementary.miscellaneous.sqrt

    """
    from .complexes import im
    from .piecewise import Piecewise
    if n is not None:
        try:
            n = as_int(n)
            arg = sympify(arg)
            if arg.is_positive or arg.is_negative:
                rv = root(arg, n)
            else:
                raise ValueError
        except ValueError:
            return root(arg, n)*Piecewise(
                (1, ~Equality(im(arg), 0)),
                (Pow(-1, Integer(1)/n)**(2*floor(n/2)), And(
                    Equality(n % 2, 1),
                    arg < 0)),
                (1, True))
    else:
        rv = sympify(arg)
    n1pow = Transform(lambda x: -(-x.base)**x.exp,
                      lambda x:
                      x.is_Pow and
                      x.base.is_negative and
                      x.exp.is_Rational and
                      x.exp.numerator == 1 and x.exp.denominator % 2)
    return rv.xreplace(n1pow)

###############################################################################
# ########################### MINIMUM and MAXIMUM ########################### #
###############################################################################


class MinMaxBase(LatticeOp):
    """Base class for Min/Max classes."""

    def __new__(cls, *args, **assumptions):
        if not args:
            raise ValueError('The Max/Min functions must have arguments.')

        args = (sympify(arg) for arg in args)

        # first standard filter, for cls.zero and cls.identity
        # also reshape Max(a, Max(b, c)) to Max(a, b, c)
        try:
            _args = frozenset(cls._new_args_filter(args))
        except ShortCircuitError:
            return cls.zero

        # second filter
        # variant I: remove ones which can be removed
        # args = cls._collapse_arguments(set(_args), **assumptions)

        # variant II: find local zeros
        args = cls._find_localzeros(set(_args), **assumptions)

        if not args:
            return cls.identity
        if len(args) == 1:
            return args.pop()
        # base creation
        # XXX should _args be made canonical with sorting?
        _args = frozenset(args)
        obj = Expr.__new__(cls, _args, **assumptions)
        obj._argset = _args
        return obj

    @classmethod
    def _new_args_filter(cls, arg_sequence):
        """
        Generator filtering args.

        first standard filter, for cls.zero and cls.identity.
        Also reshape Max(a, Max(b, c)) to Max(a, b, c),
        and check arguments for comparability

        """
        for arg in arg_sequence:

            # pre-filter, checking comparability of arguments
            if (not isinstance(arg, Expr)) or (arg.is_extended_real is False) or (arg is zoo):
                raise ValueError(f"The argument '{arg}' is not comparable.")

            if arg == cls.zero:
                raise ShortCircuitError(arg)
            if arg == cls.identity:
                continue
            if arg.func == cls:
                yield from arg.args
            else:
                yield arg

    @classmethod
    def _find_localzeros(cls, values, **options):
        """
        Sequentially allocate values to localzeros.

        When a value is identified as being more extreme than another member it
        replaces that member; if this is never true, then the value is simply
        appended to the localzeros.

        """
        localzeros = set()
        for v in values:
            is_newzero = True
            localzeros_ = list(localzeros)
            for z in localzeros_:
                assert v != z
                con = cls._is_connected(v, z)
                if con:
                    is_newzero = False
                    if con is True or con == cls:
                        localzeros.remove(z)
                        localzeros.update([v])
            if is_newzero:
                localzeros.update([v])
        return localzeros

    @classmethod
    def _is_connected(cls, x, y):
        """Check if x and y are connected somehow."""
        def hit(v, t, f):
            if not v.is_Relational:
                return t if v else f
        r = hit(Or(x >= y, y <= x), Max, Min)
        if r is not None:
            return r
        r = hit(Or(x <= y, y >= x), Min, Max)
        if r is not None:
            return r
        return False

    def _eval_derivative(self, s):
        # f(x).diff(s) -> x.diff(s) * f.fdiff(1)(s)
        i = 0
        l = []
        for a in self.args:
            i += 1
            da = a.diff(s)
            if da == 0:
                continue
            df = self.fdiff(i)
            l.append(df * da)
        return Add(*l)

    def evalf(self, dps=15, **options):
        return self.func(*[a.evalf(dps, **options) for a in self.args])

    @property
    def is_extended_real(self):
        return fuzzy_and(arg.is_extended_real for arg in self.args)

    def _eval_rewrite_as_Piecewise(self, *args):
        from .. import Heaviside, Piecewise
        return self.rewrite(Heaviside).rewrite(Piecewise)

    def _eval_rewrite_as_tractable(self, *args, wrt=None, **kwargs):
        r = args[0]
        for a in args[1:]:
            t = self.func((a - r).limit(wrt, oo), 0)
            if (self.func is Max and t != 0) or (self.func is Min and t == 0):
                r = a
        return r

    def _eval_simplify(self, ratio, measure):
        args = [arg.simplify(ratio=ratio, measure=measure) for arg in self.args]
        arg_sets = [set(Mul.make_args(arg)) for arg in args]
        common = arg_sets[0].intersection(*arg_sets[1:])
        new_args = [Mul(*[t for t in Mul.make_args(arg) if t not in common])
                    for arg in args]
        return Mul(*common)*self.func(*new_args)


class Max(MinMaxBase, Application):
    """Return, if possible, the maximum value of the list.

    When number of arguments is equal one, then
    return this argument.

    When number of arguments is equal two, then
    return, if possible, the value from (a, b) that is >= the other.

    In common case, when the length of list greater than 2, the task
    is more complicated. Return only the arguments, which are greater
    than others, if it is possible to determine directional relation.

    If is not possible to determine such a relation, return a partially
    evaluated result.

    Assumptions are used to make the decision too.

    Also, only comparable arguments are permitted.

    It is named ``Max`` and not ``max`` to avoid conflicts
    with the built-in function ``max``.

    Examples
    ========

    >>> p = Symbol('p', positive=True)
    >>> n = Symbol('n', negative=True)

    >>> Max(x, -2)
    Max(-2, x)
    >>> Max(x, -2).subs({x: 3})
    3
    >>> Max(p, -2)
    p
    >>> Max(x, y)
    Max(x, y)
    >>> Max(x, y) == Max(y, x)
    True
    >>> Max(x, Max(y, z))
    Max(x, y, z)
    >>> Max(n, 8, p, 7, -oo)
    Max(8, p)
    >>> Max(1, x, oo)
    oo

    Notes
    =====

    The task can be considered as searching of supremums in the
    directed complete partial orders.

    The source values are sequentially allocated by the isolated subsets
    in which supremums are searched and result as Max arguments.

    If the resulted supremum is single, then it is returned.

    The isolated subsets are the sets of values which are only the comparable
    with each other in the current set. E.g. natural numbers are comparable with
    each other, but not comparable with the `x` symbol. Another example: the
    symbol `x` with negative assumption is comparable with a natural number.

    Also there are "least" elements, which are comparable with all others,
    and have a zero property (maximum or minimum for all elements). E.g. `oo`.
    In case of it the allocation operation is terminated and only this value is
    returned.

    Assumption:
       - if A > B > C then A > C
       - if A == B then B can be removed

    References
    ==========

    * https://en.wikipedia.org/wiki/Directed_complete_partial_order
    * https://en.wikipedia.org/wiki/Lattice_%28order%29

    See Also
    ========

    diofant.functions.elementary.miscellaneous.Min : find minimum values

    """

    zero = oo
    identity = -oo

    def fdiff(self, argindex):
        from .. import Heaviside
        n = len(self.args)
        if 0 < argindex <= n:
            argindex -= 1
            if n == 2:
                return Heaviside(self.args[argindex] - self.args[1 - argindex])
            newargs = tuple(self.args[i] for i in range(n) if i != argindex)
            return Heaviside(self.args[argindex] - Max(*newargs))
        raise ArgumentIndexError(self, argindex)

    def _eval_rewrite_as_Heaviside(self, *args):
        from .. import Heaviside
        return Add(*[j*Mul(*[Heaviside(j - i) for i in args if i != j])
                     for j in args])


class Min(MinMaxBase, Application):
    """Return, if possible, the minimum value of the list.

    It is named ``Min`` and not ``min`` to avoid conflicts
    with the built-in function ``min``.

    Examples
    ========

    >>> p = Symbol('p', positive=True)
    >>> n = Symbol('n', negative=True)

    >>> Min(x, -2)
    Min(-2, x)
    >>> Min(x, -2).subs({x: 3})
    -2
    >>> Min(p, -3)
    -3
    >>> Min(x, y)
    Min(x, y)
    >>> Min(n, 8, p, -7, p, oo)
    Min(-7, n)

    See Also
    ========

    diofant.functions.elementary.miscellaneous.Max : find maximum values

    """

    zero = -oo
    identity = oo

    def fdiff(self, argindex):
        from .. import Heaviside
        n = len(self.args)
        if 0 < argindex <= n:
            argindex -= 1
            if n == 2:
                return Heaviside(self.args[1-argindex] - self.args[argindex])
            newargs = tuple(self.args[i] for i in range(n) if i != argindex)
            return Heaviside(Min(*newargs) - self.args[argindex])
        raise ArgumentIndexError(self, argindex)

    def _eval_rewrite_as_Heaviside(self, *args):
        from .. import Heaviside
        return Add(*[j*Mul(*[Heaviside(i-j) for i in args if i != j])
                     for j in args])
