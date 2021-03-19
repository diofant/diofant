"""
There are three types of functions implemented in Diofant:

    1) defined functions (in the sense that they can be evaluated) like
       exp or sin; they have a name and a body:
           f = exp
    2) undefined function which have a name but no body. Undefined
       functions can be defined using a Function class as follows:
           f = Function('f')
       (the result will be a Function instance)
    3) anonymous function (or lambda function) which have a body (defined
       with dummy variables) but have no name:
           f = Lambda(x, exp(x)*x)
           f = Lambda((x, y), exp(x)*y)

    Examples
    ========

    >>> f(x)
    f(x)
    >>> print(repr(f(x).func))
    Function('f')
    >>> f(x).args
    (x,)

"""

from __future__ import annotations

import collections
import inspect
import typing

import mpmath
import mpmath.libmp as mlib

from ..utilities import default_sort_key, ordered
from ..utilities.iterables import uniq
from .add import Add
from .assumptions import ManagedProperties
from .basic import Basic
from .cache import cacheit
from .compatibility import as_int, is_sequence, iterable
from .containers import Dict, Tuple
from .decorators import _sympifyit
from .evalf import PrecisionExhausted
from .evaluate import global_evaluate
from .expr import AtomicExpr, Expr
from .logic import fuzzy_and
from .numbers import Float, Integer, Rational, nan
from .operations import LatticeOp
from .rules import Transform
from .singleton import S
from .sympify import sympify


def _coeff_isneg(a):
    """Return True if the leading Number is negative.

    Examples
    ========

    >>> _coeff_isneg(-3*pi)
    True
    >>> _coeff_isneg(Integer(3))
    False
    >>> _coeff_isneg(-oo)
    True
    >>> _coeff_isneg(Symbol('n', negative=True))  # coeff is 1
    False

    """
    if a.is_Mul or a.is_MatMul:
        a = a.args[0]
    return a.is_Number and a.is_negative


class PoleError(Exception):
    """Raised when an expansion pole is encountered."""


class ArgumentIndexError(ValueError):
    """Raised when an invalid operation for positional argument happened."""

    def __str__(self):
        return ('Invalid operation with argument number %s for Function %s' %
                (self.args[1], self.args[0]))


class FunctionClass(ManagedProperties):
    """
    Base class for function classes. FunctionClass is a subclass of type.

    Use Function('<function name>' [ , signature ]) to create
    undefined function classes.

    """

    def __init__(self, *args, **kwargs):
        assert hasattr(self, 'eval')
        evalargspec = inspect.getfullargspec(self.eval)
        if evalargspec.varargs:
            evalargs = None
        else:
            evalargs = len(evalargspec.args) - 1  # subtract 1 for cls
            if evalargspec.defaults:
                # if there are default args then they are optional; the
                # fewest args will occur when all defaults are used and
                # the most when none are used (i.e. all args are given)
                evalargs = tuple(range(evalargs - len(evalargspec.defaults),
                                       evalargs + 1))
        # honor kwarg value or class-defined value before using
        # the number of arguments in the eval function (if present)
        nargs = kwargs.pop('nargs', self.__dict__.get('nargs', evalargs))
        super().__init__(args, kwargs)

        # Canonicalize nargs here; change to set in nargs.
        if is_sequence(nargs):
            if not nargs:
                raise ValueError('Incorrectly specified nargs as %s' % str(nargs))
            nargs = tuple(ordered(set(nargs)))
        elif nargs is not None:
            nargs = as_int(nargs),
        self._nargs = nargs

    @property
    def __signature__(self):
        """
        Allow inspect.signature to give a useful signature for
        Function subclasses.

        """
        # TODO: Look at nargs
        return inspect.signature(self.eval)

    @property
    def nargs(self):
        """Return a set of the allowed number of arguments for the function.

        Examples
        ========

        If the function can take any number of arguments, the set of whole
        numbers is returned:

        >>> Function('f').nargs
        Naturals0()

        If the function was initialized to accept one or more arguments, a
        corresponding set will be returned:

        >>> Function('f', nargs=1).nargs
        {1}
        >>> Function('f', nargs=(2, 1)).nargs
        {1, 2}

        The undefined function, after application, also has the nargs
        attribute; the actual number of arguments is always available by
        checking the ``args`` attribute:

        >>> f(1).nargs
        Naturals0()
        >>> len(f(1).args)
        1

        """
        from ..sets.sets import FiniteSet

        # XXX it would be nice to handle this in __init__ but there are import
        # problems with trying to import FiniteSet there
        return FiniteSet(*self._nargs) if self._nargs else S.Naturals0

    def __repr__(self):
        if issubclass(self, AppliedUndef):
            return f'Function({self.__name__!r})'
        else:
            return self.__name__

    def __str__(self):
        return self.__name__


class Application(Expr, metaclass=FunctionClass):
    """
    Base class for applied functions.

    Instances of Application represent the result of applying an application of
    any type to any object.

    """

    is_Function = True

    @cacheit
    def __new__(cls, *args, **options):
        from ..sets.fancysets import Naturals0
        from ..sets.sets import FiniteSet

        args = list(map(sympify, args))
        evaluate = options.pop('evaluate', global_evaluate[0])
        # WildFunction (and anything else like it) may have nargs defined
        # and we throw that value away here
        options.pop('nargs', None)

        if options:
            raise ValueError(f'Unknown options: {options}')

        if evaluate:
            if nan in args:
                return nan

            evaluated = cls.eval(*args)
            if evaluated is not None:
                return evaluated

        obj = super().__new__(cls, *args, **options)

        # make nargs uniform here
        try:
            # things passing through here:
            #  - functions subclassed from Function (e.g. myfunc(1).nargs)
            #  - functions like cos(1).nargs
            #  - AppliedUndef with given nargs like Function('f', nargs=1)(1).nargs
            # Canonicalize nargs here
            if is_sequence(obj.nargs):
                nargs = tuple(ordered(set(obj.nargs)))
            elif obj.nargs is not None:
                nargs = as_int(obj.nargs),
            else:
                nargs = None
        except AttributeError:
            # things passing through here:
            #  - WildFunction('f').nargs
            #  - AppliedUndef with no nargs like Function('f')(1).nargs
            nargs = obj._nargs  # note the underscore here

        obj.nargs = FiniteSet(*nargs) if nargs else Naturals0()
        return obj

    @classmethod
    def eval(cls, *args):
        """
        Returns a canonical form of cls applied to arguments args.

        The eval() method is called when the class cls is about to be
        instantiated and it should return either some simplified instance
        (possible of some other class), or if the class cls should be
        unmodified, return None.

        """
        return

    def _eval_subs(self, old, new):
        if (old.is_Function and new.is_Function and old == self.func and
                len(self.args) in new.nargs):
            return new(*self.args)


class Function(Application, Expr):
    """Base class for applied mathematical functions.

    It also serves as a constructor for undefined function classes.

    Examples
    ========

    First example shows how to use Function as a constructor for undefined
    function classes:

    >>> g = g(x)
    >>> f
    f
    >>> f(x)
    f(x)
    >>> g
    g(x)
    >>> f(x).diff(x)
    Derivative(f(x), x)
    >>> g.diff(x)
    Derivative(g(x), x)

    In the following example Function is used as a base class for
    ``MyFunc`` that represents a mathematical function *MyFunc*. Suppose
    that it is well known, that *MyFunc(0)* is *1* and *MyFunc* at infinity
    goes to *0*, so we want those two simplifications to occur automatically.
    Suppose also that *MyFunc(x)* is real exactly when *x* is real. Here is
    an implementation that honours those requirements:

    >>> class MyFunc(Function):
    ...
    ...     @classmethod
    ...     def eval(cls, x):
    ...         if x.is_Number:
    ...             if x == 0:
    ...                 return Integer(1)
    ...             elif x is oo:
    ...                 return Integer(0)
    ...
    ...     def _eval_is_real(self):
    ...         return self.args[0].is_real
    ...
    >>> MyFunc(0) + sin(0)
    1
    >>> MyFunc(oo)
    0
    >>> MyFunc(3.54).evalf()  # Not yet implemented for MyFunc.
    MyFunc(3.54)
    >>> MyFunc(I).is_real
    False

    In order for ``MyFunc`` to become useful, several other methods would
    need to be implemented. See source code of some of the already
    implemented functions for more complete examples.

    Also, if the function can take more than one argument, then ``nargs``
    must be defined, e.g. if ``MyFunc`` can take one or two arguments
    then,

    >>> class MyFunc(Function):
    ...     nargs = (1, 2)
    ...
    >>>

    """

    @property
    def _diff_wrt(self):
        """Allow derivatives wrt functions.

        Examples
        ========

        >>> f(x)._diff_wrt
        True

        """
        return True

    @cacheit
    def __new__(cls, *args, **options):
        # Handle calls like Function('f')
        if cls is Function:
            return UndefinedFunction(*args, **options)

        n = len(args)
        if n not in cls.nargs:
            # XXX: exception message must be in exactly this format to
            # make it work with NumPy's functions like vectorize(). See,
            # for example, https://github.com/numpy/numpy/issues/1697.
            # The ideal solution would be just to attach metadata to
            # the exception and change NumPy to take advantage of this.
            temp = ('%(name)s takes %(qual)s %(args)s '
                    'argument%(plural)s (%(given)s given)')
            raise TypeError(temp % {
                'name': cls,
                'qual': 'exactly' if len(cls.nargs) == 1 else 'at least',
                'args': min(cls.nargs),
                'plural': 's'*(min(cls.nargs) != 1),
                'given': n})

        evaluate = options.get('evaluate', global_evaluate[0])
        result = super().__new__(cls, *args, **options)
        if not evaluate or not isinstance(result, cls):
            return result

        pr = max(cls._should_evalf(a) for a in result.args)
        pr2 = min(cls._should_evalf(a) for a in result.args)
        if pr2 > 0:
            return result.evalf(mlib.libmpf.prec_to_dps(pr), strict=False)
        return result

    @classmethod
    def _should_evalf(cls, arg):
        """
        Decide if the function should automatically evalf().

        By default (in this implementation), this happens if (and only if) the
        ARG is a floating point number.
        This function is used by __new__.

        """
        if arg.is_Float:
            return arg._prec
        if not arg.is_Add:
            return -1
        re, im = arg.as_real_imag()
        l = [a._prec for a in [re, im] if a.is_Float]
        l.append(-1)
        return max(l)

    @classmethod
    def class_key(cls):
        """Nice order of classes."""
        from ..sets.fancysets import Naturals0
        funcs = {
            'log': 11,
            'sin': 20,
            'cos': 21,
            'tan': 22,
            'cot': 23,
            'sinh': 30,
            'cosh': 31,
            'tanh': 32,
            'coth': 33,
            'conjugate': 40,
            're': 41,
            'im': 42,
            'arg': 43,
        }
        name = cls.__name__

        try:
            i = funcs[name]
        except KeyError:
            i = 0 if isinstance(cls.nargs, Naturals0) else 10000

        return 4, i, name

    def _eval_evalf(self, prec):
        # Lookup mpmath function based on name
        try:
            if isinstance(self.func, UndefinedFunction):
                # Shouldn't lookup in mpmath but might have ._imp_
                raise AttributeError
            fname = self.func.__name__
            if not hasattr(mpmath, fname):
                from ..utilities.lambdify import MPMATH_TRANSLATIONS
                fname = MPMATH_TRANSLATIONS[fname]
            func = getattr(mpmath, fname)
        except (AttributeError, KeyError):
            try:
                return Float(self._imp_(*[i.evalf(prec) for i in self.args]), prec)
            except (AttributeError, TypeError, ValueError, PrecisionExhausted):
                return

        # Convert all args to mpf or mpc
        # Convert the arguments to *higher* precision than requested for the
        # final result.
        # XXX + 5 is a guess, it is similar to what is used in evalf.py. Should
        #     we be more intelligent about it?
        try:
            args = [arg._to_mpmath(prec + 5) for arg in self.args]
        except ValueError:
            return

        with mpmath.workprec(prec):
            v = func(*args)

        return Expr._from_mpmath(v, prec)

    def _eval_derivative(self, s):
        # f(x).diff(s) -> x.diff(s) * f.fdiff(1)(s)
        i = 0
        l = []
        for a in self.args:
            i += 1
            da = a.diff(s)
            if da == 0:
                continue
            try:
                df = self.fdiff(i)
            except ArgumentIndexError:
                df = Function.fdiff(self, i)
            l.append(df * da)
        return Add(*l)

    def _eval_is_commutative(self):
        return fuzzy_and(a.is_commutative for a in self.args)

    def as_base_exp(self):
        """Returns the method as the 2-tuple (base, exponent)."""
        return self, Integer(1)

    def _eval_aseries(self, n, args0, x, logx):
        """
        Compute an asymptotic expansion around args0, in terms of self.args.
        This function is only used internally by _eval_nseries and should not
        be called directly; derived classes can overwrite this to implement
        asymptotic expansions.

        """
        from ..utilities.misc import filldedent
        raise PoleError(filldedent("""
            Asymptotic expansion of %s around %s is
            not implemented.""" % (type(self), args0)))

    def _eval_nseries(self, x, n, logx):
        """
        This function does compute series for multivariate functions,
        but the expansion is always in terms of *one* variable.
        Examples
        ========

        >>> atan2(x, y).series(x, n=2)
        atan2(0, y) + x/y + O(x**2)
        >>> atan2(x, y).series(y, n=2)
        -y/x + atan2(x, 0) + O(y**2)

        This function also computes asymptotic expansions, if necessary
        and possible:

        >>> loggamma(1/x)._eval_nseries(x, 0, None)
        -1/x - log(x)/x + log(x)/2 + O(1)

        """
        from ..series import Order
        from ..sets.sets import FiniteSet
        from .symbol import Dummy
        args = self.args
        args0 = [t.limit(x, 0) for t in args]
        if any(isinstance(t, Expr) and t.is_finite is False for t in args0):
            from .numbers import oo, zoo

            # XXX could use t.as_leading_term(x) here but it's a little
            # slower
            a = [t.compute_leading_term(x, logx=logx) for t in args]
            a0 = [t.limit(x, 0) for t in a]
            if any(t.has(oo, -oo, zoo, nan) for t in a0):
                return self._eval_aseries(n, args0, x, logx)
            # Careful: the argument goes to oo, but only logarithmically so. We
            # are supposed to do a power series expansion "around the
            # logarithmic term". e.g.
            #      f(1+x+log(x))
            #     -> f(1+logx) + x*f'(1+logx) + O(x**2)
            # where 'logx' is given in the argument
            a = [t._eval_nseries(x, n, logx) for t in args]
            z = [r - r0 for (r, r0) in zip(a, a0)]
            p = [Dummy()]*len(z)
            q = []
            v = None
            for ai, zi, pi in zip(a0, z, p):
                if zi.has(x):
                    if v is not None:
                        raise NotImplementedError
                    q.append(ai + pi)
                    v = pi
                else:
                    q.append(ai)
            e1 = self.func(*q)
            if v is None:
                return e1
            s = e1._eval_nseries(v, n, logx)
            o = s.getO()
            s = s.removeO()
            return s.subs({v: zi}).expand() + Order(o.expr.subs({v: zi}), x)
        if (self.func.nargs is S.Naturals0
                or (self.func.nargs == FiniteSet(1) and args0[0])
                or any(c > 1 for c in self.func.nargs)):
            e = self
            e1 = e.expand()
            if e == e1:
                # for example when e = sin(x+1) or e = sin(cos(x))
                # let's try the general algorithm
                term = e.subs({x: 0})
                if term.is_finite is False:
                    raise PoleError(f'Cannot expand {self} around 0')
                series = term
                fact = Integer(1)
                _x = Dummy('x', real=True, positive=True)
                e = e.subs({x: _x})
                for i in range(n - 1):
                    i += 1
                    fact *= Rational(i)
                    e = e.diff(_x)
                    subs = e.subs({_x: 0})
                    term = subs*(x**i)/fact
                    term = term.expand()
                    series += term
                return series + Order(x**n, x)
            return e1.nseries(x, n=n, logx=logx)
        arg = self.args[0]
        f_series = order = Integer(0)
        i, terms = 0, []
        while order == 0 or i <= n:
            term = self.taylor_term(i, arg, *terms)
            term = term.nseries(x, n=n, logx=logx)
            terms.append(term)
            if term:
                f_series += term
            order = Order(term, x)
            i += 1
        return f_series + order

    def fdiff(self, argindex=1):
        """Returns the first derivative of the function."""
        from .symbol import Dummy

        if not (1 <= argindex <= len(self.args)):
            raise ArgumentIndexError(self, argindex)

        if self.args[argindex - 1].is_Symbol:
            for i in range(len(self.args)):
                if i == argindex - 1:
                    continue
                # See issue sympy/sympy#8510
                if self.args[argindex - 1] in self.args[i].free_symbols:
                    break
            else:
                return Derivative(self, self.args[argindex - 1], evaluate=False)
        # See issue sympy/sympy#4624 and issue sympy/sympy#4719
        # and issue sympy/sympy#5600
        arg_dummy = Dummy(f'xi_{argindex:d}')
        arg_dummy.dummy_index = hash(self.args[argindex - 1])
        new_args = list(self.args)
        new_args[argindex-1] = arg_dummy
        return Subs(Derivative(self.func(*new_args), arg_dummy),
                    (arg_dummy, self.args[argindex - 1]))

    def _eval_as_leading_term(self, x):
        """Stub that should be overridden by new Functions to return
        the first non-zero term in a series if ever an x-dependent
        argument whose leading term vanishes as x -> 0 might be encountered.
        See, for example, cos._eval_as_leading_term.

        """
        from ..series import Order
        args = [a.as_leading_term(x) for a in self.args]
        o = Order(1, x)
        if any(x in a.free_symbols and o.contains(a) for a in args):
            # Whereas x and any finite number are contained in O(1, x),
            # expressions like 1/x are not. If any arg simplified to a
            # vanishing expression as x -> 0 (like x or x**2, but not
            # 3, 1/x, etc...) then the _eval_as_leading_term is needed
            # to supply the first non-zero term of the series,
            #
            # e.g. expression    leading term
            #      ----------    ------------
            #      cos(1/x)      cos(1/x)
            #      cos(cos(x))   cos(1)
            #      cos(x)        1        <- _eval_as_leading_term needed
            #      sin(x)        x        <- _eval_as_leading_term needed
            #
            raise NotImplementedError(
                f'{self.func} has no _eval_as_leading_term routine')
        else:
            return self.func(*args)


class AppliedUndef(Function):
    """
    Base class for expressions resulting from the application of an undefined
    function.

    """

    def __new__(cls, *args, **options):
        args = list(map(sympify, args))
        obj = super().__new__(cls, *args, **options)
        return obj

    def _eval_as_leading_term(self, x):
        return self


class UndefinedFunction(FunctionClass):
    """The (meta)class of undefined functions."""

    def __new__(cls, name, **kwargs):
        ret = type.__new__(cls, name, (AppliedUndef,), kwargs)
        ret.__module__ = None
        return ret

    def __instancecheck__(self, instance):
        return self in type(instance).__mro__

    def __eq__(self, other):
        return (isinstance(other, self.__class__) and
                (self.class_key() == other.class_key()))

    def __hash__(self):
        return super().__hash__()


class WildFunction(Function, AtomicExpr):
    """
    A WildFunction function matches any function (with its arguments).

    Examples
    ========

    >>> F = WildFunction('F')
    >>> F.nargs
    Naturals0()
    >>> x.match(F)
    >>> F.match(F)
    {F_: F_}
    >>> f(x).match(F)
    {F_: f(x)}
    >>> cos(x).match(F)
    {F_: cos(x)}
    >>> f(x, y).match(F)
    {F_: f(x, y)}

    To match functions with a given number of arguments, set ``nargs`` to the
    desired value at instantiation:

    >>> F = WildFunction('F', nargs=2)
    >>> F.nargs
    {2}
    >>> f(x).match(F)
    >>> f(x, y).match(F)
    {F_: f(x, y)}

    To match functions with a range of arguments, set ``nargs`` to a tuple
    containing the desired number of arguments, e.g. if ``nargs = (1, 2)``
    then functions with 1 or 2 arguments will be matched.

    >>> F = WildFunction('F', nargs=(1, 2))
    >>> F.nargs
    {1, 2}
    >>> f(x).match(F)
    {F_: f(x)}
    >>> f(x, y).match(F)
    {F_: f(x, y)}
    >>> f(x, y, 1).match(F)

    """

    include: set[typing.Any] = set()

    def __init__(self, name, **assumptions):
        from ..sets.sets import FiniteSet, Set
        self.name = name
        nargs = assumptions.pop('nargs', S.Naturals0)
        if not isinstance(nargs, Set):
            # Canonicalize nargs here.  See also FunctionClass.
            if is_sequence(nargs):
                nargs = tuple(ordered(set(nargs)))
            else:
                nargs = as_int(nargs),
            nargs = FiniteSet(*nargs)
        self.nargs = nargs

    def _matches(self, expr, repl_dict={}):
        """Helper method for match()

        See Also
        ========

        diofant.core.basic.Basic.matches

        """
        if not isinstance(expr, (AppliedUndef, Function)):
            return
        if len(expr.args) not in self.nargs:
            return

        repl_dict = repl_dict.copy()
        repl_dict[self] = expr
        return repl_dict


class Derivative(Expr):
    """
    Carries out differentiation of the given expression with respect to symbols.

    expr must define ._eval_derivative(symbol) method that returns
    the differentiation result. This function only needs to consider the
    non-trivial case where expr contains symbol and it should call the diff()
    method internally (not _eval_derivative); Derivative should be the only
    one to call _eval_derivative.

    Simplification of high-order derivatives:

    Because there can be a significant amount of simplification that can be
    done when multiple differentiations are performed, results will be
    automatically simplified in a fairly conservative fashion unless the
    keyword ``simplify`` is set to False.

        >>> e = sqrt((x + 1)**2 + x)
        >>> diff(e, (x, 5), simplify=False).count_ops()
        136
        >>> diff(e, (x, 5)).count_ops()
        30

    Ordering of variables:

    If evaluate is set to True and the expression can not be evaluated, the
    list of differentiation symbols will be sorted, that is, the expression is
    assumed to have continuous derivatives up to the order asked. This sorting
    assumes that derivatives wrt Symbols commute, derivatives wrt non-Symbols
    commute, but Symbol and non-Symbol derivatives don't commute with each
    other.

    Derivative wrt non-Symbols:

    This class also allows derivatives wrt non-Symbols that have _diff_wrt
    set to True, such as Function and Derivative. When a derivative wrt a non-
    Symbol is attempted, the non-Symbol is temporarily converted to a Symbol
    while the differentiation is performed.

    Note that this may seem strange, that Derivative allows things like
    f(g(x)).diff(g(x)), or even f(cos(x)).diff(cos(x)).  The motivation for
    allowing this syntax is to make it easier to work with variational calculus
    (i.e., the Euler-Lagrange method).  The best way to understand this is that
    the action of derivative with respect to a non-Symbol is defined by the
    above description:  the object is substituted for a Symbol and the
    derivative is taken with respect to that.  This action is only allowed for
    objects for which this can be done unambiguously, for example Function and
    Derivative objects.  Note that this leads to what may appear to be
    mathematically inconsistent results.  For example::

        >>> (2*cos(x)).diff(cos(x))
        2
        >>> (2*sqrt(1 - sin(x)**2)).diff(cos(x))
        0

    This appears wrong because in fact 2*cos(x) and 2*sqrt(1 - sin(x)**2) are
    identically equal.  However this is the wrong way to think of this.  Think
    of it instead as if we have something like this::

        >>> from diofant.abc import s
        >>> def f(u):
        ...     return 2*u
        ...
        >>> def g(u):
        ...     return 2*sqrt(1 - u**2)
        ...
        >>> f(cos(x))
        2*cos(x)
        >>> g(sin(x))
        2*sqrt(-sin(x)**2 + 1)
        >>> f(c).diff(c)
        2
        >>> f(c).diff(c)
        2
        >>> g(s).diff(c)
        0
        >>> g(sin(x)).diff(cos(x))
        0

    Here, the Symbols c and s act just like the functions cos(x) and sin(x),
    respectively. Think of 2*cos(x) as f(c).subs({c: cos(x)}) (or f(c) *at*
    c = cos(x)) and 2*sqrt(1 - sin(x)**2) as g(s).subs({s: sin(x)}) (or g(s) *at*
    s = sin(x)), where f(u) == 2*u and g(u) == 2*sqrt(1 - u**2).  Here, we
    define the function first and evaluate it at the function, but we can
    actually unambiguously do this in reverse in Diofant, because
    expr.subs({Function: Symbol}) is well-defined:  just structurally replace the
    function everywhere it appears in the expression.

    This is the same notational convenience used in the Euler-Lagrange method
    when one says F(t, f(t), f'(t)).diff(f(t)).  What is actually meant is
    that the expression in question is represented by some F(t, u, v) at u =
    f(t) and v = f'(t), and F(t, f(t), f'(t)).diff(f(t)) simply means F(t, u,
    v).diff(u) at u = f(t).

    We do not allow derivatives to be taken with respect to expressions where this
    is not so well defined.  For example, we do not allow expr.diff(x*y)
    because there are multiple ways of structurally defining where x*y appears
    in an expression, some of which may surprise the reader (for example, a
    very strict definition would have that (x*y*z).diff(x*y) == 0).

        >>> (x*y*z).diff(x*y)
        Traceback (most recent call last):
        ...
        ValueError: Can't differentiate wrt the variable: x*y, 1

    Note that this definition also fits in nicely with the definition of the
    chain rule.  Note how the chain rule in Diofant is defined using unevaluated
    Subs objects::

        >>> f, g = symbols('f g', cls=Function)
        >>> f(2*g(x)).diff(x)
        2*Derivative(g(x), x)*Subs(Derivative(f(_xi_1), _xi_1), (_xi_1, 2*g(x)))
        >>> f(g(x)).diff(x)
        Derivative(g(x), x)*Subs(Derivative(f(_xi_1), _xi_1), (_xi_1, g(x)))

    Finally, note that, to be consistent with variational calculus, and to
    ensure that the definition of substituting a Function for a Symbol in an
    expression is well-defined, derivatives of functions are assumed to not be
    related to the function.  In other words, we have::

        >>> diff(f(x), x).diff(f(x))
        0

    The same is true for derivatives of different orders::

        >>> diff(f(x), (x, 2)).diff(diff(f(x), (x, 1)))
        0
        >>> diff(f(x), (x, 1)).diff(diff(f(x), (x, 2)))
        0

    Note, any class can allow derivatives to be taken with respect to itself.

    Examples
    ========

    Some basic examples:

        >>> Derivative(x**2, x, evaluate=True)
        2*x
        >>> Derivative(Derivative(f(x, y), x), y)
        Derivative(f(x, y), x, y)
        >>> Derivative(f(x), (x, 3))
        Derivative(f(x), x, x, x)
        >>> Derivative(f(x, y), y, x, evaluate=True)
        Derivative(f(x, y), x, y)

    Now some derivatives wrt functions:

        >>> Derivative(f(x)**2, f(x), evaluate=True)
        2*f(x)
        >>> Derivative(f(g(x)), x, evaluate=True)
        Derivative(g(x), x)*Subs(Derivative(f(_xi_1), _xi_1), (_xi_1, g(x)))

    """

    is_Derivative = True

    @property
    def _diff_wrt(self):
        """Allow derivatives wrt Derivatives if it contains a function.

        Examples
        ========

            >>> Derivative(f(x), x)._diff_wrt
            True
            >>> Derivative(x**2, x)._diff_wrt
            False

        """
        if self.expr.is_Function:
            return True
        else:
            return False

    def __new__(cls, expr, *args, **assumptions):
        from .symbol import Dummy

        expr = sympify(expr)

        # There are no args, we differentiate wrt all of the free symbols
        # in expr.
        if not args:
            variables = expr.free_symbols
            args = tuple(variables)
            if len(variables) != 1:
                from ..utilities.misc import filldedent
                raise ValueError(filldedent("""
                    The variable(s) of differentiation
                    must be supplied to differentiate %s""" % expr))

        # Standardize the args by sympifying them and making appending a
        # count of 1 if there is only variable: diff(e, x) -> diff(e, (x, 1)).
        args = list(sympify(args))
        for i, a in enumerate(args):
            if not isinstance(a, Tuple):
                args[i] = (a, Integer(1))

        variable_count = []
        all_zero = True
        for v, count in args:
            if not v._diff_wrt:
                from ..utilities.misc import filldedent
                ordinal = 'st' if count == 1 else 'nd' if count == 2 else 'rd' if count == 3 else 'th'
                raise ValueError(filldedent("""
                Can\'t calculate %s%s derivative wrt %s.""" % (count, ordinal, v)))
            if count:
                if all_zero:
                    all_zero = False
                variable_count.append(Tuple(v, count))

        # We make a special case for 0th derivative, because there is no
        # good way to unambiguously print this.
        if all_zero:
            return expr

        # Pop evaluate because it is not really an assumption and we will need
        # to track it carefully below.
        evaluate = assumptions.pop('evaluate', False)

        # Look for a quick exit if there are symbols that don't appear in
        # expression at all. Note, this cannnot check non-symbols like
        # functions and Derivatives as those can be created by intermediate
        # derivatives.
        if evaluate:
            symbol_set = {sc[0] for sc in variable_count if sc[0].is_Symbol}
            if symbol_set.difference(expr.free_symbols):
                return Integer(0)

        # We make a generator so as to only generate a variable when necessary.
        # If a high order of derivative is requested and the expr becomes 0
        # after a few differentiations, then we won't need the other variables.
        variablegen = (v for v, count in variable_count for i in range(count))

        # If we can't compute the derivative of expr (but we wanted to) and
        # expr is itself not a Derivative, finish building an unevaluated
        # derivative class by calling Expr.__new__.
        if (not (hasattr(expr, '_eval_derivative') and evaluate) and
                (not isinstance(expr, Derivative))):
            variables = list(variablegen)
            # If we wanted to evaluate, we sort the variables into standard
            # order for later comparisons. This is too aggressive if evaluate
            # is False, so we don't do it in that case.
            if evaluate:
                # TODO: check if assumption of discontinuous derivatives exist
                variables = cls._sort_variables(variables)
            # Here we *don't* need to reinject evaluate into assumptions
            # because we are done with it and it is not an assumption that
            # Expr knows about.
            obj = Expr.__new__(cls, expr, *variables, **assumptions)
            return obj

        # Compute the derivative now by repeatedly calling the
        # _eval_derivative method of expr for each variable. When this method
        # returns None, the derivative couldn't be computed wrt that variable
        # and we save the variable for later.
        unhandled_variables = []

        # Once we encouter a non_symbol that is unhandled, we stop taking
        # derivatives entirely. This is because derivatives wrt functions
        # don't commute with derivatives wrt symbols and we can't safely
        # continue.
        unhandled_non_symbol = False
        nderivs = 0  # how many derivatives were performed
        for v in variablegen:
            is_symbol = v.is_Symbol

            if unhandled_non_symbol:
                obj = None
            else:
                if not is_symbol:
                    new_v = Dummy(f'xi_{i:d}')
                    new_v.dummy_index = hash(v)
                    expr = expr.xreplace({v: new_v})
                    old_v = v
                    v = new_v
                obj = expr._eval_derivative(v)
                nderivs += 1
                if not is_symbol:
                    if obj is not None:
                        if obj.is_Derivative and not old_v.is_Symbol:
                            # Derivative evaluated at a generic point, i.e.
                            # that is not a symbol.
                            obj = Subs(obj, (v, old_v))
                        else:
                            obj = obj.xreplace({v: old_v})
                    v = old_v

            if obj is None:
                unhandled_variables.append(v)
                if not is_symbol:
                    unhandled_non_symbol = True
            elif obj == 0:
                return Integer(0)
            else:
                expr = obj

        if unhandled_variables:
            unhandled_variables = cls._sort_variables(unhandled_variables)
            expr = Expr.__new__(cls, expr, *unhandled_variables, **assumptions)
        else:
            # We got a Derivative at the end of it all, and we rebuild it by
            # sorting its variables.
            if isinstance(expr, Derivative):
                expr = cls(
                    expr.expr, *cls._sort_variables(expr.variables)
                )

        if nderivs > 1 and assumptions.get('simplify', True):
            from ..simplify.simplify import signsimp
            from .exprtools import factor_terms
            expr = factor_terms(signsimp(expr))
        return expr

    @classmethod
    def _sort_variables(cls, vars):
        """Sort variables, but disallow sorting of non-symbols.

        When taking derivatives, the following rules usually hold:

        * Derivative wrt different symbols commute.
        * Derivative wrt different non-symbols commute.
        * Derivatives wrt symbols and non-symbols don't commute.

        Examples
        ========

        >>> vsort = Derivative._sort_variables

        >>> vsort((x, y, z))
        [x, y, z]

        >>> vsort((h(x), g(x), f(x)))
        [f(x), g(x), h(x)]

        >>> vsort((z, y, x, h(x), g(x), f(x)))
        [x, y, z, f(x), g(x), h(x)]

        >>> vsort((x, f(x), y, f(y)))
        [x, f(x), y, f(y)]

        >>> vsort((y, x, g(x), f(x), z, h(x), y, x))
        [x, y, f(x), g(x), z, h(x), x, y]

        >>> vsort((z, y, f(x), x, f(x), g(x)))
        [y, z, f(x), x, f(x), g(x)]

        >>> vsort((z, y, f(x), x, f(x), g(x), z, z, y, x))
        [y, z, f(x), x, f(x), g(x), x, y, z, z]

        """
        sorted_vars = []
        symbol_part = []
        non_symbol_part = []
        for v in vars:
            if not v.is_Symbol:
                if len(symbol_part) > 0:
                    sorted_vars.extend(sorted(symbol_part,
                                              key=default_sort_key))
                    symbol_part = []
                non_symbol_part.append(v)
            else:
                if len(non_symbol_part) > 0:
                    sorted_vars.extend(sorted(non_symbol_part,
                                              key=default_sort_key))
                    non_symbol_part = []
                symbol_part.append(v)
        if len(non_symbol_part) > 0:
            sorted_vars.extend(sorted(non_symbol_part,
                                      key=default_sort_key))
        if len(symbol_part) > 0:
            sorted_vars.extend(sorted(symbol_part,
                                      key=default_sort_key))
        return sorted_vars

    def _eval_is_commutative(self):
        return self.expr.is_commutative

    def _eval_derivative(self, v):
        # If the variable s we are diff wrt is not in self.variables, we
        # assume that we might be able to take the derivative.
        if v not in self.variables:
            obj = self.expr.diff(v)
            if obj == 0:
                return Integer(0)
            if isinstance(obj, Derivative):
                return obj.func(obj.expr, *(self.variables + obj.variables))
            # The derivative wrt s could have simplified things such that the
            # derivative wrt things in self.variables can now be done. Thus,
            # we set evaluate=True to see if there are any other derivatives
            # that can be done. The most common case is when obj is a simple
            # number so that the derivative wrt anything else will vanish.
            return self.func(obj, *self.variables, evaluate=True)
        # In this case s was in self.variables so the derivatve wrt s has
        # already been attempted and was not computed, either because it
        # couldn't be or evaluate=False originally.
        return self.func(self.expr, *(self.variables + (v, )), evaluate=False)

    def doit(self, **hints):
        """Evaluate objects that are not evaluated by default.

        See Also
        ========

        diofant.core.basic.Basic.doit

        """
        expr = self.expr
        if hints.get('deep', True):
            expr = expr.doit(**hints)
        hints['evaluate'] = True
        return self.func(expr, *self.variables, **hints)

    @_sympifyit('z0', NotImplementedError)
    def doit_numerically(self, z0):
        """
        Evaluate the derivative at z numerically.

        When we can represent derivatives at a point, this should be folded
        into the normal evalf. For now, we need a special method.

        """
        import mpmath

        from .expr import Expr
        if len(self.free_symbols) != 1 or len(self.variables) != 1:
            raise NotImplementedError('partials and higher order derivatives')
        z = list(self.free_symbols)[0]

        def eval(x):
            f0 = self.expr.subs({z: Expr._from_mpmath(x, prec=mpmath.mp.prec)})
            f0 = f0.evalf(mlib.libmpf.prec_to_dps(mpmath.mp.prec), strict=False)
            return f0._to_mpmath(mpmath.mp.prec)
        return Expr._from_mpmath(mpmath.diff(eval,
                                             z0._to_mpmath(mpmath.mp.prec)),
                                 mpmath.mp.prec)

    @property
    def expr(self):
        """Return expression."""
        return self.args[0]

    @property
    def variables(self):
        """Return tuple of symbols, wrt derivative is taken."""
        return self.args[1:]

    @property
    def free_symbols(self):
        """Return from the atoms of self those which are free symbols.

        See Also
        ========

        diofant.core.basic.Basic.free_symbols

        """
        return self.expr.free_symbols

    def _eval_subs(self, old, new):
        if old in self.variables and not new._diff_wrt:
            # issue sympy/sympy#4719
            return Subs(self, (old, new))
        # If both are Derivatives with the same expr, check if old is
        # equivalent to self or if old is a subderivative of self.
        if old.is_Derivative and old.expr == self.expr:
            # Check if canonnical order of variables is equal.
            old_vars = collections.Counter(old.variables)
            self_vars = collections.Counter(self.variables)
            if old_vars == self_vars:
                return new

            # collections.Counter doesn't have __le__
            def _subset(a, b):
                return all(a[i] <= b[i] for i in a)

            if _subset(old_vars, self_vars):
                return Derivative(new, *(self_vars - old_vars).elements())

        return Derivative(*(x._subs(old, new) for x in self.args))

    def _eval_lseries(self, x, logx):
        for term in self.expr.lseries(x, logx=logx):
            yield self.func(term, *self.variables)

    def _eval_nseries(self, x, n, logx):
        arg = self.expr.nseries(x, n=n, logx=logx)
        o = arg.getO()
        rv = [self.func(a, *self.variables) for a in Add.make_args(arg.removeO())]
        if o:
            rv.append(o/x)
        return Add(*rv)

    def _eval_as_leading_term(self, x):
        return self.func(self.expr.as_leading_term(x), *self.variables)


class Lambda(Expr):
    """
    Lambda(x, expr) represents a lambda function similar to Python's
    'lambda x: expr'. A function of several variables is written as
    Lambda((x, y, ...), expr).

    A simple example:

    >>> f = Lambda(x, x**2)
    >>> f(4)
    16

    For multivariate functions, use:

    >>> f2 = Lambda((x, y, z, t), x + y**z + t**z)
    >>> f2(1, 2, 3, 4)
    73

    A handy shortcut for lots of arguments:

    >>> p = x, y, z
    >>> f = Lambda(p, x + y*z)
    >>> f(*p)
    x + y*z

    """

    is_Function = True

    def __new__(cls, variables, expr):
        from ..sets.sets import FiniteSet
        v = list(variables) if iterable(variables) else [variables]
        for i in v:
            if not getattr(i, 'is_Symbol', False):
                raise TypeError(f'variable is not a symbol: {i}')
        if len(v) == 1 and v[0] == expr:
            return S.IdentityFunction

        obj = Expr.__new__(cls, Tuple(*v), sympify(expr))
        obj.nargs = FiniteSet(len(v))
        return obj

    @property
    def variables(self):
        """The variables used in the internal representation of the function."""
        return self.args[0]

    @property
    def expr(self):
        """The return value of the function."""
        return self.args[1]

    @property
    def free_symbols(self):
        """Return from the atoms of self those which are free symbols.

        See Also
        ========

        diofant.core.basic.Basic.free_symbols

        """
        return self.expr.free_symbols - set(self.variables)

    def __call__(self, *args):
        n = len(args)
        if n not in self.nargs:  # Lambda only ever has 1 value in nargs
            # XXX: exception message must be in exactly this format to
            # make it work with NumPy's functions like vectorize(). See,
            # for example, https://github.com/numpy/numpy/issues/1697.
            # The ideal solution would be just to attach metadata to
            # the exception and change NumPy to take advantage of this.
            # XXX does this apply to Lambda? If not, remove this comment.
            temp = ('%(name)s takes exactly %(args)s '
                    'argument%(plural)s (%(given)s given)')
            raise TypeError(temp % {
                'name': self,
                'args': list(self.nargs)[0],
                'plural': 's'*(list(self.nargs)[0] != 1),
                'given': n})
        return self.expr.xreplace(dict(zip(self.variables, args)))

    def __eq__(self, other):
        if not isinstance(other, Lambda):
            return False
        if self.nargs != other.nargs:
            return False

        selfexpr = self.args[1]
        otherexpr = other.args[1]
        otherexpr = otherexpr.xreplace(dict(zip(other.args[0], self.args[0])))
        return selfexpr == otherexpr

    def __hash__(self):
        return super().__hash__()

    def _hashable_content(self):
        return self.expr.xreplace(self.canonical_variables),


class Subs(Expr):
    """
    Represents unevaluated substitutions of an expression.

    ``Subs`` receives at least 2 arguments: an expression, a pair of old
    and new expression to substitute or several such pairs.

    ``Subs`` objects are generally useful to represent unevaluated derivatives
    calculated at a point.

    The variables may be expressions, but they are subjected to the limitations
    of subs(), so it is usually a good practice to use only symbols for
    variables, since in that case there can be no ambiguity.

    There's no automatic expansion - use the method .doit() to effect all
    possible substitutions of the object and also of objects inside the
    expression.

    When evaluating derivatives at a point that is not a symbol, a Subs object
    is returned. One is also able to calculate derivatives of Subs objects - in
    this case the expression is always expanded (for the unevaluated form, use
    Derivative()).

    Examples
    ========

    >>> e = Subs(f(x).diff(x), (x, y))
    >>> e.subs({y: 0})
    Subs(Derivative(f(x), x), (x, 0))
    >>> e.subs({f: sin}).doit()
    cos(y)

    >>> Subs(f(x)*sin(y) + z, (x, 0), (y, 1))
    Subs(z + f(x)*sin(y), (x, 0), (y, 1))
    >>> _.doit()
    z + f(0)*sin(1)

    """

    def __new__(cls, expr, *args, **assumptions):
        from .symbol import Symbol
        args = sympify(args)
        if len(args) and all(is_sequence(_) and len(_) == 2 for _ in args):
            variables, point = zip(*args)
        else:
            raise ValueError('Subs support two or more arguments')

        if tuple(uniq(variables)) != variables:
            repeated = [ v for v in set(variables) if variables.count(v) > 1 ]
            raise ValueError('cannot substitute expressions %s more than '
                             'once.' % repeated)

        expr = sympify(expr)

        # use symbols with names equal to the point value (with preppended _)
        # to give a variable-independent expression
        pre = '_'
        pts = sorted(set(point), key=default_sort_key)
        from ..printing import StrPrinter

        class CustomStrPrinter(StrPrinter):
            def _print_Dummy(self, expr):
                return str(expr) + str(expr.dummy_index)

        def mystr(expr, **settings):
            p = CustomStrPrinter(settings)
            return p.doprint(expr)

        while 1:
            s_pts = {p: Symbol(pre + mystr(p)) for p in pts}
            reps = [(v, s_pts[p])
                    for v, p in zip(variables, point)]
            # if any underscore-preppended symbol is already a free symbol
            # and is a variable with a different point value, then there
            # is a clash, e.g. _0 clashes in Subs(_0 + _1, (_0, 1), (_1, 0))
            # because the new symbol that would be created is _1 but _1
            # is already mapped to 0 so __0 and __1 are used for the new
            # symbols
            if any(r in expr.free_symbols and
                   r in variables and
                   Symbol(pre + mystr(point[variables.index(r)])) != r
                   for _, r in reps):
                pre += '_'
                continue
            reps  # XXX "peephole" optimization, http://bugs.python.org/issue2506
            break

        obj = Expr.__new__(cls, expr, *sympify(tuple(zip(variables, point))))
        obj._expr = expr.subs(reps)
        return obj

    def _eval_is_commutative(self):
        return (self.expr.is_commutative and
                all(p.is_commutative for p in self.point))

    def doit(self, **hints):
        """Evaluate objects that are not evaluated by default.

        See Also
        ========

        diofant.core.basic.Basic.doit

        """
        return self.expr.doit(**hints).subs(list(zip(self.variables, self.point)))

    def evalf(self, dps=15, **options):
        """Evaluate the given formula to an accuracy of dps decimal digits.

        See Also
        ========

        diofant.core.evalf.EvalfMixin.evalf

        """
        return self.doit().evalf(dps, **options)

    #:
    n = evalf

    @property
    def variables(self):
        """The variables to be evaluated."""
        return Tuple(*tuple(zip(*self.args[1:])))[0]

    @property
    def expr(self):
        """The expression on which the substitution operates."""
        return self.args[0]

    @property
    def point(self):
        """The values for which the variables are to be substituted."""
        return Tuple(*tuple(zip(*self.args[1:])))[1]

    @property
    def free_symbols(self):
        """Return from the atoms of self those which are free symbols.

        See Also
        ========

        diofant.core.basic.Basic.free_symbols

        """
        return (self.expr.free_symbols - set(self.variables) |
                set(self.point.free_symbols))

    def __eq__(self, other):
        if not isinstance(other, Subs):
            return False
        return self._expr == other._expr

    def __hash__(self):
        return super().__hash__()

    def _hashable_content(self):
        return self._expr.xreplace(self.canonical_variables),

    def _eval_subs(self, old, new):
        if old in self.variables:
            return self

        if isinstance(old, Subs) and self.point == old.point:
            if self.expr.subs(zip(self.variables, old.variables)) == old.expr:
                return new

    def _eval_derivative(self, s):
        return Add((self.func(self.expr.diff(s), *self.args[1:]).doit()
                    if s not in self.variables else Integer(0)),
                   *[p.diff(s)*self.func(self.expr.diff(v), *self.args[1:]).doit()
                     for v, p in zip(self.variables, self.point)])


def diff(f, *args, **kwargs):
    """
    Differentiate f with respect to symbols.

    This is just a wrapper to unify .diff() and the Derivative class; its
    interface is similar to that of integrate().  You can use the same
    shortcuts for multiple variables as with Derivative.  For example,
    diff(f(x), x, x, x) and diff(f(x), (x, 3)) both return the third derivative
    of f(x).

    You can pass evaluate=False to get an unevaluated Derivative class.  Note
    that if there are 0 symbols (such as diff(f(x), (x, 0)), then the result will
    be the function (the zeroth derivative), even if evaluate=False.

    Examples
    ========

    >>> diff(sin(x), x)
    cos(x)
    >>> diff(f(x), x, x, x)
    Derivative(f(x), x, x, x)
    >>> diff(f(x), (x, 3))
    Derivative(f(x), x, x, x)
    >>> diff(sin(x)*cos(y), (x, 2), (y, 2))
    sin(x)*cos(y)

    >>> type(diff(sin(x), x))
    cos
    >>> type(diff(sin(x), x, evaluate=False))
    <class 'diofant.core.function.Derivative'>
    >>> type(diff(sin(x), (x, 0)))
    sin
    >>> type(diff(sin(x), (x, 0), evaluate=False))
    sin

    >>> diff(sin(x))
    cos(x)
    >>> diff(sin(x*y))
    Traceback (most recent call last):
    ...
    ValueError: specify differentiation variables to differentiate sin(x*y)

    Note that ``diff(sin(x))`` syntax is meant only for convenience
    in interactive sessions and should be avoided in library code.

    References
    ==========

    * https://reference.wolfram.com/legacy/v5_2/Built-inFunctions/AlgebraicComputation/Calculus/D.html

    See Also
    ========

    Derivative
    diofant.geometry.util.idiff: computes the derivative implicitly

    """
    kwargs.setdefault('evaluate', True)
    return Derivative(f, *args, **kwargs)


def expand(e, deep=True, modulus=None, power_base=True, power_exp=True,
           mul=True, log=True, multinomial=True, basic=True, **hints):
    r"""Expand an expression using methods given as hints.

    Hints evaluated unless explicitly set to False are:  ``basic``, ``log``,
    ``multinomial``, ``mul``, ``power_base``, and ``power_exp``.  The following
    hints are supported but not applied unless set to True:  ``complex``,
    ``func``, and ``trig``.  In addition, the following meta-hints are
    supported by some or all of the other hints:  ``frac``, ``numer``,
    ``denom``, ``modulus``, and ``force``.  ``deep`` is supported by all
    hints.  Additionally, subclasses of Expr may define their own hints or
    meta-hints.

    Parameters
    ==========

    basic : boolean, optional
        This hint is used for any special
        rewriting of an object that should be done automatically (along with
        the other hints like ``mul``) when expand is called. This is a catch-all
        hint to handle any sort of expansion that may not be described by
        the existing hint names.

    deep : boolean, optional
        If ``deep`` is set to ``True`` (the default), things like arguments of
        functions are recursively expanded.  Use ``deep=False`` to only expand on
        the top level.

    mul : boolean, optional
        Distributes multiplication over addition (``):

        >>> (y*(x + z)).expand(mul=True)
        x*y + y*z

    multinomial : boolean, optional
        Expand (x + y + ...)**n where n is a positive integer.

        >>> ((x + y + z)**2).expand(multinomial=True)
        x**2 + 2*x*y + 2*x*z + y**2 + 2*y*z + z**2

    power_exp : boolean, optional
        Expand addition in exponents into multiplied bases.

        >>> exp(x + y).expand(power_exp=True)
        E**x*E**y
        >>> (2**(x + y)).expand(power_exp=True)
        2**x*2**y

    power_base : boolean, optional
        Split powers of multiplied bases.

        This only happens by default if assumptions allow, or if the
        ``force`` meta-hint is used:

        >>> ((x*y)**z).expand(power_base=True)
        (x*y)**z
        >>> ((x*y)**z).expand(power_base=True, force=True)
        x**z*y**z
        >>> ((2*y)**z).expand(power_base=True)
        2**z*y**z

        Note that in some cases where this expansion always holds, Diofant performs
        it automatically:

        >>> (x*y)**2
        x**2*y**2

    log : boolean, optional
        Pull out power of an argument as a coefficient and split logs products
        into sums of logs.

        Note that these only work if the arguments of the log function have the
        proper assumptions--the arguments must be positive and the exponents must
        be real--or else the ``force`` hint must be True:

        >>> log(x**2*y).expand(log=True)
        log(x**2*y)
        >>> log(x**2*y).expand(log=True, force=True)
        2*log(x) + log(y)
        >>> x, y = symbols('x y', positive=True)
        >>> log(x**2*y).expand(log=True)
        2*log(x) + log(y)

    complex : boolean, optional
        Split an expression into real and imaginary parts.

        >>> x, y = symbols('x y')
        >>> (x + y).expand(complex=True)
        re(x) + re(y) + I*im(x) + I*im(y)
        >>> cos(x).expand(complex=True)
        -I*sin(re(x))*sinh(im(x)) + cos(re(x))*cosh(im(x))

        Note that this is just a wrapper around ``as_real_imag()``.  Most objects
        that wish to redefine ``_eval_expand_complex()`` should consider
        redefining ``as_real_imag()`` instead.

    func : boolean : optional
        Expand other functions.

        >>> gamma(x + 1).expand(func=True)
        x*gamma(x)

    trig : boolean, optional
        Do trigonometric expansions.

        >>> cos(x + y).expand(trig=True)
        -sin(x)*sin(y) + cos(x)*cos(y)
        >>> sin(2*x).expand(trig=True)
        2*sin(x)*cos(x)

        Note that the forms of ``sin(n*x)`` and ``cos(n*x)`` in terms of ``sin(x)``
        and ``cos(x)`` are not unique, due to the identity `\sin^2(x) + \cos^2(x)
        = 1`.  The current implementation uses the form obtained from Chebyshev
        polynomials, but this may change.

    force : boolean, optional
        If the ``force`` hint is used, assumptions about variables will be ignored
        in making the expansion.


    Notes
    =====

    - You can shut off unwanted methods::

        >>> (exp(x + y)*(x + y)).expand()
        E**x*E**y*x + E**x*E**y*y
        >>> (exp(x + y)*(x + y)).expand(power_exp=False)
        E**(x + y)*x + E**(x + y)*y
        >>> (exp(x + y)*(x + y)).expand(mul=False)
        E**x*E**y*(x + y)

    - Use deep=False to only expand on the top level::

        >>> exp(x + exp(x + y)).expand()
        E**x*E**(E**x*E**y)
        >>> exp(x + exp(x + y)).expand(deep=False)
        E**(E**(x + y))*E**x

    - Hints are applied in an arbitrary, but consistent order (in the current
      implementation, they are applied in alphabetical order, except
      multinomial comes before mul, but this may change).  Because of this,
      some hints may prevent expansion by other hints if they are applied
      first. For example, ``mul`` may distribute multiplications and prevent
      ``log`` and ``power_base`` from expanding them. Also, if ``mul`` is
      applied before ``multinomial``, the expression might not be fully
      distributed. The solution is to use the various ``expand_hint`` helper
      functions or to use ``hint=False`` to this function to finely control
      which hints are applied. Here are some examples::

        >>> x, y, z = symbols('x y z', positive=True)

        >>> expand(log(x*(y + z)))
        log(x) + log(y + z)

      Here, we see that ``log`` was applied before ``mul``.  To get the mul
      expanded form, either of the following will work::

        >>> expand_mul(log(x*(y + z)))
        log(x*y + x*z)
        >>> expand(log(x*(y + z)), log=False)
        log(x*y + x*z)

      A similar thing can happen with the ``power_base`` hint::

        >>> expand((x*(y + z))**x)
        (x*y + x*z)**x

      To get the ``power_base`` expanded form, either of the following will
      work::

        >>> expand((x*(y + z))**x, mul=False)
        x**x*(y + z)**x
        >>> expand_power_base((x*(y + z))**x)
        x**x*(y + z)**x

        >>> expand((x + y)*y/x)
        y + y**2/x

      The parts of a rational expression can be targeted::

        >>> expand((x + y)*y/x/(x + 1), frac=True)
        (x*y + y**2)/(x**2 + x)
        >>> expand((x + y)*y/x/(x + 1), numer=True)
        (x*y + y**2)/(x*(x + 1))
        >>> expand((x + y)*y/x/(x + 1), denom=True)
        y*(x + y)/(x**2 + x)

    - The ``modulus`` meta-hint can be used to reduce the coefficients of an
      expression post-expansion::

        >>> expand((3*x + 1)**2)
        9*x**2 + 6*x + 1
        >>> expand((3*x + 1)**2, modulus=5)
        4*x**2 + x + 1

    - Either ``expand()`` the function or ``.expand()`` the method can be
      used.  Both are equivalent::

        >>> expand((x + 1)**2)
        x**2 + 2*x + 1
        >>> ((x + 1)**2).expand()
        x**2 + 2*x + 1


    - Objects can define their own expand hints by defining
      ``_eval_expand_hint()``.  The function should take the form::

        def _eval_expand_hint(self, **hints):
            # Only apply the method to the top-level expression
            ...

      See also the example below.  Objects should define ``_eval_expand_hint()``
      methods only if ``hint`` applies to that specific object.  The generic
      ``_eval_expand_hint()`` method defined in Expr will handle the no-op case.

      Each hint should be responsible for expanding that hint only.
      Furthermore, the expansion should be applied to the top-level expression
      only.  ``expand()`` takes care of the recursion that happens when
      ``deep=True``.

      You should only call ``_eval_expand_hint()`` methods directly if you are
      100% sure that the object has the method, as otherwise you are liable to
      get unexpected ``AttributeError``'s.  Note, again, that you do not need to
      recursively apply the hint to args of your object: this is handled
      automatically by ``expand()``.  ``_eval_expand_hint()`` should
      generally not be used at all outside of an ``_eval_expand_hint()`` method.
      If you want to apply a specific expansion from within another method, use
      the public ``expand()`` function, method, or ``expand_hint()`` functions.

      In order for expand to work, objects must be rebuildable by their args,
      i.e., ``obj.func(*obj.args) == obj`` must hold.

      Expand methods are passed ``**hints`` so that expand hints may use
      'metahints'--hints that control how different expand methods are applied.
      For example, the ``force=True`` hint described above that causes
      ``expand(log=True)`` to ignore assumptions is such a metahint.  The
      ``deep`` meta-hint is handled exclusively by ``expand()`` and is not
      passed to ``_eval_expand_hint()`` methods.

      Note that expansion hints should generally be methods that perform some
      kind of 'expansion'.  For hints that simply rewrite an expression, use the
      .rewrite() API.

    Examples
    ========

    >>> class MyClass(Expr):
    ...     def __new__(cls, *args):
    ...         args = sympify(args)
    ...         return Expr.__new__(cls, *args)
    ...
    ...     def _eval_expand_double(self, **hints):
    ...         # Doubles the args of MyClass.
    ...         # If there more than four args, doubling is not performed,
    ...         # unless force=True is also used (False by default).
    ...         force = hints.pop('force', False)
    ...         if not force and len(self.args) > 4:
    ...             return self
    ...         return self.func(*(self.args + self.args))
    ...
    >>> a = MyClass(1, 2, MyClass(3, 4))
    >>> a
    MyClass(1, 2, MyClass(3, 4))
    >>> a.expand(double=True)
    MyClass(1, 2, MyClass(3, 4, 3, 4), 1, 2, MyClass(3, 4, 3, 4))
    >>> a.expand(double=True, deep=False)
    MyClass(1, 2, MyClass(3, 4), 1, 2, MyClass(3, 4))

    >>> b = MyClass(1, 2, 3, 4, 5)
    >>> b.expand(double=True)
    MyClass(1, 2, 3, 4, 5)
    >>> b.expand(double=True, force=True)
    MyClass(1, 2, 3, 4, 5, 1, 2, 3, 4, 5)

    See Also
    ========

    expand_log, expand_mul, expand_multinomial, expand_complex, expand_trig,
    expand_power_base, expand_power_exp, expand_func,
    diofant.simplify.hyperexpand.hyperexpand

    References
    ==========

    * https://mathworld.wolfram.com/Multiple-AngleFormulas.html

    """
    # don't modify this; modify the Expr.expand method
    hints['power_base'] = power_base
    hints['power_exp'] = power_exp
    hints['mul'] = mul
    hints['log'] = log
    hints['multinomial'] = multinomial
    hints['basic'] = basic
    return sympify(e).expand(deep=deep, modulus=modulus, **hints)

# This is a special application of two hints


def _mexpand(expr, recursive=False):
    # expand multinomials and then expand products; this may not always
    # be sufficient to give a fully expanded expression (see
    # test_sympyissue_8247_8354 in test_arit)
    was = None
    while was != expr:
        was, expr = expr, expand_mul(expand_multinomial(expr))
        if not recursive:
            break
    return expr


# These are simple wrappers around single hints.


def expand_mul(expr, deep=True):
    """
    Wrapper around expand that only uses the mul hint.  See the expand
    docstring for more information.

    Examples
    ========

    >>> x, y = symbols('x y', positive=True)
    >>> expand_mul(exp(x+y)*(x+y)*log(x*y**2))
    E**(x + y)*x*log(x*y**2) + E**(x + y)*y*log(x*y**2)

    """
    return sympify(expr).expand(deep=deep, mul=True, power_exp=False,
                                power_base=False, basic=False, multinomial=False, log=False)


def expand_multinomial(expr, deep=True):
    """
    Wrapper around expand that only uses the multinomial hint.  See the expand
    docstring for more information.

    Examples
    ========

    >>> x, y = symbols('x y', positive=True)
    >>> expand_multinomial((x + exp(x + 1))**2)
    2*E**(x + 1)*x + E**(2*x + 2) + x**2

    """
    return sympify(expr).expand(deep=deep, mul=False, power_exp=False,
                                power_base=False, basic=False, multinomial=True, log=False)


def expand_log(expr, deep=True, force=False):
    """
    Wrapper around expand that only uses the log hint.  See the expand
    docstring for more information.

    Examples
    ========

    >>> x, y = symbols('x y', positive=True)
    >>> expand_log(exp(x+y)*(x+y)*log(x*y**2))
    E**(x + y)*(x + y)*(log(x) + 2*log(y))

    """
    return sympify(expr).expand(deep=deep, log=True, mul=False,
                                power_exp=False, power_base=False, multinomial=False,
                                basic=False, force=force)


def expand_func(expr, deep=True):
    """
    Wrapper around expand that only uses the func hint.  See the expand
    docstring for more information.

    Examples
    ========

    >>> expand_func(gamma(x + 2))
    x*(x + 1)*gamma(x)

    """
    return sympify(expr).expand(deep=deep, func=True, basic=False,
                                log=False, mul=False, power_exp=False, power_base=False, multinomial=False)


def expand_trig(expr, deep=True):
    """
    Wrapper around expand that only uses the trig hint.  See the expand
    docstring for more information.

    Examples
    ========

    >>> expand_trig(sin(x+y)*(x+y))
    (x + y)*(sin(x)*cos(y) + sin(y)*cos(x))

    """
    return sympify(expr).expand(deep=deep, trig=True, basic=False,
                                log=False, mul=False, power_exp=False, power_base=False, multinomial=False)


def expand_complex(expr, deep=True):
    """
    Wrapper around expand that only uses the complex hint.  See the expand
    docstring for more information.

    Examples
    ========

    >>> expand_complex(exp(z))
    E**re(z)*I*sin(im(z)) + E**re(z)*cos(im(z))
    >>> expand_complex(sqrt(I))
    sqrt(2)/2 + sqrt(2)*I/2

    See Also
    ========

    diofant.core.expr.Expr.as_real_imag

    """
    return sympify(expr).expand(deep=deep, complex=True, basic=False,
                                log=False, mul=False, power_exp=False, power_base=False, multinomial=False)


def expand_power_base(expr, deep=True, force=False):
    """
    Wrapper around expand that only uses the power_base hint.

    A wrapper to expand(power_base=True) which separates a power with a base
    that is a Mul into a product of powers, without performing any other
    expansions, provided that assumptions about the power's base and exponent
    allow.

    deep=False (default is True) will only apply to the top-level expression.

    force=True (default is False) will cause the expansion to ignore
    assumptions about the base and exponent. When False, the expansion will
    only happen if the base is non-negative or the exponent is an integer.

    >>> (x*y)**2
    x**2*y**2

    >>> (2*x)**y
    (2*x)**y
    >>> expand_power_base(_)
    2**y*x**y

    >>> expand_power_base((x*y)**z)
    (x*y)**z
    >>> expand_power_base((x*y)**z, force=True)
    x**z*y**z
    >>> expand_power_base(sin((x*y)**z), deep=False)
    sin((x*y)**z)
    >>> expand_power_base(sin((x*y)**z), force=True)
    sin(x**z*y**z)

    >>> expand_power_base((2*sin(x))**y + (2*cos(x))**y)
    2**y*sin(x)**y + 2**y*cos(x)**y

    >>> expand_power_base((2*exp(y))**x)
    2**x*(E**y)**x

    >>> expand_power_base((2*cos(x))**y)
    2**y*cos(x)**y

    Notice that sums are left untouched. If this is not the desired behavior,
    apply full ``expand()`` to the expression:

    >>> expand_power_base(((x+y)*z)**2)
    z**2*(x + y)**2
    >>> (((x+y)*z)**2).expand()
    x**2*z**2 + 2*x*y*z**2 + y**2*z**2

    >>> expand_power_base((2*y)**(1+z))
    2**(z + 1)*y**(z + 1)
    >>> ((2*y)**(1+z)).expand()
    2*2**z*y*y**z

    See Also
    ========

    expand

    """
    return sympify(expr).expand(deep=deep, log=False, mul=False,
                                power_exp=False, power_base=True, multinomial=False,
                                basic=False, force=force)


def expand_power_exp(expr, deep=True):
    """
    Wrapper around expand that only uses the power_exp hint.

    Examples
    ========

    >>> expand_power_exp(x**(y + 2))
    x**2*x**y

    See Also
    ========

    expand

    """
    return sympify(expr).expand(deep=deep, complex=False, basic=False,
                                log=False, mul=False, power_exp=True, power_base=False, multinomial=False)


def count_ops(expr, visual=False):
    """
    Return a representation (integer or expression) of the operations in expr.

    If ``visual`` is ``False`` (default) then the sum of the coefficients of the
    visual expression will be returned.

    If ``visual`` is ``True`` then the number of each type of operation is shown
    with the core class types (or their virtual equivalent) multiplied by the
    number of times they occur.

    If expr is an iterable, the sum of the op counts of the
    items will be returned.

    Examples
    ========

    Although there isn't a SUB object, minus signs are interpreted as
    either negations or subtractions:

    >>> (x - y).count_ops(visual=True)
    SUB
    >>> (-x).count_ops(visual=True)
    NEG

    Here, there are two Adds and a Pow:

    >>> (1 + a + b**2).count_ops(visual=True)
    2*ADD + POW

    In the following, an Add, Mul, Pow and two functions:

    >>> (sin(x)*x + sin(x)**2).count_ops(visual=True)
    ADD + MUL + POW + 2*SIN

    for a total of 5:

    >>> (sin(x)*x + sin(x)**2).count_ops(visual=False)
    5

    Note that "what you type" is not always what you get. The expression
    1/x/y is translated by diofant into 1/(x*y) so it gives a DIV and MUL rather
    than two DIVs:

    >>> (1/x/y).count_ops(visual=True)
    DIV + MUL

    The visual option can be used to demonstrate the difference in
    operations for expressions in different forms. Here, the Horner
    representation is compared with the expanded form of a polynomial:

    >>> eq = x*(1 + x*(2 + x*(3 + x)))
    >>> count_ops(eq.expand(), visual=True) - count_ops(eq, visual=True)
    -MUL + 3*POW

    The count_ops function also handles iterables:

    >>> count_ops([x, sin(x), None, True, x + 2], visual=False)
    2
    >>> count_ops([x, sin(x), None, True, x + 2], visual=True)
    ADD + SIN
    >>> count_ops({x: sin(x), x + 2: y + 1}, visual=True)
    2*ADD + SIN

    """
    from ..integrals import Integral
    from ..logic.boolalg import BooleanFunction
    from ..simplify.radsimp import fraction
    from .symbol import Symbol

    expr = sympify(expr)

    if type(expr) is dict:
        ops = [count_ops(k, visual=visual) +
               count_ops(v, visual=visual) for k, v in expr.items()]
    elif iterable(expr):
        ops = [count_ops(i, visual=visual) for i in expr]
    elif isinstance(expr, Expr):

        ops = []
        args = [expr]
        NEG = Symbol('NEG')
        DIV = Symbol('DIV')
        SUB = Symbol('SUB')
        ADD = Symbol('ADD')
        while args:
            a = args.pop()

            if a.is_Rational:
                # -1/3 = NEG + DIV
                if a != 1:
                    if a.numerator < 0:
                        ops.append(NEG)
                    if a.denominator != 1:
                        ops.append(DIV)
                    # XXX "peephole" optimization, http://bugs.python.org/issue2506
                    a
                    continue
            elif a.is_Mul:
                if _coeff_isneg(a):
                    ops.append(NEG)
                    if a.args[0] == -1:
                        a = a.as_two_terms()[1]
                    else:
                        a = -a
                n, d = fraction(a)
                if n.is_Integer:
                    ops.append(DIV)
                    args.append(d)
                    continue  # won't be -Mul but could be Add
                elif d != 1:
                    if not d.is_Integer:
                        args.append(d)
                    ops.append(DIV)
                    args.append(n)
                    continue  # could be -Mul
            elif a.is_Add:
                aargs = list(a.args)
                negs = 0
                for i, ai in enumerate(aargs):
                    if _coeff_isneg(ai):
                        negs += 1
                        args.append(-ai)
                        if i > 0:
                            ops.append(SUB)
                    else:
                        args.append(ai)
                        if i > 0:
                            ops.append(ADD)
                if negs == len(aargs):  # -x - y = NEG + SUB
                    ops.append(NEG)
                elif _coeff_isneg(aargs[0]):  # -x + y = SUB, but already recorded ADD
                    ops.append(SUB - ADD)
                # XXX "peephole" optimization, http://bugs.python.org/issue2506
                a
                continue
            elif isinstance(expr, BooleanFunction):
                ops = []
                for arg in expr.args:
                    ops.append(count_ops(arg, visual=True))
                o = Symbol(expr.func.__name__.upper())
                ops.append(o)
                continue
            if a.is_Pow and a.exp == -1:
                ops.append(DIV)
                args.append(a.base)  # won't be -Mul but could be Add
                continue
            if (a.is_Mul or
                a.is_Pow or
                a.is_Function or
                isinstance(a, Derivative) or
                    isinstance(a, Integral)):

                o = Symbol(a.func.__name__.upper())
                # count the args
                if (a.is_Mul or isinstance(a, LatticeOp)):
                    ops.append(o*(len(a.args) - 1))
                else:
                    ops.append(o)
            if not a.is_Symbol:
                args.extend(a.args)

    elif not isinstance(expr, Basic):
        ops = []
    else:
        ops = []
        args = [expr]
        while args:
            a = args.pop()
            if a.args:
                o = Symbol(a.func.__name__.upper())
                ops.append(o)
                args.extend(a.args)

    if not ops:
        if visual:
            return Integer(0)
        return 0

    ops = Add(*ops)

    if visual:
        return ops

    if ops.is_Number:
        return int(ops)

    return sum(int((a.args or [1])[0]) for a in Add.make_args(ops))


def nfloat(expr, n=15, exponent=False):
    """Make all Rationals in expr Floats except those in exponents
    (unless the exponents flag is set to True).

    Examples
    ========

    >>> nfloat(x**4 + x/2 + cos(pi/3) + 1 + sqrt(y))
    x**4 + 0.5*x + sqrt(y) + 1.5
    >>> nfloat(x**4 + sqrt(y), exponent=True)
    x**4.0 + y**0.5

    """
    from ..polys.rootoftools import RootOf
    from .power import Pow
    from .symbol import Dummy

    if iterable(expr, exclude=(str,)):
        if isinstance(expr, (dict, Dict)):
            return type(expr)([(k, nfloat(v, n, exponent)) for k, v in
                               list(expr.items())])
        return type(expr)([nfloat(a, n, exponent) for a in expr])
    rv = sympify(expr)

    if rv.is_Number:
        return Float(rv, n)
    elif rv.is_number:
        # evalf doesn't always set the precision
        rv = rv.evalf(n)
        if rv.is_Number:
            rv = Float(rv, n)
        else:
            pass  # pure_complex(rv) is likely True
        return rv

    # watch out for RootOf instances that don't like to have
    # their exponents replaced with Dummies and also sometimes have
    # problems with evaluating at low precision (issue sympy/sympy#6393)
    rv = rv.xreplace({ro: ro.evalf(n) for ro in rv.atoms(RootOf)})

    if not exponent:
        reps = [(p, Pow(p.base, Dummy())) for p in rv.atoms(Pow)]
        rv = rv.xreplace(dict(reps))
    rv = rv.evalf(n, strict=False)
    if not exponent:
        rv = rv.xreplace({d.exp: p.exp for p, d in reps})
    else:
        # Pow._eval_evalf special cases Integer exponents so if
        # exponent is suppose to be handled we have to do so here
        rv = rv.xreplace(Transform(
            lambda x: Pow(x.base, Float(x.exp, n)),
            lambda x: x.is_Pow and x.exp.is_Integer))

    return rv.xreplace(Transform(
        lambda x: x.func(*nfloat(x.args, n, exponent)),
        lambda x: isinstance(x, Function)))
