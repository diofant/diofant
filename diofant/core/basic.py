"""Base class for all the objects in Diofant."""

from collections import defaultdict
from collections.abc import Mapping
from itertools import zip_longest

from ..utilities import ordered
from .cache import cacheit
from .compatibility import iterable
from .decorators import _sympifyit
from .evaluate import evaluate
from .sympify import SympifyError, sympify


class Basic:
    """
    Base class for all objects in Diofant.

    Always use ``args`` property, when accessing parameters of some instance.

    """

    # To be overridden with True in the appropriate subclasses
    is_number = False
    is_Atom: bool = False
    is_Symbol = False
    is_Dummy = False
    is_Wild = False
    is_Function = False
    is_Add = False
    is_Mul = False
    is_Pow = False
    is_Exp = False
    is_Number = False
    is_Float = False
    is_Rational = False
    is_Integer = False
    is_NumberSymbol = False
    is_Order = False
    is_Derivative = False
    is_Piecewise = False
    is_Poly = False
    is_Relational = False
    is_Equality = False
    is_Boolean = False
    is_Not = False
    is_Matrix: bool = False
    is_MatMul = False
    is_Vector = False

    def __new__(cls, *args):
        obj = object.__new__(cls)
        obj._hash = None  # will be set by __hash__ method.

        obj._args = args  # all items in args must be Basic objects
        return obj

    def copy(self):
        """Return swallow copy of self."""
        return self.func(*self.args)

    def __reduce_ex__(self, proto):
        """Pickling support."""
        return type(self), self.__getnewargs__(), self.__getstate__()

    def __getnewargs__(self):
        return self.args

    def __getstate__(self):
        return {'_hash': None}

    def __setstate__(self, state):
        for k, v in state.items():
            setattr(self, k, v)

    def __hash__(self):
        # hash cannot be cached using cache_it because infinite recurrence
        # occurs as hash is needed for setting cache dictionary keys
        h = self._hash
        if h is None:
            h = hash((type(self).__name__,) + self._hashable_content())
            self._hash = h
        return h

    def _hashable_content(self):
        """Return a tuple of information about self that can be used to
        compute the hash. If a class defines additional attributes,
        like ``name`` in Symbol, then this method should be updated
        accordingly to return such relevant attributes.

        Defining more than _hashable_content is necessary if __eq__ has
        been defined by a class. See note about this in Basic.__eq__.

        """
        return self._args

    @classmethod
    def class_key(cls):
        """Nice order of classes."""
        return 5, 0, cls.__name__

    @cacheit
    def sort_key(self, order=None):
        """Return a sort key.

        Examples
        ========

        >>> sorted([Rational(1, 2), I, -I], key=lambda x: x.sort_key())
        [1/2, -I, I]

        >>> [x, 1/x, 1/x**2, x**2, sqrt(x), root(x, 4), x**Rational(3, 2)]
        [x, 1/x, x**(-2), x**2, sqrt(x), x**(1/4), x**(3/2)]
        >>> sorted(_, key=lambda x: x.sort_key())
        [x**(-2), 1/x, x**(1/4), sqrt(x), x, x**(3/2), x**2]

        """
        from .numbers import Integer
        args = len(self.args), tuple(arg.sort_key(order)
                                     for arg in self._sorted_args)
        return self.class_key(), args, Integer(1).sort_key(), Integer(1)

    @_sympifyit('other', NotImplemented)
    def __eq__(self, other):
        """Return a boolean indicating whether a == b on the basis of
        their symbolic trees.

        Notes
        =====

        See [1]_.  If a class that overrides __eq__() needs to retain the
        implementation of __hash__() from a parent class, the
        interpreter must be told this explicitly by setting __hash__ =
        <ParentClass>.__hash__. Otherwise the inheritance of __hash__()
        will be blocked, just as if __hash__ had been explicitly set to
        None.

        References
        ==========

        * http://docs.python.org/dev/reference/datamodel.html#object.__hash__

        """
        if self is other:
            return True

        if type(self) != type(other):
            return False

        return self._hashable_content() == other._hashable_content()

    # Note, we always use the default ordering (lex) in __str__ and __repr__,
    # regardless of the global setting.  See issue sympy/sympy#5487.
    def __repr__(self):
        from ..printing import srepr
        return srepr(self, order=None)

    def __str__(self):
        from ..printing import sstr
        return sstr(self, order=None)

    def _repr_pretty_(self, p, cycle):
        from ..printing import pretty
        p.text(pretty(self))

    def _repr_latex_(self):
        from ..printing import latex
        return latex(self, mode='equation*')

    def atoms(self, *types):
        """Returns the atoms that form the current object.

        By default, only objects that are truly atomic and can't
        be divided into smaller pieces are returned: symbols, numbers,
        and number symbols like I and pi. It is possible to request
        atoms of any type, however, as demonstrated below.

        Examples
        ========

        >>> e = 1 + x + 2*sin(y + I*pi)
        >>> e.atoms()
        {1, 2, I, pi, x, y}

        If one or more types are given, the results will contain only
        those types of atoms.

        >>> e.atoms(Symbol)
        {x, y}

        >>> e.atoms(Number)
        {1, 2}

        >>> e.atoms(Number, NumberSymbol)
        {1, 2, pi}

        >>> e.atoms(Number, NumberSymbol, I)
        {1, 2, I, pi}

        Note that I (imaginary unit) and zoo (complex infinity) are special
        types of number symbols and are not part of the NumberSymbol class.

        The type can be given implicitly, too:

        >>> e.atoms(x)
        {x, y}

        Be careful to check your assumptions when using the implicit option
        since ``Integer(1).is_Integer = True`` but ``type(Integer(1))`` is
        ``One``, a special type of diofant atom, while ``type(Integer(2))``
        is type ``Integer`` and will find all integers in an expression:

        >>> e.atoms(Integer(1))
        {1}

        >>> e.atoms(Integer(2))
        {1, 2}

        Finally, arguments to atoms() can select more than atomic atoms: any
        diofant type can be listed as an argument and those types of "atoms"
        as found in scanning the arguments of the expression recursively:

        >>> from diofant.core.function import AppliedUndef

        >>> (1 + x + 2*sin(y + I*pi)).atoms(Mul)
        {I*pi, 2*sin(y + I*pi)}

        >>> f = Function('f')
        >>> e = 1 + f(x) + 2*sin(y + I*pi)
        >>> e.atoms(Function)
        {f(x), sin(y + I*pi)}
        >>> (1 + f(x) + 2*sin(y + I*pi)).atoms(AppliedUndef)
        {f(x)}

        """
        if types:
            types = tuple(t if isinstance(t, type) else type(t) for t in types)
        else:
            types = Atom,
        return set().union(*[self.find(t) for t in types])

    @property
    def free_symbols(self):
        """Return from the atoms of self those which are free symbols.

        For most expressions, all symbols are free symbols. For some classes
        this is not true. e.g. Integrals use Symbols for the dummy variables
        which are bound variables, so Integral has a method to return all
        symbols except those. Derivative keeps track of symbols with respect
        to which it will perform a derivative; those are
        bound variables, too, so it has its own free_symbols method.

        Any other method that uses bound variables should implement a
        free_symbols method.

        """
        return set().union(*[a.free_symbols for a in self.args])

    def rcall(self, *args):
        """Apply on the argument recursively through the expression tree.

        This method is used to simulate a common abuse of notation for
        operators. For instance in Diofant the the following will not work:

        ``(x+Lambda(y, 2*y))(z) == x+2*z``,

        however you can use

        >>> (x + Lambda(y, 2*y)).rcall(z)
        x + 2*z

        """
        if callable(self) and hasattr(self, '__call__'):
            return self(*args)
        if self.args:
            newargs = [sub.rcall(*args) for sub in self.args]
            return type(self)(*newargs)
        return self

    @property
    def func(self):
        """The top-level function in an expression.

        The following should hold for all objects::

            x == x.func(*x.args)

        Examples
        ========

        >>> a = 2*x
        >>> a.func
        <class 'diofant.core.mul.Mul'>
        >>> a.args
        (2, x)
        >>> a.func(*a.args)
        2*x
        >>> a == a.func(*a.args)
        True

        """
        return self.__class__

    @property
    def args(self):
        """Returns a tuple of arguments of 'self'.

        Examples
        ========

        >>> cot(x).args
        (x,)
        >>> (x*y).args
        (x, y)

        """
        return self._args

    @property
    def is_evaluated(self):
        """Test if an expession is evaluated."""
        with evaluate(True):
            expr = self.func(*self.args)
        return expr == self

    @property
    def _sorted_args(self):
        """
        The same as ``args``.  Derived classes which don't fix an
        order on their arguments should override this method to
        produce the sorted representation.

        """
        return self.args

    def subs(self, *args, **kwargs):
        """
        Substitutes old for new in an expression after sympifying args.

        `args` is either:
          - one iterable argument, e.g. foo.subs(iterable). The iterable may be
             o an iterable container with (old, new) pairs. In this case the
               replacements are processed in the order given with successive
               patterns possibly affecting replacements already made.
             o a dict or set whose key/value items correspond to old/new pairs.
               In this case the old/new pairs will be sorted by op count and in
               case of a tie, by number of args and the default_sort_key. The
               resulting sorted list is then processed as an iterable container
               (see previous).

        If the keyword ``simultaneous`` is True, the subexpressions will not be
        evaluated until all the substitutions have been made.

        Examples
        ========

        >>> (1 + x*y).subs({x: pi})
        pi*y + 1
        >>> (1 + x*y).subs({x: pi, y: 2})
        1 + 2*pi
        >>> (1 + x*y).subs([(x, pi), (y, 2)])
        1 + 2*pi
        >>> reps = [(y, x**2), (x, 2)]
        >>> (x + y).subs(reps)
        6
        >>> (x + y).subs(reversed(reps))
        x**2 + 2

        >>> (x**2 + x**4).subs({x**2: y})
        y**2 + y

        To replace only the x**2 but not the x**4, use xreplace:

        >>> (x**2 + x**4).xreplace({x**2: y})
        x**4 + y

        To delay evaluation until all substitutions have been made,
        set the keyword ``simultaneous`` to True:

        >>> (x/y).subs([(x, 0), (y, 0)])
        0
        >>> (x/y).subs([(x, 0), (y, 0)], simultaneous=True)
        nan

        This has the added feature of not allowing subsequent substitutions
        to affect those already made:

        >>> ((x + y)/y).subs({x + y: y, y: x + y})
        1
        >>> ((x + y)/y).subs({x + y: y, y: x + y}, simultaneous=True)
        y/(x + y)

        In order to obtain a canonical result, unordered iterables are
        sorted by count_op length, number of arguments and by the
        default_sort_key to break any ties. All other iterables are left
        unsorted.

        >>> from diofant.abc import e

        >>> expr = sqrt(sin(2*x))*sin(exp(x)*x)*cos(2*x) + sin(2*x)

        >>> expr.subs({sqrt(sin(2*x)): a, sin(2*x): b,
        ...            cos(2*x): c, x: d, exp(x): e})
        a*c*sin(d*e) + b

        The resulting expression represents a literal replacement of the
        old arguments with the new arguments. This may not reflect the
        limiting behavior of the expression:

        >>> (x**3 - 3*x).subs({x: oo})
        nan

        >>> limit(x**3 - 3*x, x, oo)
        oo

        If the substitution will be followed by numerical
        evaluation, it is better to pass the substitution to
        evalf as

        >>> (1/x).evalf(21, subs={x: 3.0}, strict=False)
        0.333333333333333333333

        rather than

        >>> (1/x).subs({x: 3.0}).evalf(21, strict=False)
        0.333333333333333

        as the former will ensure that the desired level of precision is
        obtained.

        See Also
        ========

        replace: replacement capable of doing wildcard-like matching,
                 parsing of match, and conditional replacements
        xreplace: exact node replacement in expr tree; also capable of
                  using matching rules
        diofant.core.evalf.EvalfMixin.evalf: calculates the given formula to
                                           a desired level of precision

        """
        from ..utilities import default_sort_key
        from .numbers import Integer
        from .symbol import Dummy

        unordered = False
        if len(args) == 1:
            sequence = args[0]
            if isinstance(sequence, set):
                unordered = True
            elif isinstance(sequence, Mapping):
                unordered = True
                sequence = sequence.items()
            elif not iterable(sequence):
                raise ValueError('Expected a mapping or iterable '
                                 'of (old, new) tuples.')
            sequence = list(sequence)
        else:
            raise ValueError('subs accepts one argument')

        sequence = [_ for _ in sympify(sequence) if not _aresame(*_)]

        if unordered:
            sequence = dict(sequence)
            if not all(k.is_Atom for k in sequence):
                d = defaultdict(list)
                for o, n in sequence.items():
                    try:
                        ops = o.count_ops(), len(o.args)
                    except TypeError:
                        ops = (0, 0)
                    d[ops].append((o, n))
                newseq = []
                for k in sorted(d, reverse=True):
                    newseq.extend(sorted((v[0] for v in d[k]),
                                         key=default_sort_key))
                sequence = [(k, sequence[k]) for k in newseq]
                del newseq, d
            else:
                sequence = sorted(((k, v) for (k, v) in sequence.items()),
                                  key=default_sort_key)

        if kwargs.pop('simultaneous', False):  # XXX should this be the default for dict subs?
            reps = {}
            rv = self
            m = Dummy()
            for old, new in sequence:
                d = Dummy(commutative=new.is_commutative)
                # using d*m so Subs will be used on dummy variables
                # in things like Derivative(f(x, y), x) in which x
                # is both free and bound
                rv = rv._subs(old, d*m, **kwargs)
                reps[d] = new
            reps[m] = Integer(1)  # get rid of m
            return rv.xreplace(reps)
        rv = self
        for old, new in sequence:
            rv = rv._subs(old, new, **kwargs)
            if not isinstance(rv, Basic):
                break
        return rv

    @cacheit
    def _subs(self, old, new, **hints):
        """Substitutes an expression old -> new.

        If self is not equal to old then _eval_subs is called.
        If _eval_subs doesn't want to make any special replacement
        then a None is received which indicates that the fallback
        should be applied wherein a search for replacements is made
        amongst the arguments of self.

        Examples
        ========

        Add's _eval_subs knows how to target x + y in the following
        so it makes the change:

            >>> (x + y + z).subs({x + y: 1})
            z + 1

        Add's _eval_subs doesn't need to know how to find x + y in
        the following:

            >>> Add._eval_subs(z*(x + y) + 3, x + y, 1) is None
            True

        The returned None will cause the fallback routine to traverse the args and
        pass the z*(x + y) arg to Mul where the change will take place and the
        substitution will succeed:

            >>> (z*(x + y) + 3).subs({x + y: 1})
            z + 3

        ** Developers Notes **

        An _eval_subs routine for a class should be written if:

            1) any arguments are not instances of Basic (e.g. bool, tuple);

            2) some arguments should not be targeted (as in integration
               variables);

            3) if there is something other than a literal replacement
               that should be attempted (as in Piecewise where the condition
               may be updated without doing a replacement).

        If it is overridden, here are some special cases that might arise:

            1) If it turns out that no special change was made and all
               the original sub-arguments should be checked for
               replacements then None should be returned.

            2) If it is necessary to do substitutions on a portion of
               the expression then _subs should be called. _subs will
               handle the case of any sub-expression being equal to old
               (which usually would not be the case) while its fallback
               will handle the recursion into the sub-arguments. For
               example, after Add's _eval_subs removes some matching terms
               it must process the remaining terms so it calls _subs
               on each of the un-matched terms and then adds them
               onto the terms previously obtained.

           3) If the initial expression should remain unchanged then
              the original expression should be returned. (Whenever an
              expression is returned, modified or not, no further
              substitution of old -> new is attempted.) Sum's _eval_subs
              routine uses this strategy when a substitution is attempted
              on any of its summation variables.

        """

        def fallback(self, old, new):
            """Try to replace old with new in any of self's arguments."""
            hit = False
            args = list(self.args)
            for i, arg in enumerate(args):
                arg = arg._subs(old, new, **hints)
                if not _aresame(arg, args[i]):
                    hit = True
                    args[i] = arg
            if hit:
                return self.func(*args)
            return self

        if _aresame(self, old):
            return new

        rv = self._eval_subs(old, new)  # pylint: disable=assignment-from-none
        if rv is None:
            rv = fallback(self, old, new)
        return rv

    def _eval_subs(self, old, new):
        """Override this stub if you want to do anything more than
        attempt a replacement of old with new in the arguments of self.

        See also
        ========

        _subs

        """
        return

    def xreplace(self, rule):
        """
        Replace occurrences of objects within the expression.

        Parameters
        ==========

        rule : dict-like
            Expresses a replacement rule

        Returns
        =======

        xreplace : the result of the replacement

        Examples
        ========

        >>> (1 + x*y).xreplace({x: pi})
        pi*y + 1
        >>> (1 + x*y).xreplace({x: pi, y: 2})
        1 + 2*pi

        Replacements occur only if an entire node in the expression tree is
        matched:

        >>> (x*y + z).xreplace({x*y: pi})
        z + pi
        >>> (x*y*z).xreplace({x*y: pi})
        x*y*z
        >>> (2*x).xreplace({2*x: y, x: z})
        y
        >>> (2*2*x).xreplace({2*x: y, x: z})
        4*z
        >>> (x + y + 2).xreplace({x + y: 2})
        x + y + 2
        >>> (x + 2 + exp(x + 2)).xreplace({x + 2: y})
        E**y + x + 2

        xreplace doesn't differentiate between free and bound symbols. In the
        following, subs(x, y) would not change x since it is a bound symbol,
        but xreplace does:

        >>> Integral(x, (x, 1, 2*x)).xreplace({x: y})
        Integral(y, (y, 1, 2*y))

        Trying to replace x with an expression raises an error:

        >>> Integral(x, (x, 1, 2*x)).xreplace({x: 2*y})
        Traceback (most recent call last):
        ...
        ValueError: Invalid limits given: ((2*y, 1, 4*y),)

        See Also
        ========

        replace: replacement capable of doing wildcard-like matching,
                 parsing of match, and conditional replacements
        subs: substitution of subexpressions as defined by the objects
              themselves.

        """
        if self in rule:
            return rule[self]
        if rule and not self.is_Atom:
            args = tuple(a.xreplace(rule) for a in self.args)
            if not _aresame(args, self.args):
                return self.func(*args)
        return self

    @cacheit
    def has(self, *patterns):
        r"""Test if any subexpression matches any of the patterns.

        Parameters
        ==========

        \*patterns : tuple of Expr
            List of expressions to search for match.

        Returns
        =======

        bool
            False if there is no match or patterns list is
            empty, else True.

        Examples
        ========

        >>> e = x**2 + sin(x*y)
        >>> e.has(z)
        False
        >>> e.has(x, y, z)
        True
        >>> x.has()
        False

        """
        from .function import Function, UndefinedFunction

        if len(patterns) != 1:
            return any(self.has(pattern) for pattern in patterns)
        pattern = sympify(patterns[0])
        if isinstance(pattern, UndefinedFunction):
            return any(pattern in (f, f.func)
                       for f in self.atoms(Function, UndefinedFunction))
        if isinstance(pattern, type):
            return any(isinstance(arg, pattern)
                       for arg in preorder_traversal(self))
        match = pattern._has_matcher()
        return any(match(arg) for arg in preorder_traversal(self))

    def _has_matcher(self):
        """Helper for .has()."""
        return lambda x: self == x

    def replace(self, query, value, exact=False):
        """Replace matching subexpressions of ``self`` with ``value``.

        Traverses an expression tree and performs replacement of matching
        subexpressions from the bottom to the top of the tree in a simultaneous
        fashion so changes made are targeted only once. In addition, if an
        expression containing more than one Wild symbol is being used to match
        subexpressions and  the ``exact`` flag is True, then the match will
        only succeed if non-zero values are received for each Wild that appears
        in the match pattern.

        The list of possible combinations of queries and replacement values
        is listed below:

        Examples
        ========

        Initial setup

            >>> f = log(sin(x)) + tan(sin(x**2))

        1.1. type -> type
            obj.replace(type, newtype)

            When object of type ``type`` is found, replace it with the
            result of passing its argument(s) to ``newtype``.

            >>> f.replace(sin, cos)
            log(cos(x)) + tan(cos(x**2))
            >>> (x*y).replace(Mul, Add)
            x + y

        1.2. type -> func
            obj.replace(type, func)

            When object of type ``type`` is found, apply ``func`` to its
            argument(s). ``func`` must be written to handle the number
            of arguments of ``type``.

            >>> f.replace(sin, lambda arg: sin(2*arg))
            log(sin(2*x)) + tan(sin(2*x**2))
            >>> (x*y).replace(Mul, lambda *args: sin(2*Mul(*args)))
            sin(2*x*y)

        2.1. pattern -> expr
            obj.replace(pattern(wild), expr(wild))

            Replace subexpressions matching ``pattern`` with the expression
            written in terms of the Wild symbols in ``pattern``.

            >>> a = Wild('a')
            >>> f.replace(sin(a), tan(a))
            log(tan(x)) + tan(tan(x**2))
            >>> f.replace(sin(a), tan(a/2))
            log(tan(x/2)) + tan(tan(x**2/2))
            >>> f.replace(sin(a), a)
            log(x) + tan(x**2)
            >>> (x*y).replace(a*x, a)
            y

            When the default value of False is used with patterns that have
            more than one Wild symbol, non-intuitive results may be obtained:

            >>> b = Wild('b')
            >>> (2*x).replace(a*x + b, b - a)
            2/x

            For this reason, the ``exact`` option can be used to make the
            replacement only when the match gives non-zero values for all
            Wild symbols:

            >>> (2*x + y).replace(a*x + b, b - a, exact=True)
            y - 2
            >>> (2*x).replace(a*x + b, b - a, exact=True)
            2*x

        2.2. pattern -> func
            obj.replace(pattern(wild), lambda wild: expr(wild))

            All behavior is the same as in 2.1 but now a function in terms of
            pattern variables is used rather than an expression:

            >>> f.replace(sin(a), lambda a: sin(2*a))
            log(sin(2*x)) + tan(sin(2*x**2))

        3.1. func -> func
            obj.replace(filter, func)

            Replace subexpression ``e`` with ``func(e)`` if ``filter(e)``
            is True.

            >>> g = 2*sin(x**3)
            >>> g.replace(lambda expr: expr.is_Number, lambda expr: expr**2)
            4*sin(x**9)

        The expression itself is also targeted by the query but is done in
        such a fashion that changes are not made twice.

            >>> e = x*(x*y + 1)
            >>> e.replace(lambda x: x.is_Mul, lambda x: 2*x)
            2*x*(2*x*y + 1)

        See Also
        ========

        subs: substitution of subexpressions as defined by the objects
              themselves.
        xreplace: exact node replacement in expr tree; also capable of
                  using matching rules

        """
        from ..simplify.simplify import bottom_up

        try:
            query = sympify(query)
        except SympifyError:
            pass
        try:
            value = sympify(value)
        except SympifyError:
            pass
        if isinstance(query, type):
            def _query(expr):
                return isinstance(expr, query)

            if isinstance(value, type) or callable(value):
                def _value(expr, result):
                    return value(*expr.args)
            else:
                raise TypeError(
                    'given a type, replace() expects another '
                    'type or a callable')
        elif isinstance(query, Basic):
            def _query(expr):
                return expr.match(query)

            # XXX remove the exact flag and make multi-symbol
            # patterns use exact=True semantics; to do this the query must
            # be tested to find out how many Wild symbols are present.
            # See https://groups.google.com/forum/
            # ?fromgroups=#!topic/sympy/zPzo5FtRiqI
            # for a method of inspecting a function to know how many
            # parameters it has.
            if isinstance(value, Basic):
                if exact:
                    def _value(expr, result):
                        return (value.subs(result)
                                if all(val for val in result.values()) else expr)
                else:
                    def _value(expr, result):
                        return value.subs(result)
            elif callable(value):
                # match dictionary keys get the trailing underscore stripped
                # from them and are then passed as keywords to the callable;
                # if ``exact`` is True, only accept match if there are no null
                # values amongst those matched.
                if exact:
                    def _value(expr, result):
                        return (value(**{str(key)[:-1]: val for key, val in result.items()})
                                if all(val for val in result.values()) else expr)
                else:
                    def _value(expr, result):
                        return value(**{str(key)[:-1]: val for key, val in result.items()})
            else:
                raise TypeError(
                    'given an expression, replace() expects '
                    'another expression or a callable')
        elif callable(query):
            _query = query

            if callable(value):
                def _value(expr, result):
                    return value(expr)
            else:
                raise TypeError(
                    'given a callable, replace() expects '
                    'another callable')
        else:
            raise TypeError(
                'first argument to replace() must be a '
                'type, an expression or a callable')

        def rec_replace(expr):
            result = _query(expr)
            if result or result == {}:
                new = _value(expr, result)
                if new is not None and new != expr:
                    expr = new
            return expr

        return bottom_up(self, rec_replace, atoms=True)

    def find(self, query):
        """Find all subexpressions matching a query."""
        try:
            query = sympify(query)
        except SympifyError:
            pass
        if isinstance(query, type):
            def _query(expr):
                return isinstance(expr, query)
        elif isinstance(query, Basic):
            def _query(expr):
                return expr.match(query) is not None
        else:
            _query = query

        groups = defaultdict(int)
        for result in filter(_query, preorder_traversal(self)):
            groups[result] += 1
        return dict(groups)

    def count(self, query):
        """Count the number of matching subexpressions."""
        return sum(self.find(query).values())

    def _matches(self, expr, repl_dict={}):
        """Helper method for match() that looks for a match between Wild
        symbols in self and expressions in expr.

        Examples
        ========

        >>> x = Wild('x')
        >>> Basic(a + x, x)._matches(Basic(a + b, c)) is None
        True
        >>> Basic(a + x, x)._matches(Basic(a + b + c, b + c))
        {x_: b + c}

        """
        expr = sympify(expr)
        if not isinstance(expr, self.func):
            return

        if self == expr:
            return repl_dict

        if self.is_Atom:
            return

        if len(self.args) != len(expr.args):
            return

        d = repl_dict.copy()
        for arg, other_arg in zip(self.args, expr.args):
            if arg == other_arg:
                continue
            d = arg.xreplace(d)._matches(other_arg, d)
            if d is None:
                return
        return d

    def match(self, pattern):
        """Pattern matching.

        Wild symbols match all.

        Parameters
        ==========

        pattern : Expr
            An expression that may contain Wild symbols.

        Returns
        =======

        dict or None
            If pattern match self, return a dictionary of
            replacement rules, such that::

                pattern.xreplace(self.match(pattern)) == self

        Examples
        ========

        >>> p = Wild('p')
        >>> q = Wild('q')
        >>> e = (x + y)**(x + y)
        >>> e.match(p**p)
        {p_: x + y}
        >>> e.match(p**q)
        {p_: x + y, q_: x + y}
        >>> (p**q).xreplace(_)
        (x + y)**(x + y)

        See Also
        ========

        xreplace
        diofant.core.symbol.Wild

        """
        from ..simplify import signsimp
        pattern = sympify(pattern)
        s = signsimp(self)
        p = signsimp(pattern)
        # if we still have the same relationship between the types of
        # input, then use the sign simplified forms
        if (pattern.func == self.func) and (s.func == p.func):
            rv = p._matches(s)
        else:
            rv = pattern._matches(self)
        return rv

    def count_ops(self, visual=None):
        """Wrapper for count_ops that returns the operation count."""
        from .function import count_ops
        return count_ops(self, visual)

    def doit(self, **hints):
        """Evaluate objects that are not evaluated by default.

        For example, limits, integrals, sums and products.  All objects of this
        kind will be evaluated recursively, unless some species were excluded
        via 'hints' or unless the 'deep' hint was set to 'False'.

        Examples
        ========

        >>> 2*Integral(x, x)
        2*Integral(x, x)

        >>> (2*Integral(x, x)).doit()
        x**2

        >>> (2*Integral(x, x)).doit(deep=False)
        2*Integral(x, x)

        """
        if hints.get('deep', True):
            terms = [term.doit(**hints) if isinstance(term, Basic) else term
                     for term in self.args]
            return self.func(*terms)
        return self

    def _eval_rewrite(self, pattern, rule, **hints):
        if self.is_Atom:
            if hasattr(self, rule):
                return getattr(self, rule)()
            return self

        if hints.get('deep', True):
            args = [a._eval_rewrite(pattern, rule, **hints)
                    if isinstance(a, Basic) else a
                    for a in self.args]
        else:
            args = self.args

        if pattern is None or isinstance(self, pattern):
            if hasattr(self, rule):
                rewritten = getattr(self, rule)(*args, **hints)
                if rewritten is not None:
                    return rewritten
        return self.func(*args)

    def rewrite(self, *args, **hints):
        """Rewrite functions in terms of other functions.

        Rewrites expression containing applications of functions
        of one kind in terms of functions of different kind. For
        example you can rewrite trigonometric functions as complex
        exponentials or combinatorial functions as gamma function.

        As a pattern this function accepts a list of functions to
        to rewrite (instances of DefinedFunction class). As rule
        you can use string or a destination function instance (in
        this case rewrite() will use the str() function).

        There is also the possibility to pass hints on how to rewrite
        the given expressions. For now there is only one such hint
        defined called 'deep'. When 'deep' is set to False it will
        forbid functions to rewrite their contents.

        Examples
        ========

        Unspecified pattern:

        >>> sin(x).rewrite(exp)
        -I*(E**(I*x) - E**(-I*x))/2

        Pattern as a single function:

        >>> sin(x).rewrite(sin, exp)
        -I*(E**(I*x) - E**(-I*x))/2

        Pattern as a list of functions:

        >>> sin(x).rewrite([sin], exp)
        -I*(E**(I*x) - E**(-I*x))/2

        """
        if not args:
            return self
        pattern = args[:-1]
        if isinstance(args[-1], str):
            rule = '_eval_rewrite_as_' + args[-1]
        else:
            rule = '_eval_rewrite_as_' + args[-1].__name__

        if not pattern:
            return self._eval_rewrite(None, rule, **hints)
        if iterable(pattern[0]):
            pattern = pattern[0]

        pattern = [p for p in pattern if self.has(p)]

        if pattern:
            return self._eval_rewrite(tuple(pattern), rule, **hints)
        return self


class Atom(Basic):
    """A parent class for atomic things.

    An atom is an expression with no subexpressions, for example Symbol,
    Number, Rational or Integer, but not Add, Mul, Pow.

    """

    is_Atom = True

    def doit(self, **hints):
        """Evaluate objects that are not evaluated by default.

        See Also
        ========

        Basic.doit

        """
        return self

    @classmethod
    def class_key(cls):
        """Nice order of classes."""
        return 2, 0, cls.__name__

    @cacheit
    def sort_key(self, order=None):
        """Return a sort key."""
        from . import Integer
        return self.class_key(), (1, (str(self),)), Integer(1).sort_key(), Integer(1)

    def _eval_simplify(self, ratio, measure):
        return self

    @property
    def _sorted_args(self):
        # this is here as a safeguard against accidentally using _sorted_args
        # on Atoms -- they cannot be rebuilt as atom.func(*atom._sorted_args)
        # since there are no args. So the calling routine should be checking
        # to see that this property is not called for Atoms.
        raise AttributeError('Atoms have no args. It might be necessary'
                             ' to make a check for Atoms in the calling code.')


def _aresame(a, b):
    """Return True if a and b are structurally the same, else False.

    Examples
    ========

    To Diofant, 2.0 == 2:

    >>> 2.0 == Integer(2)
    True

    Since a simple 'same or not' result is sometimes useful, this routine was
    written to provide that query:

    >>> _aresame(Float(2.0), Integer(2))
    False

    """
    from .function import AppliedUndef
    from .function import UndefinedFunction as UndefFunc
    for i, j in zip_longest(preorder_traversal(a), preorder_traversal(b)):
        if i != j or type(i) != type(j):
            if (isinstance(i, (UndefFunc, AppliedUndef)) and
                    isinstance(j, (UndefFunc, AppliedUndef))):
                if i.class_key() != j.class_key():
                    return False
            else:
                return False
    return True


class preorder_traversal:
    """Do a pre-order traversal of a tree.

    This iterator recursively yields nodes that it has visited in a pre-order
    fashion. That is, it yields the current node then descends through the
    tree breadth-first to yield all of a node's children's pre-order
    traversal.

    For an expression, the order of the traversal depends on the order of
    .args, which in many cases can be arbitrary.

    Parameters
    ==========

    node : diofant expression
        The expression to traverse.
    keys : (default None) sort key(s)
        The key(s) used to sort args of Basic objects. When None, args of Basic
        objects are processed in arbitrary order. If key is defined, it will
        be passed along to ordered() as the only key(s) to use to sort the
        arguments; if ``key`` is simply True then the default keys of ordered
        will be used.

    Yields
    ======

    subtree : diofant expression
        All of the subtrees in the tree.

    Examples
    ========

    The nodes are returned in the order that they are encountered unless key
    is given; simply passing key=True will guarantee that the traversal is
    unique.

    >>> list(preorder_traversal((x + y)*z, keys=True))
    [z*(x + y), z, x + y, x, y]

    """

    def __init__(self, node, keys=None):
        """Initialize self."""
        self._skip_flag = False
        self._pt = self._preorder_traversal(node, keys)

    def _preorder_traversal(self, node, keys):
        yield node
        if self._skip_flag:
            self._skip_flag = False
            return
        if isinstance(node, Basic):
            args = node.args
            if keys:
                args = ordered(args)
            for arg in args:
                yield from self._preorder_traversal(arg, keys)
        elif iterable(node):
            for item in node:
                yield from self._preorder_traversal(item, keys)

    def skip(self):
        """
        Skip yielding current node's (last yielded node's) subtrees.

        Examples
        ========

        >>> pt = preorder_traversal((x+y*z)*z)
        >>> for i in pt:
        ...     print(i)
        ...     if i == x + y*z:
        ...         pt.skip()
        z*(x + y*z)
        z
        x + y*z

        """
        self._skip_flag = True

    def __next__(self):
        return next(self._pt)

    def __iter__(self):
        return self
