import itertools
import re as _re
import string

from ..logic.boolalg import Boolean
from .assumptions import StdFactKB
from .cache import cacheit
from .expr import AtomicExpr, Expr
from .logic import fuzzy_bool
from .numbers import Integer
from .sympify import sympify


class BaseSymbol(AtomicExpr, Boolean):
    """Abstract class for Symbols.

    Do not instantiate this, use derived classes.

    Notes
    =====

    We introduce this class to prevent "flipping" of arguments
    for "rich comparison" methods [1]_.  There is no swapped-argument
    versions of these methods, like :meth:`~object.__add__` vs
    :meth:`~object.__radd__`, rather e.g. :meth:`~object.__lt__` and
    :meth:`~object.__gt__` are each other's reflection.

    According to the documentation [1]_, if the operands of such
    method are of different types, and right operand’s type is a direct or
    indirect subclass of the left operand’s type, the reflected method
    of the right operand has priority, otherwise the left operand’s
    method has priority.

    Thus, simple class hierarhy, where :class:`Symbol` is a parent
    class for both :class:`Dummy` and :class:`Wild` isn't possible,
    if we want to avoid silent switching to the reflected methods
    for rich comparisons of parent and child.  See also sympy/sympy#7951.

    Examples
    ========

    We illustrate the "flipping" problem, by using here
    BaseSymbol in place of ordinary Symbol:

    >>> p = Wild('p')
    >>> x = BaseSymbol('x')
    >>> x < p
    p_ > x

    See Also
    ========

    Symbol
    Dummy
    Wild

    References
    ==========

    * https://docs.python.org/3/reference/datamodel.html#object.__lt__

    """

    is_comparable = False

    is_Symbol = True

    @property
    def _diff_wrt(self):
        """Allow derivatives wrt Symbols.

        Examples
        ========

            >>> x._diff_wrt
            True

        """
        return True

    @staticmethod
    def _sanitize(assumptions, obj=None):
        """Remove None, covert values to bool, check commutativity *in place*."""
        # be strict about commutativity: cannot be None
        is_commutative = fuzzy_bool(assumptions.get('commutative', True))
        if is_commutative is None:
            whose = f'{obj.__name__} ' if obj else ''
            raise ValueError(
                f'{whose}commutativity must be True or False.')

        # sanitize other assumptions so 1 -> True and 0 -> False
        for key in list(assumptions):
            v = assumptions[key]
            if v is None:
                assumptions.pop(key)
                continue
            assumptions[key] = bool(v)

    def __new__(cls, name, **assumptions):
        """Symbols are identified by name and assumptions::

        >>> Symbol('x') == Symbol('x')
        True
        >>> Symbol('x', real=True) == Symbol('x', real=False)
        False

        """
        cls._sanitize(assumptions, cls)
        return BaseSymbol.__xnew_cached_(cls, name, **assumptions)

    def __new_stage2__(cls, name, **assumptions):
        if not isinstance(name, str):
            raise TypeError('name should be a string, not %s' % repr(type(name)))

        obj = Expr.__new__(cls)
        obj.name = name

        # TODO: Issue sympy/sympy#8873: Forcing the commutative assumption here means
        # later code such as ``repr()`` cannot tell whether the user
        # specified ``commutative=True`` or omitted it.  To workaround this,
        # we keep a copy of the assumptions dict, then create the StdFactKB,
        # and finally overwrite its ``._generator`` with the dict copy.  This
        # is a bit of a hack because we assume StdFactKB merely copies the
        # given dict as ``._generator``, but future modification might, e.g.,
        # compute a minimal equivalent assumption set.
        tmp_asm_copy = assumptions.copy()

        # be strict about commutativity
        is_commutative = fuzzy_bool(assumptions.get('commutative', True))
        assumptions['commutative'] = is_commutative
        obj._assumptions = StdFactKB(assumptions)
        obj._assumptions._generator = tmp_asm_copy  # Issue sympy/sympy#8873
        return obj

    __xnew__ = staticmethod(
        __new_stage2__)            # never cached (e.g. dummy)
    __xnew_cached_ = staticmethod(
        cacheit(__new_stage2__))   # symbols are always cached

    def __getnewargs__(self):
        return self.name,

    def __getstate__(self):
        state = super().__getstate__()
        state['_assumptions'] = self._assumptions
        return state

    def _hashable_content(self):
        # Note: user-specified assumptions not hashed, just derived ones
        return ((self.name,) +
                tuple(sorted((k, v) for k, v in self._assumptions.items()
                             if v is not None)))

    @classmethod
    def class_key(cls):
        """Nice order of classes."""
        return 2, 0, cls.__name__

    @cacheit
    def sort_key(self, order=None):
        """Return a sort key."""
        return self.class_key(), (1, (str(self),)), Integer(1).sort_key(), Integer(1)

    def as_dummy(self):
        """Return a Dummy having the same name and same assumptions as self."""
        return Dummy(self.name, **self._assumptions.generator)

    def is_constant(self, *wrt, **flags):
        """Test if self is constant.

        See Also
        ========

        diofant.core.expr.Expr.is_constant

        """
        if not wrt:
            return False
        return self not in wrt

    @property
    def free_symbols(self):
        """Return from the atoms of self those which are free symbols.

        See Also
        ========

        diofant.core.basic.Basic.free_symbols

        """
        return {self}


class Symbol(BaseSymbol):
    r"""Symbol is a placeholder for atomic symbolic expression.

    It has a name and a set of assumptions.

    Parameters
    ==========

    name : str
        The name for Symbol.
    \*\*assumptions : dict
        Keyword arguments to specify assumptions for
        Symbol.  Default assumption is commutative=True.

    Examples
    ========

    >>> a, b = symbols('a b')
    >>> bool(a*b == b*a)
    True

    You can override default assumptions:

    >>> A, B = symbols('A B', commutative=False)
    >>> bool(A*B != B*A)
    True
    >>> bool(A*B*2 == 2*A*B) is True  # multiplication by scalars is commutative
    True

    See Also
    ========

    :mod:`diofant.core.assumptions`
    Dummy
    Wild

    """


class Dummy(BaseSymbol):
    """Dummy symbols are each unique, identified by an internal count index:

    >>> bool(Dummy('x') == Dummy('x')) is True
    False

    If a name is not supplied then a string value of the count index will be
    used. This is useful when a temporary variable is needed and the name
    of the variable used in the expression is not important.

    >>> Dummy()  # doctest: +SKIP
    _Dummy_10

    See Also
    ========

    Symbol

    """

    _count = 0

    is_Dummy = True

    def __new__(cls, name=None, **assumptions):
        if name is None:
            name = 'Dummy_' + str(Dummy._count)

        cls._sanitize(assumptions, cls)
        obj = Symbol.__xnew__(cls, name, **assumptions)

        Dummy._count += 1
        obj.dummy_index = Dummy._count
        return obj

    def __getstate__(self):
        state = super().__getstate__()
        state['dummy_index'] = self.dummy_index
        return state

    @classmethod
    def class_key(cls):
        """Nice order of classes."""
        return 3, 0, cls.__name__

    @cacheit
    def sort_key(self, order=None):
        """Return a sort key."""
        return self.class_key(), (
            2, (str(self), self.dummy_index)), Integer(1).sort_key(), Integer(1)

    def _hashable_content(self):
        return BaseSymbol._hashable_content(self) + (self.dummy_index,)


class Wild(BaseSymbol):
    """A Wild symbol matches anything, whatever is not explicitly excluded.

    Examples
    ========

    >>> a = Wild('a')
    >>> x.match(a)
    {a_: x}
    >>> pi.match(a)
    {a_: pi}
    >>> (3*x**2).match(a*x)
    {a_: 3*x}
    >>> cos(x).match(a)
    {a_: cos(x)}
    >>> b = Wild('b', exclude=[x])
    >>> (3*x**2).match(b*x)
    >>> b.match(a)
    {a_: b_}
    >>> A = WildFunction('A')
    >>> A.match(a)
    {a_: A_}

    Notes
    =====

    When using Wild, be sure to use the exclude
    keyword to make the pattern more precise.
    Without the exclude pattern, you may get matches
    that are technically correct, but not what you
    wanted. For example, using the above without
    exclude:

    >>> a, b = symbols('a b', cls=Wild)
    >>> (2 + 3*y).match(a*x + b*y)
    {a_: 2/x, b_: 3}

    This is technically correct, because
    (2/x)*x + 3*y == 2 + 3*y, but you probably
    wanted it to not match at all. The issue is that
    you really didn't want a and b to include x and y,
    and the exclude parameter lets you specify exactly
    this.  With the exclude parameter, the pattern will
    not match.

    >>> a = Wild('a', exclude=[x, y])
    >>> b = Wild('b', exclude=[x, y])
    >>> (2 + 3*y).match(a*x + b*y)

    Exclude also helps remove ambiguity from matches.

    >>> E = 2*x**3*y*z
    >>> a, b = symbols('a b', cls=Wild)
    >>> E.match(a*b)
    {a_: 2*y*z, b_: x**3}
    >>> a = Wild('a', exclude=[x, y])
    >>> E.match(a*b)
    {a_: z, b_: 2*x**3*y}
    >>> a = Wild('a', exclude=[x, y, z])
    >>> E.match(a*b)
    {a_: 2, b_: x**3*y*z}

    See Also
    ========

    Symbol

    """

    is_Wild = True

    def __new__(cls, name, exclude=(), properties=(), **assumptions):
        exclude = tuple(sympify(x) for x in exclude)
        properties = tuple(properties)
        cls._sanitize(assumptions, cls)
        return Wild.__xnew__(cls, name, exclude, properties, **assumptions)

    def __getnewargs__(self):
        return self.name, self.exclude, self.properties

    @staticmethod
    @cacheit
    def __xnew__(cls, name, exclude, properties, **assumptions):
        obj = BaseSymbol.__xnew__(cls, name, **assumptions)
        obj.exclude = exclude
        obj.properties = properties
        return obj

    def _hashable_content(self):
        return super()._hashable_content() + (self.exclude, self.properties)

    # TODO add check against another Wild
    def _matches(self, expr, repl_dict={}):
        """Helper method for match().

        See Also
        ========

        diofant.core.basic.Basic.matches

        """
        if any(expr.has(x) for x in self.exclude):
            return
        if any(not f(expr) for f in self.properties):
            return
        repl_dict = repl_dict.copy()
        repl_dict[self] = expr
        return repl_dict


_range = _re.compile('([0-9]*:[0-9]+|[a-zA-Z]?:[a-zA-Z])')


def symbols(names, **args):
    r"""
    Transform strings into instances of :class:`Symbol` class.

    :func:`symbols` function returns a sequence of symbols with names taken
    from ``names`` argument, which can be a comma or whitespace delimited
    string, or a sequence of strings::

        >>> a, b, c = symbols('a b c')

    The type of output is dependent on the properties of input arguments::

        >>> symbols('x')
        x
        >>> symbols('x,')
        (x,)
        >>> symbols('x,y')
        (x, y)
        >>> symbols(('a', 'b', 'c'))
        (a, b, c)
        >>> symbols(['a', 'b', 'c'])
        [a, b, c]
        >>> symbols({'a', 'b', 'c'})
        {a, b, c}

    If an iterable container is needed for a single symbol, set the ``seq``
    argument to ``True`` or terminate the symbol name with a comma::

        >>> symbols('x', seq=True)
        (x,)

    To reduce typing, range syntax is supported to create indexed symbols.
    Ranges are indicated by a colon and the type of range is determined by
    the character to the right of the colon. If the character is a digit
    then all contiguous digits to the left are taken as the nonnegative
    starting value (or 0 if there is no digit left of the colon) and all
    contiguous digits to the right are taken as 1 greater than the ending
    value::

        >>> symbols('x:10')
        (x0, x1, x2, x3, x4, x5, x6, x7, x8, x9)

        >>> symbols('x5:10')
        (x5, x6, x7, x8, x9)
        >>> symbols('x5(:2)')
        (x50, x51)

        >>> symbols('x5:10 y:5')
        (x5, x6, x7, x8, x9, y0, y1, y2, y3, y4)

        >>> symbols(('x5:10', 'y:5'))
        ((x5, x6, x7, x8, x9), (y0, y1, y2, y3, y4))

    If the character to the right of the colon is a letter, then the single
    letter to the left (or 'a' if there is none) is taken as the start
    and all characters in the lexicographic range *through* the letter to
    the right are used as the range::

        >>> symbols('x:z')
        (x, y, z)
        >>> symbols('x:c')  # null range
        ()
        >>> symbols('x(:c)')
        (xa, xb, xc)

        >>> symbols(':c')
        (a, b, c)

        >>> symbols('a:d, x:z')
        (a, b, c, d, x, y, z)

        >>> symbols(('a:d', 'x:z'))
        ((a, b, c, d), (x, y, z))

    Multiple ranges are supported; contiguous numerical ranges should be
    separated by parentheses to disambiguate the ending number of one
    range from the starting number of the next::

        >>> symbols('x:2(1:3)')
        (x01, x02, x11, x12)
        >>> symbols(':3:2')  # parsing is from left to right
        (00, 01, 10, 11, 20, 21)

    Only one pair of parentheses surrounding ranges are removed, so to
    include parentheses around ranges, double them. And to include spaces,
    commas, or colons, escape them with a backslash::

        >>> symbols('x((a:b))')
        (x(a), x(b))
        >>> symbols(r'x(:1\,:2)')  # or 'x((:1)\,(:2))'
        (x(0,0), x(0,1))

    All newly created symbols have assumptions set according to ``args``::

        >>> a = symbols('a', integer=True)
        >>> a.is_integer
        True

        >>> x, y, z = symbols('x y z', real=True)
        >>> x.is_real and y.is_real and z.is_real
        True

    Despite its name, :func:`symbols` can create symbol-like objects like
    instances of Function or Wild classes. To achieve this, set ``cls``
    keyword argument to the desired type::

        >>> symbols('f g h', cls=Function)
        (f, g, h)

        >>> type(_[0])
        <class 'diofant.core.function.UndefinedFunction'>

    """
    result = []

    if isinstance(names, str):
        marker = 0
        literals = [r'\,', r'\:', r'\ ']
        for i in range(len(literals)):
            lit = literals.pop(0)
            if lit in names:
                while chr(marker) in names:
                    marker += 1
                lit_char = chr(marker)
                marker += 1
                names = names.replace(lit, lit_char)
                literals.append((lit_char, lit[1:]))

        def literal(s):
            if literals:
                for c, l in literals:
                    s = s.replace(c, l)
            return s

        names = names.strip()
        as_seq = names.endswith(',')
        if as_seq:
            names = names[:-1].rstrip()
        if not names:
            raise ValueError('no symbols given')

        # split on commas
        names = [n.strip() for n in names.split(',')]
        if not all(n for n in names):
            raise ValueError('missing symbol between commas')
        # split on spaces
        for i in range(len(names) - 1, -1, -1):
            names[i: i + 1] = names[i].split()

        cls = args.pop('cls', Symbol)
        seq = args.pop('seq', as_seq)

        for name in names:
            if ':' not in name:
                symbol = cls(literal(name), **args)
                result.append(symbol)
                continue

            split = _range.split(name)
            # remove 1 layer of bounding parentheses around ranges
            for i in range(len(split) - 1):
                if i and ':' in split[i] and split[i] != ':' and \
                        split[i - 1].endswith('(') and \
                        split[i + 1].startswith(')'):
                    split[i - 1] = split[i - 1][:-1]
                    split[i + 1] = split[i + 1][1:]
            for i, s in enumerate(split):
                if ':' in s:
                    if s[-1].endswith(':'):
                        raise ValueError('missing end range')
                    a, b = s.split(':')
                    if b[-1] in string.digits:
                        a = 0 if not a else int(a)
                        b = int(b)
                        split[i] = [str(c) for c in range(a, b)]
                    else:
                        a = a or 'a'
                        split[i] = [string.ascii_letters[c] for c in range(
                            string.ascii_letters.index(a),
                            string.ascii_letters.index(b) + 1)]  # inclusive
                    if not split[i]:
                        break
                else:
                    split[i] = [s]
            else:
                seq = True
                names = [''.join(s) for s in itertools.product(*split)]
                if literals:
                    result.extend([cls(literal(s), **args) for s in names])
                else:
                    result.extend([cls(s, **args) for s in names])

        if not seq and len(result) <= 1:
            if not result:
                return ()
            return result[0]

        return tuple(result)
    else:
        for name in names:
            result.append(symbols(name, **args))

        return type(names)(result)


def var(names, **args):
    """
    Create symbols and inject them into the global namespace.

    This calls :func:`symbols` with the same arguments and puts the results
    into the *global* namespace. It's recommended not to use :func:`var` in
    library code, where :func:`symbols` has to be used.

    Examples
    ========

    >>> var('x')
    x
    >>> x
    x

    >>> var('a ab abc')
    (a, ab, abc)
    >>> abc
    abc

    >>> var('x y', real=True)
    (x, y)
    >>> x.is_real and y.is_real
    True

    See Also
    ========

    symbols

    """
    from ..utilities import flatten
    from ..utilities.magic import pollute

    ret = symbols(names, **args)
    syms = flatten([ret])
    pollute([str(_) for _ in syms], syms)

    return ret
