"""Module for Diofant containers

(Diofant objects that store other Diofant objects)

The containers implemented in this module are subclassed to Basic.
They are supposed to work seamlessly within the Diofant framework.
"""

from collections.abc import Mapping

from .basic import Basic
from .compatibility import as_int, iterable
from .sympify import converter, sympify


class Tuple(Basic):
    """
    Wrapper around the builtin tuple object

    The Tuple is a subclass of Basic, so that it works well in the
    Diofant framework.  The wrapped tuple is available as self.args, but
    you can also access elements or slices with [:] syntax.

    >>> Tuple(a, b, c)[1:]
    (b, c)
    >>> Tuple(a, b, c).subs({a: d})
    (d, b, c)

    """

    def __new__(cls, *args):
        args = (sympify(arg) for arg in args)
        return super().__new__(cls, *args)

    def __getitem__(self, i):
        if isinstance(i, slice):
            indices = i.indices(len(self))
            return Tuple(*[self.args[j] for j in range(*indices)])
        return self.args[i]

    def __len__(self):
        return len(self.args)

    def __contains__(self, item):
        return item in self.args

    def __iter__(self):
        return iter(self.args)

    def __add__(self, other):
        if isinstance(other, Tuple):
            return Tuple(*(self.args + other.args))
        if isinstance(other, tuple):
            return Tuple(*(self.args + other))
        return NotImplemented

    def __radd__(self, other):
        if isinstance(other, Tuple):
            return Tuple(*(other.args + self.args))
        if isinstance(other, tuple):
            return Tuple(*(other + self.args))
        return NotImplemented

    def __mul__(self, other):
        try:
            n = as_int(other)
        except ValueError as exc:
            raise TypeError("Can't multiply sequence by non-integer "
                            f"of type '{type(other)!s}'") from exc
        return self.func(*(self.args*n))

    __rmul__ = __mul__

    def __eq__(self, other):
        if isinstance(other, Basic):
            return super().__eq__(other)
        return self.args == other

    def __hash__(self):
        return hash(self.args)

    def _to_mpmath(self, prec):
        return tuple(a._to_mpmath(prec) for a in self.args)

    def __lt__(self, other):
        return sympify(self.args < other.args)

    def __le__(self, other):
        return sympify(self.args <= other.args)

    # XXX: Basic defines count() as something different, so we can't
    # redefine it here. Originally this lead to cse() test failure.
    def tuple_count(self, value):
        """T.count(value) -> int -- return number of occurrences of value."""
        return self.args.count(value)

    def index(self, value, start=None, stop=None):
        """Return first index of value.

        Raises ValueError if the value is not present.

        """
        # XXX: One would expect:
        #
        # return self.args.index(value, start, stop)
        #
        # here. Any trouble with that? Yes:
        #
        # >>> [1].index(1, None, None)
        # Traceback (most recent call last):
        #   File "<stdin>", line 1, in <module>
        # TypeError: slice indices must be integers or None or have an __index__ method
        #
        # See: http://bugs.python.org/issue13340

        if start is None and stop is None:
            return self.args.index(value)
        if stop is None:
            return self.args.index(value, start)
        return self.args.index(value, start, stop)


converter[tuple] = lambda tup: Tuple(*tup)


def tuple_wrapper(method):
    """
    Decorator that converts any tuple in the function arguments into a Tuple.

    The motivation for this is to provide simple user interfaces.  The user can
    call a function with regular tuples in the argument, and the wrapper will
    convert them to Tuples before handing them to the function.

    >>> def f(*args):
    ...     return args
    >>> g = tuple_wrapper(f)

    The decorated function g sees only the Tuple argument:

    >>> g(0, (1, 2), 3)
    (0, (1, 2), 3)

    """
    def wrap_tuples(*args, **kw_args):
        newargs = []
        for arg in args:
            if type(arg) is tuple:
                newargs.append(Tuple(*arg))
            else:
                newargs.append(arg)
        return method(*newargs, **kw_args)
    return wrap_tuples


class Dict(Basic):
    """
    Wrapper around the builtin dict object

    The Dict is a subclass of Basic, so that it works well in the
    Diofant framework.  Because it is immutable, it may be included
    in sets, but its values must all be given at instantiation and
    cannot be changed afterwards.  Otherwise it behaves identically
    to the Python dict.

    >>> D = Dict({1: 'one', 2: 'two'})
    >>> for key in D:
    ...     if key == 1:
    ...         print(f'{key} {D[key]}')
    1 one

    The args are sympified so the 1 and 2 are Integers and the values
    are Symbols. Queries automatically sympify args so the following work:

    >>> 1 in D
    True
    >>> D.has('one')  # searches keys and values
    True
    >>> 'one' in D  # not in the keys
    False
    >>> D[1]
    one

    """

    def __new__(cls, *args):
        if len(args) == 1 and isinstance(args[0], (dict, Dict)):
            items = [Tuple(k, v) for k, v in args[0].items()]
        elif iterable(args) and all(len(arg) == 2 for arg in args):
            items = [Tuple(k, v) for k, v in args]
        else:
            raise TypeError('Pass Dict args as Dict((k1, v1), ...) or Dict({k1: v1, ...})')
        elements = frozenset(items)
        obj = Basic.__new__(cls, elements)
        obj.elements = elements
        obj._dict = dict(items)  # In case Tuple decides it wants to sympify
        return obj

    def __getitem__(self, key):
        """x.__getitem__(y) <==> x[y]."""
        return self._dict[sympify(key)]

    def __setitem__(self, key, value):
        raise NotImplementedError('Diofant Dicts are Immutable')

    @property
    def args(self):
        """Returns a tuple of arguments of 'self'.

        See Also
        ========

        diofant.core.basic.Basic.args

        """
        return tuple(self.elements)

    def items(self):
        """Returns a set-like object providing a view on Dict's items."""
        return self._dict.items()

    def keys(self):
        """Returns a set-like object providing a view on Dict's keys."""
        return self._dict.keys()

    def values(self):
        """Returns a set-like object providing a view on Dict's values."""
        return self._dict.values()

    def __iter__(self):
        """x.__iter__() <==> iter(x)."""
        return iter(self._dict)

    def __len__(self):
        """x.__len__() <==> len(x)."""
        return self._dict.__len__()

    def get(self, key, default=None):
        """Return the value for key if key is in the dictionary, else default."""
        return self._dict.get(sympify(key), default)

    def __contains__(self, key):
        """D.__contains__(k) -> True if D has a key k, else False."""
        return sympify(key) in self._dict

    def __lt__(self, other):
        return sympify(self.args < other.args)

    @property
    def _sorted_args(self):
        from ..utilities import default_sort_key
        return tuple(sorted(self.args, key=default_sort_key))


Mapping.register(Dict)
