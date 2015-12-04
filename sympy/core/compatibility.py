"""
Reimplementations of constructs introduced in later versions of Python than
we support. Also some functions that are needed SymPy-wide and are located
here for easy import.
"""
from __future__ import print_function, division

import operator
import os
from collections import defaultdict
from sympy.external import import_module

"""
Python 2 and Python 3 compatible imports

String and Unicode compatible changes:
    * `unicode()` removed in Python 3, import `unicode` for Python 2/3
      compatible function
    * `unichr()` removed in Python 3, import `unichr` for Python 2/3 compatible
      function
    * Use `u()` for escaped unicode sequences (e.g. u'\u2020' -> u('\u2020'))
    * Use `u_decode()` to decode utf-8 formatted unicode strings
    * `string_types` gives str in Python 3, unicode and str in Python 2,
      equivalent to basestring

Integer related changes:
    * `long()` removed in Python 3, import `long` for Python 2/3 compatible
      function
    * `integer_types` gives int in Python 3, int and long in Python 2

Types related changes:
    * `class_types` gives type in Python 3, type and ClassType in Python 2

Renamed function attributes:
    * Python 2 `.func_code`, Python 3 `.__func__`, access with
      `get_function_code()`
    * Python 2 `.func_globals`, Python 3 `.__globals__`, access with
      `get_function_globals()`
    * Python 2 `.func_name`, Python 3 `.__name__`, access with
      `get_function_name()`

Moved modules:
    * `reduce()`
    * `StringIO()`
    * `cStringIO()` (same as `StingIO()` in Python 3)
    * Python 2 `__builtins__`, access with Python 3 name, `builtins`

Iterator/list changes:
    * `xrange` removed in Python 3, import `xrange` for Python 2/3 compatible
      iterator version of range

exec:
    * Use `exec_()`, with parameters `exec_(code, globs=None, locs=None)`

Metaclasses:
    * Use `with_metaclass()`, examples below
        * Define class `Foo` with metaclass `Meta`, and no parent:
            class Foo(with_metaclass(Meta)):
                pass
        * Define class `Foo` with metaclass `Meta` and parent class `Bar`:
            class Foo(with_metaclass(Meta, Bar)):
                pass
"""

import sys
PY3 = sys.version_info[0] > 2

if PY3:
    class_types = type,
    integer_types = (int,)
    string_types = (str,)
    long = int

    # String / unicode compatibility
    unicode = str
    unichr = chr

    def u(x):
        return x

    def u_decode(x):
        return x

    Iterator = object

    # Moved definitions
    get_function_code = operator.attrgetter("__code__")
    get_function_globals = operator.attrgetter("__globals__")
    get_function_name = operator.attrgetter("__name__")

    import builtins
    from functools import reduce
    from io import StringIO
    cStringIO = StringIO

    exec_=getattr(builtins, "exec")

    range=range

    from tokenize import tokenize
else:
    import codecs
    import types

    class_types = (type, types.ClassType)
    integer_types = (int, long)
    string_types = (str, unicode)
    long = long

    # String / unicode compatibility
    unicode = unicode
    unichr = unichr

    def u(x):
        return codecs.unicode_escape_decode(x)[0]

    def u_decode(x):
        return x.decode('utf-8')

    class Iterator(object):
        def next(self):
            return type(self).__next__(self)

    # Moved definitions
    get_function_code = operator.attrgetter("func_code")
    get_function_globals = operator.attrgetter("func_globals")
    get_function_name = operator.attrgetter("func_name")

    import __builtin__ as builtins
    reduce = reduce
    from StringIO import StringIO
    from cStringIO import StringIO as cStringIO

    def exec_(_code_, _globs_=None, _locs_=None):
        """Execute code in a namespace."""
        if _globs_ is None:
            frame = sys._getframe(1)
            _globs_ = frame.f_globals
            if _locs_ is None:
                _locs_ = frame.f_locals
            del frame
        elif _locs_ is None:
            _locs_ = _globs_
        exec("exec _code_ in _globs_, _locs_")
    range=xrange

    from tokenize import generate_tokens as tokenize


def with_metaclass(meta, *bases):
    """
    Create a base class with a metaclass.

    For example, if you have the metaclass

    >>> class Meta(type):
    ...     pass

    Use this as the metaclass by doing

    >>> from sympy.core.compatibility import with_metaclass
    >>> class MyClass(with_metaclass(Meta, object)):
    ...     pass

    This is equivalent to the Python 2::

        class MyClass(object):
            __metaclass__ = Meta

    or Python 3::

        class MyClass(object, metaclass=Meta):
            pass

    That is, the first argument is the metaclass, and the remaining arguments
    are the base classes. Note that if the base class is just ``object``, you
    may omit it.

    >>> MyClass.__mro__ == (MyClass, object)
    True
    >>> type(MyClass) is Meta
    True

    """
    # This requires a bit of explanation: the basic idea is to make a dummy
    # metaclass for one level of class instantiation that replaces itself with
    # the actual metaclass.
    # Code copied from the 'six' library.
    class metaclass(meta):
        def __new__(cls, name, this_bases, d):
            return meta(name, bases, d)
    return type.__new__(metaclass, "NewBase", (), {})


# These are in here because telling if something is an iterable just by calling
# hasattr(obj, "__iter__") behaves differently in Python 2 and Python 3.  In
# particular, hasattr(str, "__iter__") is False in Python 2 and True in Python 3.
# I think putting them here also makes it easier to use them in the core.

class NotIterable:
    """
    Use this as mixin when creating a class which is not supposed to return
    true when iterable() is called on its instances. I.e. avoid infinite loop
    when calling e.g. list() on the instance
    """
    pass


def iterable(i, exclude=(string_types, dict, NotIterable)):
    """
    Return a boolean indicating whether ``i`` is SymPy iterable.
    True also indicates that the iterator is finite, i.e. you e.g.
    call list(...) on the instance.

    When SymPy is working with iterables, it is almost always assuming
    that the iterable is not a string or a mapping, so those are excluded
    by default. If you want a pure Python definition, make exclude=None. To
    exclude multiple items, pass them as a tuple.

    See Also
    ========

    is_sequence

    Examples
    ========

    >>> from sympy.utilities.iterables import iterable
    >>> from sympy import Tuple
    >>> things = [[1], (1,), {1}, Tuple(1), (j for j in [1, 2]), {1:2}, '1', 1]
    >>> for i in things:
    ...     print('%s %s' % (iterable(i), type(i)))
    True <... 'list'>
    True <... 'tuple'>
    True <... 'set'>
    True <class 'sympy.core.containers.Tuple'>
    True <... 'generator'>
    False <... 'dict'>
    False <... 'str'>
    False <... 'int'>

    >>> iterable({}, exclude=None)
    True
    >>> iterable({}, exclude=str)
    True
    >>> iterable("no", exclude=str)
    False

    """
    try:
        iter(i)
    except TypeError:
        return False
    if exclude:
        return not isinstance(i, exclude)
    return True


def is_sequence(i, include=None):
    """
    Return a boolean indicating whether ``i`` is a sequence in the SymPy
    sense. If anything that fails the test below should be included as
    being a sequence for your application, set 'include' to that object's
    type; multiple types should be passed as a tuple of types.

    Note: although generators can generate a sequence, they often need special
    handling to make sure their elements are captured before the generator is
    exhausted, so these are not included by default in the definition of a
    sequence.

    See Also
    ========

    iterable

    Examples
    ========

    >>> from sympy.utilities.iterables import is_sequence
    >>> from types import GeneratorType
    >>> is_sequence([])
    True
    >>> is_sequence(set())
    False
    >>> is_sequence('abc')
    False
    >>> is_sequence('abc', include=str)
    True
    >>> generator = (c for c in 'abc')
    >>> is_sequence(generator)
    False
    >>> is_sequence(generator, include=(str, GeneratorType))
    True
    """
    return (hasattr(i, '__getitem__') and
            iterable(i) or
            bool(include) and
            isinstance(i, include))

try:
    from itertools import zip_longest
except ImportError:  # <= Python 2.7
    from itertools import izip_longest as zip_longest


def as_int(n):
    """
    Convert the argument to a builtin integer.

    The return value is guaranteed to be equal to the input. ValueError is
    raised if the input has a non-integral value.

    Examples
    ========

    >>> from sympy.core.compatibility import as_int
    >>> from sympy import sqrt
    >>> 3.0
    3.0
    >>> as_int(3.0)  # convert to int and test for equality
    3
    >>> int(sqrt(10))
    3
    >>> as_int(sqrt(10))
    Traceback (most recent call last):
    ...
    ValueError: ... is not an integer

    """
    try:
        result = int(n)
        if result != n:
            raise TypeError
    except TypeError:
        raise ValueError('%s is not an integer' % n)
    return result


def default_sort_key(item, order=None):
    """Return a key that can be used for sorting.

    The key has the structure:

    (class_key, (len(args), args), exponent.sort_key(), coefficient)

    This key is supplied by the sort_key routine of Basic objects when
    ``item`` is a Basic object or an object (other than a string) that
    sympifies to a Basic object. Otherwise, this function produces the
    key.

    The ``order`` argument is passed along to the sort_key routine and is
    used to determine how the terms *within* an expression are ordered.
    (See examples below) ``order`` options are: 'lex', 'grlex', 'grevlex',
    and reversed values of the same (e.g. 'rev-lex'). The default order
    value is None (which translates to 'lex').

    Examples
    ========

    >>> from sympy import S, I, default_sort_key, sin, cos, sqrt
    >>> from sympy.core.function import UndefinedFunction
    >>> from sympy.abc import x

    The following are equivalent ways of getting the key for an object:

    >>> x.sort_key() == default_sort_key(x)
    True

    Here are some examples of the key that is produced:

    >>> default_sort_key(UndefinedFunction('f'))
    ((0, 0, 'UndefinedFunction'), (1, ('f',)), ((1, 0, 'Number'),
        (0, ()), (), 1), 1)
    >>> default_sort_key('1')
    ((0, 0, 'str'), (1, ('1',)), ((1, 0, 'Number'), (0, ()), (), 1), 1)
    >>> default_sort_key(S.One)
    ((1, 0, 'Number'), (0, ()), (), 1)
    >>> default_sort_key(2)
    ((1, 0, 'Number'), (0, ()), (), 2)


    While sort_key is a method only defined for SymPy objects,
    default_sort_key will accept anything as an argument so it is
    more robust as a sorting key. For the following, using key=
    lambda i: i.sort_key() would fail because 2 doesn't have a sort_key
    method; that's why default_sort_key is used. Note, that it also
    handles sympification of non-string items likes ints:

    >>> a = [2, I, -I]
    >>> sorted(a, key=default_sort_key)
    [2, -I, I]

    The returned key can be used anywhere that a key can be specified for
    a function, e.g. sort, min, max, etc...:

    >>> a.sort(key=default_sort_key); a[0]
    2
    >>> min(a, key=default_sort_key)
    2

    Note
    ----

    The key returned is useful for getting items into a canonical order
    that will be the same across platforms. It is not directly useful for
    sorting lists of expressions:

    >>> a, b = x, 1/x

    Since ``a`` has only 1 term, its value of sort_key is unaffected by
    ``order``:

    >>> a.sort_key() == a.sort_key('rev-lex')
    True

    If ``a`` and ``b`` are combined then the key will differ because there
    are terms that can be ordered:

    >>> eq = a + b
    >>> eq.sort_key() == eq.sort_key('rev-lex')
    False
    >>> eq.as_ordered_terms()
    [x, 1/x]
    >>> eq.as_ordered_terms('rev-lex')
    [1/x, x]

    But since the keys for each of these terms are independent of ``order``'s
    value, they don't sort differently when they appear separately in a list:

    >>> sorted(eq.args, key=default_sort_key)
    [1/x, x]
    >>> sorted(eq.args, key=lambda i: default_sort_key(i, order='rev-lex'))
    [1/x, x]

    The order of terms obtained when using these keys is the order that would
    be obtained if those terms were *factors* in a product.

    Although it is useful for quickly putting expressions in canonical order,
    it does not sort expressions based on their complexity defined by the
    number of operations, power of variables and others:

    >>> sorted([sin(x)*cos(x), sin(x)], key=default_sort_key)
    [sin(x)*cos(x), sin(x)]
    >>> sorted([x, x**2, sqrt(x), x**3], key=default_sort_key)
    [sqrt(x), x, x**2, x**3]

    See Also
    ========

    ordered
    sympy.core.expr.Expr.as_ordered_factors
    sympy.core.expr.Expr.as_ordered_terms
    """

    from sympy.core import S, Basic
    from sympy.core.sympify import sympify, SympifyError
    from sympy.core.compatibility import iterable

    if isinstance(item, Basic):
        return item.sort_key(order=order)

    if iterable(item, exclude=string_types):
        if isinstance(item, dict):
            args = item.items()
            unordered = True
        elif isinstance(item, set):
            args = item
            unordered = True
        else:
            # e.g. tuple, list
            args = list(item)
            unordered = False

        args = [default_sort_key(arg, order=order) for arg in args]

        if unordered:
            # e.g. dict, set
            args = sorted(args)

        cls_index, args = 10, (len(args), tuple(args))
    else:
        if not isinstance(item, string_types):
            try:
                item = sympify(item)
            except SympifyError:
                # e.g. lambda x: x
                pass
            else:
                if isinstance(item, Basic):
                    # e.g int -> Integer
                    return default_sort_key(item)
                # e.g. UndefinedFunction

        # e.g. str
        cls_index, args = 0, (1, (str(item),))

    return (cls_index, 0, item.__class__.__name__
            ), args, S.One.sort_key(), S.One


def _nodes(e):
    """
    A helper for ordered() which returns the node count of ``e`` which
    for Basic objects is the number of Basic nodes in the expression tree
    but for other objects is 1 (unless the object is an iterable or dict
    for which the sum of nodes is returned).
    """
    from .basic import Basic

    if isinstance(e, Basic):
        return e.count(Basic)
    elif iterable(e):
        return 1 + sum(_nodes(ei) for ei in e)
    elif isinstance(e, dict):
        return 1 + sum(_nodes(k) + _nodes(v) for k, v in e.items())
    else:
        return 1


def ordered(seq, keys=None, default=True, warn=False):
    """Return an iterator of the seq where keys are used to break ties in
    a conservative fashion: if, after applying a key, there are no ties
    then no other keys will be computed.

    Two default keys will be applied if 1) keys are not provided or 2) the
    given keys don't resolve all ties (but only if `default` is True). The
    two keys are `_nodes` (which places smaller expressions before large) and
    `default_sort_key` which (if the `sort_key` for an object is defined
    properly) should resolve any ties.

    If ``warn`` is True then an error will be raised if there were no
    keys remaining to break ties. This can be used if it was expected that
    there should be no ties between items that are not identical.

    Examples
    ========

    >>> from sympy.utilities.iterables import ordered
    >>> from sympy import count_ops
    >>> from sympy.abc import x, y

    The count_ops is not sufficient to break ties in this list and the first
    two items appear in their original order (i.e. the sorting is stable):

    >>> list(ordered([y + 2, x + 2, x**2 + y + 3],
    ...    count_ops, default=False, warn=False))
    ...
    [y + 2, x + 2, x**2 + y + 3]

    The default_sort_key allows the tie to be broken:

    >>> list(ordered([y + 2, x + 2, x**2 + y + 3]))
    ...
    [x + 2, y + 2, x**2 + y + 3]

    Here, sequences are sorted by length, then sum:

    >>> seq, keys = [[[1, 2, 1], [0, 3, 1], [1, 1, 3], [2], [1]], [
    ...    lambda x: len(x),
    ...    lambda x: sum(x)]]
    ...
    >>> list(ordered(seq, keys, default=False, warn=False))
    [[1], [2], [1, 2, 1], [0, 3, 1], [1, 1, 3]]

    If ``warn`` is True, an error will be raised if there were not
    enough keys to break ties:

    >>> list(ordered(seq, keys, default=False, warn=True))
    Traceback (most recent call last):
    ...
    ValueError: not enough keys to break ties


    Notes
    =====

    The decorated sort is one of the fastest ways to sort a sequence for
    which special item comparison is desired: the sequence is decorated,
    sorted on the basis of the decoration (e.g. making all letters lower
    case) and then undecorated. If one wants to break ties for items that
    have the same decorated value, a second key can be used. But if the
    second key is expensive to compute then it is inefficient to decorate
    all items with both keys: only those items having identical first key
    values need to be decorated. This function applies keys successively
    only when needed to break ties. By yielding an iterator, use of the
    tie-breaker is delayed as long as possible.

    This function is best used in cases when use of the first key is
    expected to be a good hashing function; if there are no unique hashes
    from application of a key then that key should not have been used. The
    exception, however, is that even if there are many collisions, if the
    first group is small and one does not need to process all items in the
    list then time will not be wasted sorting what one was not interested
    in. For example, if one were looking for the minimum in a list and
    there were several criteria used to define the sort order, then this
    function would be good at returning that quickly if the first group
    of candidates is small relative to the number of items being processed.

    """
    d = defaultdict(list)
    if keys:
        if not isinstance(keys, (list, tuple)):
            keys = [keys]
        keys = list(keys)
        f = keys.pop(0)
        for a in seq:
            d[f(a)].append(a)
    else:
        if not default:
            raise ValueError('if default=False then keys must be provided')
        d[None].extend(seq)

    for k in sorted(d.keys()):
        if len(d[k]) > 1:
            if keys:
                d[k] = ordered(d[k], keys, default, warn)
            elif default:
                d[k] = ordered(d[k], (_nodes, default_sort_key,),
                               default=False, warn=warn)
            elif warn:
                from sympy.utilities.iterables import uniq
                u = list(uniq(d[k]))
                if len(u) > 1:
                    raise ValueError(
                        'not enough keys to break ties: %s' % u)
        for v in d[k]:
            yield v
        d.pop(k)

# If HAS_GMPY is 0, no supported version of gmpy is available. Otherwise,
# HAS_GMPY contains the major version number of gmpy; i.e. 1 for gmpy, and
# 2 for gmpy2.

# Versions of gmpy prior to 1.03 do not work correctly with int(largempz)
# For example, int(gmpy.mpz(2**256)) would raise OverflowError.
# See issue 4980.

# Minimum version of gmpy changed to 1.13 to allow a single code base to also
# work with gmpy2.

GROUND_TYPES = os.getenv('SYMPY_GROUND_TYPES', 'auto').lower()
HAS_GMPY = 0

if GROUND_TYPES != 'python':

    # Don't try to import gmpy2 if ground types is set to gmpy1. This is
    # primarily intended for testing.

    if GROUND_TYPES != 'gmpy1':
        gmpy = import_module('gmpy2', min_module_version='2.0.0',
            module_version_attr='version', module_version_attr_call_args=())
        if gmpy:
            HAS_GMPY = 2
    else:
        GROUND_TYPES = 'gmpy'

    if not HAS_GMPY:
        gmpy = import_module('gmpy', min_module_version='1.13',
            module_version_attr='version', module_version_attr_call_args=())
        if gmpy:
            HAS_GMPY = 1

if GROUND_TYPES == 'auto':
    if HAS_GMPY:
        GROUND_TYPES = 'gmpy'
    else:
        GROUND_TYPES = 'python'

if GROUND_TYPES == 'gmpy' and not HAS_GMPY:
    from warnings import warn
    warn("gmpy library is not installed, switching to 'python' ground types")
    GROUND_TYPES = 'python'

# SYMPY_INTS is a tuple containing the base types for valid integer types.
SYMPY_INTS = integer_types

if GROUND_TYPES == 'gmpy':
    SYMPY_INTS += (type(gmpy.mpz(0)),)

# check_output() is new in Python 2.7
try:
    from subprocess import check_output
except ImportError:
    # running on platform like App Engine, no subprocess at all
    pass
