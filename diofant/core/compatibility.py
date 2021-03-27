"""
Reimplementations of constructs introduced in later versions of Python than
we support. Also some functions that are needed Diofant-wide and are located
here for easy import.
"""

import os
import typing

from ..external import import_module


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


def iterable(i, exclude=(str, dict, NotIterable)):
    """
    Return a boolean indicating whether ``i`` is Diofant iterable.
    True also indicates that the iterator is finite, i.e. you e.g.
    call list(...) on the instance.

    When Diofant is working with iterables, it is almost always assuming
    that the iterable is not a string or a mapping, so those are excluded
    by default. If you want a pure Python definition, make exclude=None. To
    exclude multiple items, pass them as a tuple.

    See Also
    ========

    is_sequence

    Examples
    ========

    >>> things = [[1], (1,), {1}, Tuple(1), (j for j in [1, 2]), {1: 2}, '1', 1]
    >>> for i in things:
    ...     print(f'{iterable(i)} {type(i)}')
    True <... 'list'>
    True <... 'tuple'>
    True <... 'set'>
    True <class 'diofant.core.containers.Tuple'>
    True <... 'generator'>
    False <... 'dict'>
    False <... 'str'>
    False <... 'int'>

    >>> iterable({}, exclude=None)
    True
    >>> iterable({}, exclude=str)
    True
    >>> iterable('no', exclude=str)
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
    Return a boolean indicating whether ``i`` is a sequence in the Diofant
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


def as_int(n):
    """
    Convert the argument to a builtin integer.

    The return value is guaranteed to be equal to the input. ValueError is
    raised if the input has a non-integral value.

    Examples
    ========

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
        raise ValueError(f'{n} is not an integer')
    return result


# If HAS_GMPY is 0, no supported version of gmpy is available. Otherwise,
# HAS_GMPY contains the major version number of gmpy.

GROUND_TYPES: str = os.getenv('DIOFANT_GROUND_TYPES', 'auto').lower()

gmpy: typing.Any = import_module('gmpy2', min_module_version='2.1.0',
                                 module_version_attr='version',
                                 module_version_attr_call_args=())
if gmpy:
    HAS_GMPY = 2
else:
    HAS_GMPY = 0

if GROUND_TYPES == 'auto':
    if HAS_GMPY:
        GROUND_TYPES = 'gmpy'
    else:
        GROUND_TYPES = 'python'

if GROUND_TYPES == 'gmpy' and not HAS_GMPY:
    from warnings import warn
    warn("gmpy library is not installed, switching to 'python' ground types")
    GROUND_TYPES = 'python'

if GROUND_TYPES == 'python':
    os.environ['MPMATH_NOGMPY'] = 'yes'
