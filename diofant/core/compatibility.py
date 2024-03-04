"""
Reimplementations of constructs introduced in later versions of Python than
we support. Also some functions that are needed Diofant-wide and are located
here for easy import.
"""

import os
import typing

from ..external import import_module


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
    except TypeError as exc:
        raise ValueError(f'{n} is not an integer') from exc
    return result


# If HAS_GMPY is 0, no supported version of gmpy is available. Otherwise,
# HAS_GMPY contains the major version number of gmpy.

GROUND_TYPES: str = os.getenv('DIOFANT_GROUND_TYPES', 'auto').lower()

gmpy: typing.Any = import_module('gmpy2')
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
