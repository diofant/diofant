import os
import typing

from .importtools import import_module


__all__ = ['GROUND_TYPES', 'HAS_GMPY', 'gmpy']


# If HAS_GMPY is 0, no supported version of gmpy is available. Otherwise,
# HAS_GMPY contains the major version number of gmpy.

GROUND_TYPES: str = os.getenv('DIOFANT_GROUND_TYPES', 'auto').lower()

gmpy: typing.Any = import_module('gmpy2')
if gmpy:
    HAS_GMPY = 2
else:
    gmpy = import_module('gmp')
    HAS_GMPY = 1 if gmpy else 0

if GROUND_TYPES == 'auto':
    if HAS_GMPY == 2:
        GROUND_TYPES = 'gmpy'
    elif HAS_GMPY == 1:
        GROUND_TYPES = 'gmp'
    else:
        GROUND_TYPES = 'python'

if GROUND_TYPES == 'gmpy' and not HAS_GMPY:
    from warnings import warn
    warn("gmpy library is not installed, switching to 'python' ground types")
    GROUND_TYPES = 'python'

if GROUND_TYPES == 'python':
    os.environ['MPMATH_NOGMPY'] = 'yes'
