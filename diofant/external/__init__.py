"""
Unified place for determining if external dependencies are installed or not.

You should import all external modules using the import_module() function.

For example

>>> numpy = import_module('numpy')

If the resulting library is not installed, the function will return None.
Otherwise, it will return the library.
"""

from .gmpy import GROUND_TYPES, HAS_GMPY, gmpy
from .importtools import import_module


__all__ = 'GROUND_TYPES', 'HAS_GMPY', 'gmpy', 'import_module'
