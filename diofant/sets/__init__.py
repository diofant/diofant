"""
Package for set theory.
"""

from .sets import (Set, Interval, Union, EmptySet, FiniteSet, ProductSet,  # noqa: F401
                   Intersection, imageset, Complement, SymmetricDifference)
from .fancysets import ImageSet, Range  # noqa: F401
from .contains import Contains  # noqa: F401
