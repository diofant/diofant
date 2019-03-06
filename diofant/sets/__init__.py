"""
Package for set theory.
"""

from .sets import (Set, Interval, Union, EmptySet, FiniteSet, ProductSet,
                   Intersection, imageset, Complement, SymmetricDifference)
from .fancysets import ImageSet, Range
from .contains import Contains

from ..core.singleton import S
Naturals = S.Naturals
Naturals0 = S.Naturals0
Integers = S.Integers
Rationals = S.Rationals
Reals = S.Reals
del S
