"""
Package for set theory.
"""

from ..core.singleton import S
from .contains import Contains
from .fancysets import ImageSet, Range
from .sets import (Complement, EmptySet, FiniteSet, Intersection, Interval,
                   ProductSet, Set, SymmetricDifference, Union, imageset)


Naturals = S.Naturals
Naturals0 = S.Naturals0
Integers = S.Integers
Rationals = S.Rationals
Reals = S.Reals
del S
