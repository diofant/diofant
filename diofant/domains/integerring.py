"""Implementation of :class:`IntegerRing` class. """

import math

from .characteristiczero import CharacteristicZero
from .ring import Ring
from .simpledomain import SimpleDomain


__all__ = ('IntegerRing',)


class IntegerRing(Ring, CharacteristicZero, SimpleDomain):
    """General class for integer rings. """

    rep = 'ZZ'

    is_IntegerRing = is_ZZ = True
    is_Numerical = True

    has_assoc_Ring = True
    has_assoc_Field = True

    @property
    def field(self):
        """Returns a field associated with ``self``. """
        from . import QQ
        return QQ

    def algebraic_field(self, *extension):
        r"""Returns an algebraic field, i.e. `\mathbb{Q}(\alpha, \ldots)`. """
        return self.field.algebraic_field(*extension)

    def from_AlgebraicField(self, a, K0):
        """Convert an algebraic number to ``dtype``. """
        if a.is_ground:
            return self.convert(a.LC(), K0.domain)

    def log(self, a, b):
        """Returns b-base logarithm of ``a``. """
        return self.dtype(math.log(int(a), b))
