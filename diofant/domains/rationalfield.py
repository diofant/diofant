"""Implementation of :class:`RationalField` class. """

from .characteristiczero import CharacteristicZero
from .field import Field
from .simpledomain import SimpleDomain


__all__ = ('RationalField',)


class RationalField(Field, CharacteristicZero, SimpleDomain):
    """General class for rational fields. """

    rep = 'QQ'

    is_RationalField = is_QQ = True
    is_Numerical = True

    has_assoc_Ring = True
    has_assoc_Field = True

    def algebraic_field(self, *extension):
        r"""Returns an algebraic field, i.e. `\mathbb{Q}(\alpha, \ldots)`. """
        from . import AlgebraicField
        return AlgebraicField(self, *extension)

    def from_AlgebraicField(self, a, K0):
        """Convert a ``ANP`` object to ``dtype``. """
        if a.is_ground:
            return self.convert(a.LC(), K0.domain)
