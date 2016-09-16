"""Implementation of :class:`RationalField` class. """

from diofant.polys.domains.field import Field
from diofant.polys.domains.simpledomain import SimpleDomain
from diofant.polys.domains.characteristiczero import CharacteristicZero
from diofant.utilities import public


@public
class RationalField(Field, CharacteristicZero, SimpleDomain):
    """General class for rational fields. """

    rep = 'QQ'

    is_RationalField = is_QQ = True
    is_Numerical = True

    has_assoc_Ring = True
    has_assoc_Field = True

    def algebraic_field(self, *extension):
        r"""Returns an algebraic field, i.e. `\mathbb{Q}(\alpha, \ldots)`. """
        from diofant.polys.domains import AlgebraicField
        return AlgebraicField(self, *extension)

    def from_AlgebraicField(self, a, K0):
        """Convert a ``ANP`` object to ``dtype``. """
        if a.is_ground:
            return self.convert(a.LC(), K0.domain)
