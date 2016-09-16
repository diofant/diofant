"""Implementation of :class:`IntegerRing` class. """

import math

from diofant.polys.domains.ring import Ring
from diofant.polys.domains.simpledomain import SimpleDomain
from diofant.polys.domains.characteristiczero import CharacteristicZero
from diofant.utilities import public


@public
class IntegerRing(Ring, CharacteristicZero, SimpleDomain):
    """General class for integer rings. """

    rep = 'ZZ'

    is_IntegerRing = is_ZZ = True
    is_Numerical = True

    has_assoc_Ring = True
    has_assoc_Field = True

    def get_field(self):
        """Returns a field associated with ``self``. """
        from diofant.polys.domains import QQ
        return QQ

    def algebraic_field(self, *extension):
        r"""Returns an algebraic field, i.e. `\mathbb{Q}(\alpha, \ldots)`. """
        return self.get_field().algebraic_field(*extension)

    def from_AlgebraicField(self, a, K0):
        """Convert a ``ANP`` object to ``dtype``. """
        if a.is_ground:
            return self.convert(a.LC(), K0.domain)

    def log(self, a, b):
        """Returns b-base logarithm of ``a``. """
        return self.dtype(math.log(int(a), b))
