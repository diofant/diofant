"""Implementaton of :class:`CharacteristicZero` class. """

from diofant.polys.domains.domain import Domain
from diofant.utilities import public


@public
class CharacteristicZero(Domain):
    """Domain that has infinite number of elements. """

    has_CharacteristicZero = True

    def characteristic(self):
        """Return the characteristic of this domain. """
        return 0
