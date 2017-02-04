"""Implementaton of :class:`CharacteristicZero` class. """

from .domain import Domain
from ...utilities import public


@public
class CharacteristicZero(Domain):
    """Domain that has infinite number of elements. """

    has_CharacteristicZero = True

    def characteristic(self):
        """Return the characteristic of this domain. """
        return 0
