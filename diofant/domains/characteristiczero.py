"""Implementaton of :class:`CharacteristicZero` class. """

from .domain import Domain


__all__ = 'CharacteristicZero',


class CharacteristicZero(Domain):
    """Domain that has infinite number of elements. """

    @property
    def characteristic(self):
        """Return the characteristic of this domain. """
        return 0
