"""Implementaton of :class:`CharacteristicZero` class."""

from .domain import Domain


class CharacteristicZero(Domain):
    """Domain that has infinite number of elements."""

    @property
    def characteristic(self):
        return 0
