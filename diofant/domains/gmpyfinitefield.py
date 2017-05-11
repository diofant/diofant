"""Implementation of :class:`GMPYFiniteField` class. """

from .finitefield import FiniteField
from .gmpyintegerring import GMPYIntegerRing


__all__ = ('GMPYFiniteField',)


class GMPYFiniteField(FiniteField):
    """Finite field based on GMPY integers. """

    alias = 'FF_gmpy'

    def __init__(self, mod, symmetric=True):
        return super(GMPYFiniteField, self).__init__(mod, GMPYIntegerRing(), symmetric)
