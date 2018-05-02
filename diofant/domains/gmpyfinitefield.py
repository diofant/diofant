"""Implementation of :class:`GMPYFiniteField` class. """

from .finitefield import FiniteField
from .gmpyintegerring import GMPYIntegerRing


__all__ = ('GMPYFiniteField',)


class GMPYFiniteField(FiniteField):
    """Finite field based on GMPY's integers. """

    def __init__(self, mod, symmetric=True):
        return super().__init__(mod, GMPYIntegerRing(), symmetric)
