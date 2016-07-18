"""Implementation of :class:`GMPYFiniteField` class. """

from diofant.polys.domains.finitefield import FiniteField
from diofant.polys.domains.gmpyintegerring import GMPYIntegerRing
from diofant.utilities import public


@public
class GMPYFiniteField(FiniteField):
    """Finite field based on GMPY integers. """

    alias = 'FF_gmpy'

    def __init__(self, mod, symmetric=True):
        return super(GMPYFiniteField, self).__init__(mod, GMPYIntegerRing(), symmetric)
