"""Implementaton of :class:`GMPYRationalField` class. """

from .groundtypes import GMPYRational, gmpy_factorial, gmpy_qdiv
from .rationalfield import RationalField


__all__ = ('GMPYRationalField',)


class GMPYRationalField(RationalField):
    """Rational field based on GMPY's rationals. """

    dtype = GMPYRational
    zero = dtype(0)
    one = dtype(1)

    def __init__(self):
        pass

    @property
    def ring(self):
        """Returns ring associated with ``self``. """
        from . import GMPYIntegerRing
        return GMPYIntegerRing()

    def exquo(self, a, b):
        """Exact quotient of `a` and `b`, implies `__truediv__`.  """
        return self.dtype(gmpy_qdiv(a, b))

    def quo(self, a, b):
        """Quotient of `a` and `b`, implies `__truediv__`. """
        return self.dtype(gmpy_qdiv(a, b))

    def factorial(self, a):
        """Returns factorial of `a`. """
        return self.dtype(gmpy_factorial(int(a)))
