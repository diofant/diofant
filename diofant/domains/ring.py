"""Implementation of :class:`Ring` class. """

from ..polys.polyerrors import ExactQuotientFailed, NotInvertible
from .domain import Domain


__all__ = 'Ring',


class Ring(Domain):
    """Represents a ring domain. """

    is_Ring = True

    @property
    def ring(self):
        """Returns a ring associated with ``self``. """
        return self

    def exquo(self, a, b):
        """Exact quotient of ``a`` and ``b``, implies ``__floordiv__``.  """
        if a % b:
            raise ExactQuotientFailed(a, b, self)
        else:
            return a // b

    def quo(self, a, b):
        """Quotient of ``a`` and ``b``, implies ``__floordiv__``. """
        return a // b

    def rem(self, a, b):
        """Remainder of ``a`` and ``b``, implies ``__mod__``.  """
        return a % b

    def div(self, a, b):
        """Division of ``a`` and ``b``, implies ``__divmod__``. """
        return divmod(a, b)

    def invert(self, a, b):
        """Returns inversion of ``a mod b``. """
        s, t, h = self.gcdex(a, b)

        if h == self.one:
            return s % b
        else:
            raise NotInvertible("zero divisor")

    def half_gcdex(self, a, b):
        """Half extended GCD of ``a`` and ``b``. """
        s, t, h = self.gcdex(a, b)
        return s, h

    def cofactors(self, a, b):
        """Returns GCD and cofactors of ``a`` and ``b``. """
        gcd = self.gcd(a, b)
        cfa = self.quo(a, gcd)
        cfb = self.quo(b, gcd)
        return gcd, cfa, cfb
