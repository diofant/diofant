"""Implementation of :class:`CommutativeRing` class."""

import abc

from ..polys.polyerrors import ExactQuotientFailed, NotInvertible
from .domain import Domain


class CommutativeRing(Domain):
    """Represents a ring domain."""

    is_Ring = True

    @property
    def ring(self):
        """Returns a ring associated with ``self``."""
        return self

    def exquo(self, a, b):
        """Exact quotient of ``a`` and ``b``, implies ``__floordiv__``."""
        if a % b:
            raise ExactQuotientFailed(a, b, self)
        else:
            return a // b

    def quo(self, a, b):
        """Quotient of ``a`` and ``b``, implies ``__floordiv__``."""
        return a // b

    def rem(self, a, b):
        """Remainder of ``a`` and ``b``, implies ``__mod__``."""
        return a % b

    def div(self, a, b):
        """Division of ``a`` and ``b``, implies ``__divmod__``."""
        return divmod(a, b)

    def invert(self, a, b):
        """Returns inversion of ``a mod b``."""
        s, h = self.half_gcdex(a, b)

        if h == 1:
            return s % b
        else:
            raise NotInvertible('zero divisor')

    def half_gcdex(self, a, b):
        """Half extended GCD of ``a`` and ``b``."""
        s, t, h = self.gcdex(a, b)
        return s, h

    def cofactors(self, a, b):
        """Returns GCD and cofactors of ``a`` and ``b``."""
        gcd, cfa, cfb = self.gcd(a, b), self.zero, self.zero
        if gcd:
            cfa = self.quo(a, gcd)
            cfb = self.quo(b, gcd)
        return gcd, cfa, cfb

    def lcm(self, a, b):
        """Returns LCM of ``a`` and ``b``."""
        return abs(a*b)//self.gcd(a, b)

    @property
    @abc.abstractmethod
    def characteristic(self):
        """Return the characteristic of this ring."""
        raise NotImplementedError

    def is_normal(self, a):
        """Returns True if ``a`` is unit normal."""
        return a >= 0
