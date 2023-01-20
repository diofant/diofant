"""Implementation of :class:`Field` class."""

from .ring import CommutativeRing


class Field(CommutativeRing):
    """Represents a field domain."""

    is_Field = True

    @property
    def ring(self):
        raise AttributeError(f'there is no ring associated with {self}')

    @property
    def field(self):
        """Return a field associated with ``self``."""
        return self

    def exquo(self, a, b):
        """Exact quotient of ``a`` and ``b``, implies ``__truediv__``."""
        return a / b

    def quo(self, a, b):
        """Quotient of ``a`` and ``b``, implies ``__truediv__``."""
        return a / b

    def rem(self, a, b):
        """Remainder of ``a`` and ``b``, implies nothing."""
        return self.zero

    def div(self, a, b):
        """Division of ``a`` and ``b``, implies ``__truediv__``."""
        return a / b, self.zero

    def gcd(self, a, b):
        """
        Return GCD of ``a`` and ``b``.

        This definition of GCD over fields allows to clear denominators
        in `primitive()`.

        >>> QQ.gcd(QQ(2, 3), QQ(4, 9))
        2/9
        >>> gcd(Rational(2, 3), Rational(4, 9))
        2/9
        >>> primitive(2*x/3 + Rational(4, 9))
        (2/9, 3*x + 2)
        """
        try:
            ring = self.ring
        except (AttributeError, NotImplementedError):
            return self.one

        p = ring.gcd(a.numerator, b.numerator)
        q = ring.lcm(a.denominator, b.denominator)

        return self.convert(p, ring)/q

    def lcm(self, a, b):
        return (a*b)/self.gcd(a, b)
