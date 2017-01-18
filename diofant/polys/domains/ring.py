"""Implementation of :class:`Ring` class. """

from diofant.polys.domains.domain import Domain
from diofant.polys.polyerrors import (ExactQuotientFailed, NotInvertible,
                                      NotReversible)
from diofant.utilities import public


@public
class Ring(Domain):
    """Represents a ring domain. """

    has_Ring = True

    def get_ring(self):
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

        if self.is_one(h):
            return s % b
        else:
            raise NotInvertible("zero divisor")

    def revert(self, a):
        """Returns ``a**(-1)`` if possible. """
        if self.is_one(a):
            return a
        else:
            raise NotReversible('only unity is reversible in a ring')

    def is_unit(self, a):
        try:
            self.revert(a)
            return True
        except NotReversible:
            return False

    def numer(self, a):
        """Returns numerator of ``a``. """
        return a

    def denom(self, a):
        """Returns denominator of `a`. """
        return self.one

    def free_module(self, rank):
        """
        Generate a free module of rank ``rank`` over self.

        >>> from diofant.abc import x
        >>> from diofant import QQ
        >>> QQ.poly_ring(x).free_module(2)
        QQ[x]**2
        """
        raise NotImplementedError

    def ideal(self, *gens):
        """
        Generate an ideal of ``self``.

        >>> from diofant.abc import x
        >>> from diofant import QQ
        >>> QQ.poly_ring(x).ideal(x**2)
        <x**2>
        """
        from diofant.polys.agca.ideals import ModuleImplementedIdeal
        return ModuleImplementedIdeal(self, self.free_module(1).submodule(
            *[[x] for x in gens]))

    def quotient_ring(self, e):
        """
        Form a quotient ring of ``self``.

        Here ``e`` can be an ideal or an iterable.

        >>> from diofant.abc import x
        >>> from diofant import QQ
        >>> QQ.poly_ring(x).quotient_ring(QQ.poly_ring(x).ideal(x**2))
        QQ[x]/<x**2>
        >>> QQ.poly_ring(x).quotient_ring([x**2])
        QQ[x]/<x**2>

        The division operator has been overloaded for this:

        >>> QQ.poly_ring(x)/[x**2]
        QQ[x]/<x**2>
        """
        from diofant.polys.agca.ideals import Ideal
        from diofant.polys.domains.quotientring import QuotientRing
        if not isinstance(e, Ideal):
            e = self.ideal(*e)
        return QuotientRing(self, e)

    def __div__(self, e):
        return self.quotient_ring(e)

    __truediv__ = __div__
