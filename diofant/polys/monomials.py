"""Tools and arithmetics for monomials of distributed polynomials."""

from collections.abc import Iterable
from math import prod
from operator import add, sub

from ..core import sympify
from ..printing.defaults import DefaultPrinting


class Monomial(tuple, DefaultPrinting):
    """Class representing a monomial, i.e. a product of powers."""

    def __new__(cls, monom, gens=()):
        """Create and return the Monomial object."""
        if not isinstance(monom, Iterable):
            monom = sympify(monom)
            poly = monom.as_poly(*gens)
            rep = poly.rep
            gens = poly.gens
            if rep and rep.is_monomial:
                monom, *_ = rep
            else:
                raise ValueError(f'Expected a monomial got {monom}')

        obj = super().__new__(cls, monom)
        obj.gens = gens

        return obj

    def as_expr(self, *gens):
        """Convert a monomial instance to a Diofant expression."""
        if gens := gens or self.gens:
            return prod(map(pow, gens, self))
        raise ValueError(f"Can't convert {self} to an expression without generators")

    def __mul__(self, other):
        """Return self*other."""
        return self.__class__(map(add, self, other), self.gens)

    def __truediv__(self, other):
        """Return self/other."""
        return self.__class__(map(sub, self, other), self.gens)

    def divides(self, other):
        """Check if self divides other."""
        return all(a <= b for a, b in zip(self, other))

    def __pow__(self, n):
        """Return pow(self, other)."""
        if isinstance(n, int) and n >= 0:
            return self.__class__((a * n for a in self), self.gens)
        raise ValueError(f'A non-negative integer expected, got {n}')

    def gcd(self, other):
        """Greatest common divisor of monomials."""
        return self.__class__(map(min, self, other), self.gens)

    def lcm(self, other):
        """Least common multiple of monomials."""
        return self.__class__(map(max, self, other), self.gens)
