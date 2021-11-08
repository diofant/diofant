"""Tools and arithmetics for monomials of distributed polynomials."""

import collections

from ..core import Integer, Mul
from ..printing.defaults import DefaultPrinting


def itermonomials(variables, degree):
    r"""
    Generate a set of monomials of the given total degree or less.

    Examples
    ========

    >>> set(itermonomials([x, y], 2))
    {1, x, x**2, y, y**2, x*y}

    """
    if not variables:
        yield Integer(1)
    else:
        x, tail = variables[0], variables[1:]

        for i in range(degree + 1):
            for m in itermonomials(tail, degree - i):
                yield m * x**i


class Monomial(tuple, DefaultPrinting):
    """Class representing a monomial, i.e. a product of powers."""

    def __new__(cls, monom, gens=None):
        """Create and return the Monomial object."""
        from .polytools import Poly

        if not isinstance(monom, collections.abc.Iterable):
            poly = Poly(monom, gens=gens)
            rep = poly.rep
            gens = poly.gens
            if rep and rep.is_monomial:
                monom = list(rep)[0]
            else:
                raise ValueError(f'Expected a monomial got {monom}')

        obj = super().__new__(cls, map(int, monom))
        obj.gens = gens

        return obj

    @classmethod
    def _get_val(cls, other):
        if not isinstance(other, cls):
            try:
                other = cls(other)
            except (ValueError, TypeError):
                return
        return other

    def as_expr(self, *gens):
        """Convert a monomial instance to a Diofant expression."""
        gens = gens or self.gens

        if not gens:
            raise ValueError(f"Can't convert {self} to an expression without generators")

        return Mul(*[gen**exp for gen, exp in zip(gens, self)])

    def __mul__(self, other):
        """Return self*other."""
        other = self._get_val(other)
        if other is not None:
            return self.__class__((a + b for a, b in zip(self, other)), self.gens)
        else:
            return NotImplemented

    def __truediv__(self, other):
        """Return self/other."""
        other = self._get_val(other)
        if other is not None:
            return self.__class__((a - b for a, b in zip(self, other)), self.gens)
        else:
            return NotImplemented

    def divides(self, other):
        """Check if self divides other."""
        other, orig = self._get_val(other), other
        if other is not None:
            return all(a <= b for a, b in zip(self, other))
        else:
            raise TypeError(f'An instance of {self.__class__.__name__} expected, got {orig}')

    def __pow__(self, n):
        """Return pow(self, other)."""
        if isinstance(n, int) and n >= 0:
            return self.__class__((a * n for a in self), self.gens)
        else:
            raise ValueError(f'A non-negative integer expected, got {n}')

    def gcd(self, other):
        """Greatest common divisor of monomials."""
        other, orig = self._get_val(other), other
        if other is not None:
            return self.__class__((min(a, b) for a, b in zip(self, other)), self.gens)
        else:
            raise TypeError(f'An instance of {self.__class__.__name__} expected, got {orig}')

    def lcm(self, other):
        """Least common multiple of monomials."""
        other, orig = self._get_val(other), other
        if other is not None:
            return self.__class__((max(a, b) for a, b in zip(self, other)), self.gens)
        else:
            raise TypeError(f'An instance of {self.__class__.__name__} expected, got {orig}')
