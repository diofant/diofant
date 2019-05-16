"""Tools and arithmetics for monomials of distributed polynomials. """

from ..core import Mul, S, Tuple, sympify
from ..core.compatibility import iterable
from .polyutils import dict_from_expr


__all__ = 'Monomial', 'itermonomials'


def itermonomials(variables, degree):
    r"""
    Generate a set of monomials of the given total degree or less.

    Examples
    ========

    >>> set(itermonomials([x, y], 2))
    {1, x, x**2, y, y**2, x*y}

    """
    if not variables:
        yield S.One
    else:
        x, tail = variables[0], variables[1:]

        for i in range(degree + 1):
            for m in itermonomials(tail, degree - i):
                yield m*x**i


def monomial_mul(A, B):
    """Multiplication of tuples representing monomials."""
    return tuple(a + b for a, b in zip(A, B))


def monomial_div(A, B):
    """Division of tuples representing monomials."""
    return tuple(a - b for a, b in zip(A, B))


def monomial_pow(A, n):
    """Return the n-th pow of the monomial. """
    return tuple(a*n for a in A)


def monomial_gcd(A, B):
    """Greatest common divisor of tuples representing monomials."""
    return tuple(min(a, b) for a, b in zip(A, B))


def monomial_lcm(A, B):
    """Least common multiple of tuples representing monomials."""
    return tuple(max(a, b) for a, b in zip(A, B))


def monomial_divides(A, B):
    """Does there exist a monomial X such that XA == B?"""
    return all(a <= b for a, b in zip(A, B))


class Monomial:
    """Class representing a monomial, i.e. a product of powers. """

    def __init__(self, monom, gens=None):
        if not iterable(monom):
            rep, gens = dict_from_expr(sympify(monom), gens=gens)
            if len(rep) == 1 and list(rep.values())[0] == 1:
                monom = list(rep)[0]
            else:
                raise ValueError("Expected a monomial got %s" % monom)

        self._exponents = tuple(map(int, monom))
        self.gens = gens

    @property
    def exponents(self):
        return self._exponents

    def __len__(self):
        return len(self._exponents)

    def __iter__(self):
        return iter(self._exponents)

    def __getitem__(self, item):
        return self._exponents[item]

    def __hash__(self):
        return hash((self.__class__.__name__, self._exponents, self.gens))

    def __str__(self):
        if self.gens:
            return "*".join(["%s**%s" % _ for _ in zip(self.gens, self)])
        else:
            return "%s(%s)" % (self.__class__.__name__, self._exponents)

    def as_expr(self, *gens):
        """Convert a monomial instance to a Diofant expression. """
        gens = gens or self.gens

        if not gens:
            raise ValueError("can't convert %s to an expression without generators" % self)

        return Mul(*[gen**exp for gen, exp in zip(gens, self)])

    def __eq__(self, other):
        if isinstance(other, Monomial):
            exponents = tuple(other)
        elif isinstance(other, (tuple, Tuple)):
            exponents = other
        else:
            return False

        return self._exponents == exponents

    def __mul__(self, other):
        if isinstance(other, Monomial):
            exponents = tuple(other)
        elif isinstance(other, (tuple, Tuple)):
            exponents = other
        else:
            return NotImplemented

        return self.__class__(monomial_mul(self, exponents), self.gens)

    def __truediv__(self, other):
        if isinstance(other, Monomial):
            exponents = tuple(other)
        elif isinstance(other, (tuple, Tuple)):
            exponents = other
        else:
            return NotImplemented

        return self.__class__(monomial_div(self, exponents), self.gens)

    def __pow__(self, other):
        n = int(other)

        if n >= 0:
            return self.__class__(monomial_pow(self, n), self.gens)
        else:
            raise ValueError("a non-negative integer expected, got %s" % other)

    def gcd(self, other):
        """Greatest common divisor of monomials. """
        if isinstance(other, Monomial):
            exponents = tuple(other)
        elif isinstance(other, (tuple, Tuple)):
            exponents = other
        else:
            raise TypeError("an instance of Monomial class expected, got %s" % other)

        return self.__class__(monomial_gcd(self, exponents), self.gens)

    def lcm(self, other):
        """Least common multiple of monomials. """
        if isinstance(other, Monomial):
            exponents = tuple(other)
        elif isinstance(other, (tuple, Tuple)):
            exponents = other
        else:
            raise TypeError("an instance of Monomial class expected, got %s" % other)

        return self.__class__(monomial_lcm(self, exponents), self.gens)
