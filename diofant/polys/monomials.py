"""Tools and arithmetics for monomials of distributed polynomials. """

from ..core import Mul, S, Tuple, sympify
from ..core.compatibility import iterable
from .polyerrors import ExactQuotientFailed
from .polyutils import dict_from_expr


__all__ = 'Monomial', 'itermonomials'


def itermonomials(variables, degree):
    r"""
    Generate a set of monomials of the given total degree or less.

    Given a set of variables `V` and a total degree `N` generate
    a set of monomials of degree at most `N`. The total number of
    monomials is huge and is given by the following formula:

    .. math::

        \frac{(\#V + N)!}{\#V! N!}

    For example if we would like to generate a dense polynomial of
    a total degree `N = 50` in 5 variables, assuming that exponents
    and all of coefficients are 32-bit long and stored in an array we
    would need almost 80 GiB of memory! Fortunately most polynomials,
    that we will encounter, are sparse.

    Examples
    ========

    >>> from diofant.polys.orderings import monomial_key

    >>> sorted(itermonomials([x, y], 2), key=monomial_key('grlex', [y, x]))
    [1, x, y, x**2, x*y, y**2]

    >>> sorted(itermonomials([x, y], 3), key=monomial_key('grlex', [y, x]))
    [1, x, y, x**2, x*y, y**2, x**3, x**2*y, x*y**2, y**3]

    """
    if not variables:
        return {S.One}
    else:
        x, tail = variables[0], variables[1:]

        monoms = itermonomials(tail, degree)

        for i in range(1, degree + 1):
            monoms |= {x**i * m for m in itermonomials(tail, degree - i)}

        return monoms


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

        self.exponents = tuple(map(int, monom))
        self.gens = gens

    def rebuild(self, exponents, gens=None):
        return self.__class__(exponents, gens or self.gens)

    def __len__(self):
        return len(self.exponents)

    def __iter__(self):
        return iter(self.exponents)

    def __getitem__(self, item):
        return self.exponents[item]

    def __hash__(self):
        return hash((self.__class__.__name__, self.exponents, self.gens))

    def __str__(self):
        if self.gens:
            return "*".join(["%s**%s" % (gen, exp) for gen, exp in zip(self.gens, self.exponents)])
        else:
            return "%s(%s)" % (self.__class__.__name__, self.exponents)

    def as_expr(self, *gens):
        """Convert a monomial instance to a Diofant expression. """
        gens = gens or self.gens

        if not gens:
            raise ValueError(
                "can't convert %s to an expression without generators" % self)

        return Mul(*[gen**exp for gen, exp in zip(gens, self.exponents)])

    def __eq__(self, other):
        if isinstance(other, Monomial):
            exponents = other.exponents
        elif isinstance(other, (tuple, Tuple)):
            exponents = other
        else:
            return False

        return self.exponents == exponents

    def __mul__(self, other):
        if isinstance(other, Monomial):
            exponents = other.exponents
        elif isinstance(other, (tuple, Tuple)):
            exponents = other
        else:
            return NotImplemented

        return self.rebuild(monomial_mul(self.exponents, exponents))

    def __truediv__(self, other):
        if isinstance(other, Monomial):
            exponents = other.exponents
        elif isinstance(other, (tuple, Tuple)):
            exponents = other
        else:
            return NotImplemented

        result = monomial_div(self.exponents, exponents)

        if all(_ >= 0 for _ in result):
            return self.rebuild(result)
        else:
            raise ExactQuotientFailed(self, Monomial(other))

    __floordiv__ = __truediv__

    def __pow__(self, other):
        n = int(other)

        if not n:
            return self.rebuild([0]*len(self))
        elif n > 0:
            exponents = self.exponents

            for i in range(1, n):
                exponents = monomial_mul(exponents, self.exponents)

            return self.rebuild(exponents)
        else:
            raise ValueError("a non-negative integer expected, got %s" % other)

    def gcd(self, other):
        """Greatest common divisor of monomials. """
        if isinstance(other, Monomial):
            exponents = other.exponents
        elif isinstance(other, (tuple, Tuple)):
            exponents = other
        else:
            raise TypeError(
                "an instance of Monomial class expected, got %s" % other)

        return self.rebuild(monomial_gcd(self.exponents, exponents))

    def lcm(self, other):
        """Least common multiple of monomials. """
        if isinstance(other, Monomial):
            exponents = other.exponents
        elif isinstance(other, (tuple, Tuple)):
            exponents = other
        else:
            raise TypeError(
                "an instance of Monomial class expected, got %s" % other)

        return self.rebuild(monomial_lcm(self.exponents, exponents))
