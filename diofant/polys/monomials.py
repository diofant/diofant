"""Tools and arithmetics for monomials of distributed polynomials. """

from ..core import Mul, S, Tuple, sympify
from ..core.compatibility import iterable
from .polyerrors import ExactQuotientFailed
from .polyutils import dict_from_expr


__all__ = ('Monomial', 'itermonomials')


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

    Consider monomials in variables `x` and `y`::

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
            monoms |= { x**i * m for m in itermonomials(tail, degree - i) }

        return monoms


def monomial_count(V, N):
    r"""
    Computes the number of monomials.

    The number of monomials is given by the following formula:

    .. math::

        \frac{(\#V + N)!}{\#V! N!}

    where `N` is a total degree and `V` is a set of variables.

    Examples
    ========

    >>> from diofant.polys.orderings import monomial_key

    >>> monomial_count(2, 2)
    6

    >>> M = itermonomials([x, y], 2)

    >>> sorted(M, key=monomial_key('grlex', [y, x]))
    [1, x, y, x**2, x*y, y**2]
    >>> len(M)
    6

    """
    from ..functions import factorial
    return factorial(V + N) / factorial(V) / factorial(N)


def monomial_mul(A, B):
    """
    Multiplication of tuples representing monomials.

    Lets multiply `x**3*y**4*z` with `x*y**2`::

        >>> monomial_mul((3, 4, 1), (1, 2, 0))
        (4, 6, 1)

    which gives `x**4*y**5*z`.

    """
    return tuple(a + b for a, b in zip(A, B))


def monomial_div(A, B):
    """
    Division of tuples representing monomials.

    Lets divide `x**3*y**4*z` by `x*y**2`::

        >>> monomial_div((3, 4, 1), (1, 2, 0))
        (2, 2, 1)

    which gives `x**2*y**2*z`. However::

        >>> monomial_div((3, 4, 1), (1, 2, 2)) is None
        True

    `x*y**2*z**2` does not divide `x**3*y**4*z`.

    """
    C = monomial_ldiv(A, B)

    if all(c >= 0 for c in C):
        return tuple(C)
    else:
        return


def monomial_ldiv(A, B):
    """
    Division of tuples representing monomials.

    Lets divide `x**3*y**4*z` by `x*y**2`::

        >>> monomial_ldiv((3, 4, 1), (1, 2, 0))
        (2, 2, 1)

    which gives `x**2*y**2*z`.

        >>> monomial_ldiv((3, 4, 1), (1, 2, 2))
        (2, 2, -1)

    which gives `x**2*y**2*z**-1`.

    """
    return tuple(a - b for a, b in zip(A, B))


def monomial_pow(A, n):
    """Return the n-th pow of the monomial. """
    return tuple(a*n for a in A)


def monomial_gcd(A, B):
    """
    Greatest common divisor of tuples representing monomials.

    Lets compute GCD of `x*y**4*z` and `x**3*y**2`::

        >>> monomial_gcd((1, 4, 1), (3, 2, 0))
        (1, 2, 0)

    which gives `x*y**2`.

    """
    return tuple(min(a, b) for a, b in zip(A, B))


def monomial_lcm(A, B):
    """
    Least common multiple of tuples representing monomials.

    Lets compute LCM of `x*y**4*z` and `x**3*y**2`::

        >>> monomial_lcm((1, 4, 1), (3, 2, 0))
        (3, 4, 1)

    which gives `x**3*y**4*z`.

    """
    return tuple(max(a, b) for a, b in zip(A, B))


def monomial_divides(A, B):
    """
    Does there exist a monomial X such that XA == B?

    >>> monomial_divides((1, 2), (3, 4))
    True
    >>> monomial_divides((1, 2), (0, 2))
    False
    """
    return all(a <= b for a, b in zip(A, B))


def monomial_max(*monoms):
    """
    Returns maximal degree for each variable in a set of monomials.

    Consider monomials `x**3*y**4*z**5`, `y**5*z` and `x**6*y**3*z**9`.
    We wish to find out what is the maximal degree for each of `x`, `y`
    and `z` variables::

        >>> monomial_max((3, 4, 5), (0, 5, 1), (6, 3, 9))
        (6, 5, 9)

    """
    M = list(monoms[0])

    for N in monoms[1:]:
        for i, n in enumerate(N):
            M[i] = max(M[i], n)

    return tuple(M)


def monomial_min(*monoms):
    """
    Returns minimal degree for each variable in a set of monomials.

    Consider monomials `x**3*y**4*z**5`, `y**5*z` and `x**6*y**3*z**9`.
    We wish to find out what is the minimal degree for each of `x`, `y`
    and `z` variables::

        >>> monomial_min((3, 4, 5), (0, 5, 1), (6, 3, 9))
        (0, 3, 1)

    """
    M = list(monoms[0])

    for N in monoms[1:]:
        for i, n in enumerate(N):
            M[i] = min(M[i], n)

    return tuple(M)


def monomial_deg(M):
    """
    Returns the total degree of a monomial.

    For example, the total degree of `xy^2` is 3:

    >>> monomial_deg((1, 2))
    3
    """
    return sum(M)


def term_div(a, b, domain):
    """Division of two terms in over a ring/field. """
    a_lm, a_lc = a
    b_lm, b_lc = b

    monom = monomial_div(a_lm, b_lm)

    if domain.is_Field:
        if monom is not None:
            return monom, domain.quo(a_lc, b_lc)
    else:
        if not (monom is None or a_lc % b_lc):
            return monom, domain.quo(a_lc, b_lc)


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
            return "*".join([ "%s**%s" % (gen, exp) for gen, exp in zip(self.gens, self.exponents) ])
        else:
            return "%s(%s)" % (self.__class__.__name__, self.exponents)

    def as_expr(self, *gens):
        """Convert a monomial instance to a Diofant expression. """
        gens = gens or self.gens

        if not gens:
            raise ValueError(
                "can't convert %s to an expression without generators" % self)

        return Mul(*[ gen**exp for gen, exp in zip(gens, self.exponents) ])

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

        if result is not None:
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
