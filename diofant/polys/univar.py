from .polyconfig import query
from .polyerrors import CoercionFailed, DomainError
from .rings import PolyElement, PolynomialRing
from .rootisolation import _FindRoot


__all__ = 'UnivarPolynomialRing', 'UnivarPolyElement'


class UnivarPolynomialRing(PolynomialRing, _FindRoot):
    """A class for representing univariate polynomial rings."""


class UnivarPolyElement(PolyElement):
    """Element of univariate distributed polynomial ring."""

    def shift(self, a):
        return self.compose(0, self.ring.gens[0] + a)

    def half_gcdex(self, other):
        """
        Half extended Euclidean algorithm in `F[x]`.

        Returns ``(s, h)`` such that ``h = gcd(self, other)``
        and ``s*self = h (mod other)``.

        Examples
        ========

        >>> R, x = ring('x', QQ)

        >>> f = x**4 - 2*x**3 - 6*x**2 + 12*x + 15
        >>> g = x**3 + x**2 - 4*x - 4

        >>> f.half_gcdex(g)
        (-1/5*x + 3/5, x + 1)

        """
        ring = self.ring
        domain = ring.domain

        if not domain.is_Field:
            raise DomainError(f"can't compute half extended GCD over {domain}")

        a, b = ring.one, ring.zero
        f, g = self, other

        while g:
            q, r = divmod(f, g)
            f, g = g, r
            a, b = b, a - q*b

        a = a.quo_ground(f.LC)
        f = f.monic()

        return a, f

    @property
    def is_cyclotomic(self):
        return self.ring._cyclotomic_p(self)

    def _right_decompose(self, s):
        ring = self.ring
        domain = ring.domain
        x = ring.gens[0]

        n = self.degree()
        lc = self.LC

        f = self.copy()
        g = x**s

        r = n // s

        for i in range(1, s):
            coeff = domain.zero

            for j in range(i):
                if not (n + j - i,) in f:
                    continue

                assert (s - j,) in g

                fc, gc = f[(n + j - i,)], g[(s - j,)]
                coeff += (i - r*j)*fc*gc

            g[(s - i,)] = domain.quo(coeff, i*r*lc)

        g._strip_zero()

        return g

    def _left_decompose(self, h):
        ring = self.ring

        g, i = ring.zero, 0
        f = self.copy()

        while f:
            q, r = divmod(f, h)

            if r.degree() > 0:
                return
            else:
                g[(i,)] = r.LC
                f, i = q, i + 1

        g._strip_zero()

        return g

    def _decompose(self):
        df = self.degree()

        for s in range(2, df):
            if df % s != 0:
                continue

            h = self._right_decompose(s)
            g = self._left_decompose(h)

            if g is not None:
                return g, h

    def decompose(self):
        """
        Compute functional decomposition of ``f`` in ``K[x]``.

        Given a univariate polynomial ``f`` with coefficients in a field of
        characteristic zero, returns list ``[f_1, f_2, ..., f_n]``, where::

                  f = f_1 o f_2 o ... f_n = f_1(f_2(... f_n))

        and ``f_2, ..., f_n`` are monic and homogeneous polynomials of at
        least second degree.

        Unlike factorization, complete functional decompositions of
        polynomials are not unique, consider examples:

        1. ``f o g = f(x + b) o (g - b)``
        2. ``x**n o x**m = x**m o x**n``
        3. ``T_n o T_m = T_m o T_n``

        where ``T_n`` and ``T_m`` are Chebyshev polynomials.

        Examples
        ========

        >>> R, x = ring('x', ZZ)

        >>> (x**4 - 2*x**3 + x**2).decompose()
        [x**2, x**2 - x]

        References
        ==========

        * :cite:`Kozen1989decomposition`

        """
        F = []
        f = self.copy()

        while True:
            result = f._decompose()

            if result is not None:
                f, h = result
                F = [h] + F
            else:
                break

        return [f] + F

    def sturm(self):
        """
        Computes the Sturm sequence of ``f`` in ``F[x]``.

        Given a univariate, square-free polynomial ``f(x)`` returns the
        associated Sturm sequence (see e.g. :cite:`Davenport1988systems`)
        ``f_0(x), ..., f_n(x)`` defined by::

           f_0(x), f_1(x) = f(x), f'(x)
           f_n = -rem(f_{n-2}(x), f_{n-1}(x))

        Examples
        ========

        >>> R, x = ring('x', QQ)

        >>> (x**3 - 2*x**2 + x - 3).sturm()
        [x**3 - 2*x**2 + x - 3, 3*x**2 - 4*x + 1, 2/9*x + 25/9, -2079/4]

        """
        return self.ring._sturm(self)

    def __mul__(self, other):
        ring = self.ring
        try:
            other = ring.convert(other)
        except CoercionFailed:
            return NotImplemented
        if max(self.degree(), other.degree()) > query('KARATSUBA_CUTOFF'):
            return self._mul_karatsuba(other)
        return super().__mul__(other)

    def _mul_karatsuba(self, other):
        """
        Multiply dense polynomials in ``K[x]`` using Karatsuba's algorithm.

        References
        ==========

        * :cite:`Hoeven02`

        """
        ring = self.ring
        domain = ring.domain

        df = self.degree()
        dg = other.degree()

        n = max(df, dg) + 1

        n2 = n//2

        fl = self.slice(0, n2)
        gl = other.slice(0, n2)

        fh = self.slice(n2, n).quo_term(((n2,), domain.one))
        gh = other.slice(n2, n).quo_term(((n2,), domain.one))

        lo = fl*gl
        hi = fh*gh

        mid = (fl + fh)*(gl + gh)
        mid -= (lo + hi)

        return lo + mid.mul_monom((n2,)) + hi.mul_monom((2*n2,))
