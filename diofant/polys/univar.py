import random

from ..config import query
from ..domains import ZZ
from .polyerrors import CoercionFailed, DomainError
from .rings import PolyElement, PolynomialRing
from .rootisolation import _FindRoot


class UnivarPolynomialRing(PolynomialRing, _FindRoot):
    """A class for representing univariate polynomial rings."""

    def __call__(self, element):
        if isinstance(element, list):
            try:
                return self.from_terms(element)
            except (TypeError, ValueError):
                return self.from_list(element)
        return super().__call__(element)

    def from_list(self, element):
        return self.from_dict({(i,): c for i, c in enumerate(element)})

    def _random(self, n, a, b, percent=None):
        domain = self.domain

        if percent is None:
            percent = 100//(b - a)
        percent = min(max(0, percent), 100)
        nz = ((n + 1)*percent)//100

        f = []
        while len(f) < n + 1:
            v = domain.convert(random.randint(a, b))
            if v:
                f.append(v)

        if nz:
            f[-nz:] = [domain.zero]*nz
            lt = f.pop(0)
            random.shuffle(f)
            f.insert(0, lt)

        return self.from_list(list(reversed(f)))

    def _gf_random(self, n, irreducible=False):
        domain = self.domain

        assert domain.is_FiniteField

        while True:
            f = [domain(random.randint(0, domain.order - 1))
                 for i in range(n)] + [domain.one]
            f = self.from_list(f)
            if not irreducible or f.is_irreducible:
                return f

    def dispersionset(self, p, q=None):
        r"""Compute the *dispersion set* of two polynomials.

        For two polynomials `f(x)` and `g(x)` with `\deg f > 0`
        and `\deg g > 0` the dispersion set `\operatorname{J}(f, g)` is defined as:

        .. math::
            \operatorname{J}(f, g)
            & := \{a \in \mathbb{N}_0 | \gcd(f(x), g(x+a)) \neq 1\} \\
            &  = \{a \in \mathbb{N}_0 | \deg \gcd(f(x), g(x+a)) \geq 1\}

        For a single polynomial one defines `\operatorname{J}(f) := \operatorname{J}(f, f)`.

        Examples
        ========

        Note that the definition of the dispersion is not symmetric:

        >>> R, x = ring('x', QQ)

        >>> fp = x**4 - 3*x**2 + 1
        >>> gp = fp.shift(-3)

        >>> R.dispersionset(fp, gp)
        {2, 3, 4}
        >>> R.dispersionset(gp, fp)
        set()

        Computing the dispersion also works over field extensions:

        >>> R, x = ring('x', QQ.algebraic_field(sqrt(5)))

        >>> fp = x**2 + sqrt(5)*x - 1
        >>> gp = x**2 + (2 + sqrt(5))*x + sqrt(5)

        >>> R.dispersionset(fp, gp)
        {2}
        >>> R.dispersionset(gp, fp)
        {1, 4}

        We can even perform the computations for polynomials
        having symbolic coefficients:

        >>> D, a = ring('a', QQ)
        >>> R, x = ring('x', D)

        >>> fp = 4*x**4 + (4*a + 8)*x**3 + (a**2 + 6*a + 4)*x**2 + (a**2 + 2*a)*x
        >>> R.dispersionset(fp)
        {0, 1}

        References
        ==========

        * :cite:`Man1994disp`
        * :cite:`Koepf98`
        * :cite:`Abramov71rat`
        * :cite:`Man1993indefsum`

        """
        # Check for valid input
        same = False if q is not None else True
        if same:
            q = p

        if p.ring is not q.ring:
            raise ValueError('Polynomials must have the same ring')

        fdomain = self.domain.field

        # We define the dispersion of constant polynomials to be zero
        if p.degree() < 1 or q.degree() < 1:
            return {0}

        # Factor p and q over the rationals
        fp = p.factor_list()
        fq = q.factor_list() if not same else fp

        # Iterate over all pairs of factors
        J = set()
        for s, unused in fp[1]:
            for t, unused in fq[1]:
                m = s.degree()
                n = t.degree()
                if n != m:
                    continue
                an = s.LC
                bn = t.LC
                if an - bn:
                    continue
                # Note that the roles of `s` and `t` below are switched
                # w.r.t. the original paper. This is for consistency
                # with the description in the book of W. Koepf.
                anm1 = s[(m - 1,)]
                bnm1 = t[(n - 1,)]
                alpha = fdomain(anm1 - bnm1)/fdomain(n*bn)
                if alpha not in ZZ:
                    continue

                alpha = ZZ.convert(alpha)

                if alpha < 0 or alpha in J:
                    continue
                if n > 1 and s - t.shift(alpha):
                    continue
                J.add(alpha)

        return J


class UnivarPolyElement(PolyElement):
    """Element of univariate distributed polynomial ring."""

    def all_coeffs(self):
        if self:
            return [self[(i,)] for i in range(self.degree() + 1)]
        else:
            return [self.parent.domain.zero]

    def shift(self, a):
        return self.compose(0, self.ring.gens[0] + a)

    def half_gcdex(self, other):
        """
        Half extended Euclidean algorithm in `F[x]`.

        Returns ``(s, h)`` such that ``h = gcd(self, other)``
        and ``s*self = h (mod other)``.

        Examples
        ========

        >>> _, x = ring('x', QQ)

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

        >>> _, x = ring('x', ZZ)

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
                result  # XXX "peephole" optimization, http://bugs.python.org/issue2506
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

        >>> _, x = ring('x', QQ)

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
