"""Implementation of :class:`PolynomialRing` class. """

from diofant.polys.domains.ring import Ring
from diofant.polys.domains.compositedomain import CompositeDomain
from diofant.polys.polyerrors import CoercionFailed, GeneratorsError
from diofant.utilities import public


@public
class PolynomialRing(Ring, CompositeDomain):
    """A class for representing multivariate polynomial rings. """

    is_PolynomialRing = is_Poly = True

    has_assoc_Ring  = True
    has_assoc_Field = True

    def __init__(self, domain_or_ring, symbols=None, order=None):
        from diofant.polys.rings import PolyRing

        if isinstance(domain_or_ring, PolyRing) and symbols is None and order is None:
            ring = domain_or_ring
        else:
            ring = PolyRing(symbols, domain_or_ring, order)

        self.ring = ring
        self.dtype = ring.dtype

        self.gens = ring.gens
        self.ngens = ring.ngens
        self.symbols = ring.symbols
        self.domain = ring.domain

    def new(self, element):
        return self.ring.ring_new(element)

    @property
    def zero(self):
        return self.ring.zero

    @property
    def one(self):
        return self.ring.one

    @property
    def order(self):
        return self.ring.order

    def __str__(self):
        return str(self.domain) + '[' + ','.join(map(str, self.symbols)) + ']'

    def __hash__(self):
        return hash((self.__class__.__name__, self.dtype, self.domain, self.symbols))

    def __eq__(self, other):
        """Returns `True` if two domains are equivalent. """
        return isinstance(other, PolynomialRing) and \
            self.dtype == other.dtype and self.ring == other.ring

    def to_diofant(self, a):
        """Convert `a` to a Diofant object. """
        return a.as_expr()

    def from_diofant(self, a):
        """Convert Diofant's expression to `dtype`. """
        return self.ring.from_expr(a)

    def from_ZZ_python(self, a, K0):
        """Convert a Python `int` object to `dtype`. """
        return self(self.domain.convert(a, K0))

    def from_QQ_python(self, a, K0):
        """Convert a Python `Fraction` object to `dtype`. """
        return self(self.domain.convert(a, K0))

    def from_ZZ_gmpy(self, a, K0):
        """Convert a GMPY `mpz` object to `dtype`. """
        return self(self.domain.convert(a, K0))

    def from_QQ_gmpy(self, a, K0):
        """Convert a GMPY `mpq` object to `dtype`. """
        return self(self.domain.convert(a, K0))

    def from_RealField(self, a, K0):
        """Convert a mpmath `mpf` object to `dtype`. """
        return self(self.domain.convert(a, K0))

    def from_AlgebraicField(self, a, K0):
        """Convert an algebraic number to ``dtype``. """
        if self.domain == K0:
            return self.new(a)

    def from_PolynomialRing(self, a, K0):
        """Convert a polynomial to ``dtype``. """
        try:
            return a.set_ring(self.ring)
        except (CoercionFailed, GeneratorsError):
            return

    def from_FractionField(self, a, K0):
        """Convert a rational function to ``dtype``. """
        denom = K0.denom(a)

        if denom.is_ground:
            return self.from_PolynomialRing(K0.numer(a)/denom, K0.field.ring.to_domain())

    def get_field(self):
        """Returns a field associated with `self`. """
        return self.ring.to_field().to_domain()

    def is_positive(self, a):
        """Returns True if `LC(a)` is positive. """
        return self.domain.is_positive(a.LC)

    def is_negative(self, a):
        """Returns True if `LC(a)` is negative. """
        return self.domain.is_negative(a.LC)

    def is_nonpositive(self, a):
        """Returns True if `LC(a)` is non-positive. """
        return self.domain.is_nonpositive(a.LC)

    def is_nonnegative(self, a):
        """Returns True if `LC(a)` is non-negative. """
        return self.domain.is_nonnegative(a.LC)

    def gcdex(self, a, b):
        """Extended GCD of `a` and `b`. """
        return a.gcdex(b)

    def gcd(self, a, b):
        """Returns GCD of `a` and `b`. """
        return a.gcd(b)

    def lcm(self, a, b):
        """Returns LCM of `a` and `b`. """
        return a.lcm(b)

    def factorial(self, a):
        """Returns factorial of `a`. """
        return self.dtype(self.domain.factorial(a))

    def free_module(self, rank):
        """
        Generate a free module of rank ``rank`` over ``self``.

        >>> from diofant.abc import x
        >>> from diofant import QQ
        >>> QQ.poly_ring(x).free_module(2)
        QQ[x]**2
        """
        from diofant.polys.agca.modules import FreeModulePolyRing
        return FreeModulePolyRing(self, rank)

    def _vector_to_sdm(self, v, order):
        """
        Turn an iterable into a sparse distributed module.

        Note that the vector is multiplied by a unit first to make all entries
        polynomials.

        >>> from diofant import ilex, QQ
        >>> from diofant.abc import x, y
        >>> from diofant.printing import pprint
        >>> R = QQ.poly_ring(x, y, order=ilex)
        >>> f = R.convert((x + 2*y) / (1 + x))
        >>> g = R.convert(x * y)
        >>> pprint(R._vector_to_sdm([f, g], ilex))
        [((0, 0, 1), 2), ((0, 1, 0), 1), ((1, 1, 1), 1), ((1,
          2, 1), 1)]
        """
        # NOTE this is quite inefficient...
        u = self.one.numer()
        for x in v:
            u *= x.denom()
        return _vector_to_sdm_helper([x.numer()*u/x.denom() for x in v], order)

    def _vector_to_sdm(self, v, order):
        """
        >>> from diofant import lex, QQ
        >>> from diofant.abc import x, y
        >>> from diofant.printing import pprint
        >>> R = QQ.poly_ring(x, y)
        >>> f = R.convert(x + 2*y)
        >>> g = R.convert(x * y)
        >>> pprint(R._vector_to_sdm([f, g], lex))
        [((1, 1, 1), 1), ((0, 1, 0), 1), ((0, 0, 1), 2)]
        """
        return _vector_to_sdm_helper(v, order)

    def _sdm_to_dics(self, s, n):
        """Helper for _sdm_to_vector."""
        from diofant.polys.distributedmodules import sdm_to_dict
        dic = sdm_to_dict(s)
        res = [{} for _ in range(n)]
        for k, v in dic.items():
            res[k[0]][k[1:]] = v
        return res

    def _sdm_to_vector(self, s, n):
        """
        For internal use by the modules class.

        Convert a sparse distributed module into a list of length ``n``.

        >>> from diofant import QQ, ilex
        >>> from diofant.abc import x, y
        >>> R = QQ.poly_ring(x, y, order=ilex)
        >>> L = [((1, 1, 1), QQ(1)), ((0, 1, 0), QQ(1)), ((0, 0, 1), QQ(2))]
        >>> R._sdm_to_vector(L, 2)
        [x + 2*y, x*y]
        """
        dics = self._sdm_to_dics(s, n)
        # NOTE this works for global and local rings!
        return [self(x) for x in dics]



def _vector_to_sdm_helper(v, order):
    """Helper method for common code in Global and Local poly rings."""
    from diofant.polys.distributedmodules import sdm_from_dict
    d = {}
    for i, e in enumerate(v):
        for key, value in e.to_dict().items():
            d[(i,) + key] = value
    return sdm_from_dict(d, order)
