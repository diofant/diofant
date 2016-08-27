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
