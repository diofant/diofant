"""Implementation of :class:`PolynomialRing` class. """

from diofant.polys.domains.ring import Ring
from diofant.polys.domains.compositedomain import CompositeDomain
from diofant.polys.domains.characteristiczero import CharacteristicZero
from diofant.polys.domains.old_fractionfield import FractionField
from diofant.polys.polyclasses import DMP, DMF
from diofant.polys.polyerrors import (GeneratorsNeeded, PolynomialError,
                                      CoercionFailed, ExactQuotientFailed,
                                      NotReversible)
from diofant.polys.polyutils import (dict_from_basic, basic_from_dict,
                                     _dict_reorder)
from diofant.polys.orderings import monomial_key, build_product_order
from diofant.polys.agca.modules import FreeModulePolyRing
from diofant.core.compatibility import iterable
from diofant.utilities import public

# XXX why does this derive from CharacteristicZero???


@public
class PolynomialRingBase(Ring, CharacteristicZero, CompositeDomain):
    """
    Base class for generalized polynomial rings.

    This base class should be used for uniform access to generalized polynomial
    rings. Subclasses only supply information about the element storage etc.

    Do not instantiate.
    """

    has_assoc_Ring = True
    has_assoc_Field = True

    default_order = "grevlex"

    def __init__(self, dom, *gens, **opts):
        if not gens:
            raise GeneratorsNeeded("generators not specified")

        lev = len(gens) - 1
        self.ngens = len(gens)

        self.zero = self.dtype.zero(lev, dom, ring=self)
        self.one = self.dtype.one(lev, dom, ring=self)

        self.domain = dom
        self.symbols = self.gens = gens
        # NOTE 'order' may not be set if inject was called through CompositeDomain
        self.order = opts.get('order', monomial_key(self.default_order))

    def new(self, element):
        return self.dtype(element, self.domain, len(self.gens) - 1, ring=self)

    def __str__(self):
        s_order = str(self.order)
        orderstr = (
            " order=" + s_order) if s_order != self.default_order else ""
        return str(self.domain) + '[' + ','.join(map(str, self.gens)) + orderstr + ']'

    def __hash__(self):
        return hash((self.__class__.__name__, self.dtype, self.domain,
                     self.gens, self.order))

    def __eq__(self, other):
        """Returns `True` if two domains are equivalent. """
        return isinstance(other, PolynomialRingBase) and \
            self.dtype == other.dtype and self.domain == other.domain and \
            self.gens == other.gens and self.order == other.order

    def from_ZZ_python(self, a, K0):
        """Convert a Python `int` object to `dtype`. """
        return self(self.domain.convert(a, K0))

    def from_ZZ_gmpy(self, a, K0):
        """Convert a GMPY `mpz` object to `dtype`. """
        return self(self.domain.convert(a, K0))

    def from_GlobalPolynomialRing(self, a, K0):
        """Convert a `DMP` object to `dtype`. """
        if self.gens == K0.gens:
            if self.domain == K0.domain:
                return self(a.rep)  # set the correct ring
            else:
                return self(a.convert(self.domain).rep)
        else:
            monoms, coeffs = _dict_reorder(a.to_dict(), K0.gens, self.gens)

            if self.domain != K0.domain:
                coeffs = [self.domain.convert(c, K0.domain) for c in coeffs]

            return self(dict(zip(monoms, coeffs)))

    def get_field(self):
        """Returns a field associated with `self`. """
        return FractionField(self.domain, *self.gens)

    def revert(self, a):
        try:
            return 1/a
        except (ExactQuotientFailed, ZeroDivisionError):
            raise NotReversible('%s is not a unit' % a)

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

        >>> R = QQ.old_poly_ring(x, y, order=ilex)
        >>> L = [((1, 1, 1), QQ(1)), ((0, 1, 0), QQ(1)), ((0, 0, 1), QQ(2))]
        >>> R._sdm_to_vector(L, 2)
        [x + 2*y, x*y]
        """
        dics = self._sdm_to_dics(s, n)
        # NOTE this works for global and local rings!
        return [self(x) for x in dics]

    def free_module(self, rank):
        """
        Generate a free module of rank ``rank`` over ``self``.

        >>> from diofant.abc import x
        >>> from diofant import QQ
        >>> QQ.old_poly_ring(x).free_module(2)
        QQ[x]**2
        """
        return FreeModulePolyRing(self, rank)


def _vector_to_sdm_helper(v, order):
    """Helper method for common code in Global and Local poly rings."""
    from diofant.polys.distributedmodules import sdm_from_dict
    d = {}
    for i, e in enumerate(v):
        for key, value in e.to_dict().items():
            d[(i,) + key] = value
    return sdm_from_dict(d, order)


@public
class GlobalPolynomialRing(PolynomialRingBase):
    """A true polynomial ring, with objects DMP. """

    is_PolynomialRing = is_Poly = True
    dtype = DMP

    def from_FractionField(self, a, K0):
        """
        Convert a ``DMF`` object to ``DMP``.

        Examples
        ========

        >>> from diofant.polys.polyclasses import DMP, DMF
        >>> from diofant.polys.domains import ZZ
        >>> from diofant.abc import x

        >>> f = DMF(([ZZ(1), ZZ(1)], [ZZ(1)]), ZZ)
        >>> K = ZZ.old_frac_field(x)

        >>> F = ZZ.old_poly_ring(x).from_FractionField(f, K)

        >>> F == DMP([ZZ(1), ZZ(1)], ZZ)
        True
        >>> type(F)
        <class 'diofant.polys.polyclasses.DMP'>

        """
        if a.denom().is_one:
            return self.from_GlobalPolynomialRing(a.numer(), K0)

    def to_diofant(self, a):
        """Convert `a` to a Diofant object. """
        return basic_from_dict(a.to_diofant_dict(), *self.gens)

    def from_diofant(self, a):
        """Convert Diofant's expression to `dtype`. """
        try:
            rep, _ = dict_from_basic(a, gens=self.gens)
        except PolynomialError:
            raise CoercionFailed("can't convert %s to type %s" % (a, self))

        for k, v in rep.items():
            rep[k] = self.domain.from_diofant(v)

        return self(rep)

    def _vector_to_sdm(self, v, order):
        """
        >>> from diofant import lex, QQ
        >>> from diofant.abc import x, y

        >>> R = QQ.old_poly_ring(x, y)
        >>> f = R.convert(x + 2*y)
        >>> g = R.convert(x * y)
        >>> R._vector_to_sdm([f, g], lex)
        [((1, 1, 1), 1), ((0, 1, 0), 1), ((0, 0, 1), 2)]
        """
        return _vector_to_sdm_helper(v, order)


class GeneralizedPolynomialRing(PolynomialRingBase):
    """A generalized polynomial ring, with objects DMF. """

    dtype = DMF

    def new(self, a):
        """Construct an element of `self` domain from `a`. """
        res = self.dtype(a, self.domain, len(self.gens) - 1, ring=self)

        # make sure res is actually in our ring
        if res.denom().terms(order=self.order)[0][0] != (0,)*len(self.gens):
            from diofant.printing.str import sstr
            raise CoercionFailed("denominator %s not allowed in %s"
                                 % (sstr(res), self))
        return res

    def __contains__(self, a):
        try:
            a = self.convert(a)
        except CoercionFailed:
            return False
        return a.denom().terms(order=self.order)[0][0] == (0,)*len(self.gens)

    def from_FractionField(self, a, K0):
        dmf = self.get_field().from_FractionField(a, K0)
        return self((dmf.num, dmf.den))

    def to_diofant(self, a):
        """Convert `a` to a Diofant object. """
        return (basic_from_dict(a.numer().to_diofant_dict(), *self.gens) /
                basic_from_dict(a.denom().to_diofant_dict(), *self.gens))

    def from_diofant(self, a):
        """Convert Diofant's expression to `dtype`. """
        p, q = a.as_numer_denom()

        num, _ = dict_from_basic(p, gens=self.gens)
        den, _ = dict_from_basic(q, gens=self.gens)

        for k, v in num.items():
            num[k] = self.domain.from_diofant(v)

        for k, v in den.items():
            den[k] = self.domain.from_diofant(v)

        return self((num, den)).cancel()

    def _vector_to_sdm(self, v, order):
        """
        Turn an iterable into a sparse distributed module.

        Note that the vector is multiplied by a unit first to make all entries
        polynomials.

        >>> from diofant import ilex, QQ
        >>> from diofant.abc import x, y

        >>> R = QQ.old_poly_ring(x, y, order=ilex)
        >>> f = R.convert((x + 2*y) / (1 + x))
        >>> g = R.convert(x * y)
        >>> R._vector_to_sdm([f, g], ilex)
        [((0, 0, 1), 2), ((0, 1, 0), 1), ((1, 1, 1), 1), ((1,
          2, 1), 1)]
        """
        # NOTE this is quite inefficient...
        u = self.one.numer()
        for x in v:
            u *= x.denom()
        return _vector_to_sdm_helper([x.numer()*u/x.denom() for x in v], order)


@public
def PolynomialRing(dom, *gens, **opts):
    r"""
    Create a generalized multivariate polynomial ring.

    A generalized polynomial ring is defined by a ground field `K`, a set
    of generators (typically `x_1, \ldots, x_n`) and a monomial order `<`.
    The monomial order can be global, local or mixed. In any case it induces
    a total ordering on the monomials, and there exists for every (non-zero)
    polynomial `f \in K[x_1, \ldots, x_n]` a well-defined "leading monomial"
    `LM(f) = LM(f, >)`. One can then define a multiplicative subset
    `S = S_> = \{f \in K[x_1, \ldots, x_n] | LM(f) = 1\}`. The generalized
    polynomial ring corresponding to the monomial order is
    `R = S^{-1}K[x_1, \ldots, x_n]`.

    If `>` is a so-called global order, that is `1` is the smallest monomial,
    then we just have `S = K` and `R = K[x_1, \ldots, x_n]`.

    Examples
    ========

    A few examples may make this clearer.

    >>> from diofant.abc import x, y
    >>> from diofant import QQ

    Our first ring uses global lexicographic order.

    >>> R1 = QQ.old_poly_ring(x, y, order=(("lex", x, y),))

    The second ring uses local lexicographic order. Note that when using a
    single (non-product) order, you can just specify the name and omit the
    variables:

    >>> R2 = QQ.old_poly_ring(x, y, order="ilex")

    The third and fourth rings use a mixed orders:

    >>> o1 = (("ilex", x), ("lex", y))
    >>> o2 = (("lex", x), ("ilex", y))
    >>> R3 = QQ.old_poly_ring(x, y, order=o1)
    >>> R4 = QQ.old_poly_ring(x, y, order=o2)

    We will investigate what elements of `K(x, y)` are contained in the various
    rings.

    >>> L = [x, 1/x, y/(1 + x), 1/(1 + y), 1/(1 + x*y)]
    >>> test = lambda R: [f in R for f in L]

    The first ring is just `K[x, y]`:

    >>> test(R1)
    [True, False, False, False, False]

    The second ring is R1 localised at the maximal ideal (x, y):

    >>> test(R2)
    [True, False, True, True, True]

    The third ring is R1 localised at the prime ideal (x):

    >>> test(R3)
    [True, False, True, False, True]

    Finally the fourth ring is R1 localised at `S = K[x, y] \setminus yK[y]`:

    >>> test(R4)
    [True, False, False, True, False]
    """

    order = opts.get("order", GeneralizedPolynomialRing.default_order)
    if iterable(order):
        order = build_product_order(order, gens)
    order = monomial_key(order)
    opts['order'] = order

    if order.is_global:
        return GlobalPolynomialRing(dom, *gens, **opts)
    else:
        return GeneralizedPolynomialRing(dom, *gens, **opts)
