"""Implementation of :class:`AlgebraicField` class. """

import functools

from ..core import Integer, sympify
from ..core.sympify import CantSympify
from ..polys.densearith import (dmp_neg, dmp_pow, dmp_rem, dup_add, dup_mul,
                                dup_sub)
from ..polys.densebasic import dmp_LC, dmp_strip, dmp_to_dict, dmp_to_tuple
from ..polys.densetools import dmp_compose
from ..polys.euclidtools import dup_invert
from ..polys.polyerrors import CoercionFailed, DomainError, NotAlgebraic
from ..printing.defaults import DefaultPrinting
from .characteristiczero import CharacteristicZero
from .domainelement import DomainElement
from .field import Field
from .simpledomain import SimpleDomain


__all__ = ('AlgebraicField',)


_algebraic_numbers_cache = {}


class AlgebraicField(Field, CharacteristicZero, SimpleDomain):
    """A class for representing algebraic number fields. """

    is_AlgebraicField = is_Algebraic = True
    is_Numerical = True

    has_assoc_Ring = False
    has_assoc_Field = True

    def __init__(self, dom, *ext):
        if not dom.is_RationalField:
            raise DomainError("ground domain must be a rational field")

        ext = [sympify(_).as_expr() for _ in ext]

        from ..polys.numberfields import primitive_element

        minpoly, coeffs, H = primitive_element(ext)

        self.ext = sum(c*e for c, e in zip(coeffs, ext))
        self.minpoly = minpoly
        self.mod = minpoly.rep
        self.domain = dom

        self.ngens = 1
        self.symbols = self.gens = (self.ext.as_expr(),)

        try:
            self.dtype = _algebraic_numbers_cache[(self.domain, self.ext)]
        except KeyError:
            self.dtype = type("AlgebraicElement", (AlgebraicElement,), {"_parent": self})
            _algebraic_numbers_cache[(self.domain, self.ext)] = self.dtype

        self.root = sum(self.dtype(h) for h in H)
        self.unit = self.dtype([dom(1), dom(0)])

        self.zero = self.dtype([dom(0)])
        self.one = self.dtype([dom(1)])

        self.rep = str(self.domain) + '<' + str(self.ext) + '>'

    def new(self, element):
        if isinstance(element, list):
            return self.dtype(element)
        else:
            return self.convert(element)

    def __hash__(self):
        return hash((self.__class__.__name__, self.domain, self.ext))

    def __eq__(self, other):
        """Returns ``True`` if two domains are equivalent. """
        return isinstance(other, AlgebraicField) and self.ext == other.ext

    def algebraic_field(self, *extension):
        r"""Returns an algebraic field, i.e. `\mathbb{Q}(\alpha, \ldots)`. """
        return self.domain.algebraic_field(*((self.ext.as_expr(),) + extension))

    def to_expr(self, a):
        """Convert ``a`` to a Diofant object. """
        return sum((c*self.ext**n for n, c in enumerate(reversed(a.rep))), Integer(0))

    def from_expr(self, a):
        """Convert Diofant's expression to ``dtype``. """
        try:
            K0 = self.domain.algebraic_field(a)
        except NotAlgebraic:
            raise CoercionFailed("%s is not a valid algebraic number in %s" % (a, self))
        return self.convert(K0.root, K0)

    def _from_PythonIntegerRing(self, a, K0):
        return self([self.domain.convert(a, K0)])

    def _from_PythonRationalField(self, a, K0):
        return self([self.domain.convert(a, K0)])

    def _from_GMPYIntegerRing(self, a, K0):
        return self([self.domain.convert(a, K0)])

    def _from_GMPYRationalField(self, a, K0):
        return self([self.domain.convert(a, K0)])

    def _from_RealField(self, a, K0):
        return self([self.domain.convert(a, K0)])

    def _from_AlgebraicField(self, a, K0):
        from ..polys import field_isomorphism

        coeffs = field_isomorphism(K0, self)

        if coeffs is not None:
            a_coeffs = dmp_compose(a.rep, coeffs, 0, self.domain)
            return self(a_coeffs)
        else:
            raise CoercionFailed("%s is not in a subfield of %s" % (K0, self))

    @property
    def ring(self):
        """Returns a ring associated with ``self``. """
        raise AttributeError('there is no ring associated with %s' % self)

    def is_positive(self, a):
        """Returns True if ``a`` is positive. """
        return self.domain.is_positive(a.LC())

    def is_negative(self, a):
        """Returns True if ``a`` is negative. """
        return self.domain.is_negative(a.LC())

    def is_nonpositive(self, a):
        """Returns True if ``a`` is non-positive. """
        return self.domain.is_nonpositive(a.LC())

    def is_nonnegative(self, a):
        """Returns True if ``a`` is non-negative. """
        return self.domain.is_nonnegative(a.LC())


class AlgebraicElement(DomainElement, CantSympify, DefaultPrinting):
    """Dense Algebraic Number Polynomials over a field. """

    def __init__(self, rep):
        dom = self.domain

        if type(rep) is not list:
            rep = [dom.convert(rep)]
        else:
            rep = [dom.convert(_) for _ in rep]

        self.rep = dmp_rem(dmp_strip(rep, 0), self.mod, 0, self.domain)

    @property
    def parent(self):
        return self._parent

    @property
    def mod(self):
        return self._parent.mod.rep

    @property
    def domain(self):
        return self._parent.domain

    def __hash__(self):
        return hash((self.__class__.__name__, dmp_to_tuple(self.rep, 0),
                     self._parent.domain, self._parent.ext))

    def per(self, rep):
        return type(self)(rep)

    def to_dict(self):
        """Convert ``self`` to a dict representation with native coefficients. """
        return dmp_to_dict(self.rep, 0, self.domain)

    def to_diofant_list(self):
        """Convert ``self`` to a list representation with Diofant coefficients. """
        return [self.domain.to_expr(c) for c in self.rep]

    @classmethod
    def from_list(cls, rep):
        return cls(dmp_strip(list(map(cls._parent.domain.convert, rep)), 0))

    def LC(self):
        """Returns the leading coefficient of ``self``. """
        return dmp_LC(self.rep, self.domain)

    @property
    def is_ground(self):
        """Returns ``True`` if ``self`` is an element of the ground domain. """
        return len(self.rep) <= 1

    def __pos__(self):
        return self

    def __neg__(self):
        return self.per(dmp_neg(self.rep, 0, self.domain))

    def __add__(self, other):
        if not isinstance(other, self.parent.dtype):
            try:
                other = self.per(other)
            except CoercionFailed:
                return NotImplemented

        return self.per(dup_add(self.rep, other.rep, self.domain))

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        if not isinstance(other, self.parent.dtype):
            try:
                other = self.per(other)
            except CoercionFailed:
                return NotImplemented

        return self.per(dup_sub(self.rep, other.rep, self.domain))

    def __rsub__(self, other):
        return (-self).__add__(other)

    def __mul__(self, other):
        if not isinstance(other, self.parent.dtype):
            try:
                other = self.per(other)
            except CoercionFailed:
                return NotImplemented

        return self.per(dmp_rem(dup_mul(self.rep, other.rep, self.domain),
                                self.mod, 0, self.domain))

    def __rmul__(self, other):
        return self.__mul__(other)

    def __pow__(self, n):
        if isinstance(n, int):
            if n < 0:
                F, n = dup_invert(self.rep, self.mod, self.domain), -n
            else:
                F = self.rep

            return self.per(dmp_rem(dmp_pow(F, n, 0, self.domain),
                                    self.mod, 0, self.domain))
        else:
            raise TypeError("``int`` expected, got %s" % type(n))

    def __truediv__(self, other):
        if not isinstance(other, self.parent.dtype):
            try:
                other = self.per(other)
            except CoercionFailed:
                return NotImplemented

        return self.per(dmp_rem(dup_mul(self.rep, dup_invert(other.rep, self.mod, self.domain), self.domain),
                                self.mod, 0, self.domain))

    def __eq__(self, other):
        if not isinstance(other, self.parent.dtype):
            try:
                other = self.per(other)
            except CoercionFailed:
                return False

        return self.rep == other.rep

    def __lt__(self, other):
        return self.rep.__lt__(other.rep)

    def __bool__(self):
        return bool(self.rep)

    @property
    def numerator(self):
        return self*self.denominator

    @property
    def denominator(self):
        return self.per(functools.reduce(self.domain.ring.lcm,
                                         (_.denominator for _ in self.rep),
                                         self.domain.ring.one))
