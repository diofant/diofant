"""Implementation of :class:`AlgebraicField` class. """

import functools
import operator

from ..core import I, Integer, sympify
from ..core.sympify import CantSympify
from ..polys.densearith import (dmp_neg, dmp_pow, dmp_rem, dup_add, dup_mul,
                                dup_sub)
from ..polys.densebasic import dmp_LC, dmp_strip, dmp_to_dict, dmp_to_tuple
from ..polys.densetools import dmp_compose, dmp_eval_in
from ..polys.euclidtools import dup_invert
from ..polys.polyerrors import CoercionFailed, DomainError, NotAlgebraic
from ..printing.defaults import DefaultPrinting
from .characteristiczero import CharacteristicZero
from .domainelement import DomainElement
from .field import Field
from .rationalfield import RationalField
from .simpledomain import SimpleDomain


__all__ = ('AlgebraicField',)


_algebraic_numbers_cache = {}


class AlgebraicField(Field, CharacteristicZero, SimpleDomain):
    """A class for representing algebraic number fields. """

    is_AlgebraicField = is_Algebraic = True
    is_Numerical = True

    has_assoc_Ring = False
    has_assoc_Field = True

    def __new__(cls, dom, *ext):
        if not (dom.is_RationalField or dom.is_AlgebraicField):
            raise DomainError("ground domain must be a rational "
                              "or an algebraic field")

        ext = [sympify(_).as_expr() for _ in ext]
        ext = [_ for _ in ext if _ not in dom]

        if not ext:
            return dom

        from ..polys.numberfields import primitive_element

        minpoly, coeffs, H = primitive_element(ext, domain=dom)
        ext = sum(c*e for c, e in zip(coeffs, ext))

        is_real = ext.is_real
        if is_real is None:
            ext_root = cls._compute_ext_root(ext, minpoly)
            is_real = ext_root[1].is_real
        else:
            ext_root = None

        if is_real:
            new_cls = RealAlgebraicField
        elif isinstance(dom, (RationalField, RealAlgebraicField)) and ext == I:
            new_cls = ComplexAlgebraicField
        else:
            new_cls = cls
        obj = super().__new__(new_cls)

        obj.ext = ext
        obj._ext_root = ext_root
        obj.minpoly = minpoly
        obj.mod = minpoly.rep
        obj.domain = dom

        obj.ngens = 1
        obj.symbols = obj.gens = (obj.ext.as_expr(),)

        try:
            obj.dtype = _algebraic_numbers_cache[(obj.domain, obj.ext)]
        except KeyError:
            if new_cls == RealAlgebraicField:
                dtype_cls = RealAlgebraicElement
            elif new_cls == ComplexAlgebraicField:
                dtype_cls = ComplexAlgebraicElement
            else:
                dtype_cls = AlgebraicElement
            obj.dtype = type(dtype_cls.__name__, (dtype_cls,), {"_parent": obj})
            _algebraic_numbers_cache[(obj.domain, obj.ext)] = obj.dtype

        obj.root = sum(obj.dtype(h) for h in H)
        obj.unit = obj.dtype([dom(1), dom(0)])

        obj.zero = obj.dtype([dom(0)])
        obj.one = obj.dtype([dom(1)])

        obj.rep = str(obj.domain) + '<' + str(obj.ext) + '>'

        return obj

    def new(self, element):
        if isinstance(element, list):
            return self.dtype(element)
        else:
            return self.convert(element)

    def __hash__(self):
        return hash((self.__class__.__name__, self.domain, self.ext))

    def __eq__(self, other):
        """Returns ``True`` if two domains are equivalent. """
        return isinstance(other, AlgebraicField) and self.domain == other.domain and self.ext == other.ext

    def algebraic_field(self, *extension):
        r"""Returns an algebraic field, i.e. `\mathbb{Q}(\alpha, \ldots)`. """
        return AlgebraicField(self, *extension)

    def to_expr(self, a):
        """Convert ``a`` to a Diofant object. """
        return sum(((a.domain.to_expr(c)*self.ext**n).expand() for n, c in enumerate(reversed(a.rep))), Integer(0))

    def from_expr(self, a):
        """Convert Diofant's expression to ``dtype``. """
        try:
            K0 = self.domain.algebraic_field(a)
        except NotAlgebraic:
            raise CoercionFailed("%s is not a valid algebraic number in %s" % (a, self))
        if a in self.domain:
            return self.new([a])
        else:
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
        if K0 == self.domain:
            return self([a])
        elif self == K0.domain and len(a.rep) <= 1:
            return a.rep[0] if a else self.zero

        from ..polys import field_isomorphism

        coeffs = field_isomorphism(K0, self)

        if coeffs is not None:
            if K0.domain == self.domain:
                return self(dmp_compose(a.rep, coeffs, 0, self.domain))
            else:
                return self.from_expr(K0.to_expr(a))
        else:
            raise CoercionFailed("%s is not in a subfield of %s" % (K0, self))

    def _from_ExpressionDomain(self, a, K0):
        expr = K0.to_expr(a)
        return self.from_expr(expr)

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

    @staticmethod
    def _compute_ext_root(ext, minpoly):
        from ..polys import minimal_polynomial

        for r in minpoly.all_roots(radicals=False):  # pragma: no branch
            if not minimal_polynomial(ext - r)(0):
                return r.as_content_primitive()

    @property
    def ext_root(self):
        if self._ext_root is None:
            self._ext_root = self._compute_ext_root(self.ext, self.minpoly)
        return self._ext_root


class RealAlgebraicField(AlgebraicField):
    """A class for representing real algebraic number fields. """

    is_RealAlgebraicField = True

    def is_positive(self, a):
        """Returns True if ``a`` is positive. """
        return a > 0

    def is_negative(self, a):
        """Returns True if ``a`` is negative. """
        return a < 0


class ComplexAlgebraicField(AlgebraicField):
    """A class for representing complex algebraic number fields. """

    is_ComplexAlgebraicField = True


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

    def __bool__(self):
        return bool(self.rep)

    def __int__(self):
        try:
            from . import QQ
            return int(QQ.convert(self))
        except CoercionFailed:
            raise TypeError("Can't convert algebraic number to int")

    @property
    def numerator(self):
        return self*self.denominator

    @property
    def denominator(self):
        from . import ZZ
        return self.per(functools.reduce(ZZ.lcm, (ZZ.convert(_.denominator)
                                                  for _ in self.rep), ZZ.one))


class ComplexAlgebraicElement(AlgebraicElement):
    """Elements of complex algebraic numbers field. """

    @property
    def real(self):
        """Returns real part of ``self``. """
        return self.domain.convert(self.rep[-1]) if self else self.domain.zero

    @property
    def imag(self):
        """Returns imaginary part of ``self``. """
        return self.domain.convert((self - self.real)/self.parent.unit)

    def conjugate(self):
        """Returns the complex conjugate of ``self``. """
        return self.real - self.parent.unit*self.imag


class RealAlgebraicElement(ComplexAlgebraicElement):
    """Elements of real algebraic numbers field. """

    def __abs__(self):
        return self if self >= 0 else -self

    def _cmp(self, other, op):
        from ..polys.rootisolation import dup_count_real_roots

        if not isinstance(other, self.parent.dtype):
            try:
                other = self.per(other)
            except CoercionFailed:
                return NotImplemented

        diff = self - other
        rep = dmp_compose(diff.rep,
                          [self.domain.from_expr(self.parent.ext_root[0]), 0],
                          0, self.domain)
        while dup_count_real_roots(rep, self.domain,
                                   inf=self.parent.ext_root[1].interval.a,
                                   sup=self.parent.ext_root[1].interval.b):
            self.parent.ext_root[1].refine()
        v = dmp_eval_in(rep, diff.parent.ext_root[1].interval.center,
                        0, 0, diff.domain)
        return bool(op(v, 0))

    def __lt__(self, other):
        return self._cmp(other, operator.lt)

    def __le__(self, other):
        return self._cmp(other, operator.le)

    def __gt__(self, other):
        return self._cmp(other, operator.gt)

    def __ge__(self, other):
        return self._cmp(other, operator.ge)

    @property
    def real(self):
        """Returns real part of ``self``. """
        return self

    @property
    def imag(self):
        """Returns imaginary part of ``self``. """
        return self.parent.zero
