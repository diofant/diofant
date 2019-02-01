"""Implementation of :class:`AlgebraicField` class. """

import functools
import numbers

from ..core import I, Integer, sympify
from ..core.sympify import CantSympify
from ..polys.densetools import dmp_compose, dmp_eval_in
from ..polys.polyerrors import CoercionFailed, DomainError, NotAlgebraic
from .characteristiczero import CharacteristicZero
from .domainelement import DomainElement
from .field import Field
from .rationalfield import RationalField
from .simpledomain import SimpleDomain


__all__ = 'AlgebraicField',


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

        minpoly, coeffs, _ = primitive_element(ext, domain=dom)
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
        obj.domain = dom

        obj.ngens = 1
        obj.symbols = obj.gens = obj.ext.as_expr(),

        rep_ring = dom.poly_ring(obj.ext)
        obj.mod = mod = rep_ring.from_dense(minpoly.rep.rep)

        try:
            obj.dtype = _algebraic_numbers_cache[(obj.domain, obj.ext)]
        except KeyError:
            if new_cls == RealAlgebraicField:
                dtype_cls = RealAlgebraicElement
            elif new_cls == ComplexAlgebraicField:
                dtype_cls = ComplexAlgebraicElement
            else:
                dtype_cls = AlgebraicElement
            obj.dtype = type(dtype_cls.__name__, (dtype_cls,),
                             {"mod": mod, "domain": rep_ring, "_parent": obj})
            _algebraic_numbers_cache[(obj.domain, obj.ext)] = obj.dtype

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
        return sum(((self.domain.to_expr(c)*self.ext**n).expand()
                    for n, c in enumerate(reversed(a.rep.to_dense()))), Integer(0))

    def from_expr(self, a):
        """Convert Diofant's expression to ``dtype``. """
        try:
            K0 = self.domain.algebraic_field(a)
        except NotAlgebraic:
            raise CoercionFailed("%s is not a valid algebraic number in %s" % (a, self))
        if a in self.domain:
            return self.new([a])
        else:
            from ..polys import field_isomorphism

            coeffs = field_isomorphism(K0, self)
            factor = Integer((K0.to_expr(K0.unit)/a).simplify())

            return self.dtype(coeffs)/factor

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
            return a.rep.coeff(1) if a else self.zero

        from ..polys import field_isomorphism

        coeffs = field_isomorphism(K0, self)

        if coeffs is not None:
            if K0.domain == self.domain:
                return self(dmp_compose(a.rep.to_dense(), coeffs, 0, self.domain))
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


class ComplexAlgebraicField(AlgebraicField):
    """A class for representing complex algebraic number fields. """

    is_ComplexAlgebraicField = True


class RealAlgebraicField(ComplexAlgebraicField):
    """A class for representing real algebraic number fields. """

    is_RealAlgebraicField = True

    def is_positive(self, a):
        """Returns True if ``a`` is positive. """
        return a > 0

    def is_negative(self, a):
        """Returns True if ``a`` is negative. """
        return a < 0


class AlgebraicElement(DomainElement, CantSympify):
    """Dense Algebraic Number Polynomials over a field. """

    def __init__(self, rep):
        dom = self.domain

        if isinstance(rep, dict):
            rep = dom.from_dict(rep)
        else:
            if type(rep) is not list:
                rep = [dom.domain.convert(rep)]
            else:
                rep = [dom.domain.convert(_) for _ in rep]
            rep = dom.from_dense(rep)

        self.rep = rep % self.mod

    @property
    def parent(self):
        return self._parent

    def __hash__(self):
        return hash((self.__class__.__name__, self.rep, self.domain.domain,
                     self.parent.ext))

    def to_dict(self):
        """Convert ``self`` to a dict representation with native coefficients. """
        return self.rep.to_dict()

    def LC(self):
        """Returns the leading coefficient of ``self``. """
        return self.rep.LC

    @property
    def is_ground(self):
        """Returns ``True`` if ``self`` is an element of the ground domain. """
        return self.rep.is_ground

    def __pos__(self):
        return self

    def __neg__(self):
        return self.__class__(-self.rep)

    def __add__(self, other):
        try:
            other = self.parent.convert(other)
        except CoercionFailed:
            return NotImplemented
        return self.__class__(self.rep + other.rep)

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        try:
            other = self.parent.convert(other)
        except CoercionFailed:
            return NotImplemented
        return self.__class__(self.rep - other.rep)

    def __rsub__(self, other):
        return (-self).__add__(other)

    def __mul__(self, other):
        try:
            other = self.parent.convert(other)
        except CoercionFailed:
            return NotImplemented
        return self.__class__(self.rep * other.rep)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __pow__(self, exp):
        if not isinstance(exp, numbers.Integral):
            raise TypeError("Integer exponent expected, got %s" % type(exp))
        if exp < 0:
            rep, exp = self.domain.invert(self.rep, self.mod), -exp
        else:
            rep = self.rep
        return self.__class__(rep**exp)

    def __truediv__(self, other):
        try:
            other = self.parent.convert(other)
        except CoercionFailed:
            return NotImplemented
        return self * other**-1

    def __eq__(self, other):
        try:
            other = self.parent.convert(other)
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
        return self.__class__(functools.reduce(ZZ.lcm, (ZZ.convert(_.denominator)
                                                        for _ in self.rep.to_dense()), ZZ.one))


class ComplexAlgebraicElement(AlgebraicElement):
    """Elements of complex algebraic numbers field. """

    @property
    def real(self):
        """Returns real part of ``self``. """
        return self.domain.domain.convert(self.rep.coeff(1)) if self else self.domain.domain.zero

    @property
    def imag(self):
        """Returns imaginary part of ``self``. """
        return self.domain.domain.convert((self - self.real)/self.parent.unit)

    def conjugate(self):
        """Returns the complex conjugate of ``self``. """
        return self.parent.one*self.real - self.parent.unit*self.imag


@functools.total_ordering
class RealAlgebraicElement(ComplexAlgebraicElement):
    """Elements of real algebraic numbers field. """

    def __abs__(self):
        return self if self >= 0 else -self

    def __lt__(self, other):
        from ..polys.rootisolation import dup_count_real_roots

        try:
            other = self.parent.convert(other)
        except CoercionFailed:
            return NotImplemented

        parent = self.parent
        dom = parent.domain

        if parent._ext_root is None:
            parent._ext_root = parent._compute_ext_root(parent.ext,
                                                        parent.minpoly)
        coeff, root = parent._ext_root

        rep = dmp_compose((self - other).rep.to_dense(),
                          (parent.unit.rep*coeff).to_dense(), 0, dom)

        while dup_count_real_roots(rep, dom, root.interval.a, root.interval.b):
            root.refine()

        self.parent._ext_root = (coeff, root)
        return dmp_eval_in(rep, root.interval.center, 0, 0, dom) < 0

    @property
    def real(self):
        """Returns real part of ``self``. """
        return self

    @property
    def imag(self):
        """Returns imaginary part of ``self``. """
        return self.parent.zero
