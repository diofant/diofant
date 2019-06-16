"""Implementation of :class:`AlgebraicField` class. """

import functools

from ..core import I, Integer, cacheit, sympify
from ..core.sympify import CantSympify
from ..polys.densetools import dmp_compose, dmp_diff_in, dmp_eval_in
from ..polys.polyerrors import CoercionFailed, DomainError, NotAlgebraic
from .characteristiczero import CharacteristicZero
from .field import Field
from .quotientring import QuotientRingElement
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
        if is_real is not False:
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
        obj.mod = mod = rep_ring.from_dense(minpoly.rep.to_dense())

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
            return self([a])
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


class AlgebraicElement(QuotientRingElement, CantSympify):
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

    def to_dict(self):
        """Convert ``self`` to a dict representation with native coefficients. """
        return dict(self.rep)

    def LC(self):
        """Returns the leading coefficient of ``self``. """
        return self.rep.LC

    @property
    def is_ground(self):
        """Returns ``True`` if ``self`` is an element of the ground domain. """
        return self.rep.is_ground

    @property
    def numerator(self):
        return self*self.denominator

    @property
    def denominator(self):
        from . import ZZ
        return ZZ.convert(self.rep.content().denominator)


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

    @cacheit
    def __lt__(self, other):
        from ..polys.rootisolation import dup_count_real_roots

        try:
            other = self.parent.convert(other)
        except CoercionFailed:
            return NotImplemented

        dom = self.parent.domain
        coeff, root = self.parent._ext_root

        rep = dmp_compose((self - other).rep.to_dense(),
                          (self.parent.unit.rep*coeff).to_dense(), 0, dom)

        while dup_count_real_roots(rep, dom, root.interval.a, root.interval.b):
            root.refine()

        self.parent._ext_root = coeff, root
        return dmp_eval_in(rep, root.interval.center, 0, 0, dom) < 0

    @cacheit
    def __int__(self):
        from ..polys.rootisolation import dup_count_real_roots

        dom = self.parent.domain
        coeff, root = self.parent._ext_root

        rep = dmp_compose(self.rep.to_dense(),
                          (self.parent.unit.rep*coeff).to_dense(), 0, dom)
        df = dmp_diff_in(rep, 1, 0, 0, dom)

        while (dup_count_real_roots(df, dom, root.interval.a, root.interval.b) or
               int(dmp_eval_in(rep, root.interval.b, 0, 0, dom)) !=
               int(dmp_eval_in(rep, root.interval.a, 0, 0, dom))):
            root.refine()

        self.parent._ext_root = coeff, root
        return int(dmp_eval_in(rep, root.interval.a, 0, 0, dom))

    @property
    def real(self):
        """Returns real part of ``self``. """
        return self

    @property
    def imag(self):
        """Returns imaginary part of ``self``. """
        return self.parent.zero
