"""Implementation of :class:`AlgebraicField` class."""

from __future__ import annotations

import functools

from ..core import I, cacheit
from ..core.sympify import CantSympify, sympify
from ..polys.polyerrors import CoercionFailed, DomainError, NotAlgebraic
from .characteristiczero import CharacteristicZero
from .field import Field
from .quotientring import QuotientRingElement
from .rationalfield import RationalField
from .simpledomain import SimpleDomain


class AlgebraicField(CharacteristicZero, SimpleDomain, Field):
    """A class for representing algebraic number fields."""

    is_AlgebraicField = True
    is_Numerical = True

    def __new__(cls, dom, *ext):
        if not (dom.is_RationalField or dom.is_AlgebraicField):
            raise DomainError('ground domain must be a rational '
                              'or an algebraic field')

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

        rep_ring = dom.inject(obj.ext)
        obj.mod = mod = rep_ring.from_dict(minpoly.rep)

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
                             {'mod': mod, 'domain': rep_ring, '_parent': obj})
            _algebraic_numbers_cache[(obj.domain, obj.ext)] = obj.dtype

        obj.unit = obj.dtype([dom(0), dom(1)])

        obj.zero = obj.dtype([dom(0)])
        obj.one = obj.dtype([dom(1)])

        obj.rep = str(obj.domain) + '<' + str(obj.ext) + '>'

        return obj

    def __getnewargs_ex__(self):
        return (self.domain, self.ext), {}

    def __hash__(self):
        return hash((self.__class__.__name__, self.domain, self.ext))

    def __eq__(self, other):
        return isinstance(other, AlgebraicField) and self.domain == other.domain and self.ext == other.ext

    def algebraic_field(self, *extension):
        r"""Returns an algebraic field, i.e. `\mathbb{Q}(\alpha, \ldots)`."""
        return AlgebraicField(self, *extension)

    def to_expr(self, element):
        rep = element.rep
        return rep.ring.to_expr(rep)

    def from_expr(self, expr):
        from ..polys import primitive_element

        if expr in self.domain:
            return self([expr])

        try:
            _, (c,), (rep,) = primitive_element([expr], domain=self.domain)
        except NotAlgebraic:
            raise CoercionFailed(f'{expr} is not an algebraic number')

        K0 = self.domain.algebraic_field(c*expr)
        return self.convert(K0(rep), K0)

    def _from_PythonIntegerRing(self, a, K0):
        return self([self.domain.convert(a, K0)])
    _from_PythonRationalField = _from_PythonIntegerRing
    _from_GMPYIntegerRing = _from_PythonIntegerRing
    _from_GMPYRationalField = _from_PythonIntegerRing
    _from_RealField = _from_PythonIntegerRing

    def _from_ComplexField(self, a, K0):
        if self.ext == I:
            return self.from_expr(K0.to_expr(a))

    def _from_AlgebraicField(self, a, K0):
        if K0 == self.domain:
            return self([a])
        elif self == K0.domain and len(a.rep) <= 1:
            return a.rep[1] if a else self.zero

        from ..polys import field_isomorphism

        coeffs = field_isomorphism(K0, self)

        if coeffs is not None:
            if K0.domain == self.domain:
                return self(a.rep.compose(0, a.rep.ring.from_list(coeffs)))
            else:
                return self.from_expr(K0.to_expr(a))
        else:
            raise CoercionFailed(f'{K0} is not in a subfield of {self}')

    def _from_ExpressionDomain(self, a, K0):
        return self.from_expr(K0.to_expr(a))

    @property
    def ring(self):
        raise NotImplementedError(f'ring of integers of {self} is not implemented yet')

    def is_normal(self, a):
        return self.domain.is_normal(a.rep.LC)

    @staticmethod
    def _compute_ext_root(ext, minpoly):
        from ..polys import minimal_polynomial

        for r in minpoly.all_roots(radicals=False):  # pragma: no branch
            if not minimal_polynomial(ext - r)(0):
                return r.as_content_primitive()


_algebraic_numbers_cache: dict[tuple, AlgebraicField] = {}


class ComplexAlgebraicField(AlgebraicField):
    """A class for representing complex algebraic number fields."""

    is_ComplexAlgebraicField = True


class RealAlgebraicField(ComplexAlgebraicField):
    """A class for representing real algebraic number fields."""

    is_RealAlgebraicField = True

    def is_normal(self, a):
        return a >= 0


class AlgebraicElement(QuotientRingElement, CantSympify):
    """Dense Algebraic Number Polynomials over a field."""

    def __init__(self, rep):
        dom = self.domain

        if isinstance(rep, dict):
            rep = dom.from_dict(rep)
        else:
            if type(rep) is not list:
                rep = [dom.domain.convert(rep)]
            else:
                rep = [dom.domain.convert(_) for _ in rep]
            rep = dom.from_list(rep)

        self.rep = rep % self.mod

    def to_dict(self):
        """Convert ``self`` to a dict representation with native coefficients."""
        return dict(self.rep)

    @property
    def is_ground(self):
        """Returns ``True`` if ``self`` is an element of the ground domain."""
        return self.rep.is_ground

    @property
    def numerator(self):
        return self*self.denominator

    @property
    def denominator(self):
        from . import ZZ
        return ZZ.convert(self.rep.content().denominator)


class ComplexAlgebraicElement(AlgebraicElement):
    """Elements of complex algebraic numbers field."""

    @property
    def real(self):
        """Returns real part of ``self``."""
        return self.domain.domain.convert(self.rep[1]) if self else self.domain.domain.zero

    @property
    def imag(self):
        """Returns imaginary part of ``self``."""
        return self.domain.domain.convert((self - self.real)/self.parent.unit)

    def conjugate(self):
        """Returns the complex conjugate of ``self``."""
        return self.parent.one*self.real - self.parent.unit*self.imag


@functools.total_ordering
class RealAlgebraicElement(ComplexAlgebraicElement):
    """Elements of real algebraic numbers field."""

    def __abs__(self):
        return self if self >= 0 else -self

    @cacheit
    def __lt__(self, other):
        try:
            other = self.parent.convert(other)
        except CoercionFailed:
            return NotImplemented

        coeff, root = self.parent._ext_root

        ring = self.rep.ring
        rep = (self - other).rep.compose(0, self.parent.unit.rep*coeff)

        while ring._count_real_roots(rep, root.interval.a, root.interval.b):
            root.refine()

        self.parent._ext_root = coeff, root
        return rep(root.interval.center) < 0

    @cacheit
    def __int__(self):
        coeff, root = self.parent._ext_root

        ring = self.rep.ring
        rep = self.rep.compose(0, self.parent.unit.rep*coeff)
        df = rep.diff()

        while (ring._count_real_roots(df, root.interval.a, root.interval.b) or
               int(rep(root.interval.b)) != int(rep(root.interval.a))):
            root.refine()

        self.parent._ext_root = coeff, root
        return int(rep(root.interval.a))

    @property
    def real(self):
        """Returns real part of ``self``."""
        return self

    @property
    def imag(self):
        """Returns imaginary part of ``self``."""
        return self.parent.zero
