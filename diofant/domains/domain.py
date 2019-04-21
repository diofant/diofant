"""Implementation of :class:`Domain` class. """

import abc
import inspect

from ..core import Expr
from ..core.compatibility import HAS_GMPY
from ..polys.orderings import lex
from ..polys.polyerrors import CoercionFailed, UnificationFailed
from ..polys.polyutils import _unify_gens
from ..printing.defaults import DefaultPrinting
from .domainelement import DomainElement


__all__ = 'Domain',


class Domain(DefaultPrinting, abc.ABC):
    """Represents an abstract domain. """

    dtype = None
    zero = None
    one = None

    is_Ring = False
    is_Field = False

    has_assoc_Ring = False
    has_assoc_Field = False

    is_FiniteField = is_FF = False
    is_IntegerRing = is_ZZ = False
    is_RationalField = is_QQ = False
    is_RealField = is_RR = False
    is_ComplexField = is_CC = False
    is_AlgebraicField = is_Algebraic = False
    is_RealAlgebraicField = False
    is_ComplexAlgebraicField = False
    is_PolynomialRing = is_Poly = False
    is_FractionField = is_Frac = False
    is_SymbolicDomain = is_EX = False

    is_Exact = True
    is_Numerical = False

    is_Simple = False
    is_Composite = False

    rep = None

    def __hash__(self):
        return hash((self.__class__.__name__, self.dtype))

    def __call__(self, *args):
        """Construct an element of ``self`` domain from ``args``. """
        return self.dtype(*args)

    @abc.abstractmethod
    def from_expr(self, element):
        """Convert Diofant's expression to ``dtype``. """
        raise NotImplementedError

    @abc.abstractmethod
    def to_expr(self, element):
        """Convert ``element`` to Diofant expression. """
        raise NotImplementedError

    def convert_from(self, element, base):
        """Convert ``element`` to ``self.dtype`` given the base domain. """
        for superclass in inspect.getmro(base.__class__):
            method = "_from_" + superclass.__name__

            convert = getattr(self, method, None)

            if convert:
                result = convert(element, base)

                if result is not None:
                    return result

        raise CoercionFailed("can't convert %s of type %s from %s "
                             "to %s" % (element, type(element), base, self))

    def convert(self, element, base=None):
        """Convert ``element`` to ``self.dtype``. """
        if base is not None:
            return self.convert_from(element, base)

        if isinstance(element, self.dtype):
            return element

        from .integerring import PythonIntegerRing, GMPYIntegerRing
        from .rationalfield import PythonRationalField, GMPYRationalField
        from . import RealField, ComplexField, PythonRational

        if isinstance(element, int):
            return self.convert_from(element, PythonIntegerRing())

        if isinstance(element, PythonRational):
            return self.convert_from(element, PythonRationalField())

        if HAS_GMPY:
            integers = GMPYIntegerRing()
            if isinstance(element, integers.dtype):
                return self.convert_from(element, integers)

            rationals = GMPYRationalField()
            if isinstance(element, rationals.dtype):
                return self.convert_from(element, rationals)

        if isinstance(element, float):
            parent = RealField(tol=False)
            return self.convert_from(parent(element), parent)

        if isinstance(element, complex):
            parent = ComplexField(tol=False)
            return self.convert_from(parent(element), parent)

        if isinstance(element, DomainElement):
            return self.convert_from(element, element.parent)

        if isinstance(element, Expr):
            try:
                return self.from_expr(element)
            except (TypeError, ValueError):
                pass

        raise CoercionFailed("can't convert %s of type %s to %s" % (element, type(element), self))

    def __contains__(self, a):
        """Check if ``a`` belongs to this domain. """
        try:
            self.convert(a)
            return True
        except CoercionFailed:
            return False

    def _from_PolynomialRing(self, a, K0):
        if a.is_ground:
            return self.convert(a.LC, K0.domain)

    def unify(self, K1, symbols=()):
        """
        Construct a minimal domain that contains elements of ``self`` and ``K1``.

        Known domains (from smallest to largest):

        - ``GF(p)``
        - ``ZZ``
        - ``QQ``
        - ``RR(prec, tol)``
        - ``CC(prec, tol)``
        - ``ALG(a, b, c)``
        - ``K[x, y, z]``
        - ``K(x, y, z)``
        - ``EX``
        """
        if symbols:
            if any(d.is_Composite and (set(d.symbols) & set(symbols))
                   for d in [self, K1]):
                raise UnificationFailed("Can't unify %s with %s, given %s"
                                        " generators" % (self, K1, tuple(symbols)))

            return self.unify(K1)

        if self == K1:
            return self

        if self.is_SymbolicDomain:
            return self
        if K1.is_SymbolicDomain:
            return K1

        if self.is_Composite or K1.is_Composite:
            self_ground = self.domain if self.is_Composite else self
            K1_ground = K1.domain if K1.is_Composite else K1

            self_symbols = self.symbols if self.is_Composite else ()
            K1_symbols = K1.symbols if K1.is_Composite else ()

            domain = self_ground.unify(K1_ground)
            symbols = _unify_gens(self_symbols, K1_symbols)
            order = self.order if self.is_Composite else K1.order

            if ((self.is_FractionField and K1.is_PolynomialRing or
                 K1.is_FractionField and self.is_PolynomialRing) and
                    (not self_ground.is_Field or not K1_ground.is_Field) and domain.has_assoc_Ring):
                domain = domain.ring

            if self.is_Composite and (not K1.is_Composite or self.is_FractionField or K1.is_PolynomialRing):
                cls = self.__class__
            else:
                cls = K1.__class__

            return cls(domain, symbols, order)

        def mkinexact(cls, K0, K1):
            prec = max(K0.precision, K1.precision)
            tol = max(K0.tolerance, K1.tolerance)
            return cls(prec=prec, tol=tol)

        if self.is_ComplexField and K1.is_ComplexField:
            return mkinexact(self.__class__, self, K1)
        if self.is_ComplexField and K1.is_RealField:
            return mkinexact(self.__class__, self, K1)
        if self.is_RealField and K1.is_ComplexField:
            return mkinexact(K1.__class__, K1, self)
        if self.is_RealField and K1.is_RealField:
            return mkinexact(self.__class__, self, K1)
        if self.is_ComplexField or self.is_RealField:
            return self
        if K1.is_ComplexField or K1.is_RealField:
            return K1

        if self.is_AlgebraicField and K1.is_AlgebraicField:
            return self.__class__(self.domain.unify(K1.domain), *_unify_gens(self.gens, K1.gens))
        elif self.is_AlgebraicField:
            return self
        elif K1.is_AlgebraicField:
            return K1

        if self.is_RationalField:
            return self
        if K1.is_RationalField:
            return K1

        if self.is_FiniteField and self.domain == K1:
            return self
        if K1.is_FiniteField and K1.domain == self:
            return K1

        raise NotImplementedError

    def __eq__(self, other):
        """Returns ``True`` if two domains are equivalent. """
        return isinstance(other, Domain) and self.dtype == other.dtype

    def get_exact(self):
        """Returns an exact domain associated with ``self``. """
        return self

    def poly_ring(self, *symbols, **kwargs):
        """Returns a polynomial ring, i.e. `K[X]`. """
        from ..polys import PolynomialRing
        return PolynomialRing(self, symbols, kwargs.get("order", lex))

    def frac_field(self, *symbols, **kwargs):
        """Returns a fraction field, i.e. `K(X)`. """
        from ..polys import FractionField
        return FractionField(self, symbols, kwargs.get("order", lex))

    def is_positive(self, a):
        """Returns True if ``a`` is positive. """
        return a > 0

    def is_negative(self, a):
        """Returns True if ``a`` is negative. """
        return a < 0
