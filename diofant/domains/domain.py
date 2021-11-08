"""Implementation of :class:`Domain` class."""

import abc
import inspect

from ..core import Expr
from ..core.compatibility import HAS_GMPY
from ..polys.orderings import lex
from ..polys.polyerrors import CoercionFailed, UnificationFailed
from ..polys.polyutils import _unify_gens
from ..printing.defaults import DefaultPrinting
from .domainelement import DomainElement


class Domain(DefaultPrinting, abc.ABC):
    """Represents an abstract domain."""

    is_Ring = False
    is_Field = False

    has_assoc_Ring = False

    is_FiniteField = False
    is_IntegerRing = False
    is_RationalField = False
    is_RealField = False
    is_ComplexField = False
    is_AlgebraicField = False
    is_RealAlgebraicField = False
    is_ComplexAlgebraicField = False
    is_PolynomialRing = False
    is_FractionField = False
    is_ExpressionDomain = False

    is_Exact = True
    is_Numerical = False

    def __hash__(self):
        return hash((self.__class__.__name__, self.dtype))

    def __call__(self, *args):
        """Construct an element of ``self`` domain from ``args``."""
        return self.dtype(*args)

    def __getstate__(self):
        return {}

    @abc.abstractmethod
    def from_expr(self, expr):
        """Convert Diofant's expression ``expr`` to ``dtype``."""
        raise NotImplementedError

    @abc.abstractmethod
    def to_expr(self, element):
        """Convert domain ``element`` to Diofant expression."""
        raise NotImplementedError

    def convert_from(self, element, base):
        """Convert ``element`` to ``self.dtype`` given the base domain."""
        for superclass in inspect.getmro(base.__class__):
            method = '_from_' + superclass.__name__

            convert = getattr(self, method, None)

            if convert:
                result = convert(element, base)

                if result is not None:
                    return result

        raise CoercionFailed(f"can't convert {element} of type {type(element)} "
                             f'from {base} to {self}')

    def convert(self, element, base=None):
        """Convert ``element`` to ``self.dtype``."""
        if base is not None:
            return self.convert_from(element, base)

        if isinstance(element, self.dtype):
            return element

        from . import ComplexField, PythonRational, RealField
        from .expressiondomain import ExpressionDomain
        from .integerring import GMPYIntegerRing, PythonIntegerRing
        from .rationalfield import GMPYRationalField, PythonRationalField

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

        if isinstance(element, ExpressionDomain.Expression):
            return self.convert_from(element, ExpressionDomain())

        if isinstance(element, Expr):
            try:
                return self.from_expr(element)
            except (TypeError, ValueError):
                pass

        raise CoercionFailed(f"can't convert {element} of type {type(element)} to {self}")

    def __contains__(self, a):
        """Check if ``a`` belongs to this domain."""
        try:
            self.convert(a)
            return True
        except CoercionFailed:
            return False

    def _from_PolynomialRing(self, a, K0):
        if a.is_ground:
            return self.convert(a.LC, K0.domain)

    def _from_FractionField(self, a, K0):
        if a.numerator.is_ground and a.denominator == 1:
            return self.convert(a.numerator.LC, K0.domain.ring)

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
        from .compositedomain import CompositeDomain

        if symbols:
            if any(isinstance(d, CompositeDomain) and (set(d.symbols) & set(symbols))
                   for d in [self, K1]):
                raise UnificationFailed(f"Can't unify {self} with {K1}, "
                                        f'given {symbols} generators')

            return self.unify(K1)

        if self == K1:
            return self

        if self.is_ExpressionDomain:
            return self
        if K1.is_ExpressionDomain:
            return K1

        if any(isinstance(d, CompositeDomain) for d in (self, K1)):
            if isinstance(self, CompositeDomain):
                self_ground = self.domain
                self_symbols = self.symbols
                order = self.order
            else:
                self_ground = self
                self_symbols = ()
                order = K1.order

            if isinstance(K1, CompositeDomain):
                K1_ground = K1.domain
                K1_symbols = K1.symbols
            else:
                K1_ground = K1
                K1_symbols = ()

            domain = self_ground.unify(K1_ground)
            symbols = _unify_gens(self_symbols, K1_symbols)

            if ((self.is_FractionField and K1.is_PolynomialRing or
                 K1.is_FractionField and self.is_PolynomialRing) and
                    (not self_ground.is_Field or not K1_ground.is_Field) and domain.has_assoc_Ring):
                domain = domain.ring

            if isinstance(self, CompositeDomain) and (not isinstance(K1, CompositeDomain) or self.is_FractionField or K1.is_PolynomialRing):
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
        """Returns ``True`` if two domains are equivalent."""
        return isinstance(other, Domain) and self.dtype == other.dtype

    def get_exact(self):
        return self

    def poly_ring(self, *symbols, **kwargs):
        """Returns a polynomial ring, i.e. `K[X]`."""
        from ..polys import PolynomialRing
        return PolynomialRing(self, symbols, kwargs.get('order', lex))

    def frac_field(self, *symbols, **kwargs):
        """Returns a fraction field, i.e. `K(X)`."""
        from ..polys import FractionField
        return FractionField(self, symbols, kwargs.get('order', lex))
