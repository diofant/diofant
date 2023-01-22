"""Definitions of common exceptions for :mod:`~diofant.polys` module."""


class BasePolynomialError(Exception):
    """Base class for polynomial related exceptions."""


class ExactQuotientFailedError(BasePolynomialError):
    """Raised when exact quotient is failed."""

    def __init__(self, f, g, dom=None):
        """Initialize self."""
        super().__init__()

        self.f = f
        self.g = g
        self.domain = dom

    def __str__(self):
        from ..printing import sstr

        if self.domain is None:
            return f'{sstr(self.g)} does not divide {sstr(self.f)}'
        return f'{sstr(self.g)} does not divide {sstr(self.f)} in {sstr(self.domain)}'


class PolynomialDivisionFailedError(BasePolynomialError):
    """Raised when polynomial division is failed."""

    def __init__(self, f, g, domain):
        """Initialize self."""
        super().__init__()

        self.f = f
        self.g = g
        self.domain = domain

    def __str__(self):
        if self.domain.is_ExpressionDomain:
            msg = 'You may want to use a different simplification algorithm. Note ' \
                  "that in general it's not possible to guarantee to detect zero "  \
                  'in this domain.'
        elif not self.domain.is_Exact:
            msg = 'Your working precision or tolerance of computations may be set ' \
                  'improperly. Adjust those parameters of the coefficient domain '  \
                  'and try again.'
        else:
            msg = 'Zero detection is guaranteed in this coefficient domain. This '  \
                  'may indicate a bug in Diofant or the domain is user defined and '  \
                  "doesn't implement zero detection properly."

        return "couldn't reduce degree in a polynomial division algorithm when "    \
               f"dividing {self.f} by {self.g}. This can happen when it's not possible to "      \
               'detect zero in the coefficient domain. The domain of computation '  \
               f'is {self.domain}. {msg}'


class OperationNotSupportedError(BasePolynomialError):
    """Raised when an operation is not supported."""

    def __init__(self, poly, func):
        """Initialize self."""
        super().__init__()

        self.poly = poly
        self.func = func

    def __str__(self):
        return f'`{self.func}` operation not supported by {self.poly.rep.__class__.__name__} representation'


class HeuristicGCDFailedError(BasePolynomialError):
    """Raised when a heuristic GCD is failed."""


class ModularGCDFailedError(BasePolynomialError):
    """Raised when a modular GCD is failed."""


class UnluckyLeadingCoefficientError(BasePolynomialError):
    """Raised when there are unlucky LC."""


class HomomorphismFailedError(BasePolynomialError):
    """Raised when a homomorphism is failed."""


class IsomorphismFailedError(BasePolynomialError):
    """Raised when an isomprphism is failed."""


class ExtraneousFactorsError(BasePolynomialError):
    """Raised when there are extraneous factors."""


class EvaluationFailedError(BasePolynomialError):
    """Raised when a polynomial evaluation is failed."""


class RefinementFailedError(BasePolynomialError):
    """Raised when a root refinement is failed."""


class CoercionFailedError(BasePolynomialError):
    """Raised when a coercion is failed."""


class NotInvertibleError(BasePolynomialError):
    """Raised when a element is not invertible."""


class NotReversibleError(BasePolynomialError):
    """Raised when a element is not reversible."""


class NotAlgebraicError(BasePolynomialError):
    """Raised when a non algebraic element occurred."""


class DomainError(BasePolynomialError):
    """Generic domain error."""


class PolynomialError(BasePolynomialError):
    """Generic polynomial error."""


class UnificationFailedError(BasePolynomialError):
    """Raised when domains unification failed."""


class GeneratorsError(BasePolynomialError):
    """Raised when polynomial generators are unsuitable."""


class GeneratorsNeededError(GeneratorsError):
    """Raised when more generators needed."""


class ComputationFailedError(BasePolynomialError):
    """Raised when polynomial computation failed."""

    def __init__(self, func, nargs, exc):
        """Initialize self."""
        super().__init__()

        self.func = func
        self.nargs = nargs
        self.exc = exc

    def __str__(self):
        return f'{self.func}({", ".join(map(str, self.exc.exprs[:self.nargs]))}) failed without generators'


class UnivariatePolynomialError(PolynomialError):
    """Generic univariate polynomial error."""


class MultivariatePolynomialError(PolynomialError):
    """Generic multivariate polynomial error."""


class PolificationFailedError(PolynomialError):
    """Raised if polunomial construction is failed."""

    def __init__(self, opt, origs, exprs, seq=False):
        """Initialize self."""
        super().__init__()

        if not seq:
            self.orig = origs
            self.expr = exprs
            self.origs = [origs]
            self.exprs = [exprs]
        else:
            self.origs = origs
            self.exprs = exprs

        self.opt = opt
        self.seq = seq

    def __str__(self):
        if not self.seq:
            return f"can't construct a polynomial from {self.orig}"
        return f"can't construct polynomials from {', '.join(map(str, self.origs))}"


class OptionError(BasePolynomialError):
    """Generic option error."""


class FlagError(OptionError):
    """Generic flag error."""
