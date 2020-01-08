"""Definitions of common exceptions for :mod:`~diofant.polys` module."""


__all__ = ('BasePolynomialError', 'ExactQuotientFailed',
           'PolynomialDivisionFailed', 'OperationNotSupported',
           'HeuristicGCDFailed', 'HomomorphismFailed',
           'IsomorphismFailed', 'ExtraneousFactors',
           'EvaluationFailed', 'RefinementFailed', 'CoercionFailed',
           'NotInvertible', 'NotReversible', 'NotAlgebraic',
           'DomainError', 'PolynomialError', 'UnificationFailed',
           'GeneratorsError', 'GeneratorsNeeded', 'ComputationFailed',
           'UnivariatePolynomialError', 'MultivariatePolynomialError',
           'PolificationFailed', 'OptionError', 'FlagError')


class BasePolynomialError(Exception):
    """Base class for polynomial related exceptions."""

    def new(self, *args):
        raise NotImplementedError("abstract base class")


class ExactQuotientFailed(BasePolynomialError):
    """Raised when exact quotient is failed."""

    def __init__(self, f, g, dom=None):
        self.f, self.g, self.domain = f, g, dom

    def __str__(self):
        from ..printing import sstr

        if self.domain is None:
            return "%s does not divide %s" % (sstr(self.g), sstr(self.f))
        else:
            return "%s does not divide %s in %s" % (sstr(self.g), sstr(self.f), sstr(self.domain))

    def new(self, f, g):
        return self.__class__(f, g, self.domain)


class PolynomialDivisionFailed(BasePolynomialError):
    """Raised when polynomial division is failed."""

    def __init__(self, f, g, domain):
        self.f = f
        self.g = g
        self.domain = domain

    def __str__(self):
        if self.domain.is_SymbolicDomain:
            msg = "You may want to use a different simplification algorithm. Note " \
                  "that in general it's not possible to guarantee to detect zero "  \
                  "in this domain."
        elif not self.domain.is_Exact:
            msg = "Your working precision or tolerance of computations may be set " \
                  "improperly. Adjust those parameters of the coefficient domain "  \
                  "and try again."
        else:
            msg = "Zero detection is guaranteed in this coefficient domain. This "  \
                  "may indicate a bug in Diofant or the domain is user defined and "  \
                  "doesn't implement zero detection properly."

        return "couldn't reduce degree in a polynomial division algorithm when "    \
               "dividing %s by %s. This can happen when it's not possible to "      \
               "detect zero in the coefficient domain. The domain of computation "  \
               "is %s. %s" % (self.f, self.g, self.domain, msg)


class OperationNotSupported(BasePolynomialError):
    """Raised when an operation is not supported."""

    def __init__(self, poly, func):
        self.poly = poly
        self.func = func

    def __str__(self):
        return "`%s` operation not supported by %s representation" % (self.func, self.poly.rep.__class__.__name__)


class HeuristicGCDFailed(BasePolynomialError):
    """Raised when a heuristic GCD is failed."""

    pass


class ModularGCDFailed(BasePolynomialError):
    """Raised when a modular GCD is failed."""

    pass


class HomomorphismFailed(BasePolynomialError):
    """Raised when a homomorphism is failed."""

    pass


class IsomorphismFailed(BasePolynomialError):
    """Raised when an isomprphism is failed."""

    pass


class ExtraneousFactors(BasePolynomialError):
    """Raised when there are extraneous factors."""

    pass


class EvaluationFailed(BasePolynomialError):
    """Raised when a polynomial evaluation is failed."""

    pass


class RefinementFailed(BasePolynomialError):
    """Raised when a root refinement is failed."""

    pass


class CoercionFailed(BasePolynomialError):
    """Raised when a coercion is failed."""

    pass


class NotInvertible(BasePolynomialError):
    """Raised when a element is not invertible."""

    pass


class NotReversible(BasePolynomialError):
    """Raised when a element is not reversible."""

    pass


class NotAlgebraic(BasePolynomialError):
    """Raised when a non algebraic element occurred."""

    pass


class DomainError(BasePolynomialError):
    """Generic domain error."""

    pass


class PolynomialError(BasePolynomialError):
    """Generic polynomial error."""

    pass


class UnificationFailed(BasePolynomialError):
    """Raised when domains unification failed."""

    pass


class GeneratorsError(BasePolynomialError):
    """Raised when polynomial generators are unsuitable."""

    pass


class GeneratorsNeeded(GeneratorsError):
    """Raised when more generators needed."""

    pass


class ComputationFailed(BasePolynomialError):
    """Raised when polynomial computation failed."""

    def __init__(self, func, nargs, exc):
        self.func = func
        self.nargs = nargs
        self.exc = exc

    def __str__(self):
        return "%s(%s) failed without generators" % (self.func, ', '.join(map(str, self.exc.exprs[:self.nargs])))


class UnivariatePolynomialError(PolynomialError):
    """Generic univariate polynomial error."""

    pass


class MultivariatePolynomialError(PolynomialError):
    """Generic multivariate polynomial error."""

    pass


class PolificationFailed(PolynomialError):
    """Raised if polunomial construction is failed."""

    def __init__(self, opt, origs, exprs, seq=False):
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
            return "can't construct a polynomial from %s" % str(self.orig)
        else:
            return "can't construct polynomials from %s" % ', '.join(map(str, self.origs))


class OptionError(BasePolynomialError):
    """Generic option error."""

    pass


class FlagError(OptionError):
    """Generic flag error."""

    pass
