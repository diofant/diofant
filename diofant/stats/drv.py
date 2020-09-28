from ..concrete import Sum, summation
from ..core import Dummy, Expr, Lambda, cacheit, symbols, sympify
from ..functions import Piecewise
from ..sets import Integers
from .rv import NamedArgsMixin, SingleDomain, SinglePSpace


class SingleDiscreteDistribution(Expr, NamedArgsMixin):
    """Discrete distribution of a single variable.

    Serves as superclass for PoissonDistribution etc....

    Provides methods for pdf, cdf, and sampling

    See Also:
        diofant.stats.crv_types.*

    """

    set = Integers

    def __new__(cls, *args):
        args = list(map(sympify, args))
        return Expr.__new__(cls, *args)

    @cacheit
    def compute_cdf(self, **kwargs):
        """Compute the CDF from the PDF.

        Returns a Lambda

        """
        x, z = symbols('x, z', integer=True, finite=True, cls=Dummy)
        left_bound = self.set.inf

        # CDF is integral of PDF from left bound to z
        pdf = self.pdf(x)
        cdf = summation(pdf, (x, left_bound, z), **kwargs)
        # CDF Ensure that CDF left of left_bound is zero
        cdf = Piecewise((cdf, z >= left_bound), (0, True))
        return Lambda(z, cdf)

    def cdf(self, x, **kwargs):
        """Cumulative density function."""
        return self.compute_cdf(**kwargs)(x)

    def expectation(self, expr, var, evaluate=True, **kwargs):
        """Expectation of expression over distribution."""
        # TODO: support discrete sets with non integer stepsizes
        if evaluate:
            return summation(expr * self.pdf(var),
                             (var, self.set.inf, self.set.sup), **kwargs)
        else:
            return Sum(expr * self.pdf(var),
                       (var, self.set.inf, self.set.sup), **kwargs)

    def __call__(self, *args):
        return self.pdf(*args)


class SingleDiscreteDomain(SingleDomain):
    """Base class for a discrete domain."""


class SingleDiscretePSpace(SinglePSpace):
    """Discrete probability space over a single univariate variable."""

    @property
    def set(self):
        return self.distribution.set

    @property
    def domain(self):
        return SingleDiscreteDomain(self.symbol, self.set)

    def integrate(self, expr, rvs=None, **kwargs):
        rvs = rvs or (self.value,)
        expr = expr.xreplace({rv: rv.symbol for rv in rvs})
        x = self.value.symbol
        return self.distribution.expectation(expr, x, evaluate=False, **kwargs)
