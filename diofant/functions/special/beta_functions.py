from ...core.function import ArgumentIndexError, Function
from .gamma_functions import digamma, gamma


###############################################################################
# ########################## COMPLETE BETA  FUNCTION ######################## #
###############################################################################


class beta(Function):
    r"""
    The beta integral is called the Eulerian integral of the first kind by
    Legendre:

    .. math::
        \mathrm{B}(x,y) := \int^{1}_{0} t^{x-1} (1-t)^{y-1} \mathrm{d}t.

    Beta function or Euler's first integral is closely associated with gamma function.
    The Beta function often used in probability theory and mathematical statistics.
    It satisfies properties like:

    .. math::
        \mathrm{B}(a,1) = \frac{1}{a} \\
        \mathrm{B}(a,b) = \mathrm{B}(b,a)  \\
        \mathrm{B}(a,b) = \frac{\Gamma(a) \Gamma(b)}{\Gamma(a+b)}

    Therefore for integral values of a and b:

    .. math::
        \mathrm{B} = \frac{(a-1)! (b-1)!}{(a+b-1)!}

    Examples
    ========

    The Beta function obeys the mirror symmetry:

    >>> conjugate(beta(x, y))
    beta(conjugate(x), conjugate(y))

    Differentiation with respect to both x and y is supported:

    >>> diff(beta(x, y), x)
    (polygamma(0, x) - polygamma(0, x + y))*beta(x, y)

    >>> diff(beta(x, y), y)
    (polygamma(0, y) - polygamma(0, x + y))*beta(x, y)

    We can numerically evaluate the gamma function to arbitrary precision
    on the whole complex plane:

    >>> beta(pi, pi).evalf(40)
    0.02671848900111377452242355235388489324562

    >>> beta(1 + I, 1 + I).evalf(20)
    -0.2112723729365330143 - 0.7655283165378005676*I

    See Also
    ========

    diofant.functions.special.gamma_functions.gamma: Gamma function.
    diofant.functions.special.gamma_functions.uppergamma: Upper incomplete gamma function.
    diofant.functions.special.gamma_functions.lowergamma: Lower incomplete gamma function.
    diofant.functions.special.gamma_functions.polygamma: Polygamma function.
    diofant.functions.special.gamma_functions.loggamma: Log Gamma function.
    diofant.functions.special.gamma_functions.digamma: Digamma function.
    diofant.functions.special.gamma_functions.trigamma: Trigamma function.

    References
    ==========

    * https://en.wikipedia.org/wiki/Beta_function
    * https://mathworld.wolfram.com/BetaFunction.html
    * https://dlmf.nist.gov/5.12

    """

    unbranched = True

    def fdiff(self, argindex=1):
        x, y = self.args
        if argindex == 1:
            return beta(x, y)*(digamma(x) - digamma(x + y))
        if argindex == 2:
            return beta(x, y)*(digamma(y) - digamma(x + y))
        raise ArgumentIndexError(self, argindex)

    @classmethod
    def eval(cls, x, y):
        pass

    def _eval_expand_func(self, **hints):
        x, y = self.args
        return gamma(x)*gamma(y) / gamma(x + y)

    def _eval_conjugate(self):
        return self.func(self.args[0].conjugate(), self.args[1].conjugate())
