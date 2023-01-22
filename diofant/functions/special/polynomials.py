"""
This module mainly implements special orthogonal polynomials.

See also functions.combinatorial.numbers which contains some
combinatorial polynomials.

"""

from ...core import Dummy, Function, Integer, Rational, oo, pi, zoo
from ...core.function import ArgumentIndexError
from ...polys.orthopolys import (chebyshevt_poly, chebyshevu_poly,
                                 gegenbauer_poly, hermite_poly, jacobi_poly,
                                 laguerre_poly, legendre_poly)
from ..combinatorial.factorials import RisingFactorial, binomial, factorial
from ..elementary.complexes import re
from ..elementary.exponential import exp
from ..elementary.miscellaneous import sqrt
from ..elementary.trigonometric import cos
from .gamma_functions import gamma
from .hyper import hyper


_x = Dummy('dummy_for_special_polynomials')


class OrthogonalPolynomial(Function):
    """Base class for orthogonal polynomials."""

    @classmethod
    def _eval_at_order(cls, n, x):
        if n.is_integer and n >= 0:
            return cls._ortho_poly(int(n), _x).subs({_x: x})

    def _eval_conjugate(self):
        return self.func(self.args[0], self.args[1].conjugate())

############################################################################
# Jacobi polynomials
#


class jacobi(OrthogonalPolynomial):
    r"""
    Jacobi polynomial `P_n^{\left(\alpha, \beta\right)}(x)`

    jacobi(n, alpha, beta, x) gives the nth Jacobi polynomial
    in x, `P_n^{\left(\alpha, \beta\right)}(x)`.

    The Jacobi polynomials are orthogonal on `[-1, 1]` with respect
    to the weight `\left(1-x\right)^\alpha \left(1+x\right)^\beta`.

    Examples
    ========

    >>> jacobi(0, a, b, x)
    1
    >>> jacobi(1, a, b, x)
    a/2 - b/2 + x*(a/2 + b/2 + 1)

    >>> jacobi(n, a, b, x)
    jacobi(n, a, b, x)

    >>> jacobi(n, a, a, x)
    RisingFactorial(a + 1, n)*gegenbauer(n,
        a + 1/2, x)/RisingFactorial(2*a + 1, n)

    >>> jacobi(n, 0, 0, x)
    legendre(n, x)

    >>> jacobi(n, Rational(1, 2), Rational(1, 2), x)
    RisingFactorial(3/2, n)*chebyshevu(n, x)/factorial(n + 1)

    >>> jacobi(n, -Rational(1, 2), -Rational(1, 2), x)
    RisingFactorial(1/2, n)*chebyshevt(n, x)/factorial(n)

    >>> jacobi(n, a, b, -x)
    (-1)**n*jacobi(n, b, a, x)

    >>> jacobi(n, a, b, 0)
    2**(-n)*gamma(a + n + 1)*hyper((-b - n, -n), (a + 1,), -1)/(factorial(n)*gamma(a + 1))
    >>> jacobi(n, a, b, 1)
    RisingFactorial(a + 1, n)/factorial(n)

    >>> conjugate(jacobi(n, a, b, x))
    jacobi(n, conjugate(a), conjugate(b), conjugate(x))

    >>> diff(jacobi(n, a, b, x), x)
    (a + b + n + 1)*jacobi(n - 1, a + 1, b + 1, x)/2

    See Also
    ========

    gegenbauer,
    chebyshevt_root, chebyshevu, chebyshevu_root,
    legendre, assoc_legendre,
    hermite,
    laguerre, assoc_laguerre,
    diofant.polys.orthopolys.jacobi_poly,
    diofant.polys.orthopolys.gegenbauer_poly
    diofant.polys.orthopolys.chebyshevt_poly
    diofant.polys.orthopolys.chebyshevu_poly
    diofant.polys.orthopolys.hermite_poly
    diofant.polys.orthopolys.legendre_poly
    diofant.polys.orthopolys.laguerre_poly

    References
    ==========

    * https://en.wikipedia.org/wiki/Jacobi_polynomials
    * https://mathworld.wolfram.com/JacobiPolynomial.html
    * http://functions.wolfram.com/Polynomials/JacobiP/

    """

    @classmethod
    def eval(cls, n, a, b, x):
        # Simplify to other polynomials
        # P^{a, a}_n(x)
        if a == b:
            if a == -Rational(1, 2):
                return RisingFactorial(Rational(1, 2), n) / factorial(n) * chebyshevt(n, x)
            if a == 0:
                return legendre(n, x)
            if a == Rational(1, 2):
                return RisingFactorial(Rational(3, 2), n) / factorial(n + 1) * chebyshevu(n, x)
            return RisingFactorial(a + 1, n) / RisingFactorial(2*a + 1, n) * gegenbauer(n, a + Rational(1, 2), x)
        if b == -a:
            # P^{a, -a}_n(x)
            return gamma(n + a + 1) / gamma(n + 1) * (1 + x)**(a/2) / (1 - x)**(a/2) * assoc_legendre(n, -a, x)

        if not n.is_Number:
            # Symbolic result P^{a,b}_n(x)
            # P^{a,b}_n(-x)  --->  (-1)**n * P^{b,a}_n(-x)
            if x.could_extract_minus_sign():
                return (-1)**n * jacobi(n, b, a, -x)
            # We can evaluate for some special values of x
            if x == 0:
                return (2**(-n) * gamma(a + n + 1) / (gamma(a + 1) * factorial(n)) *
                        hyper([-b - n, -n], [a + 1], -1))
            if x == 1:
                return RisingFactorial(a + 1, n) / factorial(n)
            if x == oo:
                if n.is_positive:
                    # Make sure a+b+2*n \notin Z
                    if (a + b + 2*n).is_integer:
                        raise ValueError('a + b + 2*n should not be an integer.')
                    return RisingFactorial(a + b + n + 1, n) * oo
        else:
            # n is a given fixed integer, evaluate into polynomial
            return jacobi_poly(n, a, b, x)

    def fdiff(self, argindex=4):
        from ...concrete import Sum
        if argindex == 1:
            # Diff wrt n
            raise ArgumentIndexError(self, argindex)
        if argindex == 2:
            # Diff wrt a
            n, a, b, x = self.args
            k = Dummy('k')
            f1 = 1 / (a + b + n + k + 1)
            f2 = ((a + b + 2*k + 1) * RisingFactorial(b + k + 1, n - k) /
                  ((n - k) * RisingFactorial(a + b + k + 1, n - k)))
            return Sum(f1 * (jacobi(n, a, b, x) + f2*jacobi(k, a, b, x)), (k, 0, n - 1))
        if argindex == 3:
            # Diff wrt b
            n, a, b, x = self.args
            k = Dummy('k')
            f1 = 1 / (a + b + n + k + 1)
            f2 = (-1)**(n - k) * ((a + b + 2*k + 1) * RisingFactorial(a + k + 1, n - k) /
                                  ((n - k) * RisingFactorial(a + b + k + 1, n - k)))
            return Sum(f1 * (jacobi(n, a, b, x) + f2*jacobi(k, a, b, x)), (k, 0, n - 1))
        if argindex == 4:
            # Diff wrt x
            n, a, b, x = self.args
            return (a + b + n + 1) * jacobi(n - 1, a + 1, b + 1, x) / 2
        raise ArgumentIndexError(self, argindex)

    def _eval_conjugate(self):
        n, a, b, x = self.args
        return self.func(n, a.conjugate(), b.conjugate(), x.conjugate())


def jacobi_normalized(n, a, b, x):
    r"""
    Jacobi polynomial `P_n^{\left(\alpha, \beta\right)}(x)`

    jacobi_normalized(n, alpha, beta, x) gives the nth Jacobi polynomial
    in x, `P_n^{\left(\alpha, \beta\right)}(x)`.

    The Jacobi polynomials are orthogonal on `[-1, 1]` with respect
    to the weight `\left(1-x\right)^\alpha \left(1+x\right)^\beta`.

    This functions returns the polynomials normilzed:

    .. math::

        \int_{-1}^{1}
          P_m^{\left(\alpha, \beta\right)}(x)
          P_n^{\left(\alpha, \beta\right)}(x)
          (1-x)^{\alpha} (1+x)^{\beta} \mathrm{d}x
        = \delta_{m,n}

    Examples
    ========

    >>> jacobi_normalized(n, a, b, x)
    jacobi(n, a, b, x)/sqrt(2**(a + b + 1)*gamma(a + n + 1)*gamma(b + n + 1)/((a + b + 2*n + 1)*factorial(n)*gamma(a + b + n + 1)))

    See Also
    ========

    gegenbauer,
    chebyshevt_root, chebyshevu, chebyshevu_root,
    legendre, assoc_legendre,
    hermite,
    laguerre, assoc_laguerre,
    diofant.polys.orthopolys.jacobi_poly,
    diofant.polys.orthopolys.gegenbauer_poly
    diofant.polys.orthopolys.chebyshevt_poly
    diofant.polys.orthopolys.chebyshevu_poly
    diofant.polys.orthopolys.hermite_poly
    diofant.polys.orthopolys.legendre_poly
    diofant.polys.orthopolys.laguerre_poly

    References
    ==========

    * https://en.wikipedia.org/wiki/Jacobi_polynomials
    * https://mathworld.wolfram.com/JacobiPolynomial.html
    * http://functions.wolfram.com/Polynomials/JacobiP/

    """
    nfactor = (Integer(2)**(a + b + 1) * (gamma(n + a + 1) * gamma(n + b + 1))
               / (2*n + a + b + 1) / (factorial(n) * gamma(n + a + b + 1)))

    return jacobi(n, a, b, x) / sqrt(nfactor)


############################################################################
# Gegenbauer polynomials
#


class gegenbauer(OrthogonalPolynomial):
    r"""
    Gegenbauer polynomial `C_n^{\left(\alpha\right)}(x)`

    gegenbauer(n, alpha, x) gives the nth Gegenbauer polynomial
    in x, `C_n^{\left(\alpha\right)}(x)`.

    The Gegenbauer polynomials are orthogonal on `[-1, 1]` with
    respect to the weight `\left(1-x^2\right)^{\alpha-\frac{1}{2}}`.

    Examples
    ========

    >>> gegenbauer(0, a, x)
    1
    >>> gegenbauer(1, a, x)
    2*a*x
    >>> gegenbauer(2, a, x)
    -a + x**2*(2*a**2 + 2*a)
    >>> gegenbauer(3, a, x)
    x**3*(4*a**3/3 + 4*a**2 + 8*a/3) + x*(-2*a**2 - 2*a)

    >>> gegenbauer(n, a, x)
    gegenbauer(n, a, x)
    >>> gegenbauer(n, a, -x)
    (-1)**n*gegenbauer(n, a, x)

    >>> gegenbauer(n, a, 0)
    2**n*sqrt(pi)*gamma(a + n/2)/(gamma(a)*gamma(-n/2 + 1/2)*gamma(n + 1))
    >>> gegenbauer(n, a, 1)
    gamma(2*a + n)/(gamma(2*a)*gamma(n + 1))

    >>> conjugate(gegenbauer(n, a, x))
    gegenbauer(n, conjugate(a), conjugate(x))

    >>> diff(gegenbauer(n, a, x), x)
    2*a*gegenbauer(n - 1, a + 1, x)

    See Also
    ========

    jacobi,
    chebyshevt_root, chebyshevu, chebyshevu_root,
    legendre, assoc_legendre,
    hermite,
    laguerre, assoc_laguerre,
    diofant.polys.orthopolys.jacobi_poly
    diofant.polys.orthopolys.gegenbauer_poly
    diofant.polys.orthopolys.chebyshevt_poly
    diofant.polys.orthopolys.chebyshevu_poly
    diofant.polys.orthopolys.hermite_poly
    diofant.polys.orthopolys.legendre_poly
    diofant.polys.orthopolys.laguerre_poly

    References
    ==========

    * https://en.wikipedia.org/wiki/Gegenbauer_polynomials
    * https://mathworld.wolfram.com/GegenbauerPolynomial.html
    * http://functions.wolfram.com/Polynomials/GegenbauerC3/

    """

    @classmethod
    def eval(cls, n, a, x):
        # For negative n the polynomials vanish
        # See http://functions.wolfram.com/Polynomials/GegenbauerC3/03/01/03/0012/
        if n.is_negative:
            return Integer(0)

        # Some special values for fixed a
        if a == Rational(1, 2):
            return legendre(n, x)
        if a == 1:
            return chebyshevu(n, x)
        if a == -1:
            return Integer(0)

        if not n.is_Number:
            # Handle this before the general sign extraction rule
            if x == -1:
                if (re(a) - Rational(1, 2)).is_positive:
                    return zoo
                return

            # Symbolic result C^a_n(x)
            # C^a_n(-x)  --->  (-1)**n * C^a_n(x)
            if x.could_extract_minus_sign():
                return (-1)**n * gegenbauer(n, a, -x)
            # We can evaluate for some special values of x
            if x == 0:
                return (2**n * sqrt(pi) * gamma(a + n/2) /
                        (gamma((1 - n)/2) * gamma(n + 1) * gamma(a)))
            if x == 1:
                return gamma(2*a + n) / (gamma(2*a) * gamma(n + 1))
            if x == oo:
                if n.is_positive:
                    return RisingFactorial(a, n) * oo
        else:
            # n is a given fixed integer, evaluate into polynomial
            return gegenbauer_poly(n, a, x)

    def fdiff(self, argindex=3):
        from ...concrete import Sum
        if argindex == 1:
            # Diff wrt n
            raise ArgumentIndexError(self, argindex)
        if argindex == 2:
            # Diff wrt a
            n, a, x = self.args
            k = Dummy('k')
            factor1 = 2 * (1 + (-1)**(n - k)) * (k + a) / ((k +
                                                            n + 2*a) * (n - k))
            factor2 = 2*(k + 1) / ((k + 2*a) * (2*k + 2*a + 1)) + \
                2 / (k + n + 2*a)
            kern = factor1*gegenbauer(k, a, x) + factor2*gegenbauer(n, a, x)
            return Sum(kern, (k, 0, n - 1))
        if argindex == 3:
            # Diff wrt x
            n, a, x = self.args
            return 2*a*gegenbauer(n - 1, a + 1, x)
        raise ArgumentIndexError(self, argindex)

    def _eval_conjugate(self):
        n, a, x = self.args
        return self.func(n, a.conjugate(), x.conjugate())

############################################################################
# Chebyshev polynomials of first and second kind
#


class chebyshevt(OrthogonalPolynomial):
    r"""
    Chebyshev polynomial of the first kind, `T_n(x)`

    chebyshevt(n, x) gives the nth Chebyshev polynomial (of the first
    kind) in x, `T_n(x)`.

    The Chebyshev polynomials of the first kind are orthogonal on
    `[-1, 1]` with respect to the weight `\frac{1}{\sqrt{1-x^2}}`.

    Examples
    ========

    >>> chebyshevt(0, x)
    1
    >>> chebyshevt(1, x)
    x
    >>> chebyshevt(2, x)
    2*x**2 - 1

    >>> chebyshevt(n, x)
    chebyshevt(n, x)
    >>> chebyshevt(n, -x)
    (-1)**n*chebyshevt(n, x)
    >>> chebyshevt(-n, x)
    chebyshevt(n, x)

    >>> chebyshevt(n, 0)
    cos(pi*n/2)
    >>> chebyshevt(n, -1)
    (-1)**n

    >>> diff(chebyshevt(n, x), x)
    n*chebyshevu(n - 1, x)

    See Also
    ========

    jacobi, gegenbauer,
    chebyshevt_root, chebyshevu, chebyshevu_root,
    legendre, assoc_legendre,
    hermite,
    laguerre, assoc_laguerre,
    diofant.polys.orthopolys.jacobi_poly
    diofant.polys.orthopolys.gegenbauer_poly
    diofant.polys.orthopolys.chebyshevt_poly
    diofant.polys.orthopolys.chebyshevu_poly
    diofant.polys.orthopolys.hermite_poly
    diofant.polys.orthopolys.legendre_poly
    diofant.polys.orthopolys.laguerre_poly

    References
    ==========

    * https://en.wikipedia.org/wiki/Chebyshev_polynomial
    * https://mathworld.wolfram.com/ChebyshevPolynomialoftheFirstKind.html
    * https://mathworld.wolfram.com/ChebyshevPolynomialoftheSecondKind.html
    * http://functions.wolfram.com/Polynomials/ChebyshevT/
    * http://functions.wolfram.com/Polynomials/ChebyshevU/

    """

    _ortho_poly = staticmethod(chebyshevt_poly)

    @classmethod
    def eval(cls, n, x):
        if not n.is_Number:
            # Symbolic result T_n(x)
            # T_n(-x)  --->  (-1)**n * T_n(x)
            if x.could_extract_minus_sign():
                return (-1)**n * chebyshevt(n, -x)
            # T_{-n}(x)  --->  T_n(x)
            if n.could_extract_minus_sign():
                return chebyshevt(-n, x)
            # We can evaluate for some special values of x
            if x == 0:
                return cos(pi*n/2)
            if x == 1:
                return Integer(1)
            if x == oo:
                return oo
        else:
            # n is a given fixed integer, evaluate into polynomial
            if n.is_negative:
                # T_{-n}(x) == T_n(x)
                return cls._eval_at_order(-n, x)
            return cls._eval_at_order(n, x)

    def fdiff(self, argindex=2):
        if argindex == 2:
            # Diff wrt x
            n, x = self.args
            return n * chebyshevu(n - 1, x)
        # wrt n
        raise ArgumentIndexError(self, argindex)


class chebyshevu(OrthogonalPolynomial):
    r"""
    Chebyshev polynomial of the second kind, `U_n(x)`

    chebyshevu(n, x) gives the nth Chebyshev polynomial of the second
    kind in x, `U_n(x)`.

    The Chebyshev polynomials of the second kind are orthogonal on
    `[-1, 1]` with respect to the weight `\sqrt{1-x^2}`.

    Examples
    ========

    >>> chebyshevu(0, x)
    1
    >>> chebyshevu(1, x)
    2*x
    >>> chebyshevu(2, x)
    4*x**2 - 1

    >>> chebyshevu(n, x)
    chebyshevu(n, x)
    >>> chebyshevu(n, -x)
    (-1)**n*chebyshevu(n, x)
    >>> chebyshevu(-n, x)
    -chebyshevu(n - 2, x)

    >>> chebyshevu(n, 0)
    cos(pi*n/2)
    >>> chebyshevu(n, 1)
    n + 1

    >>> diff(chebyshevu(n, x), x)
    (-x*chebyshevu(n, x) + (n + 1)*chebyshevt(n + 1, x))/(x**2 - 1)

    See Also
    ========

    jacobi, gegenbauer,
    chebyshevt, chebyshevt_root, chebyshevu_root,
    legendre, assoc_legendre,
    hermite,
    laguerre, assoc_laguerre,
    diofant.polys.orthopolys.jacobi_poly
    diofant.polys.orthopolys.gegenbauer_poly
    diofant.polys.orthopolys.chebyshevt_poly
    diofant.polys.orthopolys.chebyshevu_poly
    diofant.polys.orthopolys.hermite_poly
    diofant.polys.orthopolys.legendre_poly
    diofant.polys.orthopolys.laguerre_poly

    References
    ==========

    * https://en.wikipedia.org/wiki/Chebyshev_polynomial
    * https://mathworld.wolfram.com/ChebyshevPolynomialoftheFirstKind.html
    * https://mathworld.wolfram.com/ChebyshevPolynomialoftheSecondKind.html
    * http://functions.wolfram.com/Polynomials/ChebyshevT/
    * http://functions.wolfram.com/Polynomials/ChebyshevU/

    """

    _ortho_poly = staticmethod(chebyshevu_poly)

    @classmethod
    def eval(cls, n, x):
        if not n.is_Number:
            # Symbolic result U_n(x)
            # U_n(-x)  --->  (-1)**n * U_n(x)
            if x.could_extract_minus_sign():
                return (-1)**n * chebyshevu(n, -x)
            # U_{-n}(x)  --->  -U_{n-2}(x)
            if n.could_extract_minus_sign():
                return -chebyshevu(-n - 2, x)
            # We can evaluate for some special values of x
            if x == 0:
                return cos(pi*n/2)
            if x == 1:
                return 1 + n
            if x == oo:
                return oo
        else:
            # n is a given fixed integer, evaluate into polynomial
            if n.is_negative:
                # U_{-n}(x)  --->  -U_{n-2}(x)
                if n == -1:
                    return Integer(0)
                return -cls._eval_at_order(-n - 2, x)
            return cls._eval_at_order(n, x)

    def fdiff(self, argindex=2):
        if argindex == 2:
            # Diff wrt x
            n, x = self.args
            return ((n + 1) * chebyshevt(n + 1, x) - x * chebyshevu(n, x)) / (x**2 - 1)
        # wrt n
        raise ArgumentIndexError(self, argindex)


class chebyshevt_root(Function):
    r"""
    chebyshev_root(n, k) returns the kth root (indexed from zero) of
    the nth Chebyshev polynomial of the first kind; that is, if
    0 <= k < n, chebyshevt(n, chebyshevt_root(n, k)) == 0.

    Examples
    ========

    >>> chebyshevt_root(3, 2)
    -sqrt(3)/2
    >>> chebyshevt(3, chebyshevt_root(3, 2))
    0

    See Also
    ========

    jacobi, gegenbauer,
    chebyshevt, chebyshevu, chebyshevu_root,
    legendre, assoc_legendre,
    hermite,
    laguerre, assoc_laguerre,
    diofant.polys.orthopolys.jacobi_poly
    diofant.polys.orthopolys.gegenbauer_poly
    diofant.polys.orthopolys.chebyshevt_poly
    diofant.polys.orthopolys.chebyshevu_poly
    diofant.polys.orthopolys.hermite_poly
    diofant.polys.orthopolys.legendre_poly
    diofant.polys.orthopolys.laguerre_poly

    """

    @classmethod
    def eval(cls, n, k):
        if not 0 <= k < n:
            raise ValueError(f'must have 0 <= k < n, got k = {k} and n = {n}')
        return cos(pi*(2*k + 1)/(2*n))


class chebyshevu_root(Function):
    r"""
    chebyshevu_root(n, k) returns the kth root (indexed from zero) of the
    nth Chebyshev polynomial of the second kind; that is, if 0 <= k < n,
    chebyshevu(n, chebyshevu_root(n, k)) == 0.

    Examples
    ========

    >>> chebyshevu_root(3, 2)
    -sqrt(2)/2
    >>> chebyshevu(3, chebyshevu_root(3, 2))
    0

    See Also
    ========

    chebyshevt, chebyshevt_root, chebyshevu,
    legendre, assoc_legendre,
    hermite,
    laguerre, assoc_laguerre,
    diofant.polys.orthopolys.jacobi_poly
    diofant.polys.orthopolys.gegenbauer_poly
    diofant.polys.orthopolys.chebyshevt_poly
    diofant.polys.orthopolys.chebyshevu_poly
    diofant.polys.orthopolys.hermite_poly
    diofant.polys.orthopolys.legendre_poly
    diofant.polys.orthopolys.laguerre_poly

    """

    @classmethod
    def eval(cls, n, k):
        if not 0 <= k < n:
            raise ValueError(f'must have 0 <= k < n, got k = {k} and n = {n}')
        return cos(pi*(k + 1)/(n + 1))

############################################################################
# Legendre polynomials and Associated Legendre polynomials
#


class legendre(OrthogonalPolynomial):
    r"""
    legendre(n, x) gives the nth Legendre polynomial of x, `P_n(x)`

    The Legendre polynomials are orthogonal on [-1, 1] with respect to
    the constant weight 1. They satisfy `P_n(1) = 1` for all n; further,
    `P_n` is odd for odd n and even for even n.

    Examples
    ========

    >>> legendre(0, x)
    1
    >>> legendre(1, x)
    x
    >>> legendre(2, x)
    3*x**2/2 - 1/2
    >>> legendre(n, x)
    legendre(n, x)
    >>> diff(legendre(n, x), x)
    n*(x*legendre(n, x) - legendre(n - 1, x))/(x**2 - 1)

    See Also
    ========

    jacobi, gegenbauer,
    chebyshevt, chebyshevt_root, chebyshevu, chebyshevu_root,
    assoc_legendre,
    hermite,
    laguerre, assoc_laguerre,
    diofant.polys.orthopolys.jacobi_poly
    diofant.polys.orthopolys.gegenbauer_poly
    diofant.polys.orthopolys.chebyshevt_poly
    diofant.polys.orthopolys.chebyshevu_poly
    diofant.polys.orthopolys.hermite_poly
    diofant.polys.orthopolys.legendre_poly
    diofant.polys.orthopolys.laguerre_poly

    References
    ==========

    * https://en.wikipedia.org/wiki/Legendre_polynomial
    * https://mathworld.wolfram.com/LegendrePolynomial.html
    * http://functions.wolfram.com/Polynomials/LegendreP/
    * http://functions.wolfram.com/Polynomials/LegendreP2/

    """

    _ortho_poly = staticmethod(legendre_poly)

    @classmethod
    def eval(cls, n, x):
        if not n.is_Number:
            # Symbolic result L_n(x)
            # L_n(-x)  --->  (-1)**n * L_n(x)
            if x.could_extract_minus_sign():
                return (-1)**n * legendre(n, -x)
            # L_{-n}(x)  --->  L_{n-1}(x)
            if n.could_extract_minus_sign():
                return legendre(-n - Integer(1), x)
            # We can evaluate for some special values of x
            if x == 0:
                return sqrt(pi)/(gamma((1 - n)/2)*gamma(Integer(1) + n/2))
            if x == 1:
                return Integer(1)
            if x == oo:
                return oo
        else:
            # n is a given fixed integer, evaluate into polynomial
            if n.is_negative:
                raise ValueError('The index n must be nonnegative '
                                 f'integer (got {n!r})')
            return cls._eval_at_order(n, x)

    def fdiff(self, argindex=2):
        if argindex == 2:
            # Diff wrt x
            # Find better formula, this is unsuitable for x = 1
            n, x = self.args
            return n/(x**2 - 1)*(x*legendre(n, x) - legendre(n - 1, x))
        # wrt n
        raise ArgumentIndexError(self, argindex)


class assoc_legendre(Function):
    r"""
    assoc_legendre(n,m, x) gives `P_n^m(x)`, where n and m are
    the degree and order or an expression which is related to the nth
    order Legendre polynomial, `P_n(x)` in the following manner:

    .. math::
        P_n^m(x) = (-1)^m (1 - x^2)^{\frac{m}{2}}
                   \frac{\mathrm{d}^m P_n(x)}{\mathrm{d} x^m}

    Associated Legendre polynomial are orthogonal on [-1, 1] with:

    - weight = 1            for the same m, and different n.
    - weight = 1/(1-x**2)   for the same n, and different m.

    Examples
    ========

    >>> assoc_legendre(0, 0, x)
    1
    >>> assoc_legendre(1, 0, x)
    x
    >>> assoc_legendre(1, 1, x)
    -sqrt(-x**2 + 1)
    >>> assoc_legendre(n, m, x)
    assoc_legendre(n, m, x)

    See Also
    ========

    jacobi, gegenbauer,
    chebyshevt, chebyshevt_root, chebyshevu, chebyshevu_root,
    legendre,
    hermite,
    laguerre, assoc_laguerre,
    diofant.polys.orthopolys.jacobi_poly
    diofant.polys.orthopolys.gegenbauer_poly
    diofant.polys.orthopolys.chebyshevt_poly
    diofant.polys.orthopolys.chebyshevu_poly
    diofant.polys.orthopolys.hermite_poly
    diofant.polys.orthopolys.legendre_poly
    diofant.polys.orthopolys.laguerre_poly

    References
    ==========

    * https://en.wikipedia.org/wiki/Associated_Legendre_polynomials
    * https://mathworld.wolfram.com/LegendrePolynomial.html
    * http://functions.wolfram.com/Polynomials/LegendreP/
    * http://functions.wolfram.com/Polynomials/LegendreP2/

    """

    @classmethod
    def _eval_at_order(cls, n, m):
        P = legendre_poly(n, _x, polys=True).diff((_x, m))
        return (-1)**m * (1 - _x**2)**Rational(m, 2) * P.as_expr()

    @classmethod
    def eval(cls, n, m, x):
        if m.could_extract_minus_sign():
            # P^{-m}_n  --->  F * P^m_n
            return (-1)**(-m) * (factorial(m + n)/factorial(n - m)) * assoc_legendre(n, -m, x)
        if m == 0:
            # P^0_n  --->  L_n
            return legendre(n, x)
        if x == 0:
            return 2**m*sqrt(pi) / (gamma((1 - m - n)/2)*gamma(1 - (m - n)/2))
        if n.is_Number and m.is_Number and n.is_integer and m.is_integer:
            if n.is_negative:
                raise ValueError(f'{cls} : 1st index must be nonnegative integer (got {n!r})')
            if abs(m) > n:
                raise ValueError(f"{cls} : abs('2nd index') must be <= '1st index' (got {n!r}, {m!r})")
            return cls._eval_at_order(int(n), abs(int(m))).subs({_x: x})

    def fdiff(self, argindex=3):
        if argindex == 3:
            # Diff wrt x
            # Find better formula, this is unsuitable for x = 1
            n, m, x = self.args
            return 1/(x**2 - 1)*(x*n*assoc_legendre(n, m, x) - (m + n)*assoc_legendre(n - 1, m, x))
        # wrt n, m
        raise ArgumentIndexError(self, argindex)

    def _eval_conjugate(self):
        n, m, x = self.args
        return self.func(n, m.conjugate(), x.conjugate())

############################################################################
# Hermite polynomials
#


class hermite(OrthogonalPolynomial):
    r"""
    hermite(n, x) gives the nth Hermite polynomial in x, `H_n(x)`

    The Hermite polynomials are orthogonal on `(-\infty, \infty)`
    with respect to the weight `\exp\left(-x^2\right)`.

    Examples
    ========

    >>> hermite(0, x)
    1
    >>> hermite(1, x)
    2*x
    >>> hermite(2, x)
    4*x**2 - 2
    >>> hermite(n, x)
    hermite(n, x)
    >>> diff(hermite(n, x), x)
    2*n*hermite(n - 1, x)
    >>> hermite(n, -x)
    (-1)**n*hermite(n, x)

    See Also
    ========

    jacobi, gegenbauer,
    chebyshevt, chebyshevt_root, chebyshevu, chebyshevu_root,
    legendre, assoc_legendre,
    laguerre, assoc_laguerre,
    diofant.polys.orthopolys.jacobi_poly
    diofant.polys.orthopolys.gegenbauer_poly
    diofant.polys.orthopolys.chebyshevt_poly
    diofant.polys.orthopolys.chebyshevu_poly
    diofant.polys.orthopolys.hermite_poly
    diofant.polys.orthopolys.legendre_poly
    diofant.polys.orthopolys.laguerre_poly

    References
    ==========

    * https://en.wikipedia.org/wiki/Hermite_polynomial
    * https://mathworld.wolfram.com/HermitePolynomial.html
    * http://functions.wolfram.com/Polynomials/HermiteH/

    """

    _ortho_poly = staticmethod(hermite_poly)

    @classmethod
    def eval(cls, n, x):
        if not n.is_Number:
            # Symbolic result H_n(x)
            # H_n(-x)  --->  (-1)**n * H_n(x)
            if x.could_extract_minus_sign():
                return (-1)**n * hermite(n, -x)
            # We can evaluate for some special values of x
            if x == 0:
                return 2**n * sqrt(pi) / gamma((1 - n)/2)
            if x == oo:
                return oo
        else:
            # n is a given fixed integer, evaluate into polynomial
            if n.is_negative:
                raise ValueError('The index n must be nonnegative '
                                 f'integer (got {n!r})')
            return cls._eval_at_order(n, x)

    def fdiff(self, argindex=2):
        if argindex == 2:
            # Diff wrt x
            n, x = self.args
            return 2*n*hermite(n - 1, x)
        # wrt n
        raise ArgumentIndexError(self, argindex)

############################################################################
# Laguerre polynomials
#


class laguerre(OrthogonalPolynomial):
    r"""
    Returns the nth Laguerre polynomial in x, `L_n(x)`.

    Parameters
    ==========

    n : int
        Degree of Laguerre polynomial. Must be ``n >= 0``.

    Examples
    ========

    >>> laguerre(0, x)
    1
    >>> laguerre(1, x)
    -x + 1
    >>> laguerre(2, x)
    x**2/2 - 2*x + 1
    >>> laguerre(3, x)
    -x**3/6 + 3*x**2/2 - 3*x + 1

    >>> laguerre(n, x)
    laguerre(n, x)

    >>> diff(laguerre(n, x), x)
    -assoc_laguerre(n - 1, 1, x)

    See Also
    ========

    jacobi, gegenbauer,
    chebyshevt, chebyshevt_root, chebyshevu, chebyshevu_root,
    legendre, assoc_legendre,
    hermite,
    assoc_laguerre,
    diofant.polys.orthopolys.jacobi_poly
    diofant.polys.orthopolys.gegenbauer_poly
    diofant.polys.orthopolys.chebyshevt_poly
    diofant.polys.orthopolys.chebyshevu_poly
    diofant.polys.orthopolys.hermite_poly
    diofant.polys.orthopolys.legendre_poly
    diofant.polys.orthopolys.laguerre_poly

    References
    ==========

    * https://en.wikipedia.org/wiki/Laguerre_polynomial
    * https://mathworld.wolfram.com/LaguerrePolynomial.html
    * http://functions.wolfram.com/Polynomials/LaguerreL/
    * http://functions.wolfram.com/Polynomials/LaguerreL3/

    """

    _ortho_poly = staticmethod(laguerre_poly)

    @classmethod
    def eval(cls, n, x):
        if not n.is_Number:
            # Symbolic result L_n(x)
            # L_{n}(-x)  --->  exp(-x) * L_{-n-1}(x)
            # L_{-n}(x)  --->  exp(x) * L_{n-1}(-x)
            if n.could_extract_minus_sign():
                return exp(x) * laguerre(-n - 1, -x)
            # We can evaluate for some special values of x
            if x == 0:
                return Integer(1)
            if x == -oo:
                return oo
            if x == oo:
                return (-1)**n * oo
        else:
            # n is a given fixed integer, evaluate into polynomial
            if n.is_negative:
                raise ValueError('The index n must be nonnegative '
                                 f'integer (got {n!r})')
            return cls._eval_at_order(n, x)

    def fdiff(self, argindex=2):
        if argindex == 2:
            # Diff wrt x
            n, x = self.args
            return -assoc_laguerre(n - 1, 1, x)
        # wrt n
        raise ArgumentIndexError(self, argindex)


class assoc_laguerre(OrthogonalPolynomial):
    r"""
    Returns the nth generalized Laguerre polynomial in x, `L_n(x)`.

    Parameters
    ==========

    n : int
        Degree of Laguerre polynomial. Must be ``n >= 0``.

    alpha : Expr
        Arbitrary expression. For ``alpha=0`` regular Laguerre
        polynomials will be generated.

    Examples
    ========

    >>> assoc_laguerre(0, a, x)
    1
    >>> assoc_laguerre(1, a, x)
    a - x + 1
    >>> assoc_laguerre(2, a, x)
    a**2/2 + 3*a/2 + x**2/2 + x*(-a - 2) + 1
    >>> assoc_laguerre(3, a, x)
    a**3/6 + a**2 + 11*a/6 - x**3/6 + x**2*(a/2 + 3/2) +
        x*(-a**2/2 - 5*a/2 - 3) + 1

    >>> assoc_laguerre(n, a, 0)
    binomial(a + n, a)

    >>> assoc_laguerre(n, a, x)
    assoc_laguerre(n, a, x)

    >>> assoc_laguerre(n, 0, x)
    laguerre(n, x)

    >>> diff(assoc_laguerre(n, a, x), x)
    -assoc_laguerre(n - 1, a + 1, x)

    >>> diff(assoc_laguerre(n, a, x), a)
    Sum(assoc_laguerre(_k, a, x)/(-a + n), (_k, 0, n - 1))

    See Also
    ========

    jacobi, gegenbauer,
    chebyshevt, chebyshevt_root, chebyshevu, chebyshevu_root,
    legendre, assoc_legendre,
    hermite,
    laguerre,
    diofant.polys.orthopolys.jacobi_poly
    diofant.polys.orthopolys.gegenbauer_poly
    diofant.polys.orthopolys.chebyshevt_poly
    diofant.polys.orthopolys.chebyshevu_poly
    diofant.polys.orthopolys.hermite_poly
    diofant.polys.orthopolys.legendre_poly
    diofant.polys.orthopolys.laguerre_poly

    References
    ==========

    * https://mathworld.wolfram.com/AssociatedLaguerrePolynomial.html
    * http://functions.wolfram.com/Polynomials/LaguerreL/
    * http://functions.wolfram.com/Polynomials/LaguerreL3/

    """

    @classmethod
    def _eval_at_order(cls, n, x, alpha):
        if n.is_integer and n >= 0:
            return laguerre_poly(int(n), _x, alpha).subs({_x: x})

    @classmethod
    def eval(cls, n, alpha, x):
        # L_{n}^{0}(x)  --->  L_{n}(x)
        if alpha == 0:
            return laguerre(n, x)

        if not n.is_Number:
            # We can evaluate for some special values of x
            if x == 0:
                return binomial(n + alpha, alpha)
            if x == oo and n.is_positive:
                return (-1)**n * oo
            if x == -oo and n.is_positive:
                return oo
        else:
            # n is a given fixed integer, evaluate into polynomial
            if n.is_negative:
                raise ValueError('The index n must be nonnegative'
                                 f' integer (got {n!r})')
            return cls._eval_at_order(n, x, alpha)

    def fdiff(self, argindex=3):
        from ...concrete import Sum
        if argindex == 2:
            # Diff wrt alpha
            n, alpha, x = self.args
            k = Dummy('k')
            return Sum(assoc_laguerre(k, alpha, x) / (n - alpha), (k, 0, n - 1))
        if argindex == 3:
            # Diff wrt x
            n, alpha, x = self.args
            return -assoc_laguerre(n - 1, alpha + 1, x)
        # wrt n
        raise ArgumentIndexError(self, argindex)

    def _eval_conjugate(self):
        n, alpha, x = self.args
        return self.func(n, alpha.conjugate(), x.conjugate())
