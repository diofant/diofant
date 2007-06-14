from sympy.core.symbol import Symbol
from sympy.core.numbers import Rational, Real
from sympy.modules.simplify import simplify
from factorials import Function2
import decimal


# Simple implementation of Newton's method for root-finding
def _newton(h, x0, eps):
    x = x0
    prevdiff = 1
    while 1:
        new = x - h(x)
        diff = abs(x - new)
        if diff <= eps:
            break
        prevdiff = diff
        x = new
    return new


class legendre(Function2):
    """
    Usage
    =====
        legendre(n, x) - nth Legendre polynomial of x, P_n(x)

    Notes
    =====
        The Legendre polynomials are orthogonal on [-1, 1]
        For all n, P_n(1) = 1
        P_n is odd for odd n and even for even n

    Examples
    ========

        >>> x = Symbol('x')
        >>> legendre(3, x)
        -3/2*x+5/2*x**3


    See also
    ========

       External links
       --------------
         U{Wikipedia: Legendre polynomial<http://en.wikipedia.org/wiki/Legendre_polynomial>}
    
    """

    _x = Symbol('x')
    _memo = {0:Rational(1), 1:_x}

    def poly(self):
        n, x = self._args
        assert n.is_integer and n >= 0
        n = int(n)
        m = len(self._memo)
        if n < m:
            return self._memo[n]
        else:
            for i in range(m, n+1):
                L = ((2*i-1)*self._x*self._memo[i-1] - (i-1)*self._memo[i-2])/i
                L = simplify(L)
                self._memo[i] = (L)
            return self._memo[n]

    def eval(self):
        n, x = self._args
        if isinstance(x, legendre_zero) and x._args[0] == n:
            return 0
        elif n.is_integer and n >= 0:
            return self.poly().subs(self._x, x)
        return self


class legendre_zero(Function2):
    """
    Usage
    =====

        legendre_zero(n, k) represents the kth zero (counting from zero)
        of the nth Legendre polynomial; that is, if 0 <= k < n,
        legendre(n, legendre_zero(n, k)) == 0.

        All zeros for a given Legendre polynomial are located symmetrically
        around 0 in the open interval (-1, 1). The zeros are indexed from
        left to right.

    Examples
    ========

        >>> legendre(5, legendre_zero(5, 3)) == 0
        True

    """

    def evalf(self, prec=10):

        # Increasing the precision is really just a matter of using
        # a lower epsilon; the problem is that numerical evaluation of
        # polynomials currently doesn't work as it should
        if prec > 10:
            raise NotImplementedError
        eps = 1e-10

        n, k = self._args
        assert 0 <= k < n

        L = lambda x: legendre(n, x)

        t = Symbol('t')
        Ldpol = legendre(n, t).diff(t)
        Ld = lambda x: Ldpol.subs(t, x)

        # Good initial estimate for use with Newton's method
        import math
        x = -math.cos(math.pi*(k+1-0.25)/(n+0.5))

        return _newton(lambda t: L(t)/Ld(t), x, eps)
