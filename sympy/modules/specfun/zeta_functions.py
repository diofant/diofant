from sympy.core.functions import Function, log
from sympy.core.numbers import Rational, pi, oo
from factorials import factorial, binomial, Function2

# first some utilities for calculating Bernoulli numbers

class _memo:
    def __init__(self, f):
        self.func = f
        self.cache = {}
        self.highest = 0
    def __call__(self, m):
        if m < self.highest:
            return self.cache[m]
        else:
            x = self.func(m)
            self.cache[m] = x
            self.highest = m
            return x

def _product(low, high):
    p = 1
    for k in range(low, high+1):
        p *= k
    return p

def _bernoulli_sum(m, mtop):
    s = 0
    a = int(binomial(m+3, m-6))
    for j in range(1, mtop+1):
        s += a*bernoulli(m-6*j)
        a *= _product(m-6 - 6*j + 1, m-6*j)
        a //= _product(6*j+4, 6*j+9)
    return s

@_memo
def _b0mod6(m):
    return (Rational(m+3, 3) - _bernoulli_sum(m, m//6)) / binomial(m+3, m)

@_memo
def _b2mod6(m):
    return (Rational(m+3, 3) - _bernoulli_sum(m, (m-2)//6)) / binomial(m+3, m)

@_memo
def _b4mod6(m):
    return (-Rational(m+3, 6) - _bernoulli_sum(m, (m-4)//6)) / binomial(m+3, m)

class bernoulli(Function):
    """
    Usage
    =====
        bernoulli(n) -> nth Bernoulli number, B_n

    Notes
    =====
        * Bernoulli numbers are rational numbers
        * For odd integers n > 1, bernoulli(n) = 0

    Examples
    ========
        >>> from sympy.modules.specfun.zeta_functions import *
        >>> [bernoulli(n) for n in range(11)]
        [1, -1/2, 1/6, 0, -1/30, 0, 1/42, 0, -1/30, 0, 5/66]
        >>> bernoulli(1000001)
        0

    """
    def eval(self):
        m = self._args
        if not (m.is_integer and m >= 0):
            return self
        m = int(m)
        if m == 0: return 1
        if m == 1: return -Rational(1,2)
        if m % 6 == 0: return _b0mod6(m)
        if m % 6 == 2: return _b2mod6(m)
        if m % 6 == 4: return _b4mod6(m)
        return 0


# TODO: speed up
class bernoulli_poly(Function2):
    """
    bernoulli_poly(n, x) - nth Bernoulli polynomial of x
    """
    def eval(self):
        n, x = self._args

        if isinstance(n, Rational) and n.is_integer:
            s = 0
            for k in range(n+1):
                s += binomial(n,k)*bernoulli(k)*x**(n-k)
            return s

        return self


class zeta(Function):
    """
    Usage
    =====
        zeta(s) -> Riemann zeta function of s

    Notes
    =====
        * The zeta function has a pole at s = 1
        * For even positive integers n, zeta(n) is a rational number
          times a power of pi
        * For nonpositive integers n, zeta(n) is a rational number

    Examples
    ========
        >>> from sympy.modules.specfun.zeta_functions import *
        >>> zeta(4)
        1/90*pi**4
        >>> zeta(5)  # no closed form known
        zeta(5)

    """
    def eval(self):
        s = self._args
        if s.is_integer:
            if s == 0:
                return Rational(-1,2)
            if s == 1:
                return oo
            if s > 1 and int(s) % 2 == 0:
                return abs(bernoulli(s)) * 2**(s-1) / factorial(s) * pi**s
            if s < 1:
                return -bernoulli(-s+1)/(-s+1)
        return self


class dirichlet_eta(Function):
    """
    Dirichlet eta function
    """
    def eval(self):
        s = self._args
        if s == 1:
            return log(2)
        else:
            return (1-2**(1-s)) * zeta(s)
