from sympy.core.functions import Function, exp, sqrt
from sympy.core.numbers import Number, Real, Rational, pi, I, oo
from sympy.core import Symbol, Add, Mul, Basic
from sympy.modules.simplify import simplify
from sympy import Order
from sympy.modules.trigonometric import sin

# Factorial and gamma related functions


class factorial(Function):
    """
    Usage
    =====
        factorial(x) -> Returns the factorial of x, defined as
        x! = 1*2*3*...*x if x is a positive integer

    Notes
    =====
        factorial(x) is evaluated explicitly if x is an integer or half
        an integer. If x is a negative integer, the value is infinite.

    Examples
    ========
        >>> from sympy import *
        >>> from sympy.modules.specfun.factorials import *
        >>> factorial(5)
        120
        >>> factorial(0)
        1
        >>> factorial(Rational(5,2))
        15/8*pi**(1/2)

    """
    def eval(self):
        x = self._args
        if x.is_integer:
            if x < 0:
                return oo
            y = Rational(1)
            for m in range(1, x+1):
                y *= m
            return y
        elif isinstance(x, Rational) and x.q == 2:
            n = (x.p + 1) / 2
            if n < 0:
                return (-1)**(-n+1) * pi * x / factorial(-x)
            return sqrt(pi) * Rational(1, 2**n) * factorial2(2*n-1)
        return self

    # This should give a series expansion around x = oo. Needs fixing
    def _series(self, x, n):
        return sqrt(2*pi*x) * x**x * exp(-x) * (1 + Order(1/x))

    def __latex__(self):
        x = self._args
        if (isinstance(x, Rational) and x.is_integer and x >= 0) or \
            isinstance(x, Symbol):
            s = x.__latex__()
        else:
            s = "(" + x.__latex__() + ")"
        return s + "!"


class factorial2(Function):
    """
    Usage
    =====
        factorial2(x) -> Returns the double factorial of x, defined as
        x!! = 2*4*6*...*x if x is a positive even integer and as
        x!! = 1*3*5*...*x if x is a positive odd integer.

    Notes
    =====
        Also defined for negative odd integers, but infinite for
        negative even integers.

    Examples
    ========
        >>> from sympy import *
        >>> from sympy.modules.specfun.factorials import *
        >>> factorial2(5)
        15
        >>> factorial2(6)
        48

    """
    def eval(self):
        x = self._args
        if x.is_integer:
            if int(x) % 2 == 0:
                if x < 0:
                    return oo
                else:
                    return 2**(x/2) * factorial(x/2)
            else:
                if x < 0:
                    return factorial2(x+2) / (x+2)
                else:
                    return factorial(x) / 2**((x-1)/2) / factorial((x-1)/2)
        return self

    def __latex__(self):
        x = self._args
        if (isinstance(x, Rational) and x.is_integer and x >= 0) or \
            isinstance(x, Symbol):
            s = x.__latex__()
        else:
            s = "(" + x.__latex__() + ")"
        return s + "!!"


def factorial_quotient(p, q):
    """
    Usage
    =====
        Calculate the quotient p!/q!, simplifying symbolically if possible.
        If both p! and q! are infinite, the correct limit is returned

    Examples
    ========
        >>> from sympy import *
        >>> from sympy.modules.specfun.factorials import *
        >>> factorial_quotient(pi+1, pi)
        1+pi

    """
    p = Basic.sympify(p)
    q = Basic.sympify(q)
    delta = simplify(p-q)
    if delta == 0:
        return 1
    if delta.is_integer:
        t = Rational(1)
        if delta > 0:
            for k in range(1, delta+1):
                t *= (q+k)
        else:
            for k in range(1, -delta+1):
                t /= (p+k)
        return t
    return factorial(p) / factorial(q)

# This class is a temporary solution
class Function2(Function):
    def __init__(self, x, y):
        Basic.__init__(self, is_commutative=True)
        self._args = self.sympify(x), self.sympify(y)


class rising_factorial(Function2):
    """
    Usage
    =====
        Calculate the rising factorial (x)^(n) = x(x+1)...(x+n-1)
        as a quotient of factorials

    Examples
    ========
        >>> from sympy.modules.specfun.factorials import *
        >>> rising_factorial(3, 2)
        12

    """

    def eval(self):
        x, n = self._args
        return factorial_quotient(x+n-1, x-1)

    def __latex__(self):
        x, n = self._args
        return "{(%s)}^{(%s)}" % (x.__latex__(), n.__latex__())


class falling_factorial(Function2):
    """
    Usage
    =====
        Calculate the falling factorial (x)_(n) = x(x-1)...(x-n+1)
        as a quotient of factorials

    Examples
    ========
        >>> from sympy.modules.specfun.factorials import *
        >>> falling_factorial(5, 3)
        60

    """
    def eval(self):
        x, n = self._args
        return factorial_quotient(x, x-n)

    def __latex__(self):
        x, n = self._args
        return "{(%s)}_{(%s)}" % (x.__latex__(), n.__latex__())


class binomial(Function2):
    """
    Usage
    =====
        Calculate the binomial coefficient C(n,k) = n!/(k!(n-k)!)

    Notes
    =====
        When n and k are positive integers, the result is always
        a positive integer

        If k is a positive integer, the result is a polynomial in n
        that is evaluated explicitly.

    Examples
    ========
        >>> from sympy import *
        >>> from sympy.modules.specfun.factorials import *
        >>> binomial(15,8)
        6435
        >>> # Building Pascal's triangle
        >>> [binomial(0,k) for k in range(1)]
        [1]
        >>> [binomial(1,k) for k in range(2)]
        [1, 1]
        >>> [binomial(2,k) for k in range(3)]
        [1, 2, 1]
        >>> [binomial(3,k) for k in range(4)]
        [1, 3, 3, 1]
        >>> # n can be arbitrary if k is a positive integer
        >>> binomial(Rational(5,4), 3)
        -5/128
        >>> x = Symbol('x')
        >>> binomial(x, 3)
        1/6*x*(-2+x)*(-1+x)

    """
    def eval(self):
        n, k = self._args
        if k == 0 or n == k:
            return Rational(1)
        if n.is_integer and k.is_integer:
            if n >= 0 and (k < 0 or (n-k) < 0):
                return Rational(0)
            # Todo: better support for negative integer arguments:
            # handle factorial poles that cancel
            pass
        if n.is_integer and k.is_integer and n >= 0 and k >= 0:
            # Choose the faster way to do the calculation
            if k > n-k:
                return factorial_quotient(n, k) / factorial(n-k)
            else:
                return factorial_quotient(n, n-k) / factorial(k)
        if not n.is_integer and k.is_integer and k >= 0:
            return factorial_quotient(n, n-k) / factorial(k)
        if n == 0:
            return sin(pi*k)/(pi*k)
        # Probably no simplification possible
        return self

    def __latex__(self):
        n, k = self._args
        return r"{{%s}\choose{%s}}" % (n.__latex__(), k.__latex__())

class gamma(Function):
    """
    Usage
    =====
        gamma(x) -> calculate the gamma function of x

    Notes
    =====
        gamma(x) = (x-1)!

        When x is an integer or half-integer, the result is automatically
        simplified to the corresponding factorial

    Examples
    ========
        >>> from sympy import *
        >>> from sympy.modules.specfun.factorials import *
        >>> gamma(3)
        2
        >>> gamma(Rational(1,2))
        pi**(1/2)

    """

    def eval(self):
        x = self._args
        y = factorial(x-1)
        if isinstance(y, factorial):
            return self
        else:
            return y

    def __latex__(self):
        return "\Gamma(" + self._args.__latex__() + ")"

