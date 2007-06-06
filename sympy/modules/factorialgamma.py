from sympy.core.functions import Function, exp, sqrt
from sympy.core.numbers import Number, Real, Rational, pi, I, oo
from sympy.core import Symbol, Add, Mul
from sympy.modules.simplify import simplify
from sympy import Order

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
        >>> from sympy.modules.factorialgamma import *
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
        >>> from sympy.modules.factorialgamma import *
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
        Calculate the quotient p!/q!, simplifying symbolically if possible

    Examples
    ========
        >>> from sympy import *
        >>> from sympy.modules.factorialgamma import *
        >>> factorial_quotient(pi+1, pi)
        1+pi

    """
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
        >>> from sympy.modules.factorialgamma import *
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
