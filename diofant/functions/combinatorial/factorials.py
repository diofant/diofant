from functools import reduce
from math import sqrt as _sqrt

from ...core import Dummy, E, Function, Integer, S, cacheit, oo, sympify, zoo
from ...core.function import ArgumentIndexError
from ...ntheory import sieve


class CombinatorialFunction(Function):
    """Base class for combinatorial functions. """

    def _eval_simplify(self, ratio, measure):
        from ...simplify import combsimp
        expr = combsimp(self)
        if measure(expr) <= ratio*measure(self):
            return expr
        return self

###############################################################################
# ###################### FACTORIAL and MULTI-FACTORIAL ###################### #
###############################################################################


class factorial(CombinatorialFunction):
    """Implementation of factorial function over nonnegative integers.

    By convention (consistent with the gamma function and the binomial
    coefficients), factorial of a negative integer is complex infinity.

    The factorial is very important in combinatorics where it gives
    the number of ways in which `n` objects can be permuted. It also
    arises in calculus, probability, number theory, etc.

    There is strict relation of factorial with gamma function. In
    fact n! = gamma(n+1) for nonnegative integers. Rewrite of this
    kind is very useful in case of combinatorial simplification.

    Computation of the factorial is done using two algorithms. For
    small arguments naive product is evaluated. However for bigger
    input algorithm Prime-Swing is used. It is the fastest algorithm
    known and computes n! via prime factorization of special class
    of numbers, called here the 'Swing Numbers'.

    Examples
    ========

    >>> factorial(0)
    1

    >>> factorial(7)
    5040

    >>> factorial(-2)
    zoo

    >>> factorial(n)
    factorial(n)

    >>> factorial(2*n)
    factorial(2*n)

    >>> factorial(Rational(1, 2))
    factorial(1/2)

    See Also
    ========

    diofant.functions.combinatorial.factorials.factorial2
    diofant.functions.combinatorial.factorials.RisingFactorial
    diofant.functions.combinatorial.factorials.FallingFactorial
    """

    def fdiff(self, argindex=1):
        from .. import gamma, polygamma
        if argindex == 1:
            return gamma(self.args[0] + 1)*polygamma(0, self.args[0] + 1)
        else:
            raise ArgumentIndexError(self, argindex)

    _small_swing = [
        1, 1, 1, 3, 3, 15, 5, 35, 35, 315, 63, 693, 231, 3003, 429, 6435, 6435, 109395,
        12155, 230945, 46189, 969969, 88179, 2028117, 676039, 16900975, 1300075,
        35102025, 5014575, 145422675, 9694845, 300540195, 300540195
    ]

    @classmethod
    def _swing(cls, n):
        if n < 33:
            return cls._small_swing[n]
        else:
            N, primes = int(_sqrt(n)), []

            for prime in sieve.primerange(3, N + 1):
                p, q = 1, n

                while True:
                    q //= prime

                    if q > 0:
                        if q & 1 == 1:
                            p *= prime
                    else:
                        break

                if p > 1:
                    primes.append(p)

            for prime in sieve.primerange(N + 1, n//3 + 1):
                if (n // prime) & 1 == 1:
                    primes.append(prime)

            L_product = R_product = 1

            for prime in sieve.primerange(n//2 + 1, n + 1):
                L_product *= prime

            for prime in primes:
                R_product *= prime

            return L_product*R_product

    @classmethod
    def _recursive(cls, n):
        if n < 2:
            return 1
        else:
            return (cls._recursive(n//2)**2)*cls._swing(n)

    @classmethod
    def eval(cls, n):
        n = sympify(n)

        if n.is_Number:
            if n is S.Zero:
                return S.One
            elif n is oo:
                return oo
            elif n.is_Integer:
                if n.is_negative:
                    return zoo
                else:
                    n, result = n.numerator, 1

                    if n < 20:
                        for i in range(2, n + 1):
                            result *= i
                    else:
                        N, bits = n, 0

                        while N != 0:
                            if N & 1 == 1:
                                bits += 1

                            N = N >> 1

                        result = cls._recursive(n)*2**(n - bits)

                    return Integer(result)

    def _eval_rewrite_as_gamma(self, n):
        from .. import gamma
        return gamma(n + 1)

    def _eval_rewrite_as_tractable(self, n):
        from .. import loggamma, exp
        return exp(loggamma(n + 1))

    def _eval_rewrite_as_Product(self, n):
        from ...concrete import Product
        if n.is_nonnegative and n.is_integer:
            i = Dummy('i', integer=True)
            return Product(i, (i, 1, n))

    def _eval_is_integer(self):
        n = self.args[0]
        if n.is_integer and n.is_nonnegative:
            return True

    def _eval_is_positive(self):
        n = self.args[0]
        if n.is_integer and n.is_nonnegative:
            return True

    def _eval_is_composite(self):
        n = self.args[0]
        if n.is_integer:
            return (n - 3).is_nonnegative

    def _eval_is_extended_real(self):
        n = self.args[0]
        if n.is_nonnegative or n.is_noninteger:
            return True


class MultiFactorial(CombinatorialFunction):
    pass


class subfactorial(CombinatorialFunction):
    r"""The subfactorial counts the derangements of n items and is
    defined for non-negative integers as::

              ,
             |  1                             for n = 0
        !n = {  0                             for n = 1
             |  (n - 1)*(!(n - 1) + !(n - 2)) for n > 1
              `

    It can also be written as int(round(n!/exp(1))) but the recursive
    definition with caching is implemented for this function.

    An interesting analytic expression is the following

    .. math:: !x = \Gamma(x + 1, -1)/e

    which is valid for non-negative integers x. The above formula
    is not very useful in case of non-integers. `\Gamma(x + 1, -1)` is
    single-valued only for integral arguments x, elsewhere on the positive real
    axis it has an infinite number of branches none of which are real.

    References
    ==========

    * https//en.wikipedia.org/wiki/Subfactorial
    * http://mathworld.wolfram.com/Subfactorial.html

    Examples
    ========

    >>> subfactorial(n + 1)
    subfactorial(n + 1)
    >>> subfactorial(5)
    44

    See Also
    ========

    diofant.functions.combinatorial.factorials.factorial,
    diofant.utilities.iterables.generate_derangements,
    diofant.functions.special.gamma_functions.uppergamma
    """

    @classmethod
    @cacheit
    def _eval(cls, n):
        if not n:
            return S.One
        elif n == 1:
            return S.Zero
        return (n - 1)*(cls._eval(n - 1) + cls._eval(n - 2))

    @classmethod
    def eval(cls, arg):
        if arg.is_Number:
            if arg.is_Integer and arg.is_nonnegative:
                return cls._eval(arg)
            elif arg is oo:
                return oo

    def _eval_is_even(self):
        n = self.args[0]
        if n.is_odd and n.is_nonnegative:
            return True

    def _eval_is_integer(self):
        n = self.args[0]
        if n.is_integer and n.is_nonnegative:
            return True

    def _eval_rewrite_as_uppergamma(self, n):
        from .. import uppergamma
        return uppergamma(n + 1, -1)/E

    def _eval_is_nonnegative(self):
        n = self.args[0]
        if n.is_integer and n.is_nonnegative:
            return True

    def _eval_is_odd(self):
        n = self.args[0]
        if n.is_even and n.is_nonnegative:
            return True


class factorial2(CombinatorialFunction):
    """The double factorial n!!, not to be confused with (n!)!

    The double factorial is defined for nonnegative integers and for odd
    negative integers as::

               ,
              |  n*(n - 2)*(n - 4)* ... * 1    for n positive odd
        n!! = {  n*(n - 2)*(n - 4)* ... * 2    for n positive even
              |  1                             for n = 0
              |  (n+2)!! / (n+2)               for n negative odd
               `

    References
    ==========

    * https://en.wikipedia.org/wiki/Double_factorial

    Examples
    ========

    >>> factorial2(n + 1)
    factorial2(n + 1)
    >>> factorial2(5)
    15
    >>> factorial2(-1)
    1
    >>> factorial2(-5)
    1/3

    See Also
    ========

    diofant.functions.combinatorial.factorials.factorial
    diofant.functions.combinatorial.factorials.RisingFactorial
    diofant.functions.combinatorial.factorials.FallingFactorial
    """

    @classmethod
    def eval(cls, n):
        # TODO: extend this to complex numbers?
        if n.is_Number:
            if n.is_infinite:
                return

            if n.is_negative:
                if n.is_odd:
                    return n * (S.NegativeOne) ** ((1 - n) / 2) / factorial2(-n)
                elif n.is_even:
                    raise ValueError("argument must be nonnegative or odd")
            else:
                if n.is_even:
                    k = n / 2
                    return 2 ** k * factorial(k)
                elif n.is_integer:
                    return factorial(n) / factorial2(n - 1)

    def _eval_is_even(self):
        n = self.args[0]
        if n.is_integer:
            if n.is_odd:
                return False
            if n.is_even:
                if n.is_positive:
                    return True
                if n.is_zero:
                    return False

    def _eval_is_integer(self):
        n = self.args[0]
        if n.is_integer:
            if (n + 1).is_nonnegative:
                return True
            if n.is_odd:
                return (n + 3).is_nonnegative

    def _eval_is_positive(self):
        n = self.args[0]
        if n.is_integer:
            if (n + 1).is_nonnegative:
                return True
            if n.is_odd:
                return ((n + 1)/2).is_even


###############################################################################
# ###################### RISING and FALLING FACTORIALS ###################### #
###############################################################################


class RisingFactorial(CombinatorialFunction):
    """Rising factorial (also called Pochhammer symbol) is a double valued
    function arising in concrete mathematics, hypergeometric functions
    and series expansions.

    It is defined by::

                   rf(x, k) = x * (x+1) * ... * (x + k-1)

    where 'x' can be arbitrary expression and 'k' is an integer. For
    more information check "Concrete mathematics" by Graham, pp. 66
    or visit http://mathworld.wolfram.com/RisingFactorial.html page.

    Examples
    ========

    >>> rf(x, 0)
    1

    >>> rf(1, 5)
    120

    >>> rf(x, 5) == x*(1 + x)*(2 + x)*(3 + x)*(4 + x)
    True

    See Also
    ========

    diofant.functions.combinatorial.factorials.factorial
    diofant.functions.combinatorial.factorials.factorial2
    diofant.functions.combinatorial.factorials.FallingFactorial
    """

    @classmethod
    def eval(cls, x, k):
        x = sympify(x)
        k = sympify(k)

        if x is S.One:
            return factorial(k)
        elif k.is_Integer:
            if k is S.Zero:
                return S.One
            else:
                if k.is_positive:
                    if x is oo:
                        return oo
                    elif x == -oo:
                        if k.is_odd:
                            return -oo
                        else:
                            return oo
                    else:
                        return reduce(lambda r, i: r*(x + i), range(int(k)), 1)
                else:
                    if x in (oo, -oo):
                        return oo
                    else:
                        return 1/reduce(lambda r, i: r*(x - i), range(1, abs(int(k)) + 1), 1)

    def _eval_rewrite_as_gamma(self, x, k):
        from .. import gamma
        return gamma(x + k) / gamma(x)

    def _eval_rewrite_as_tractable(self, x, k):
        return self._eval_rewrite_as_gamma(x, k).rewrite('tractable')

    def _eval_is_integer(self):
        x, k = self.args
        if x.is_integer and k.is_integer:
            if k.is_nonnegative:
                return True


class FallingFactorial(CombinatorialFunction):
    """Falling factorial (related to rising factorial) is a double valued
    function arising in concrete mathematics, hypergeometric functions
    and series expansions.

    It is defined by::

                   ff(x, k) = x * (x-1) * ... * (x - k+1)

    where 'x' can be arbitrary expression and 'k' is an integer. For
    more information check "Concrete mathematics" by Graham, pp. 66
    or visit http://mathworld.wolfram.com/FallingFactorial.html page.

    >>> ff(x, 0)
    1

    >>> ff(5, 5)
    120

    >>> ff(x, 5) == x*(x-1)*(x-2)*(x-3)*(x-4)
    True

    See Also
    ========

    diofant.functions.combinatorial.factorials.factorial
    diofant.functions.combinatorial.factorials.factorial2
    diofant.functions.combinatorial.factorials.RisingFactorial
    """

    @classmethod
    def eval(cls, x, k):
        x = sympify(x)
        k = sympify(k)

        if k.is_Integer:
            if k is S.Zero:
                return S.One
            else:
                if k.is_positive:
                    if x is oo:
                        return oo
                    elif x == -oo:
                        if k.is_odd:
                            return -oo
                        else:
                            return oo
                    else:
                        return reduce(lambda r, i: r*(x - i), range(int(k)), 1)
                else:
                    if x in (oo, -oo):
                        return oo
                    else:
                        return 1/reduce(lambda r, i: r*(x + i), range(1, abs(int(k)) + 1), 1)

    def _eval_rewrite_as_gamma(self, x, k):
        from .. import gamma
        return (-1)**k * gamma(-x + k) / gamma(-x)

    def _eval_is_integer(self):
        x, k = self.args
        if x.is_integer and k.is_integer:
            if k.is_nonnegative:
                return True


rf = RisingFactorial
ff = FallingFactorial

###############################################################################
# ######################### BINOMIAL COEFFICIENTS ########################### #
###############################################################################


class binomial(CombinatorialFunction):
    """Implementation of the binomial coefficient. It can be defined
    in two ways depending on its desired interpretation:

        C(n,k) = n!/(k!(n-k)!)   or   C(n, k) = ff(n, k)/k!

    First, in a strict combinatorial sense it defines the
    number of ways we can choose 'k' elements from a set of
    'n' elements. In this case both arguments are nonnegative
    integers and binomial is computed using an efficient
    algorithm based on prime factorization.

    The other definition is generalization for arbitrary 'n',
    however 'k' must also be nonnegative. This case is very
    useful when evaluating summations.

    For the sake of convenience for negative 'k' this function
    will return zero no matter what valued is the other argument.

    To expand the binomial when n is a symbol, use either
    expand_func() or expand(func=True). The former will keep the
    polynomial in factored form while the latter will expand the
    polynomial itself. See examples for details.

    Examples
    ========

    >>> n = Symbol('n', integer=True, positive=True)

    >>> binomial(15, 8)
    6435

    >>> binomial(n, -1)
    0

    Rows of Pascal's triangle can be generated with the binomial function:

    >>> for N in range(8):
    ...     [ binomial(N, i) for i in range(N + 1)]
    ...
    [1]
    [1, 1]
    [1, 2, 1]
    [1, 3, 3, 1]
    [1, 4, 6, 4, 1]
    [1, 5, 10, 10, 5, 1]
    [1, 6, 15, 20, 15, 6, 1]
    [1, 7, 21, 35, 35, 21, 7, 1]

    As can a given diagonal, e.g. the 4th diagonal:

    >>> N = -4
    >>> [binomial(N, i) for i in range(1 - N)]
    [1, -4, 10, -20, 35]

    >>> binomial(Rational(5, 4), 3)
    -5/128
    >>> binomial(Rational(-5, 4), 3)
    -195/128

    >>> binomial(n, 3)
    binomial(n, 3)

    >>> binomial(n, 3).expand(func=True)
    n**3/6 - n**2/2 + n/3

    >>> expand_func(binomial(n, 3))
    n*(n - 2)*(n - 1)/6

    """

    def fdiff(self, argindex=1):
        from .. import polygamma
        if argindex == 1:
            # http://functions.wolfram.com/GammaBetaErf/Binomial/20/01/01/
            n, k = self.args
            return binomial(n, k)*(polygamma(0, n + 1) -
                                   polygamma(0, n - k + 1))
        elif argindex == 2:
            # http://functions.wolfram.com/GammaBetaErf/Binomial/20/01/02/
            n, k = self.args
            return binomial(n, k)*(polygamma(0, n - k + 1) -
                                   polygamma(0, k + 1))
        else:
            raise ArgumentIndexError(self, argindex)

    @classmethod
    def _eval(cls, n, k):
        assert n.is_Number and k.is_Integer and k != 1 and n != k
        if n.is_Integer and n >= 0:
            n, k = int(n), int(k)

            if k > n // 2:
                k = n - k

            M, result = int(_sqrt(n)), 1

            for prime in sieve.primerange(2, n + 1):
                if prime > n - k:
                    result *= prime
                elif prime > n // 2:
                    continue
                elif prime > M:
                    if n % prime < k % prime:
                        result *= prime
                else:
                    N, K = n, k
                    exp = a = 0

                    while N > 0:
                        a = int((N % prime) < (K % prime + a))
                        N, K = N // prime, K // prime
                        exp = a + exp

                    if exp > 0:
                        result *= prime**exp
            return Integer(result)
        else:
            d = result = n - k + 1
            for i in range(2, k + 1):
                d += 1
                result *= d
                result /= i
            return result

    @classmethod
    def eval(cls, n, k):
        n, k = map(sympify, (n, k))
        d = n - k
        if d.is_zero or k.is_zero:
            return S.One
        elif d.is_zero is False:
            if (k - 1).is_zero:
                return n
            elif k.is_negative:
                return S.Zero
            elif n.is_integer and n.is_nonnegative and d.is_negative:
                return S.Zero
        if k.is_Integer and k > 0 and n.is_Number:
            return cls._eval(n, k)

    def _eval_expand_func(self, **hints):
        """
        Function to expand binomial(n,k) when n is positive integer
        """
        n = self.args[0]
        if n.is_Number:
            return self.func(*self.args)

        k = self.args[1]
        if k.is_Add and n in k.args:
            k = n - k

        if k.is_Integer:
            if k == 0:
                return S.One
            elif k > 0:
                n = self.args[0]
                result = n - k + 1
                for i in range(2, k + 1):
                    result *= n - k + i
                    result /= i
                return result

        return self.func(*self.args)

    def _eval_rewrite_as_factorial(self, n, k):
        return factorial(n)/(factorial(k)*factorial(n - k))

    def _eval_rewrite_as_gamma(self, n, k):
        from .. import gamma
        return gamma(n + 1)/(gamma(k + 1)*gamma(n - k + 1))

    def _eval_rewrite_as_tractable(self, n, k):
        return self._eval_rewrite_as_gamma(n, k).rewrite('tractable')

    def _eval_is_integer(self):
        n, k = self.args
        if n.is_integer and k.is_integer:
            return True
